# ------------------------------------------------------------
# Ne simulation plots: one PDF, one page per species,
# each page has 4 panels for t = 4, 8, 12, 16 generations
# ------------------------------------------------------------

library(tidyverse)
library(patchwork)

rm(list = ls())
# ---- full names --------------------------------------------------
species_lookup <- tibble::tribble(
    ~species,                 ~abbreviation,
    "Ajuga reptans",           "AjuRep",
    "Alopecurus pratensis",    "AloPra",
    "Arrhenatherum elatius",   "ArrEla",
    "Bellis perennis",         "BelPer",
    "Crepis biennis",          "CreBie",
    "Daucus carota",           "DauCar",
    "Galium mollugo",          "GalMol",
    "Geranium pratense",       "GerPra",
    "Knautia arvensis",        "KnaArv",
    "Lathyrus pratensis",      "LatPra",
    "Leucanthemum vulgare",    "LeuVul",
    "Lotus corniculatus",      "LotCor",
    "Luzula campestris",       "LuzCam",
    "Medicago x varia",        "MedVar",
    "Plantago lanceolata",     "PlaLan",
    "Plantago media",          "PlaMed",
    "Primula veris",           "PriVer",
    "Ranunculus acris",        "RanAcr",
    "Trisetum flavescens",     "TriFla",
    "Trifolium pratense",      "TriPra",
    "Veronica chamaedrys",     "VerCha",
    "Vicia cracca",            "VicCra"
)

# ---- Inputs --------------------------------------------------
ploidy_df <- read_delim("sp.ploiy.txt", delim = "\t")

df <- read_delim(file = "all.ExpHet.csv", delim = ",", col_names = FALSE)
names(df) <- c("sp", "group_name", "criteria", "cutoff", "Expected_H", "diversity")

df0 <- df %>%
    filter(criteria == "DP>=30", cutoff == 80) %>%
    mutate(species = str_replace(sp, "_[0-9]+x$", "")) %>%
    mutate(sp_plot = paste0(species, "_", group_name)) 

df1 <- df0 %>%
    left_join(ploidy_df, by = "species") %>%
    mutate(
        used_ploidy = case_when(
            str_detect(sp, "_") ~ str_sub(sp, 8, 8),
            TRUE ~ str_sub(ploidy, 1, 1)
        ),
        ploidy_num = as.numeric(str_extract(used_ploidy, "\\d+"))
    ) %>%
    filter(!sp %in% c("AjuRep", "", "GerPra", "TriFla")) %>% 
    filter(!str_detect(sp, "_2x"))

# ---- Drift curves settings -----------------------------------
Ne_grid <- c(20, 50, 100, 200, 500, 1000)
t_vals  <- c(4, 8, 12, 16)

# ---- Reference heterozygosity per species --------------------
ref_tbl <- df1 %>%
    group_by(sp) %>%
    summarise(
        ref_diversity = if (all(is.na(diversity))) NA_real_ else min(diversity, na.rm = TRUE),
        ref_Expected_H = if (all(is.na(diversity))) NA_real_ else
            mean(Expected_H[diversity == min(diversity, na.rm = TRUE)], na.rm = TRUE),
        .groups = "drop"
    ) %>%
    select(sp, ref_diversity, ref_Expected_H)

# ---- Observed x,y per row ------------------------------------
df_obs <- df1 %>%
    left_join(ref_tbl, by = "sp") %>%
    mutate(
        # c_denom: 2 for diploid, 4 for tetraploid (fallback 2)
        c_denom = case_when(
            ploidy_num == 4 ~ 4,
            ploidy_num == 2 ~ 2,
            TRUE ~ 2
        ),
        # heterozygosity ratio
        D = if_else(!is.na(ref_Expected_H) & ref_Expected_H > 0 & !is.na(Expected_H),
                    Expected_H / ref_Expected_H, NA_real_),
        # richness-based proxy for census ratio (as in your script)
        r_proxy = if_else(!is.na(diversity) & !is.na(ref_diversity) & diversity > 0,
                          ref_diversity / diversity, NA_real_),
        x = if_else(!is.na(r_proxy) & r_proxy > 0, r_proxy, NA_real_),
        y = if_else(!is.na(D) & D > 0, D, NA_real_),
        is_ref = !is.na(ref_diversity) & (diversity == ref_diversity)
    )

# ---- Output PDF ----------------------------------------------
pdf_out <- "heterozygosity_ratio_by_species_t4_8_12_16.pdf"
pdf(pdf_out, width = 15, height = 6)  # landscape
on.exit(dev.off(), add = TRUE)
unique(df_obs$sp) |> purrr::walk(function(sp_i) {
    
    d <- df_obs %>% filter(sp == sp_i)
    if (nrow(d) == 0) return(invisible(NULL))
    
    # --- NEW TITLE LOGIC ---
    # 1. Get the ploidy string (e.g., "2x", "4x")
    curr_ploidy <- as.character(unique(d$c_denom)[1])
    if(is.na(curr_ploidy)) curr_ploidy <- "?"
    
    # 2. Get the bare abbreviation (remove any suffixes if they exist)
    curr_abbr <- str_replace(sp_i, "_.*", "") 
    
    # 3. Look up full name
    full_name <- species_lookup$species[species_lookup$abbreviation == curr_abbr]
    
    # 4. Safety: If abbreviation not found in lookup, revert to abbreviation
    if(length(full_name) == 0) full_name <- curr_abbr
    
    # 5. Construct Final Title: "Ajuga reptans (2x)"
    final_title <- paste0(full_name, " (", curr_ploidy, "x)")
    # -----------------------
    
    # c_denom for this species
    c_i <- d %>% filter(!is.na(c_denom)) %>% slice(1) %>% pull(c_denom)
    if (length(c_i) == 0 || is.na(c_i)) c_i <- 2
    
    # r-range covering points
    r_min <- if (all(is.na(d$x))) 1/60 else max(1/60, min(d$x, na.rm = TRUE) * 0.8)
    r_max <- if (all(is.na(d$x))) 2.0  else max(d$x, na.rm = TRUE) * 1.2
    r_seq <- seq(r_min, r_max, length.out = 400)
    
    # reference vertical line
    x_ref <- d %>% filter(is_ref) %>% distinct(x) %>% pull()
    if (length(x_ref) == 0 || is.na(x_ref[1])) x_ref <- 1
    
    # mean points
    mean_pts <- d %>%
        filter(!is.na(diversity)) %>%
        group_by(diversity) %>%
        summarise(
            x = first(na.omit(x)),
            y = mean(y, na.rm = TRUE),
            .groups = "drop"
        ) %>%
        filter(!is.na(x), !is.na(y))
    
    # --- CALCULATE GLOBAL Y LIMITS FOR THIS SPECIES ---
    all_lines_range <- tidyr::crossing(Ne = Ne_grid, r = r_seq, t_val = t_vals) %>%
        mutate(
            slope = t_val / (c_i * Ne),
            D_th  = exp(slope * (r - 1) / r)
        ) %>%
        summarise(min_y = min(D_th, na.rm = TRUE), max_y = max(D_th, na.rm = TRUE))
    
    obs_range <- range(d$y, na.rm = TRUE)
    
    final_ymin <- min(all_lines_range$min_y, obs_range[1], na.rm = TRUE)
    final_ymax <- max(all_lines_range$max_y, obs_range[2], na.rm = TRUE)
    
    # Add 5% padding
    final_ylim <- c(final_ymin * 0.95, final_ymax * 1.05)
    # --------------------------------------------------
    # Build one panel for a given t
    make_plot_t <- function(t_val) {
        
        lines_df <- tidyr::crossing(Ne = Ne_grid, r = r_seq) %>%
            mutate(
                slope = t_val / (c_i * Ne),
                D_th  = exp(slope * (r - 1) / r),
                Ne_lab = factor(as.character(Ne), levels = as.character(Ne_grid))
            )
        
        ggplot() +
            geom_line(
                data = lines_df,
                aes(x = r, y = D_th, color = Ne_lab),
                linewidth = 0.9
            ) +
            geom_point(
                data = d,
                aes(x = x, y = y, shape = is_ref),
                size = 2.4
            ) +
            geom_point(
                data = mean_pts,
                aes(x = x, y = y),
                inherit.aes = FALSE,
                color = "red",
                size = 2.2,
                shape = 17
            ) +
            geom_vline(xintercept = x_ref[1], linetype = "dashed") +
            scale_shape_manual(values = c(`TRUE` = 17, `FALSE` = 16), guide = "none") +
            coord_cartesian(ylim = final_ylim) + 
            labs(
                x = "Expected population size ratio",
                y = "Heterozygosity ratio",
                color = "Ne",
                title = paste0("t = ", t_val, " generations")
            ) +
            theme_minimal(base_size = 12) +
            theme(
                plot.title = element_text(size = 12, face = "bold"),
                legend.position = "bottom"
            )
    }
    
    p_list <- lapply(t_vals, make_plot_t)
    
    combined <- wrap_plots(p_list, nrow = 1, guides = "collect") +
        plot_annotation(
            title = final_title,  # <--- Uses the new name
            theme = theme(
                plot.title = element_text(size = 16, face = "bold.italic"), # Italic for scientific names
                legend.position = "bottom"
            )
        )
    
    print(combined)
})
dev.off()

message("Done. Wrote: ", pdf_out)


