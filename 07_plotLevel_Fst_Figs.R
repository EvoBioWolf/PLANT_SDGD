library(dplyr)
library(stringr)
library(tidyr)
library(ggplot2)
library(readr)

# species_factor <- plot_df_div$species.y
# saveRDS(species_factor, file = "Figure2_species_order_factor.rds")

species_order <- readRDS("Figure2_species_order_factor.rds")

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

df <- read_delim("Fst.all.txt", delim = ",", col_names = TRUE) %>%
    rename(sp = base_name) %>%
    filter(!str_detect(sp, "_2x")) %>%
    filter(!sp %in% c("AjuRep", "", "GerPra", "TriFla")) %>% 
    filter(criteria == "DP>=30") %>%
    mutate(species = str_replace(sp, "_[0-9]+x$", ""))

df_keep <- df %>%
    mutate(
        g1 = str_sub(group1, 1, 1),
        g2 = str_sub(group2, 1, 1),
        
        # unordered "pair type" ignoring direction (so LM and ML are the same)
        pair_type = paste0(pmin(g1, g2), pmax(g1, g2)),
        
        class = case_when(
            pair_type == "HH" ~ "H",
            pair_type %in% c("LL", "MM", "LM") ~ "L",
            TRUE ~ NA_character_
        )
    ) %>%
    filter(!is.na(class)) %>%
    mutate(class = factor(class, levels = c("L", "H")))

df_keep_nodup <- df_keep %>%
    mutate(
        gmin = pmin(group1, group2),
        gmax = pmax(group1, group2)
    ) %>%
    distinct(sp, gmin, gmax, cutoff, criteria, .keep_all = TRUE)

ggplot(df_keep_nodup, aes(x = class, y = Fst)) +
    # 1. ORIGINAL RAW DATA DOTS
    geom_point(
        position = position_jitter(width = 0.18, height = 0),
        alpha = 0.6,
        size = 2
    ) +
    # 2. NEW RED MEAN DOTS
    stat_summary(
        fun = mean,
        geom = "point",
        shape = 16,      # Shape 18 is a diamond, or use 16 for a circle
        size = 3,        # Make it larger than the regular points to stand out
        color = "red"
    ) +
    facet_wrap(~ species, scales = "free_y") +
    theme_bw() +
    labs(x = "Class (L vs H)", y = "Fst")

ggsave("Fst_dots_figureS3.png", height = 7, width = 10)


means_LH <- df_keep_nodup %>%
    group_by(species, class) %>%
    summarise(Fst_mean = mean(Fst, na.rm = TRUE), .groups = "drop") %>%
    complete(species, class, fill = list(Fst_mean = NA_real_)) %>%
    mutate(
        class = factor(class, levels = c("L", "H"))) %>%
    left_join(species_lookup, by = c("species" = "abbreviation"))

means_LH <- means_LH %>%
    mutate(species.y = factor(species.y, levels = levels(species_order)))

wide_LH <- means_LH %>%
    select(species.y, class, Fst_mean) %>%
    pivot_wider(names_from = class, values_from = Fst_mean)

pts_df <- wide_LH %>%
    pivot_longer(cols = c(L, H), names_to = "class", values_to = "Fst_mean") %>%
    filter(!is.na(Fst_mean)) %>%
    mutate(class = factor(class, levels = c("L", "H")))

# ---- plot ----
p_dumbbell <- ggplot() +
    geom_segment(
        data = wide_LH,
        aes(y = species.y, yend = species.y, x = L, xend = H),
        linewidth = 0.6
    ) +
    geom_point(
        data = pts_df,
        aes(y = species.y, x = Fst_mean, color = class),
        size = 3, alpha = 0.6
    ) +
    scale_color_manual(
        values = c(L = "#1f77b4", H = "#ff7f0e")
    ) +
    theme_bw() +
    labs(x = "Mean Fst", y = NULL, color = "Species diversity") +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          axis.text.y = element_text(face="italic"),
          legend.position = c(0.98, 0.6),          # inside: (x, y) in [0,1]
          legend.justification = c(1, 1),           # anchor top-right of legend box
          legend.direction = "horizontal",
          
          legend.background = element_rect(fill = "white", color = NA),
          legend.key = element_rect(fill = "white", color = NA),
          
          legend.title = element_text(size = 12),
          legend.text  = element_text(size = 11)) +
    guides(color = guide_legend(nrow = 1, 
                                title.position = "top", 
                                title.hjust = 0.5))

ggsave("Fst_dumbbell.png", height = 4, width = 5)

