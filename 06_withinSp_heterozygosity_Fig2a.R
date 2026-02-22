library(tidyverse)
library(forcats)
library(ggplot2)

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


grp <- read_delim("forestPlot/sp_functionalGroup.txt", delim = "\t")
grp2 <- grp %>%
    rename(Functional_group = `Funtional group`) %>%  # fix typo if present
    mutate(
        species_base = species
    ) %>%
    select(species_base, Functional_group)

df <- read_delim("forestPlot/div_biomassResid_withinDiv.effect.txt", delim = "\t")
# columns: species, variable, Estimate, Std_Error, z_value, `Pr(>|z|)`

# 1) Normalize species names (drop trailing _2x/_4x/etc.)
df2 <- df %>%
    mutate(
        species_base = str_remove(species, "_\\d+x$")  # removes _2x, _4x, _6x, etc.
    )

# 2) Join functional group onto effects
plot_df <- df2 %>%
    left_join(grp2, by = "species_base") %>%
    mutate(
        Functional_group = coalesce(Functional_group, "Unknown"),
        lower = Estimate - Std_Error,            # ± SE bars
        upper = Estimate + Std_Error,
        species_label = species, 
        species_label = fct_reorder(species_label, Estimate)
    ) %>% 
    filter(!str_detect(species, "_2x")) %>% 
    filter(!species %in% c('AjuRep', '', 'GerPra', 'TriFla')) # species replicates failed

order_tbl <- plot_df %>%
    filter(variable == "diversity_log2") %>%
    arrange(Estimate) %>%
    distinct(species_label, .keep_all = TRUE)

plot_df <- plot_df %>%
    mutate(
        species_label = factor(species_label, levels = order_tbl$species_label),
        variable = fct_recode(variable, "biomass_residual" = "resid_species_div")
    )

ggplot(plot_df, aes(x = Estimate, y = species_label, color = Functional_group)) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0, na.rm = TRUE) +
    geom_point(size = 2) +
    facet_wrap(~ variable, scales = "free_x") +
    labs(
        x = "Estimate (± Std. Error)",
        y = NULL,
        color = "Functional group"
    ) +
    xlim(-0.05, 0.05) +
    theme_classic() +
    theme(
        legend.position = c(0.99, 0.45),
        legend.justification = c("right", "bottom"),
        legend.title = element_text(size = 15),
        legend.text  = element_text(size = 14),
        axis.title.x = element_text(size = 12),
        strip.text.x = element_text(size = 14),
        axis.text    = element_text(size = 10)
    )

###================###
### only diversity ### 
###================###
df  <- read_delim("forestPlot/div_biomassResid_withinDiv.effect.txt", delim = "\t")

plot_df_div <- df %>%
    mutate(species_base = str_remove(species, "_\\d+x$")) %>%
    filter(!str_detect(species, "_2x"),
           !species %in% c("AjuRep", "", "GerPra", "TriFla"),
           variable == "diversity_log2") %>%
    left_join(
        species_lookup,
        by = c("species_base" = "abbreviation")
    ) %>%
    mutate(
        lower = Estimate - Std_Error,
        upper = Estimate + Std_Error,
        species.y = fct_reorder(species.y, Estimate)
    )

ggplot(plot_df_div, aes(x = Estimate, y = species.y)) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0, na.rm = TRUE) +
    geom_point(size = 2) +
    labs(
        x = "Estimate (± Std. Error)",
        y = NULL
    ) +
    xlim(-0.05, 0.05) +
    theme_bw() +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          axis.text.y = element_text(face="italic"))

ggsave("forestPlot/diversity_log2_forest.png", height = 4, width = 5)


#####==============================###
##            only biomass         ###
#####==============================###
df  <- read_delim("forestPlot/div_biomassResid_withinDiv.effect.txt", delim = "\t")

plot_df_biomass <- df %>%
    mutate(species_base = str_remove(species, "_\\d+x$")) %>%
    filter(!str_detect(species, "_2x"),
           !species %in% c("AjuRep", "", "GerPra", "TriFla"),
           variable == "resid_species_div") %>%
    left_join(
        species_lookup,
        by = c("species_base" = "abbreviation")
    ) %>%
    mutate(
        lower = Estimate - Std_Error,
        upper = Estimate + Std_Error,
        species.y = factor(species.y, levels = plot_df_div$species.y)
    )

ggplot(plot_df_biomass, aes(x = Estimate, y = species.y)) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0, na.rm = TRUE) +
    geom_point(size = 1) +
    labs(
        x = "Estimate (± Std. Error)",
        y = NULL
    ) +
    xlim(-0.03, 0.15) +
    theme_bw() +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          axis.text.y = element_text(face="italic"))

ggsave("forestPlot/biomass_forest.png", height = 4, width = 5)

#####==============================###
###         only expression        ###
#####==============================###
df  <- read_delim("forestPlot/div_biomassResid_withinDiv.effect.txt", delim = "\t")

plot_df_expression <- df %>%
    mutate(species_base = str_remove(species, "_\\d+x$")) %>%
    filter(!str_detect(species, "_2x"),
           !species %in% c("AjuRep", "", "GerPra", "TriFla"),
           variable == "Expr_ln") %>%
    left_join(
        species_lookup,
        by = c("species_base" = "abbreviation")
    ) %>%
    mutate(
        lower = Estimate - Std_Error,
        upper = Estimate + Std_Error,
        species.y = factor(species.y, levels = plot_df_div$species.y)
    )

ggplot(plot_df_expression, aes(x = Estimate, y = species.y)) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0, na.rm = TRUE) +
    geom_point(size = 1) +
    labs(
        x = "Estimate (± Std. Error)",
        y = NULL
    ) +
    # xlim(-0.03, 0.15) +
    theme_bw() +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          axis.text.y = element_text(face="italic"))
ggsave("forestPlot/expression_forest.png", height = 4, width = 5)

species_factor <- plot_df_div$species.y
saveRDS(species_factor, file = "Figure2_species_order_factor.rds")
