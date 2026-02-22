library(metafor)
library(tidyverse)
library(forcats)
library(ggplot2)

# columns: species, variable, Estimate, Std_Error, z_value, `Pr(>|z|)`
df <- read_delim("forestPlot/div_biomassResid_withinDiv.effect.txt", delim = "\t")

meta_df <- df %>%
    filter(!str_detect(species, "_2x")) %>% 
    filter(!species %in% c('AjuRep', '', 'GerPra', 'TriFla'))
    mutate(
        species_base = str_remove(species, "_\\d+x$")  # removes _2x, _4x, _6x, etc.
    ) %>%
    filter(variable == "diversity_log2") %>% 
    mutate(
        yi = Estimate,
        sei = Std_Error
    )

# Random-effects meta-analysis
res <- rma(yi = yi, sei = sei, method = "REML", data = meta_df)

summary(res)
