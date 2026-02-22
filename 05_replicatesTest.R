library(tidyverse)
library(ggrepel)

rm(list=ls())

ploidy_df <- read_delim("sp.ploiy.txt", delim = "\t")

df <- read_delim(file = "replicates_randoms/replicates.genotypeCounts.tsv", delim = "\t") %>% 
    filter(criteria != "GQ20") %>% 
    mutate(non_ref = identical_nonref_exact + different_exact) %>% 
    mutate(non_ref_identical_rate = identical_nonref_exact / non_ref) %>% 
    mutate(type = case_when(
        grepl("rep", pairID) ~ "replicates",
        grepl("random", pairID) ~ "randoms",
        TRUE ~ "other"
    )) %>% 
    left_join(ploidy_df) %>% 
    mutate(sp_ploidy = paste0(species, "_", ploidy)) %>% 
    select(pairID, sample1, sample2, criteria, 
           identical_nonref_exact, different_exact, non_ref,
           non_ref_identical_rate, species, type, sp_ploidy)
df$sp_ploidy <- gsub("_NA$", "", df$sp_ploidy)

df2 <- df %>% 
    group_by(sp_ploidy, type, criteria) %>% 
    mutate(
        q1    = quantile(non_ref_identical_rate, 0.25, na.rm = TRUE),
        q3    = quantile(non_ref_identical_rate, 0.75, na.rm = TRUE),
        iqr   = q3 - q1,
        lower = q1 - 1.5 * iqr,
        upper = q3 + 1.5 * iqr,
        is_outlier = non_ref_identical_rate < lower | non_ref_identical_rate > upper,
        pair_label = paste(sample1, sample2, sep = " vs ")
    ) %>% 
    ungroup()

pdf_out <- "replicates_randoms/All_sp.pdf"
pdf(pdf_out)

species <- unique(df$sp_ploidy)
for (sp in species){
    
    sp_df <- df2 %>% filter(sp_ploidy == sp)
    p <- ggplot(sp_df, aes(x = type, y = non_ref_identical_rate, color = criteria)) +
        geom_boxplot(size = 1, outlier.shape = NA) +   # let *us* control outliers
        geom_point(position = position_jitter(width = 0.1), alpha = 0.4) +
        geom_text_repel(
            data = subset(sp_df, is_outlier),
            aes(label = pair_label),
            position = position_jitter(width = 0.1),
            show.legend = FALSE,
            size = 3
        ) +
        labs(
            x = "Pair Type",
            y = "Non-reference Identical Rate",
            title = sp
        ) +
        theme(
            axis.text  = element_text(size = 14),
            axis.title = element_text(size = 15),
            title      = element_text(size = 16)
        )
    
    print(p)
}
dev.off()
