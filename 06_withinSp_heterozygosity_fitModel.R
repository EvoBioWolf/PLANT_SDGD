library(DHARMa)
library(glmmTMB)
library(tidyverse)

# Read command-line arguments
args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
output_summary_file <- args[2]
output_RData <- args[3]
# Load data
df <- read.table(input_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

sp <- str_remove(basename(input_file), "\\..*$")
sp <- str_remove(sp, "_\\d+x$")
biomass_file <- "06_withinSp_heterozygosity_BiomassResid_withinDiv.txt"
biomass <- read_delim(biomass_file, delim = "\t") %>% filter(species == sp) %>% select(plotID, resid_species_div)
# species	ID	Plot	biomass_across_years	div	sp_plot	species_div_mean	resid_species_div	plotID
# AjuRep	28-H1	B1A22	0	60	AjuRep_H1	0.140972222222222	-0.140972222222222	H1

criterias <- unique(df$criteria)
cat("", file = output_summary_file)
models <- list()

for (i in criterias) {
    df_fit <- df %>% filter(criteria == i) %>% left_join(biomass, by = "plotID")
    model_nbinom2 <- glmmTMB(
      polymorphic_site ~ diversity_log2 + Expr_ln + resid_species_div +
        (1 | geneID) + (1 | sampleID) + (1 | plotID) +
        offset(genotyped_site_ln),
      ziformula = ~1,   # keep if many extra zeros
      family = nbinom2,
      data = df_fit
    )
    a_nbinom2 <- AIC(model_nbinom2)


    set.seed(1)
    res <- simulateResiduals(model_nbinom2, n = 1000)
    disp <- tryCatch(testDispersion(res), error = function(e) NA)
    zi   <- tryCatch(testZeroInflation(res), error = function(e) NA)
    uni  <- tryCatch(testUniformity(res), error = function(e) NA)

    # Append summary
    summ <- capture.output(summary(model_nbinom2))
    cat(paste0("### Criterion: ", i, " — Model_nbinom2 Summary\n"), file = output_summary_file, append = TRUE)
    cat(paste(summ, collapse = "\n"), "\n\n", file = output_summary_file, append = TRUE)
    cat(paste("AIC:", a_nbinom2, collapse = "\n"), "\n\n", file = output_summary_file, append = TRUE)
    if (!is.na(disp)[1]) cat(sprintf("DHARMa testDispersion:  stat=%0.3f, p=%0.3g\n",
                                     disp$statistic, disp$p.value),
                             file = output_summary_file, append = TRUE)
    if (!is.na(zi)[1])   cat(sprintf("DHARMa testZeroInflation:  stat=%0.3f, p=%0.3g\n",
                                     zi$statistic, zi$p.value),
                             file = output_summary_file, append = TRUE)
    if (!is.na(uni)[1])  cat(sprintf("DHARMa testUniformity:  KS=%0.3f, p=%0.3g\n",
                                     uni$statistic, uni$p.value),
                             file = output_summary_file, append = TRUE)

  
    # Store
    tag = paste0(as.character(i), '_nbinom2') 
    models[[tag]] <- model_nbinom2
}

saveRDS(models, file = output_RData)
