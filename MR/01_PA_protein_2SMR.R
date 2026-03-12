## MR script - physical activity (exposure) and circulating proteins (outcome)
## 5/2/25
## Lucy Goudswaard
## Using R version 4.4.1 on BC4

## Install and load libraries
#install.packages("remotes")
#remotes::install_github("MRCIEU/TwoSampleMR")
library(TwoSampleMR)
library(ieugwasr)
library(data.table)
library(dplyr)
library(readxl)
library(R.utils)
library(ggplot2)

## Check token and user
ieugwasr::get_opengwas_jwt()
user()

## Parameter file with paths - need to sort setwd to scripts
source("~/parameters/parameters_for_MR.R")

## File to map protein rsids - from Lisa Hobson
protein_map <- fread(paste0(working_dir, "olink_rsid_map.txt"))
protein_map <- as.data.frame(protein_map)
protein_map_subset <- protein_map[,c("ID", "rsid", "POS19", "POS38")] 

## Read exposure data - Klimentidis overall activity SNPs (HG19/GRCh37)
PA_exposure <- read_xlsx(paste0(working_dir, "PA_sumstats/Klimentidis_PA_SNPs.xlsx"))
PA_exposure <- as.data.frame(PA_exposure)
colnames(PA_exposure) <- c("SNP", "effect_allele", "other_allele", "chr", "pos",
                           "gene", "eaf", "beta", "se" ,"samplesize", "R2", "F", "units")
PA_exp_dat <- format_data(PA_exposure, type = "exposure")
PA_exp_dat$exposure <- "overall_physical_activity_klimentidis_2_SNPs"

## overall PA SNPs
PA_exposure_10 <- read_xlsx(paste0(working_dir, "PA_sumstats/Klimentidis_PA_SNPs_10.xlsx"))
PA_exposure_10 <- as.data.frame(PA_exposure_10)
colnames(PA_exposure_10) <- c("SNP", "effect_allele", "other_allele", "chr", "pos",
                           "gene", "eaf", "beta", "se" ,"samplesize", "R2", "F", "units")
PA_exp_dat_10 <- format_data(PA_exposure_10, type = "exposure")
PA_exp_dat_10$exposure <- "overall_physical_activity_klimentidis_10_SNPs"

## Also read in sedentary PA SNPs from Doherty et al - pos38 (Grch38/hg38) derived and needs to be kept
PA_sedentary_exposure <- read_xlsx(paste0(working_dir, "PA_sumstats/Doherty_sedentary_PA_SNPs.xlsx"), col_names = TRUE)
PA_sedentary_exposure <- as.data.frame(PA_sedentary_exposure)
PA_sedentary_exposure <- PA_sedentary_exposure %>%
  select(-pos19)
colnames(PA_sedentary_exposure) <- c('SNP', "chr", 'effect_allele', "other_allele", "beta", "se", "eaf", "eaf2", "pos")
PA_sedentary_exp_dat <- format_data(PA_sedentary_exposure, type = "exposure")
PA_sedentary_exp_dat$exposure <- "sedentary_doherty"

## overall PA SNPs from Doherty et al
PA_doherty_exposure <- read_xlsx(paste0(working_dir, "PA_sumstats/Doherty_overall_activity_SNPs.xlsx"))
PA_doherty_exposure <- as.data.frame(PA_doherty_exposure)
colnames(PA_doherty_exposure) <- c("SNP", "effect_allele", "other_allele", "chr", "pos",
                           "gene", "eaf", "beta", "se" ,"samplesize", "R2", "F")
PA_exp_dat_doherty <- format_data(PA_doherty_exposure, type = "exposure")
PA_exp_dat_doherty$exposure <- "overall_physical_activity_doherty"

## Combine the physical activity exposure data
# Find common columns across all 4 exposure datasets
common_cols <- Reduce(intersect, list(
  names(PA_exp_dat),
  names(PA_exp_dat_10),
  names(PA_sedentary_exp_dat),
  names(PA_exp_dat_doherty)
))

# Select only the common columns from each dataset
PA_exp_dat_filtered <- PA_exp_dat %>% select(all_of(common_cols))
PA_exp_dat_10_filtered <- PA_exp_dat_10 %>% select(all_of(common_cols))
PA_sedentary_exp_dat_filtered <- PA_sedentary_exp_dat %>% select(all_of(common_cols))
PA_exp_dat_doherty_filtered <- PA_exp_dat_doherty %>% select(all_of(common_cols))

# Combine all filtered datasets
PA_combined_exp_dat <- bind_rows(
  PA_exp_dat_filtered,
  PA_exp_dat_10_filtered,
  PA_sedentary_exp_dat_filtered,
  PA_exp_dat_doherty_filtered
)

## Clump exposure data with default settings (r2 = 0.001)
PA_clump <- data.frame(ld_clump(
  dplyr::tibble(rsid=PA_combined_exp_dat$SNP, pval=PA_combined_exp_dat$pval.exposure, id=PA_combined_exp_dat$exposure),
  clump_kb = 10000,
  clump_r2 = 0.001,
  clump_p= 5e-8,
  bfile = paste0(working_dir, "GRCh37/EUR"),
  plink_bin = "/plink"
))

## Keep the SNPs in PA_clump
PA_combined_exp_dat <- PA_combined_exp_dat %>%
  filter(SNP %in% PA_clump$rsid)

## Ran script 03_extract_outcome_SNPs_v2.sh on BC4 to extract the PA SNPs
## from the protein summary stats
protein_outcome <- fread(paste0(working_dir, "UKB_PPP_sumstats/concatenated_output_overall_sedentary.txt"), sep = " ")
protein_outcome <- as.data.table(protein_outcome)
colnames(protein_outcome)[colnames(protein_outcome) == "EXTRA\tgene"] <- "gene"
protein_outcome <- protein_outcome %>%
  mutate(gene = sub(".*\t", "", gene)) 

protein_outcome <- protein_outcome %>%
  group_by(gene) %>%
  filter(!(gene == "TNF" & row_number() > 2)) %>%
  ungroup()
protein_outcome <- as.data.frame(protein_outcome)

## Merge in the rsID (the build is 38 looking at position, and so are the PA SNPs)
protein_SNPs_with_rsid <- merge(protein_outcome, protein_map_subset, all.x = T, all.y = F, by = "ID")
protein_SNPs_with_rsid$p <- 10^-(protein_SNPs_with_rsid$LOG10P)
outcome_protein_dat <- protein_SNPs_with_rsid[,c("CHROM", "GENPOS", "ALLELE0", "ALLELE1", "A1FREQ",
                                                 'N', "BETA", "SE", "gene", "gene", "rsid", "p")]
colnames(outcome_protein_dat) <- c("chr.outcome", "pos.outcome", "other_allele.outcome", "effect_allele.outcome",
                                      "eaf.outcome", "samplesize", "beta.outcome", "se.outcome", "outcome", "id.outcome", 
                                      "SNP", "pval.outcome")

## Harmonise data 
dat <- harmonise_data(
  exposure_dat = PA_combined_exp_dat,
  outcome_dat = outcome_protein_dat
)
write.table(dat, file = paste0(working_dir, "results/tables/01_PA_protein_harmonised_data.txt"), 
            sep = "\t", row.names = FALSE, col.names = TRUE)

## Function to choose MR method based on number of SNPs
run_mr <- function(dat_subset) {
  nsnps <- length(unique(dat_subset$SNP))
  
  if (nsnps == 2) {
    method_list <- c("mr_ivw_fe", "mr_weighted_median", "mr_weighted_mode")
  } else {
    method_list <- c("mr_ivw_mre", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode")
  }
  
  mr(dat_subset, method_list = method_list)
}

## Run MR for each exposure
exposures <- unique(dat$exposure)
mr_results_all <- lapply(exposures, function(e) {
  dat_sub <- subset(dat, exposure == e)
  res <- run_mr(dat_sub)
  res$exposure <- e
  return(res)
})

mr_results <- do.call(rbind, mr_results_all)
write.table(mr_results, file = paste0(working_dir, "results/tables/01_PA_protein_mr_results.txt"), 
            sep = "\t", row.names = FALSE, col.names = TRUE)

## Sensitivity analyses
## Single SNP MR
res_single <- mr_singlesnp(dat, all_method = "mr_two_sample_ml")
write.table(res_single, file = paste0(working_dir, "results/tables/01_PA_protein_single_SNP_mr_results.txt"), 
            sep = "\t", col.names = TRUE, row.names = FALSE)

p2 <- mr_forest_plot(res_single)
pdf(paste0(working_dir, "results/figures/01_PA_protein_single_SNP_mr_results.pdf"), width = 7, height = 7)  
for (i in seq_along(p2)) {
  print(p2[[i]])
}
dev.off()  

## Copy files across from BC4 to RDSF
