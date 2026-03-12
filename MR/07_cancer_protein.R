#### Cancer and protein MR ####
### Run on BC4 with script 

library(TwoSampleMR)
library(ieugwasr)
library(data.table)
library(dplyr)
library(readxl)
library(R.utils)
library(ggplot2)

source("~/parameters/parameters_for_MR.R")

#ieugwasr::get_opengwas_jwt()
#user()
#Sys.setenv(OPENGWAS_JWT = "") ## add in my token

## Exposure SNPs for cancer
cancers <- extract_instruments(outcomes = c('ieu-b-85', 'ieu-a-1126', 'ebi-a-GCST006464'))
crc <- fread(file.path(working_dir, "cancer_sumstats/overall_CRC_GWAS_noUKBio_summary_stats_annotated.txt"))
crc <- as.data.frame(crc)
crc_exposure <- format_data(
  crc,
  type = "exposure",
  snp_col = "SNP",
  beta_col = "Effect",
  se_col = "StdErr",
  effect_allele_col = "Allele1",
  other_allele_col = "Allele2",
  pval_col = "P.value",
  chr_col = "chr_name",
  pos_col = "chrom_start",
  eaf_col = "Freq1"
)
crc_exposure$exposure <- "colorectal cancer (GECCO)"
crc_exposure <- crc_exposure[crc_exposure$pval.exposure < 5e-8, ]

all_cols <- union(names(crc_exposure), names(cancers))

crc_exposure_full <- crc_exposure
for (col in setdiff(all_cols, names(crc_exposure_full))) {
  crc_exposure_full[[col]] <- NA
}

cancers_full <- cancers
for (col in setdiff(all_cols, names(cancers_full))) {
  cancers_full[[col]] <- NA
}

crc_exposure_full <- crc_exposure_full[, all_cols]
cancers_full <- cancers_full[, all_cols]

cancer_exposure <- rbind(crc_exposure_full, cancers_full)

## write the full list of exposure SNPs to look up in protein outcome GWAS summary stats
#write.table(cancer_exposure, file = paste0(working_dir, "cancer_sumstats/cancer_exposure.txt"), sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

## Clump exposure SNPs - default is r2 0.001
cancer_clump <- clump_data(cancer_exposure)

## write.table 
#write.table(cancer_clump, file = paste0(working_dir, "cancer_sumstats/cancer_clump.txt"), sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

## Use script 06_extract_cancer_SNPs_from_protein.sh to extract these clumped SNPs from the proteins
## File to map protein rsids - from Lisa Hobson
protein_map <- fread(paste0(working_dir, "olink_rsid_map.txt"))
protein_map <- as.data.frame(protein_map)
protein_map_subset <- protein_map[,c("ID", "rsid", "POS19", "POS38")] 

## from the protein summary stats
protein_outcome <- fread(paste0(working_dir, "UKB_PPP_sumstats/concatenated_output_cancer.txt"), sep = " ")
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
  exposure_dat = cancer_clump,
  outcome_dat = outcome_protein_dat,
  action = 2 
)
write.table(dat, file = paste0(working_dir, "results/tables/07_cancer_protein_harmonised_data.txt"), sep = "\t")

## MR of cancer on protein
mr_results <- mr(dat, method = c("mr_ivw_mre", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode"))
mr_results$b <- as.numeric(mr_results$b)
mr_results$se <- as.numeric(mr_results$se)
mr_results$lower_CI <- mr_results$b - 1.96 * mr_results$se
mr_results$upper_CI <- mr_results$b + 1.96 * mr_results$se

## Derive odds ratios and 95% CIs
mr_results$OR <- exp(mr_results$b)  # Calculate OR from beta
mr_results$OR_lower_CI <- exp(mr_results$b - 1.96 * mr_results$se)  # Lower CI for OR
mr_results$OR_upper_CI <- exp(mr_results$b + 1.96 * mr_results$se)  # Upper CI for OR
write.table(mr_results, file = paste0(working_dir, "results/tables/07_cancer_protein_mr_results.txt"), sep = "\t", col.names = T, row.names = F)

## single SNP MR
res_single <- mr_singlesnp(dat, all_method = "mr_two_sample_ml")
write.table(res_single, file = paste0(working_dir, "results/tables/07_cancer_protein_single_SNP_mr.txt"), sep = "\t", col.names = T, row.names = F)

p2 <- mr_forest_plot(res_single)
pdf(paste0(working_dir, "results/figures/07_cancer_protein_single_SNP_mr.pdf"), width = 7, height = 7)  
for (i in seq_along(p2)) {
  print(p2[[i]])  # Print each plot into the PDF
}
dev.off()  

## Cochran's Q statistic
het <- mr_heterogeneity(dat)
write.table(het, file = paste0(working_dir, "results/tables/07_cancer_protein_heterogeneity.txt"), sep = "\t", col.names = T, row.names = F)

