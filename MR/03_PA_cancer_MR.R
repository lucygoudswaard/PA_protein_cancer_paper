#### PA - cancer script #####
## submitted with 03c_submit_PA_cancer_MR_script.sh

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

## Read exposure data - Klimentidis overall activity SNPs (from ebi-a-GCST006099), build is GRCh38
PA_exposure <- read_xlsx(paste0(working_dir, "PA_sumstats/Klimentidis_PA_SNPs.xlsx"))
PA_exposure <- as.data.frame(PA_exposure)
colnames(PA_exposure) <- c("SNP", "effect_allele", "other_allele", "chr", "pos",
                           "gene", "eaf", "beta", "se" ,"samplesize", "R2", "F", "units")
PA_exp_dat <- format_data(PA_exposure, type = "exposure")
PA_exp_dat$exposure <- "overall_physical_activity_klimentidis"

## overall PA SNPs
PA_exposure_10 <- read_xlsx(paste0(working_dir, "PA_sumstats/Klimentidis_PA_SNPs_10.xlsx"))
PA_exposure_10 <- as.data.frame(PA_exposure_10)
colnames(PA_exposure_10) <- c("SNP", "effect_allele", "other_allele", "chr", "pos",
                              "gene", "eaf", "beta", "se" ,"samplesize", "R2", "F", "units")
PA_exp_dat_10 <- format_data(PA_exposure_10, type = "exposure")
PA_exp_dat_10$exposure <- "overall_physical_activity_klimentidis_10_SNPs"

## Also read in sedentary PA SNPs from Doherty et al, build is GRCh38
PA_sedentary_exposure <- read_xlsx(paste0(working_dir, "PA_sumstats/Doherty_sedentary_PA_SNPs.xlsx"), col_names = TRUE)
PA_sedentary_exposure <- as.data.frame(PA_sedentary_exposure)
colnames(PA_sedentary_exposure) <- c('SNP', "chr", "pos", 'effect_allele', "other_allele", "beta", "se", "eaf", "eaf2", "pos19")
PA_sedentary_exp_dat <- format_data(PA_sedentary_exposure, type = "exposure")
PA_sedentary_exp_dat$exposure <- "sedentary_doherty"

## Doherty overall PA SNPs, build is GRCh38/hg38
PA_doherty_exposure <- read_xlsx(paste0(working_dir, "PA_sumstats/Doherty_overall_activity_SNPs.xlsx"))
PA_doherty_exposure <- as.data.frame(PA_doherty_exposure)
colnames(PA_doherty_exposure) <- c("SNP", "effect_allele", "other_allele", "chr", "pos",
                                   "gene", "eaf", "beta", "se" ,"samplesize", "R2", "F")
PA_exp_dat_doherty <- format_data(PA_doherty_exposure, type = "exposure")
PA_exp_dat_doherty$exposure <- "overall_physical_activity_doherty"

## Combine the physical activity exposure data
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

## Clump exposure data with default settings (r2 = 0.001) - one SNP dropped for doherty overall physical activity (rs2696625)
PA_clump <- clump_data(PA_combined_exp_dat)

## Cancer outcomes from IEU open GWAS
ao <- available_outcomes()

## Extract cancer SNPs.
## Prostate cancer is ieu-b-85, build hg19/grch37
## Breast cancer ieu-a-1126, build hg19/grch37
## Endometrial cancer ebi-a-GCST006464, build hg19/grch37
## Colorectal cancer from GECCO (without UKB)
out_ieu_b_85 <- extract_outcome_data(snps = PA_clump$SNP, outcomes = "ieu-b-85", proxies = TRUE, rsq = 0.8)
out_ieu_a_1126 <- extract_outcome_data(snps = PA_clump$SNP, outcomes = "ieu-a-1126", proxies = TRUE, rsq = 0.8)
out_ebi <- extract_outcome_data(snps = PA_clump$SNP, outcomes = "ebi-a-GCST006464", proxies = TRUE, rsq = 0.8)
common_cols <- c(
  "SNP", "chr", "pos", "beta.outcome", "se.outcome",
  "samplesize.outcome", "pval.outcome", "eaf.outcome",
  "effect_allele.outcome", "other_allele.outcome",
  "outcome", "id.outcome", "originalname.outcome",
  "outcome.deprecated", "mr_keep.outcome", "data_source.outcome"
)
out_ieu_b_85_subset <- out_ieu_b_85[, common_cols]
out_ieu_a_1126_subset <- out_ieu_a_1126[, common_cols]
out_ebi_subset <- out_ebi[, common_cols]
cancer_out_data <- rbind(out_ieu_b_85_subset, out_ieu_a_1126_subset, out_ebi_subset)
crc <- fread(file.path(working_dir, "cancer_sumstats/overall_CRC_GWAS_noUKBio_summary_stats_annotated.txt"))
crc <- as.data.frame(crc)
crc_outcome <- format_data(
  crc,
  type = "outcome",
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
crc_outcome$outcome <- "colorectal cancer (GECCO)"

## Harmonise the data - rs25981 and rs26579 drop from sedentary SNPs because mr_keep = FALSE - potentially palindromic SNPs)
dat_cancer <- harmonise_data(exposure_dat = PA_clump, outcome_dat = cancer_out_data)
dat_crc <- harmonise_data(
  exposure_dat = PA_clump,
  outcome_dat = crc_outcome,
  action = 2
)
common_cols <- intersect(names(dat_cancer), names(dat_crc))

dat <- rbind(
  dat_cancer[, common_cols],
  dat_crc[, common_cols]
)
duplicated_snps_crc <- crc$SNP[duplicated(crc$SNP)]
length(which(PA_combined_exp_dat$SNP %in% duplicated_snps_crc)) ## none of the duplicated crc rsids are in the PA SNP lists
write.table(dat, file = paste0(working_dir, "results/tables/03_PA_cancer_harmonised_data.txt"), sep = "\t")

## MR of PA on cancer - 2 sedentary doherty SNPs not retaining
dat_two <- subset(dat, exposure == "overall_physical_activity_klimentidis")
dat_rest <- subset(dat, exposure != "overall_physical_activity_klimentidis")

# Run MR separately
mr_two <- mr(dat_two, method_list = c("mr_ivw_fe", "mr_weighted_median", "mr_weighted_mode"))
mr_rest <- mr(dat_rest, method_list = c("mr_ivw_mre", "mr_egger_regression",
                                        "mr_weighted_median", "mr_weighted_mode"))

# Combine results
mr_results <- rbind(mr_two, mr_rest)
mr_results$b <- as.numeric(mr_results$b)
mr_results$se <- as.numeric(mr_results$se)
mr_results$lower_CI <- mr_results$b - 1.96 * mr_results$se
mr_results$upper_CI <- mr_results$b + 1.96 * mr_results$se

## Derive odds ratios and 95% CIs
mr_results$OR <- exp(mr_results$b)  # Calculate OR from beta
mr_results$OR_lower_CI <- exp(mr_results$b - 1.96 * mr_results$se)  # Lower CI for OR
mr_results$OR_upper_CI <- exp(mr_results$b + 1.96 * mr_results$se)  # Upper CI for OR
write.table(mr_results, file = paste0(working_dir, "results/tables/03_PA_cancer_mr_results.txt"), sep = "\t", col.names = T, row.names = F)

## Plot the MR estimates
pdf(file = paste0(working_dir, "results/figures/03_PA_cancer_mr_results.pdf"), height = 6, width = 12)
ggplot(mr_results, aes(x = b, y = outcome, shape = exposure, color = method)) +
  geom_pointrange(aes(xmin = lower_CI, xmax = upper_CI), size = 1, position = position_dodge(width = 0.5)) +  # Apply dodge
  geom_vline(xintercept = 0, linetype = "dashed") +  # Add vertical line at 0 (null effect)
  labs(x = "Log odds cancer risk per unit higher PA (+/- 95% CI)", y = "Cancer Type", title = "MR estimate for effect of physical activity on cancer outcomes") +
  theme_minimal() +
  theme(legend.position = "top", legend.title = element_blank(), plot.title = element_text(hjust = 0.5)) +
  theme(axis.text.y = element_text(size = 10)) +
  theme(axis.title.y = element_text(size = 12)) +  # Adjust axis title size
  guides(shape = guide_legend(ncol = 2), color = guide_legend(ncol = 2))  # Split legend into two columns
dev.off()

## single SNP MR
res_single <- mr_singlesnp(dat, all_method = "mr_two_sample_ml")
write.table(res_single, file = paste0(working_dir, "results/tables/03_PA_cancer_single_SNP_mr_results.txt"), sep = "\t", col.names = T, row.names = F)

p2 <- mr_forest_plot(res_single)
pdf(paste0(working_dir, "results/figures/03_PA_cancer_single_SNP_mr_results.pdf"), width = 7, height = 7)  
for (i in seq_along(p2)) {
  print(p2[[i]])  # Print each plot into the PDF
}
dev.off()  

## Cochran's Q statistic
het <- mr_heterogeneity(dat)
write.table(het, file = paste0(working_dir, "results/tables/03_PA_cancer_heterogeneity.txt"), sep = "\t", col.names = T, row.names = F)

# Perform Steiger filtering across all harmonised exposure-outcome pairs
dat$samplesize.exposure <- NA  # start with NA
dat$samplesize.exposure[grepl("klimentidis", dat$exposure)] <- 91084
dat$samplesize.exposure[grepl("doherty", dat$exposure)] <- 91105
steiger_res <- steiger_filtering(dat)

# Save full results (including steiger_dir, variance explained, etc.)
write.table(steiger_res, 
            file = paste0(working_dir, "results/tables/03_PA_cancer_steiger_filtering_all.txt"), 
            sep = "\t", 
            quote = FALSE, 
            row.names = FALSE)


