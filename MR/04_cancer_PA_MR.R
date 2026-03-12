#### cancer-PA script #####
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

## Clump exposure SNPs - default is r2 0.001
cancer_clump <- clump_data(cancer_exposure)
duplicated_snps_crc <- crc$SNP[duplicated(crc$SNP)]
length(which(cancer_clump$SNP %in% duplicated_snps_crc)) ## 0

## Outcome SNPs for PA
## Doherty overall summary statistics 
gwas_overall_activity_doherty <- fread(
  cmd = paste("zcat", paste0(working_dir, "PA_sumstats/Doherty-2018-NatureComms-overall-activity.csv.gz")), 
  sep = ",", 
  header = TRUE
)
gwas_overall_activity_doherty <- gwas_overall_activity_doherty[, .(SNP, CHR, BP, BETA, SE, ALLELE1, ALLELE0, P_BOLT_LMM, A1FREQ)]

outcome_data_overall_activity_doherty <- gwas_overall_activity_doherty %>%
  rename(
    SNP = SNP,                  # SNP ID
    beta.outcome = BETA,        # Effect size (beta) for the outcome
    se.outcome = SE,            # Standard error for the outcome
    pval.outcome = P_BOLT_LMM,  # P-value for the outcome
    eaf.outcome = A1FREQ,       # Effect allele frequency for the outcome
    chr = CHR,          # Chromosome
    pos = BP,            # Base pair position
    effect_allele.outcome = ALLELE1,  # Effect allele
    other_allele.outcome = ALLELE0   # Other allele (reference allele)
  ) %>%
  select(SNP, beta.outcome, se.outcome, pval.outcome, eaf.outcome, chr, pos, effect_allele.outcome, other_allele.outcome)
outcome_data_overall_activity_doherty$id.outcome <- "overall_activity_doherty_SD"

## Doherty sedentary summary statistics
gwas_sedentary_doherty <- fread(cmd = paste("zcat", paste0(working_dir, "PA_sumstats/Doherty-2018-NatureComms-sedentary.csv.gz")), 
                                sep = ",", 
                                header = TRUE)

gwas_sedentary_doherty <- gwas_sedentary_doherty[, .(SNP, CHR, BP, BETA, SE, ALLELE1, ALLELE0, P_BOLT_LMM, A1FREQ)]

outcome_data_sedentary_doherty <- gwas_sedentary_doherty %>%
  select(
    SNP,                  # SNP ID
    beta.outcome = BETA,  # Effect size
    se.outcome = SE,      # Standard error
    pval.outcome = P_BOLT_LMM,  # P-value
    eaf.outcome = A1FREQ, # Effect allele frequency
    chr = CHR,    # Chromosome
    pos = BP,      # Base pair position
    effect_allele.outcome = ALLELE1,  # Effect allele
    other_allele.outcome = ALLELE0    # Reference allele
  )
outcome_data_sedentary_doherty$id.outcome <- "sedentary_doherty_SD"

## Read the outcome PA summary statistics - klimentidis overall acceleration
snps <- cancer_clump$SNP
chunk_size <- 50
snps_chunks <- split(snps, ceiling(seq_along(snps) / chunk_size))

outcome_data_list <- lapply(snps_chunks, function(snp_chunk) {
  extract_outcome_data(snps = snp_chunk, outcomes = "ebi-a-GCST006099")
})

# Combine all chunks into one dataframe
outcome_data_klimentidis <- bind_rows(outcome_data_list)
outcome_data_klimentidis$id.outcome <- "overall_activity_klimentidis_mg"

## Combine all the PA outcome GWAS
outcome_data_sedentary_doherty <- outcome_data_sedentary_doherty %>%
  mutate(outcome = "sedentary_doherty_SD")

outcome_data_overall_activity_doherty <- outcome_data_overall_activity_doherty %>%
  mutate(outcome = "overall_activity_doherty_SD")

outcome_data_klimentidis <- outcome_data_klimentidis %>%
  mutate(outcome = "overall_activity_klimentidis_mg")

## Harmonise
harmonised_data_1 <- harmonise_data(cancer_clump, outcome_data_sedentary_doherty)
harmonised_data_2 <- harmonise_data(cancer_clump, outcome_data_overall_activity_doherty)
harmonised_data_3 <- harmonise_data(cancer_clump, outcome_data_klimentidis)

common_columns <- intersect(intersect(colnames(harmonised_data_1), colnames(harmonised_data_2)), colnames(harmonised_data_3))

harmonised_data_1 <- harmonised_data_1 %>%
  select(all_of(common_columns)) %>%
  mutate(chr = as.character(chr))

harmonised_data_2 <- harmonised_data_2 %>%
  select(all_of(common_columns)) %>%
  mutate(chr = as.character(chr))

harmonised_data_3 <- harmonised_data_3 %>%
  select(all_of(common_columns)) %>%
  mutate(chr = as.character(chr))

harmonised_data_combined <- rbind(harmonised_data_1, harmonised_data_2, harmonised_data_3)

## Write harmonised data to file
write.table(harmonised_data_combined, file = paste0(working_dir, "results/tables/04_cancer_PA_harmonised_data.txt"), sep = "\t")

## Run the MR
mr_results_cancer_PA <- mr(harmonised_data_combined, method = c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode"))
write.table(mr_results_cancer_PA, file = paste0(working_dir, "results/tables/04_cancer_PA_mr_results.txt"), sep = "\t")

## heterogeneity
het <- mr_heterogeneity(harmonised_data_combined)
write.table(het, file = paste0(working_dir, "results/tables/04_cancer_PA_heterogeneity.txt"), sep = "\t", col.names = T, row.names = F)

# Perform Steiger filtering across all harmonised exposure-outcome pairs
harmonised_data_combined$samplesize.outcome <- NA  # start with NA
harmonised_data_combined$samplesize.outcome[grepl("klimentidis", harmonised_data_combined$outcome)] <- 91084
harmonised_data_combined$samplesize.outcome[grepl("doherty", harmonised_data_combined$outcome)] <- 91105
prevalence_map <- c(
  "Breast cancer (Combined Oncoarray; iCOGS; GWAS meta analysis) || id:ieu-a-1126" = 0.12,
  "colectal cancer (GECCO)" = 0.045,
  "Prostate cancer || id:ieu-b-85" = 0.125,
  "Endometrial cancer || id:ebi-a-GCST006464" = 0.03
)

# Add prevalence column based on exposure
harmonised_data_combined$prevalence.exposure <- prevalence_map[harmonised_data_combined$exposure]
steiger_ready <- harmonised_data_combined[
  complete.cases(harmonised_data_combined[, c(
    "beta.exposure", "se.exposure", "samplesize.exposure",
    "beta.outcome", "se.outcome", "samplesize.outcome",
    "eaf.exposure", "eaf.outcome", "prevalence.exposure"
  )]),
]

# Now run
steiger_res <- steiger_filtering(steiger_ready)

# Save full results (including steiger_dir, variance explained, etc.)
write.table(steiger_res, 
            file = paste0(working_dir, "results/tables/04_cancer_PA_steiger_filtering_all.txt"), 
            sep = "\t", 
            quote = FALSE, 
            row.names = FALSE)


