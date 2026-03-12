## MR script - proteins (exposure) and circulating proteins (outcome)
## 5/3/25
## Lucy Goudswaard
## Using R version 4.4.1 on BC4 - 04b_protein_PA_MR.sh

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
library(vcfR)
library(tidyr)
library(stringr)

## Check token and user
ieugwasr::get_opengwas_jwt()
user()
#Sys.setenv(OPENGWAS_JWT = "") ## add in my token

## Parameter file with paths - need to sort setwd to scripts
source("~/parameters/parameters_for_MR.R")

## Read exposure data - these SNPs are from Sun 
mr_proteins <- read.table(paste0(working_dir, "protein_names_MR.txt"), header = F)

## Protein map
protein_map <- fread(paste0(working_dir, "olink_protein_map_3k_v1.tsv"))
protein_map <- as.data.frame(protein_map)

## Replace with finemapped genetic variants
finemap_instruments <- read.table(paste0(working_dir, "UKB_PPP_sumstats/instruments.txt"), header = T)

## Keep relevant proteins
finemap_instruments_edited <- finemap_instruments %>%
  separate(exposure, into = c("exposure_study", "exposure_population", "exposure_sex", "exposure", "exposure_cohort", "exposure_finemap"), sep = ";", remove = FALSE) %>%
  mutate(
    SeqId = case_when(
      exposure_study == "ferkingstad_2021_PMID34857953" ~ sub("^([A-Za-z0-9]+_[A-Za-z0-9]+).*", "\\1", exposure),
      exposure_study == "zhang_2022_PMID35501419" ~ sub("^SeqId_(.*)", "\\1", exposure),
      exposure_study == "pietzner_2021_PMID34648354" ~ sub(".*_([^_]+_[^_]+)$", "\\1", exposure),
      TRUE ~ NA_character_
    ),
    UniProt = case_when(
      exposure_study == "sun_2023_PMID37794186" ~ sub("^[^_]+_([^_]+)_.*", "\\1", exposure),
      TRUE ~ NA_character_
    )
  ) %>%
  mutate(SeqId = str_replace_all(SeqId, "_", "-")) %>%
  filter(
    exposure_finemap != "Finimom"
  )

## Only keep the proteins interested in the current analysis 
finemap_instruments_ukb <- finemap_instruments_edited %>%
  filter(exposure_study == "sun_2023_PMID37794186") %>%
  filter(exposure_population == "EUR")  %>%
  mutate(gene = sub("_.*", "", exposure)) %>%
  filter(gene %in% mr_proteins$V1)

## Merge finemapped instruments with map
finemap_instruments_ukb <- finemap_instruments_ukb %>%
  mutate(olink_id = str_extract(exposure, "OID\\d+"))
finemap_instruments_ukb_prelim <- finemap_instruments_ukb %>%
  left_join(protein_map, by = c("olink_id" = "OlinkID"))

## Prioritise by ld 0.001 and p value priority, then by coloc
finemap_instruments_ukb_filtered <- finemap_instruments_ukb_prelim %>%
  # Step 1: Filter out palindromic SNPs with intermediate allele frequency
  filter(
    !((
      (effect_allele.exposure %in% c("C", "G") & other_allele.exposure %in% c("C", "G")) |
        (effect_allele.exposure %in% c("A", "T") & other_allele.exposure %in% c("A", "T"))
    ) & eaf.exposure > 0.4 & eaf.exposure < 0.6)
  ) %>%
  # Step 2: Filter based on finemapping method & significance threshold
  group_by(exposure) %>%
  filter(
    (any(exposure_finemap == "p_ld-0.001" & pval.exposure < 5e-8) & 
       exposure_finemap == "p_ld-0.001" & pval.exposure < 5e-8) |
      (!any(exposure_finemap == "p_ld-0.001" & pval.exposure < 5e-8) & 
         exposure_finemap == "coloc.abf" & pval.exposure < 5e-8)
  ) %>%
  ungroup() %>%
  # Step 3: Filter SNPs within ±500kb window of gene boundaries, same chromosome
  filter(
    !is.na(gene_start) & !is.na(gene_end) & !is.na(chr) & !is.na(pos.exposure),
    chr == chr.exposure,
    pos.exposure >= (gene_start - 500000),
    pos.exposure <= (gene_end + 500000)
  )

## keep the SNP for every protein with the smallest p value
finemap_instruments_ukb_filtered <- finemap_instruments_ukb_filtered %>%
  # Step 2: Keep only the SNP with the smallest p-value for each exposure
  group_by(exposure) %>%
  slice_min(pval.exposure, n = 1, with_ties = FALSE) %>%
  ungroup()
finemap_instruments_ukb_filtered <- as.data.frame(finemap_instruments_ukb_filtered)

## Outcome data
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
    chr.outcome = CHR,          # Chromosome
    bp.outcome = BP,            # Base pair position
    effect_allele.outcome = ALLELE1,  # Effect allele
    other_allele.outcome = ALLELE0   # Other allele (reference allele)
    ) %>%
  select(SNP, beta.outcome, se.outcome, pval.outcome, eaf.outcome, chr.outcome, bp.outcome, effect_allele.outcome, other_allele.outcome)
outcome_data_overall_activity_doherty$id.outcome <- "overall_activity_doherty_SD"

## Doherty sedentary summary statistics
gwas_sedentary_doherty <- fread(cmd = paste("zcat", paste0(working_dir, "PA_sumstats/Doherty-2018-NatureComms-sedentary.csv.gz")), 
                                sep = ",", 
                                header = TRUE)

gwas_sedentary_doherty <- gwas_sedentary_doherty[, .(SNP, CHR, BP, BETA, SE, ALLELE1, ALLELE0, P_BOLT_LMM, A1FREQ)]

outcome_data_sedentary_doherty <- gwas_sedentary_doherty %>%
  rename(
    SNP = SNP,                  # SNP ID
    beta.outcome = BETA,        # Effect size (beta) for the outcome
    se.outcome = SE,            # Standard error for the outcome
    pval.outcome = P_BOLT_LMM,  # P-value for the outcome
    eaf.outcome = A1FREQ,       # Effect allele frequency for the outcome
    chr.outcome = CHR,          # Chromosome
    bp.outcome = BP,            # Base pair position
    effect_allele.outcome = ALLELE1,  # Effect allele
    other_allele.outcome = ALLELE0   # Other allele (reference allele) 
    ) %>%
  select(SNP, beta.outcome, se.outcome, pval.outcome, eaf.outcome, chr.outcome, bp.outcome, effect_allele.outcome, other_allele.outcome)
outcome_data_sedentary_doherty$id.outcome <- "sedentary_doherty_SD"

## Read the outcome PA summary statistics - klimentidis overall acceleration
# Your list of SNPs
snp_list <- finemap_instruments_ukb_filtered$SNP

# Batch size
batch_size <- 30

# Split into chunks
batches <- split(snp_list, ceiling(seq_along(snp_list) / batch_size))

# Loop through batches
all_results <- list()

for (i in seq_along(batches)) {
  cat("Processing batch", i, "of", length(batches), "\n")
  try({
    res <- extract_outcome_data(
      snps = batches[[i]],
      outcomes = "ebi-a-GCST006099",
      proxies = TRUE
    )
    all_results[[i]] <- res
    Sys.sleep(2)  # avoid hammering the API
  })
}

# Combine all results
# Find common columns across all batches
common_cols <- Reduce(intersect, lapply(all_results, colnames))

# Subset each batch to common columns
all_results_aligned <- lapply(all_results, function(x) x[, common_cols, drop = FALSE])

# Bind safely
outcome_data_klimentidis <- do.call(rbind, all_results_aligned)
outcome_data_klimentidis$id.outcome <- "overall_activity_klimentidis_mg"

## Add outcome column
outcome_data_sedentary_doherty <- outcome_data_sedentary_doherty %>%
  mutate(outcome = "sedentary_doherty_SD")

outcome_data_overall_activity_doherty <- outcome_data_overall_activity_doherty %>%
  mutate(outcome = "overall_activity_doherty_SD")

outcome_data_klimentidis <- outcome_data_klimentidis %>%
  mutate(outcome = "overall_activity_klimentidis_mg")

## Harmonise
harmonised_data_1 <- harmonise_data(finemap_instruments_ukb_filtered, outcome_data_sedentary_doherty)
harmonised_data_2 <- harmonise_data(finemap_instruments_ukb_filtered, outcome_data_overall_activity_doherty)
harmonised_data_3 <- harmonise_data(finemap_instruments_ukb_filtered, outcome_data_klimentidis)

common_columns <- intersect(intersect(colnames(harmonised_data_1), colnames(harmonised_data_2)), colnames(harmonised_data_3))

harmonised_data_1 <- harmonised_data_1 %>%
  select(all_of(common_columns))

harmonised_data_2 <- harmonised_data_2 %>%
  select(all_of(common_columns)) 

harmonised_data_3 <- harmonised_data_3 %>%
  select(all_of(common_columns)) 

harmonised_data_combined <- rbind(harmonised_data_1, harmonised_data_2, harmonised_data_3)
write.table(harmonised_data_combined, file = paste0(working_dir, "results/tables/02_protein_PA_harmonised_data.txt"), sep = "\t")

## Run the MR
mr_results <- mr(harmonised_data_combined, method_list = "mr_wald_ratio")
write.table(mr_results, file = paste0(working_dir, "results/tables/02_protein_PA_mr_results.txt"), sep = "\t")

# Perform Steiger filtering across all harmonised exposure-outcome pairs
harmonised_data_combined$samplesize.outcome <- NA  # start with NA
harmonised_data_combined$samplesize.outcome[grepl("klimentidis", harmonised_data_combined$outcome)] <- 91084
harmonised_data_combined$samplesize.outcome[grepl("doherty", harmonised_data_combined$outcome)] <- 91105
steiger_res <- steiger_filtering(harmonised_data_combined)

# Save full results (including steiger_dir, variance explained, etc.)
write.table(steiger_res, 
            file = paste0(working_dir, "results/tables/02_protein_PA_steiger_filtering_all.txt"), 
            sep = "\t", 
            quote = FALSE, 
            row.names = FALSE)