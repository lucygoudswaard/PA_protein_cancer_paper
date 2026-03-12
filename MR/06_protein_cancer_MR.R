#### Protein and cancer MR ###
## Ran with 05_protein_cancer_MR.sh

## Load required packages
library(TwoSampleMR)
library(ieugwasr)
library(data.table)
library(dplyr)
library(readxl)
library(stringr)
library(tidyr)
library(ggplot2)

## Source parameters
source("~/parameters/parameters_for_MR.R")

## Read exposure data - these SNPs are from Sun 
## Just do it for the 40 proteins that have an association with PA and cancer observationally
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

## Cancer outcomes: Prostate (ieu-b-85), Breast (ieu-a-1126) and endometrial (ebi-a-GCST006464)
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
      outcomes = c("ieu-b-85", "ieu-a-1126", "ebi-a-GCST006464"),
      proxies = TRUE
    )
    all_results[[i]] <- res
    Sys.sleep(2)  # avoid hammering the API
  })
}


# Combine all successfully retrieved results
all_results_clean <- Filter(function(x) {
  !is.null(x) && is.data.frame(x) && nrow(x) > 0 && all(c("SNP", "beta.outcome") %in% names(x))
}, all_results)

cancer_outcome_dat <- do.call(rbind, all_results_clean)
print(paste0("dim cancer_outcome_dat = ", paste(dim(cancer_outcome_dat), collapse = " x ")))

## Also extract SNPs from this CRC GWAS
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

## Harmonise
cancer_dat <- harmonise_data(exposure_dat = finemap_instruments_ukb_filtered, outcome_dat = cancer_outcome_dat, action = 2)
dat_crc <- harmonise_data(
  exposure_dat = finemap_instruments_ukb_filtered,
  outcome_dat = crc_outcome,
  action = 2
)
common_cols <- intersect(names(cancer_dat), names(dat_crc))

harmonised_dat <- rbind(
  cancer_dat[, common_cols],
  dat_crc[, common_cols]
)
write.table(harmonised_dat, file = paste0(working_dir, "results/tables/06_protein_cancer_harmonised.txt"), sep = "\t")

## Run MR
mr_results <- mr(harmonised_dat, method_list = "mr_wald_ratio")
mr_results <- mr_results %>%
  mutate(
    OR = exp(b),
    OR_lower_CI = exp(b - 1.96 * se),
    OR_upper_CI = exp(b + 1.96 * se)
  )

write.table(mr_results, file = paste0(working_dir, "results/tables/06_protein_cancer_mr_results.txt"), sep = "\t", col.names = TRUE, row.names = FALSE)

## Forest plot
# Split MR results by protein (exposure) for manageable plots
mr_split <- split(mr_results, mr_results$exposure)

pdf(file = paste0(working_dir, "results/figures/06_protein_cancer_mr_forest_split.pdf"), width = 10, height = 6)

for (protein in names(mr_split)) {
  plot_dat <- mr_split[[protein]]
  
  p <- ggplot(plot_dat, aes(x = b, y = outcome, color = method, shape = method)) +
    geom_pointrange(aes(xmin = b - 1.96 * se, xmax = b + 1.96 * se),
                    position = position_dodge(width = 0.6)) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    labs(title = paste0("MR Results for Protein: ", protein),
         x = "Effect estimate (log OR per SD higher protein)", y = "Cancer Outcome") +
    theme_minimal(base_size = 12) +
    theme(legend.position = "top")
  
  print(p)
}

dev.off()

# Perform Steiger filtering across all harmonised exposure-outcome pairs
steiger_res <- steiger_filtering(harmonised_dat)

# Save full results (including steiger_dir, variance explained, etc.)
write.table(steiger_res, 
            file = paste0(working_dir, "results/tables/06_protein_cancer_steiger_filtering_all.txt"), 
            sep = "\t", 
            quote = FALSE, 
            row.names = FALSE)


