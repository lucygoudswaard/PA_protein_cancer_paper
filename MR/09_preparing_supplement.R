### Prepare supplementary tables ###

## libraries
library(dplyr)
library(tidyr)
library(devtools)
library(stringr)
library(ggplot2)
library(patchwork)
library(cowplot)

## clear environment
rm(list=ls())

## connect to project and setwd to script folder

## read in parameter file (specified on command line)
source("parameter_file/parameters_for_r.R")

## Read in all the MR results files
PA_protein_MR <- read.table(paste0(results_dir, "tables/MR/01_PA_protein_mr_results.txt"), header = T, sep = "\t")
protein_PA_MR <- read.table(paste0(results_dir, "tables/MR/02_protein_PA_mr_results.txt"), header = T, sep = "\t")
PA_cancer_MR <- read.table(paste0(results_dir, "tables/MR/03_PA_cancer_mr_results.txt"), header = T, sep = "\t")
cancer_PA_MR <- read.table(paste0(results_dir, "tables/MR/04_cancer_PA_mr_results.txt"), header = T, sep = "\t")
protein_cancer_MR <- read.table(paste0(results_dir, "tables/MR/06_protein_cancer_mr_results.txt"), header = T, sep = "\t")
cancer_protein_MR <- read.table(paste0(results_dir, "tables/MR/07_cancer_protein_mr_results.txt"), header = T, sep = "\t")

## Restrict to proteins associated with both overall physical activity and cancer
## Observational overall physical activity and protein results
overall_PA_protein_regression_bmi <- read.table(paste0(results_dir, "tables/overall_PA_protein_regression_fully_adjusted_bmi.txt"), header = T, sep = "\t")
overall_PA_protein_regression_bmi$activity_type <- "overall_acceleration_average"
overall_PA_protein_regression_bmi$category <- ifelse(overall_PA_protein_regression_bmi$P_lm < 0.05/1885, "associated", "not associated")
overall_PA_protein_sig <- overall_PA_protein_regression_bmi %>% 
  filter(category == "associated") %>% 
  dplyr::select(protein)


## protein-cancer tables
protein_cancer_cox_PH_bmi <- read.table(paste0(results_dir, "tables/protein_cancer_incidence_cox_age_2years.txt"), header = T, sep = "\t")
pa_prot_cancer_list <- protein_cancer_cox_PH_bmi %>%
  filter(p.value < 0.001) %>%
  mutate(protein = str_remove(protein, "^rnt_"),
         protein = toupper(protein)) %>%
  distinct(protein) %>%
  pull(protein)


## Subset PA_protein_MR
PA_protein_MR_overall <- PA_protein_MR %>%
  filter(grepl("overall_physical_activity", exposure)) %>%
  filter(outcome %in% pa_prot_cancer_list)

PA_protein_MR_overall_2SNP <- PA_protein_MR_overall %>%
  filter(grepl("overall_physical_activity_klimentidis_2_SNPs", exposure))

## Visualise the PA_protein estimates with the different instruments
# Reshape to wide format (1 row per protein, columns = exposures)
PA_protein_MR_overall_unique <- PA_protein_MR_overall %>%
  filter(str_detect(method, paste("Inverse variance weighted", collapse = "|"))) %>%
  group_by(outcome, exposure) %>%
  summarise(
    b = mean(b, na.rm = TRUE),
    se = mean(se, na.rm = TRUE),
    pval = mean(pval, na.rm = TRUE),
    .groups = "drop"
  )

PA_protein_MR_wide <- PA_protein_MR_overall_unique %>%
  pivot_wider(
    names_from = exposure,
    values_from = c(b, se, pval)
  )

## Subset protein_PA_MR
protein_PA_MR_overall <- protein_PA_MR %>%
  filter(grepl("overall_activity", outcome)) %>%
  filter(str_detect(exposure, paste(pa_prot_cancer_list, collapse = "|")))
protein_PA_MR_overall <- protein_PA_MR_overall %>%
  mutate(
    b = ifelse(outcome == "overall_activity_klimentidis_mg", b / 8.14, b),
    se = ifelse(outcome == "overall_activity_klimentidis_mg", se / 8.14, se),
  )
protein_PA_MR_overall$lower_CI <- protein_PA_MR_overall$b - 1.96 * protein_PA_MR_overall$se
protein_PA_MR_overall$upper_CI <- protein_PA_MR_overall$b + 1.96 * protein_PA_MR_overall$se
write.table(protein_PA_MR_overall, file = paste0(results_dir, "tables/MR/09_protein_PA_supp_table.txt"), col.names = T, row.names = F, sep = "\t")

## Subset PA_cancer_MR
PA_cancer_MR_overall <- PA_cancer_MR %>%
  filter(exposure %in% c("overall_physical_activity_klimentidis", 
                         "overall_physical_activity_doherty"))
PA_cancer_MR_overall <- PA_cancer_MR_overall %>%
  mutate(
    b = ifelse(exposure == "overall_physical_activity_klimentidis", b * 8.14, b),
    se = ifelse(exposure == "overall_physical_activity_klimentidis", se * 8.14, se),
    lower_CI = ifelse(exposure == "overall_physical_activity_klimentidis", b - 1.96 * se, lower_CI),
    upper_CI = ifelse(exposure == "overall_physical_activity_klimentidis", b + 1.96 * se, upper_CI)
  )
PA_cancer_MR_overall$lower_CI <- PA_cancer_MR_overall$b - 1.96 * PA_cancer_MR_overall$se
PA_cancer_MR_overall$upper_CI <- PA_cancer_MR_overall$b + 1.96 * PA_cancer_MR_overall$se
PA_cancer_MR_overall <- PA_cancer_MR_overall %>%
  mutate(
    lower_CI = b - 1.96 * se,
    upper_CI = b + 1.96 * se,
    OR = exp(b),
    OR_lower_CI = exp(lower_CI),
    OR_upper_CI = exp(upper_CI)
  )
write.table(PA_cancer_MR_overall, file = paste0(results_dir, "tables/MR/09_PA_cancer_supp_table.txt"), col.names = T, row.names = F, sep = "\t")

## Subset cancer_PA_MR
cancer_PA_MR_overall <- cancer_PA_MR %>%
  filter(grepl("overall_activity", outcome)) 
cancer_PA_MR_overall <- cancer_PA_MR_overall %>%
  mutate(
    b = ifelse(outcome == "overall_activity_klimentidis_mg", b / 8.14, b),
    se = ifelse(outcome == "overall_activity_klimentidis_mg", se / 8.14, se),
  )
cancer_PA_MR_overall$lower_CI <- cancer_PA_MR_overall$b - 1.96 * cancer_PA_MR_overall$se
cancer_PA_MR_overall$upper_CI <- cancer_PA_MR_overall$b + 1.96 * cancer_PA_MR_overall$se
write.table(cancer_PA_MR_overall, file = paste0(results_dir, "tables/MR/09_cancer_PA_supp_table.txt"), col.names = T, row.names = F, sep = "\t")

## Subset protein_cancer_MR
protein_cancer_MR_overall <- protein_cancer_MR %>%
  filter(str_detect(exposure, paste(pa_prot_cancer_list, collapse = "|")))
write.table(protein_cancer_MR_overall, file = paste0(results_dir, "tables/MR/09_protein_cancer_supp_table.txt"), col.names = T, row.names = F, sep = "\t")

## Subset cancer_protein_MR
cancer_protein_MR_overall <- cancer_protein_MR %>%
  filter(str_detect(outcome, paste(pa_prot_cancer_list, collapse = "|")))
write.table(cancer_protein_MR_overall, file = paste0(results_dir, "tables/MR/09_cancer_protein_supp_table.txt"), col.names = T, row.names = F, sep = "\t")
