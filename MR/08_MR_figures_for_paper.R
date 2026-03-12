##### Results and figures for paper #####

## libraries
library(dplyr)
library(ggplot2)
library(ggalluvial)
library(tidyr)
library(EpiViz)
library(gridExtra)
library(networkD3)
library(igraph)
library(htmlwidgets)
library(webshot)
library(ggpubr)
library(cowplot)
library(patchwork)
library(cowplot)
library(devtools)
library(ggforce) 
library(forcats)
library(stringr)

## clear environment
rm(list=ls())

## connect to project and setwd to script folder

# read in parameter file (specified on command line)
source("parameter_file/parameters_for_r.R")

## Observational overall physical activity and protein results
overall_PA_protein_regression_bmi <- read.table(paste0(results_dir, "tables/overall_PA_protein_regression_fully_adjusted_bmi.txt"), header = T, sep = "\t")
overall_PA_protein_regression_bmi$activity_type <- "overall_acceleration_average"
overall_PA_protein_regression_bmi$category <- ifelse(overall_PA_protein_regression_bmi$P_lm < 0.05/1885, "associated", "not associated")
overall_PA_protein_sig <- overall_PA_protein_regression_bmi %>% 
  filter(category == "associated") %>% 
  dplyr::select(protein)

## PA-cancer tables observational
PA_cancer_regression_bmi <- read.table(paste0(results_dir, "tables/all_cancer_PA_regression_results_filtered_fully_adjusted_bmi.txt"), header = T, sep = "\t")

## Protein and cancer Cox PH
protein_cancer_cox_PH_bmi <- read.table(paste0(results_dir, "tables/protein_cancer_incidence_cox_age_2years.txt"), header = T, sep = "\t")
pa_prot_cancer_list <- protein_cancer_cox_PH_bmi %>%
  filter(p.value < 0.001) %>%
  mutate(protein = str_remove(protein, "^rnt_"),
         protein = toupper(protein)) %>%
  distinct(protein) %>%
  pull(protein)

## read in all the MR results
PA_protein_MR <- read.table(paste0(results_dir, "tables/MR/01_PA_protein_mr_results.txt"), header = T, sep = "\t")
protein_PA_MR <- read.table(paste0(results_dir, "tables/MR/02_protein_PA_mr_results.txt"), header = T, sep = "\t")
PA_cancer_MR <- read.table(paste0(results_dir, "tables/MR/03_PA_cancer_mr_results.txt"), header = T, sep = "\t")
cancer_PA_MR <- read.table(paste0(results_dir, "tables/MR/04_cancer_PA_mr_results.txt"), header = T, sep = "\t")
protein_cancer_MR <- read.table(paste0(results_dir, "tables/MR/06_protein_cancer_mr_results.txt"), header = T, sep = "\t")
cancer_protein_MR <- read.table(paste0(results_dir, "tables/MR/07_cancer_protein_mr_results.txt"), header = T, sep = "\t")

## plot overall PA-cancer results and restrict to IVW
PA_cancer_overall <- PA_cancer_MR %>%
  filter(grepl("overall_physical_activity", exposure)) %>%
  filter(grepl("Inverse variance", method)) %>%
  arrange(OR)

PA_cancer_overall <- PA_cancer_overall %>% 
  mutate(cancer_outcome_simple = sub(" cancer.*", " cancer", outcome),
         exposure_name = exposure)

PA_cancer_MR_overall_sig <-PA_cancer_overall  %>% 
  filter(pval < 0.05)

max_CI <- max(c(abs(log(PA_cancer_overall$OR_lower_CI)),
                abs(log(PA_cancer_overall$OR_upper_CI))), na.rm = TRUE)

PA_cancer_MR_overall_2SNP <- PA_cancer_overall %>%
  filter(exposure == "overall_physical_activity_klimentidis")
PA_cancer_MR_overall_2SNP <- PA_cancer_MR_overall_2SNP %>%
  mutate(exposure = "Overall physical activity (2 SNP instrument)")

# Scale b and se
PA_cancer_MR_overall_2SNP$b_SD <- PA_cancer_MR_overall_2SNP$b * 8.14
PA_cancer_MR_overall_2SNP$se_SD <- PA_cancer_MR_overall_2SNP$se * 8.14

# Recalculate OR and 95% CI
PA_cancer_MR_overall_2SNP$OR_SD <- exp(PA_cancer_MR_overall_2SNP$b_SD)
PA_cancer_MR_overall_2SNP$OR_lower_CI_SD <- exp(PA_cancer_MR_overall_2SNP$b_SD - 1.96 * PA_cancer_MR_overall_2SNP$se_SD)
PA_cancer_MR_overall_2SNP$OR_upper_CI_SD <- exp(PA_cancer_MR_overall_2SNP$b_SD + 1.96 * PA_cancer_MR_overall_2SNP$se_SD)

pdf(paste0(results_dir, "figures/MR_figures/08_PA_cancer_IVW_2_SNP.pdf"), height = 4, width = 8)
ggplot(PA_cancer_MR_overall_2SNP, aes(x = fct_rev(factor(exposure)), 
                               y = OR_SD, 
                               ymin = OR_lower_CI_SD, 
                               ymax = OR_upper_CI_SD,
                               color = cancer_outcome_simple)) +
  geom_pointrange(position = position_dodge(width = 0.6)) +
  coord_flip() +
  geom_hline(yintercept = 1, linetype = "dashed") +
  labs(x = "Exposure", y = "Odds Ratio per SD higher overall acceleration average (95% CI)", color = "Cancer Type") +
  theme_minimal(base_size = 14)
dev.off()
