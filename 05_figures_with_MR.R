### Figures plus MR results ####
## 12/2/25

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
library(EpiViz)
library(ggrepel)
library(scales)

## clear environment
rm(list=ls())

## connect to project and setwd to script folder

# read in parameter file (specified on command line)
source("parameter_file/parameters_for_r.R")

## PA-protein tables logistic regression
overall_PA_protein_regression_bmi <- read.table(paste0(results_dir, "tables/overall_PA_protein_regression_fully_adjusted_bmi.txt"), header = T, sep = "\t")
overall_PA_protein_regression_bmi$activity_type <- "overall_acceleration_average"

## protein-cancer table cox PH
protein_cancer_cox_bmi <- read.table(paste0(results_dir, "tables/protein_cancer_incidence_cox_age_2years.txt"), header = T, sep = "\t")

## Set threshold for associations in each dataframe
overall_PA_protein_regression_bmi$category <- ifelse(overall_PA_protein_regression_bmi$P_lm < 0.05/1885, "associated", "not associated")

## Forest plot for overall physical activity
overall_PA_protein_regression_bmi_plot <- overall_PA_protein_regression_bmi
overall_PA_protein_regression_bmi_plot$protein <- gsub("rnt_", "", overall_PA_protein_regression_bmi_plot$protein)
pdf(file = paste0(results_dir, "figures/volcano_overall_PA_protein_bmi.pdf"), height = 6, width = 6)
overall_PA_protein_regression_bmi_plot %>%
  filter(n_lm != 34) %>%
  mutate(
    neg_log10_p = -log10(P_lm)
  ) %>%
  ggplot(aes(x = beta_lm, y = neg_log10_p, colour = category)) +
  geom_point(alpha = 0.7, size = 1.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "grey50") +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_text_repel(
    data = function(d) d %>% filter(neg_log10_p > 30),
    aes(label = protein),
    size = 3,
    alpha = 0.7,
    nudge_y = 1,           # move labels slightly above points
    box.padding = 0.3,     # spacing between labels
    point.padding = 0.2,   # spacing from points
    segment.color = "grey50",
    max.overlaps = Inf,
    show.legend = FALSE
  ) +
  scale_colour_brewer(palette = "Set1") +
  labs(
    x = "Normalised SD difference in protein levels\nper normalised SD higher overall acceleration average physical activity",
    y = expression(-log[10](P)),
    colour = "Category",
    title = ""
  ) +
  coord_cartesian(xlim = c(-0.3, NA)) +
  theme_minimal()
dev.off()


## plot ALL proteins associated with overall physical activity if they have an association with cancer
## Use the Cox PH estimate
# --- Start with PA → Protein data ---
overall_PA_protein <- bind_rows(
  overall_PA_protein_regression_bmi %>%
    filter(P_lm < 0.05 / 1885) %>%
    dplyr::select(protein, activity_type, beta_lm, P_lm, se_lm)
) %>%
  rename(
    PA = activity_type,
    Protein = protein,
    Beta = beta_lm,
    P_value = P_lm
  )

overall_PA_protein$Protein <- gsub("rnt_", "", overall_PA_protein$Protein)
overall_PA_protein$PA <- recode(overall_PA_protein$PA, "overall_acceleration_average" = "overall")

# --- Prepare cancer dataset ---
protein_cancer_cox_bmi <- protein_cancer_cox_bmi %>%
  mutate(
    Protein = sub("^rnt_", "", protein)
  )
protein_cancer_PA_associated_proteins <- protein_cancer_cox_bmi %>%
  filter(p.value < 0.001) 
PA_cancer_proteins <- protein_cancer_PA_associated_proteins$Protein

# --- Merge PA-associated proteins with all cancer estimates ---
data_for_plot_overall <- merge(
  overall_PA_protein %>% filter(Protein %in% PA_cancer_proteins),
  protein_cancer_cox_bmi,  # merged Cox PH results
  by = "Protein"
)

# This ensures we include incident estimates even if prevalence is the only significant
data_long_overall <- bind_rows(
  
  # ---- PA → Protein ----
  data_for_plot_overall %>%
    distinct(Protein, PA, Beta, se_lm) %>%   # avoid duplication across cancers
    transmute(
      Protein,
      PA,
      cancer_type = NA,
      Type = "PA → Protein",
      Estimate = Beta,
      CI_lower = Beta - 1.96 * se_lm,
      CI_upper = Beta + 1.96 * se_lm,
      p.signif = NA
    ),
  
  # ---- Protein → Cancer (incident Cox PH) ----
  data_for_plot_overall %>%
    transmute(
      Protein,
      PA = NA,
      cancer_type,
      Type = "Protein → Cancer",
      Estimate = HR,      # log(HR)
      CI_lower = HR_conf_low,
      CI_upper = HR_conf_high,
      p.signif = p.value < 0.001
    )
)

sig_cancer_proteins <- data_long_overall %>%
  filter(
    Type == "Protein → Cancer",
    p.signif = TRUE,
  ) %>%
  pull(Protein) %>%
  unique()

data_long_overall <- data_long_overall %>%
  filter(
    Protein %in% sig_cancer_proteins,
    # keep only significant cancer rows
    Type == "PA → Protein" |
      (Type == "Protein → Cancer" & p.signif == TRUE)
  )

# --- Ensure cancer labels match your forest plot ---
# Define cancer colors to match your desired scheme
gg_colors <- hue_pal()(4)  # gives 4 default ggplot colors

# Assign to cancers: Endometrial, Colorectal, Prostate, Breast
cancer_colors <- c(
  "Endometrial" = gg_colors[1],  # default blue
  "Colorectal"  = gg_colors[2],  # default green
  "Prostate"    = gg_colors[3],  # default pink/purple
  "Breast"      = gg_colors[4]   # default red
)

# Make sure cancer_label matches exactly
data_long_overall <- data_long_overall %>%
  mutate(
    cancer_label = recode(
      cancer_type,
      PC  = "Prostate",
      BC  = "Breast",
      EC  = "Endometrial",
      CRC = "Colorectal"
    )
  )

# --- PA → Protein plot ---
plot_pa_protein_overall <- ggplot(
  data_long_overall %>% filter(Type == "PA → Protein"),
  aes(x = Estimate, y = Protein, color = PA)
) +
  geom_point(size = 3, position = position_dodge(width = 0.7)) +
  geom_errorbarh(aes(xmin = CI_lower, xmax = CI_upper), height = 0.2, position = position_dodge(width = 0.7)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  labs(
    x = "Difference in protein (SDs) per SD unit higher physical activity (±95% CI)",
    y = "Proteins",
    title = ""
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 12, face = "bold"),
    axis.title.y = element_blank(),
    legend.position = "none"  # <-- remove legend
  )

# --- Protein → Cancer plot (HR scale) ---
plot_protein_cancer_overall <- ggplot(
  data_long_overall %>% filter(Type == "Protein → Cancer"),
  aes(
    x = Estimate,         # HR
    y = Protein,
    color = cancer_label
  )
) +
  geom_point(size = 3, position = position_dodge(width = 0.7)) +
  geom_errorbarh(
    aes(xmin = CI_lower, xmax = CI_upper),
    height = 0.2,
    position = position_dodge(width = 0.7)
  ) +
  scale_color_manual(
    values = cancer_colors,
    name = "Cancer site",
    guide = guide_legend(nrow = 2, byrow = TRUE)
  ) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "black") +
  labs(
    x = "Hazard ratio per SD higher protein (±95% CI)",
    y = ""
  ) +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.y = element_text(size = 12, face = "bold"),
    axis.title.y = element_blank(),
    legend.position = "bottom",
    legend.title = element_text(size = 10, face = "bold")
  )

# --- Combine side by side ---
plot_overall <- plot_pa_protein_overall | plot_protein_cancer_overall
plot_overall

ggsave(
  paste0(results_dir, "figures/forest_overal_PA_protein_and_protein_cancer_v2.pdf"), 
  plot = plot_overall,  # use the correct object
  width = 14, 
  height = 8
)

## Do any of these have MR evidence for the PA -> protein step? Using proteins only associated with overall physical activity.
PA_protein_MR <- read.table(paste0(results_dir, "tables/MR/01_PA_protein_mr_results.txt"), header = T, sep = "\t")
PA_protein_MR$outcome <- tolower(PA_protein_MR$outcome)
PA_protein_MR$pval[PA_protein_MR$pval == 0] <- .Machine$double.xmin
PA_protein_MR <- PA_protein_MR %>%
  filter(exposure %in% c("overall_physical_activity_klimentidis_2_SNPs", "overall_physical_activity_doherty"))
prots <- unique(data_long_overall$Protein)
proteins_MR_rows <- which(PA_protein_MR$outcome %in% prots)
proteins_MR <- PA_protein_MR[proteins_MR_rows,]
proteins_MR$b <- ifelse(
  proteins_MR$exposure == "overall_physical_activity_klimentidis_2_SNPs",
  proteins_MR$b * 8.14,
  proteins_MR$b
)

proteins_MR$se <- ifelse(
  proteins_MR$exposure == "overall_physical_activity_klimentidis_2_SNPs",
  proteins_MR$se * 8.14,
  proteins_MR$se
)

## Write this to file 
write.table(proteins_MR, file = paste0(results_dir, "tables/MR/PA_protein_MR_final_table.txt"),
            row.names = FALSE, col.names = TRUE, sep = "\t")

## Subset to method for plot
proteins_MR <- proteins_MR %>%
  filter(method %in% c("Inverse variance weighted (fixed effects)"))

# **PA → Protein MR plot**
protein_levels <- data_long_overall %>%
  filter(Type == "Protein → Cancer") %>%
  distinct(Protein) %>%
  pull(Protein) %>%
  sort(decreasing=TRUE) %>%
  rev()

# PA → Protein MR plot
plot_pa_protein_mr_overall <- ggplot(proteins_MR, 
                                     aes(x = b, y = outcome, color = exposure, shape = exposure)) +  
  geom_point(size = 3, position = position_dodge(width = 0.7)) +  
  geom_errorbarh(aes(xmin = b - 1.96*se, xmax = b + 1.96*se), height = 0.2, position = position_dodge(width = 0.7)) +  
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +   
  scale_y_discrete(limits = protein_levels, drop = FALSE) +
  labs(
    x = "Difference in protein (SD ±95% CI) per SD higher overall physical activity",
    y = "",
    title = ""
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 12, face = "bold"),
    axis.title.y = element_blank(),
    legend.position = "none",
    legend.title = element_text(size = 10, face = "bold")
  )


plot_pa_protein_mr_overall 

# Save plot
ggsave(paste0(results_dir, "figures/MR_figures/PA_protein_MR_forest.pdf"), 
       plot = plot_pa_protein_mr_overall, width = 15, height = 8)

## Save these gene names
genes_pa_protein <- unique(data_long_overall %>% 
                             filter(Type == "PA → Protein") %>% 
                             pull(Protein))

genes_protein_cancer <- unique(data_long_overall %>% 
                                 filter(Type == "Protein → Cancer") %>% 
                                 pull(Protein))
genes_pa_protein_mr <- unique(proteins_MR$outcome)
all_genes <- toupper(unique(c(genes_pa_protein, genes_protein_cancer, genes_pa_protein_mr)))
write.table(all_genes, file = paste0(results_dir, "tables/genes_concordant_PA_protein_protein_cancer.txt"),
            row.names = FALSE, col.names = FALSE, quote = FALSE)

## Protein->cancer MR results
# Read MR results
protein_cancer_MR <- read.table(
  paste0(results_dir, "tables/MR/06_protein_cancer_mr_results.txt"),
  header = TRUE,
  sep = "\t"
)

# Clean MR results 
protein_cancer_MR_clean <- protein_cancer_MR %>% 
  mutate( Gene = tolower(sub("_.*$", "", exposure)), 
          cancer_label = case_when( grepl("colorectal", outcome, ignore.case = TRUE) ~ "Colorectal", grepl("breast", outcome, ignore.case = TRUE) ~ "Breast", grepl("prostate", outcome, ignore.case = TRUE) ~ "Prostate", grepl("endometrial", outcome, ignore.case = TRUE)~ "Endometrial" ) )

# Valid protein–cancer pairs from downstream analysis 
valid_pairs <- data_long_overall %>% 
  filter(Type == "Protein → Cancer") %>% 
  mutate(Protein = tolower(Protein)) %>% 
  dplyr::select(Protein, cancer_label) %>% 
  distinct() 

# Restrict to proteins of interest AND valid cancer outcomes 
protein_cancer_MR_filtered <- protein_cancer_MR_clean %>% 
  filter(Gene %in% tolower(all_genes)) %>% 
  inner_join( valid_pairs, by = c("Gene" = "Protein", "cancer_label") )

## Write this table to file
write.table(protein_cancer_MR_filtered, file = paste0(results_dir, "tables/MR/protein_cancer_MR_final_table.txt"),
            row.names = FALSE, col.names = TRUE, sep = "\t")

## Subset cancer->protein 
cancer_protein_MR <- read.table(paste0(results_dir, "tables/MR/07_cancer_protein_mr_results.txt"), header = T, sep = "\t")

cancer_protein_MR_clean <- cancer_protein_MR %>%
  mutate(
    Protein = tolower(outcome),
    cancer_label = case_when(
      grepl("colorectal", exposure, ignore.case = TRUE) ~ "Colorectal",
      grepl("breast", exposure, ignore.case = TRUE)     ~ "Breast",
      grepl("prostate", exposure, ignore.case = TRUE)   ~ "Prostate",
      grepl("endometrial", exposure, ignore.case = TRUE)~ "Endometrial"
    )
  )

cancer_protein_MR_subset <- cancer_protein_MR_clean %>%
  inner_join(
    valid_pairs,
    by = c("Protein", "cancer_label")
  )

write.table(cancer_protein_MR_subset, file = paste0(results_dir, "tables/MR/cancer_protein_MR_final_table.txt"),
            row.names = FALSE, col.names = TRUE, sep = "\t")

# Protein → Cancer MR plot
protein_cancer_plot <- protein_cancer_MR_filtered %>%
  mutate(
    Gene = factor(Gene, levels = protein_levels),
    cancer = case_when(
      grepl("colorectal", outcome, ignore.case = TRUE) ~ "Colorectal",
      grepl("breast", outcome, ignore.case = TRUE)     ~ "Breast",
      grepl("prostate", outcome, ignore.case = TRUE)   ~ "Prostate",
      grepl("endometrial", outcome, ignore.case = TRUE)~ "Endometrial",
      TRUE ~ "Other"
    )
  )

plot_protein_cancer_mr <- ggplot(
  protein_cancer_plot,
  aes(x = OR, y = Gene, color = cancer)
) +
  geom_point(size = 3, position = position_dodge(width = 0.6)) +
  geom_errorbarh(
    aes(xmin = OR_lower_CI, xmax = OR_upper_CI),
    height = 0.2,
    position = position_dodge(width = 0.6)
  ) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "black") +
  scale_x_log10() +
  scale_y_discrete(limits = protein_levels, drop = FALSE) +
  scale_color_manual(values = cancer_colors) +
  labs(
    x = "Odds ratio for cancer per SD higher protein (±95% CI)",
    y = "Protein",
    color = "Cancer site",
    title = ""
  ) +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.y = element_text(size = 11, face = "bold"),
    axis.title.y = element_blank(),
    legend.position = "bottom"
  )

# Combine plots side by side
plot_MR <- plot_pa_protein_mr_overall | plot_protein_cancer_mr

# Save combined figure
ggsave(
  paste0(results_dir, "figures/forest_MR_overall_PA_protein_and_protein_cancer.pdf"), 
  plot = plot_MR,
  width = 14, 
  height = 8
)

## contatenate all heterogeneity tables into one
# List files with "het" in the filename
# Define file paths
files <- list(
  PA_protein = paste0(results_dir, "tables/MR/01_PA_protein_heterogeneity.txt"),
  PA_cancer = paste0(results_dir, "tables/MR/03_PA_cancer_heterogeneity.txt"),
  cancer_PA = paste0(results_dir, "tables/MR/04_cancer_PA_heterogeneity.txt"),
  cancer_protein = paste0(results_dir, "tables/MR/07_cancer_protein_heterogeneity.txt")
)

# PA -> Protein
PA_protein_results <- read.table(files$PA_protein, header = TRUE, sep = "\t", stringsAsFactors = FALSE) %>%
  filter(tolower(outcome) %in% valid_pairs$Protein) %>%
  filter(exposure %in% c(
    "overall_physical_activity_klimentidis_2_SNPs",
    "overall_physical_activity_doherty"
  ))


# PA -> Cancer
PA_cancer_results <- read.table(files$PA_cancer, header = TRUE, sep = "\t", stringsAsFactors = FALSE) %>%
  filter(exposure %in% c(
    "overall_physical_activity_klimentidis",
    "overall_physical_activity_doherty"
  ))

# Cancer -> PA
cancer_PA_results <- read.table(files$cancer_PA, header = TRUE, sep = "\t", stringsAsFactors = FALSE) %>%
  filter(outcome %in% c(
    "overall_activity_doherty_SD",
    "overall_activity_klimentidis_mg"
  ))

# Cancer -> Protein
cancer_protein_results <- read.table(files$cancer_protein, header = TRUE, sep = "\t", stringsAsFactors = FALSE) %>%
  mutate(
    Protein = tolower(outcome),
    cancer_label = case_when(
      grepl("Prostate cancer", exposure) ~ "Prostate",
      grepl("Breast cancer", exposure) ~ "Breast",
      grepl("colorectal cancer", exposure, ignore.case = TRUE) ~ "Colorectal",
      grepl("Endometrial cancer", exposure) ~ "Endometrial",
      TRUE ~ NA_character_
    )
  ) %>%
  semi_join(valid_pairs, by = c("Protein", "cancer_label"))

# Optional: concatenate all four into one dataframe if needed
all_het_results <- bind_rows(
  PA_protein_results,
  PA_cancer_results,
  cancer_PA_results,
  cancer_protein_results
)

write.table(all_het_results, file = paste0(results_dir, "tables/MR/MR_het_final.txt"),
                        row.names = FALSE, col.names = TRUE, sep = "\t")
