### Figures from result files ###
### Lucy Goudswaard - 30th July 2024 ###

## load packages
library(data.table)
library(psych)
library(dplyr)
library(tidyr)
library(ggplot2)
library(reshape2)
library(moosefun)
library(broom)
library(readxl)
library(ggupset)
library(tibble)
library(purrr)
library(RColorBrewer) 

## clear environment
rm(list=ls())

## connect to project and setwd to script folder

# read in parameter file (specified on command line)
source("parameter_file/parameters_for_r.R")

## Read in regression results
mod_vig_protein_regression <- read.table(paste0(results_dir, "tables/mod_vig_PA_quantile_numeric_protein_regression_fully_adjusted.txt"), header = T, sep = "\t")
overall_PA_protein_regression <- read.table(paste0(results_dir, "tables/overall_PA_protein_regression_fully_adjusted.txt"), header = T, sep = "\t")
sedentary_PA_protein_regression <- read.table(paste0(results_dir, "tables/sedentary_PA_protein_regression_full_adjusted.txt"), header = T, sep = "\t")
olink_info <- read.csv(paste0(input_dir, "UKB_PPP_olink_info.csv"), header = T)
olink_info$Assay.Target <- tolower(olink_info$Assay.Target)
olink_gene_panel <- olink_info[, c("Protein.panel", "Assay.Target")]

# Create dataframe for upset plot
mod_vig_protein_regression$mod_vig_PA <- ifelse(mod_vig_protein_regression$P_lm < 0.05/1885, "mod_vig_PA", NA)
mod_vig_plot_results <- mod_vig_protein_regression[, c("protein", "mod_vig_PA")]

overall_PA_protein_regression$overall_PA <- ifelse(overall_PA_protein_regression$P_lm < 0.05/1885, "overall_PA", NA)
overall_plot_results <- overall_PA_protein_regression[, c("protein", "overall_PA")]

sedentary_PA_protein_regression$sedentary_PA <- ifelse(sedentary_PA_protein_regression$P_lm < 0.05/1885, "sedentary_PA", NA)
sedentary_plot_results <- sedentary_PA_protein_regression[, c("protein", "sedentary_PA")]

# Merge the dataframes
PA_protein_results_1 <- merge(mod_vig_plot_results, overall_plot_results, by = "protein")
PA_protein_results <- merge(PA_protein_results_1, sedentary_plot_results, by = "protein")

# Replace NA with NA_character_
PA_protein_results$mod_vig_PA <- ifelse(is.na(PA_protein_results$mod_vig_PA), NA_character_, PA_protein_results$mod_vig_PA)
PA_protein_results$overall_PA <- ifelse(is.na(PA_protein_results$overall_PA), NA_character_, PA_protein_results$overall_PA)
PA_protein_results$sedentary_PA <- ifelse(is.na(PA_protein_results$sedentary_PA), NA_character_, PA_protein_results$sedentary_PA)

# Combine changes into a single column and count the occurrences of each combination
PA_protein_results_combined <- PA_protein_results %>%
  mutate(changes = pmap(list(mod_vig_PA, overall_PA, sedentary_PA), 
                        ~paste(sort(na.omit(c(..1, ..2, ..3))), collapse = ", "))) %>%
  count(changes) %>%
  mutate(changes = factor(changes, levels = unique(changes))) # Ensure the factor levels are in order
PA_protein_results_combined <- PA_protein_results_combined %>%
  arrange(desc(n))

# Print the combined dataframe to check
print(PA_protein_results_combined)
str(PA_protein_results_combined)

PA_protein_results_combined <- PA_protein_results_combined %>%
  mutate(changes = as.character(changes)) %>%  # Convert factor to character
  mutate(changes = ifelse(changes == "", NA, changes)) %>%  # Handle empty string
  mutate(changes = strsplit(changes, ", "))  # Convert to list of vectors

# Plot for combined changes
combined_plot <- ggplot(PA_protein_results_combined, aes(x = changes, y = n)) +
  geom_bar(stat = 'identity', fill = "#8278BF") +
  scale_x_upset() +
  geom_text(aes(label = n), vjust = -0.5) +
  scale_y_continuous(name = "Protein Count") +
  labs(title = "Proteins associated with physical activity",
       x = "Physical activity measurements",
       y = "Protein Count") +
  theme(plot.margin = margin(0.4, 0.4, 0.4, 1.5, "cm"),
        axis.text.x = element_text(angle = 45, hjust = 1))  # Adjust text angle for readability

# Print the plot
print(combined_plot)

# save plot to pdf
pdf(file = paste0(results_dir, "figures/upset_plot_PA_results.pdf"), height = 6, width = 6)
combined_plot
dev.off()

## Remake plots adjusting for BMI
## Read in regression results with "_bmi" added
mod_vig_protein_regression_bmi <- read.table(paste0(results_dir, "tables/mod_vig_PA_quantile_numeric_protein_regression_fully_adjusted_bmi.txt"), header = T, sep = "\t")
overall_PA_protein_regression_bmi <- read.table(paste0(results_dir, "tables/overall_PA_protein_regression_fully_adjusted_bmi.txt"), header = T, sep = "\t")
sedentary_PA_protein_regression_bmi <- read.table(paste0(results_dir, "tables/sedentary_PA_protein_regression_full_adjusted_bmi.txt"), header = T, sep = "\t")

# Create dataframe for upset plot with "_bmi" results
mod_vig_protein_regression_bmi$mod_vig_PA <- ifelse(mod_vig_protein_regression_bmi$P_lm < 0.05/1885, "mod_vig_PA", NA)
mod_vig_plot_results_bmi <- mod_vig_protein_regression_bmi[, c("protein", "mod_vig_PA")]

overall_PA_protein_regression_bmi$overall_PA <- ifelse(overall_PA_protein_regression_bmi$P_lm < 0.05/1885, "overall_PA", NA)
overall_plot_results_bmi <- overall_PA_protein_regression_bmi[, c("protein", "overall_PA")]

sedentary_PA_protein_regression_bmi$sedentary_PA <- ifelse(sedentary_PA_protein_regression_bmi$P_lm < 0.05/1885, "sedentary_PA", NA)
sedentary_plot_results_bmi <- sedentary_PA_protein_regression_bmi[, c("protein", "sedentary_PA")]

# Merge the dataframes
PA_protein_results_1_bmi <- merge(mod_vig_plot_results_bmi, overall_plot_results_bmi, by = "protein")
PA_protein_results_bmi <- merge(PA_protein_results_1_bmi, sedentary_plot_results_bmi, by = "protein")

## make upset plot with PA and protein results, adjusted for BMI but colour coded
PA_protein_results_bmi_upset <- PA_protein_results_bmi
PA_protein_results_bmi_upset$protein <- gsub("^rnt_", "", PA_protein_results_bmi_upset$protein)
PA_protein_results_bmi_upset <- merge(PA_protein_results_bmi_upset, olink_gene_panel, by.x = "protein", by.y = "Assay.Target", all.x = TRUE)

#### Select other rows to drop 
PA_protein_results_bmi_upset <- PA_protein_results_bmi_upset %>%
  filter(!(protein == "tnf" & Protein.panel %in% c("Oncology", "Inflammation", "Neurology")) & 
           !(protein == "il6" & Protein.panel %in% c("Cardiometabolic", "Neurology", "Inflammation")) & 
           !(protein == "ido1" & Protein.panel %in% c("Oncology_II", "Neurology_II", "Inflammation_II")) & 
           !(protein == "lmod1" & Protein.panel %in% c("Oncology_II", "Cardiometabolic_II", "Inflammation_II")) & 
           !(protein == "scrib" & Protein.panel %in% c("Oncology_II", "Neurology_II", "Inflammation_II")) &
           !(protein == "cxcl8" & Protein.panel %in% c("Cardiometabolic", "Neurology", "Inflammation")))


# Handle NA for physical activity columns
PA_protein_results_bmi_upset$mod_vig_PA <- ifelse(is.na(PA_protein_results_bmi_upset$mod_vig_PA), NA_character_, PA_protein_results_bmi_upset$mod_vig_PA)
PA_protein_results_bmi_upset$overall_PA <- ifelse(is.na(PA_protein_results_bmi_upset$overall_PA), NA_character_, PA_protein_results_bmi_upset$overall_PA)
PA_protein_results_bmi_upset$sedentary_PA <- ifelse(is.na(PA_protein_results_bmi_upset$sedentary_PA), NA_character_, PA_protein_results_bmi_upset$sedentary_PA)

PA_protein_results_with_panel_upset <- PA_protein_results_bmi_upset %>%
  mutate(changes = pmap_chr(list(mod_vig_PA, overall_PA, sedentary_PA), function(...){
    activity_list <- na.omit(c(...))
    if (length(activity_list) > 0) {
      paste(sort(activity_list), collapse = ", ")
    } else {
      NA_character_ # Set as NA if all values are NA
    }
  }))

print(PA_protein_results_with_panel_upset)
str(PA_protein_results_with_panel_upset)

PA_protein_results_bmi_upset <- PA_protein_results_with_panel_upset %>%
  # First, calculate `n` for each `changes` value
  count(changes, Protein.panel, name = "panel_count") %>%
  group_by(changes) %>%
  mutate(n = sum(panel_count)) %>%              # Sum panel counts to get total `n` for each `changes`
  mutate(proportion = panel_count / n) %>%      # Calculate proportion for each panel
  dplyr::select(-panel_count) %>%               # Drop intermediate `panel_count` column
  pivot_wider(names_from = Protein.panel,       # Pivot to create columns for each panel with proportions
              values_from = proportion, 
              values_fill = list(proportion = 0)) %>%  # Fill missing values with 0
  ungroup()

# Print the resulting dataframe to verify
print(PA_protein_results_bmi_upset)
str(PA_protein_results_bmi_upset)
PA_protein_results_bmi_upset <- PA_protein_results_bmi_upset %>%
  mutate(changes = strsplit(as.character(changes), ", "))

PA_protein_results_long_upset <- PA_protein_results_bmi_upset %>%
  pivot_longer(cols = c(Cardiometabolic, Cardiometabolic_II, Inflammation, Inflammation_II, Neurology, Neurology_II, Oncology, Oncology_II),  
               names_to = "Protein.panel",
               values_to = "proportion") %>%
  filter(proportion > 0) %>%  # Keep only the rows where the proportion is greater than 0
  mutate(count = proportion * n)  # Calculate the actual count for each panel

# Create a summary dataframe for total counts
total_counts_upset <- PA_protein_results_long_upset %>%
  group_by(changes) %>%
  summarise(total = sum(count), n = first(n), .groups = "drop")  # Get total and original count for each changes group
#total_counts_upset <- total_counts$n

# Create the plot
combined_plot_bmi_coloured_2 <- ggplot(PA_protein_results_long_upset, aes(x = changes, y = count, fill = Protein.panel)) +
  geom_bar(stat = 'identity', position = "stack") +  # Stacked bar chart to show panel contributions
  # geom_text(aes(label = count), position = position_stack(vjust = 0.5), color = "white") +  # Center the text on the bars
  scale_x_upset() +
  scale_y_continuous(name = "Protein Count") +
  labs(title = "Proteins associated with physical activity (BMI adjusted)",
       x = "Physical activity measurements",
       y = "Protein Count") +
  theme(plot.margin = margin(0.4, 0.4, 0.4, 1.5, "cm"),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_brewer(palette = "Set3")  # Adjust the color palette as desired

# Print the plot
print(combined_plot_bmi_coloured_2)

# Save plot to pdf
pdf(file = paste0(results_dir, "figures/upset_plot_PA_results_bmi_coloured.pdf"), height = 6, width = 6)
combined_plot_bmi_coloured_2
dev.off()

# Replace NA with NA_character_
PA_protein_results_bmi$mod_vig_PA <- ifelse(is.na(PA_protein_results_bmi$mod_vig_PA), NA_character_, PA_protein_results_bmi$mod_vig_PA)
PA_protein_results_bmi$overall_PA <- ifelse(is.na(PA_protein_results_bmi$overall_PA), NA_character_, PA_protein_results_bmi$overall_PA)
PA_protein_results_bmi$sedentary_PA <- ifelse(is.na(PA_protein_results_bmi$sedentary_PA), NA_character_, PA_protein_results_bmi$sedentary_PA)

# Combine changes into a single column and count the occurrences of each combination
PA_protein_results_combined_bmi <- PA_protein_results_bmi %>%
  mutate(changes = pmap(list(mod_vig_PA, overall_PA, sedentary_PA), 
                        ~paste(sort(na.omit(c(..1, ..2, ..3))), collapse = ", "))) %>%
  count(changes) %>%
  mutate(changes = factor(changes, levels = unique(changes))) # Ensure the factor levels are in order
PA_protein_results_combined_bmi <- PA_protein_results_combined_bmi %>%
  arrange(desc(n))

# Print the combined dataframe to check
print(PA_protein_results_combined_bmi)
str(PA_protein_results_combined_bmi)

PA_protein_results_combined_bmi <- PA_protein_results_combined_bmi %>%
  mutate(changes = as.character(changes)) %>%  # Convert factor to character
  mutate(changes = ifelse(changes == "", NA, changes)) %>%  # Handle empty string
  mutate(changes = strsplit(changes, ", "))  # Convert to list of vectors

# Plot for combined changes for "_bmi"
combined_plot_bmi <- ggplot(PA_protein_results_combined_bmi, aes(x = changes, y = n)) +
  geom_bar(stat = 'identity', fill = "#8278BF") +
  scale_x_upset() +
  geom_text(aes(label = n), vjust = -0.5) +
  scale_y_continuous(name = "Protein Count") +
  labs(title = "Proteins associated with physical activity (BMI adjusted)",
       x = "Physical activity measurements",
       y = "Protein Count") +
  theme(plot.margin = margin(0.4, 0.4, 0.4, 1.5, "cm"),
        axis.text.x = element_text(angle = 45, hjust = 1))  # Adjust text angle for readability

# Print the plot
print(combined_plot_bmi)

# Save plot to pdf
pdf(file = paste0(results_dir, "figures/upset_plot_PA_results_bmi.pdf"), height = 6, width = 6)
combined_plot_bmi
dev.off()

## Save a list of proteins that are altered by physical activity adjusted for BMI and the direction of effect
# Filter proteins altered by physical activity adjusted for BMI
proteins_all_PA_bmi <- PA_protein_results_bmi %>%
  filter(mod_vig_PA == "mod_vig_PA" & 
           overall_PA == "overall_PA" & 
           sedentary_PA == "sedentary_PA")

# Extract the protein names
names_PA_prot_bmi <- proteins_all_PA_bmi$protein

# Start building the data frame with the protein names
PA_prot_direction_of_effect <- data.frame(Protein = names_PA_prot_bmi)

# Add beta and standard error from mod_vig_protein_regression_bmi
PA_prot_direction_of_effect <- PA_prot_direction_of_effect %>%
  left_join(
    mod_vig_protein_regression_bmi %>%
      filter(protein %in% names_PA_prot_bmi) %>%
      dplyr::select(protein, beta_lm, se_lm) %>%
      rename(mod_vig_beta = beta_lm, mod_vig_se = se_lm),
    by = c("Protein" = "protein")
  )

# Add beta and standard error from overall_PA_protein_regression_bmi
PA_prot_direction_of_effect <- PA_prot_direction_of_effect %>%
  left_join(
    overall_PA_protein_regression_bmi %>%
      filter(protein %in% names_PA_prot_bmi) %>%
      dplyr::select(protein, beta_lm, se_lm) %>%
      rename(overall_PA_beta = beta_lm, overall_PA_se = se_lm),
    by = c("Protein" = "protein")
  )

# Add beta and standard error from sedentary_PA_protein_regression_bmi
PA_prot_direction_of_effect <- PA_prot_direction_of_effect %>%
  left_join(
    sedentary_PA_protein_regression_bmi %>%
      filter(protein %in% names_PA_prot_bmi) %>%
      dplyr::select(protein, beta_lm, se_lm) %>%
      rename(sedentary_PA_beta = beta_lm, sedentary_PA_se = se_lm),
    by = c("Protein" = "protein")
  )

# Print the resulting data frame
print(PA_prot_direction_of_effect)

# Print the resulting data frame
print(PA_prot_direction_of_effect)
PA_prot_direction_of_effect$estimates_agreed <- ifelse(PA_prot_direction_of_effect$mod_vig_beta > 0 & 
                                                         PA_prot_direction_of_effect$overall_PA_beta > 0 & 
                                                         PA_prot_direction_of_effect$sedentary_PA_beta < 0, "Yes", 
                                                       ifelse(PA_prot_direction_of_effect$mod_vig_beta < 0 &
                                                                PA_prot_direction_of_effect$overall_PA_beta < 0 &
                                                                PA_prot_direction_of_effect$sedentary_PA_beta > 0, "Yes", "No"))
write.table(PA_prot_direction_of_effect, file = paste0(results_dir, "tables/PA_prot_consistent_bmi_adjusted.txt"), col.names = T, row.names = F, sep = "\t")

## Plot these results
# Calculate 95% Confidence Intervals
PA_prot_direction_of_effect <- PA_prot_direction_of_effect %>%
  mutate(
    mod_vig_ci = 1.96 * mod_vig_se,
    overall_PA_ci = 1.96 * overall_PA_se,
    sedentary_PA_ci =  1.96 * sedentary_PA_se,
  )

# Reshape the data for plotting and combine lower and upper CI into one column
PA_prot_long <- PA_prot_direction_of_effect %>%
  dplyr::select(
    Protein,
    mod_vig_beta, mod_vig_se, mod_vig_ci,
    overall_PA_beta, overall_PA_se, overall_PA_ci,
    sedentary_PA_beta, sedentary_PA_se, sedentary_PA_ci
  ) %>%
  pivot_longer(
    cols = -Protein,
    names_to = c("measure", "stat"),
    names_pattern = "(.*)_(.*)",
    values_to = "value"
  ) %>%
  # Reshape further to separate beta, se, lower_ci, upper_ci
  pivot_wider(
    names_from = stat,
    values_from = value,
    names_prefix = ""
  )

# Check the structure of the reshaped data
str(PA_prot_long)  # Verify the names and structure
PA_prot_long$upper_ci <- PA_prot_long$beta + PA_prot_long$ci
PA_prot_long$lower_ci <- PA_prot_long$beta - PA_prot_long$ci
PA_prot_long$Protein <- gsub("^rnt_", "", PA_prot_long$Protein)

# Create the forest plot
pdf(paste0(results_dir, "figures/forest_plot_PA_protein_bmi_consistent.pdf"), height = 12, width = 6)
ggplot(PA_prot_long, aes(x = beta, y = Protein, color = measure)) +
  geom_point(size = 3) +  # Points for beta coefficients
  geom_errorbarh(aes(xmin = lower_ci, xmax = upper_ci), height = 0.2) +  # Error bars for CI
  geom_vline(xintercept = 0, linetype = "dotted", color = "black") +  # Dotted vertical line at x = 0
  facet_wrap(~ measure, scales = "free_x") +  # Facet by measure type
  labs(x = "Difference in protein (SDs) ± 95% CI per unit higher physical activity", y = "Protein",
       title = "Forest Plot of Physical Activity Association with Olink Proteins") +
  theme_minimal() +
  theme(legend.position = "bottom")
dev.off()

