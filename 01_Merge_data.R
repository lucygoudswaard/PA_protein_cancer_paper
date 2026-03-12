### Physical activity and protein analysis ###
### Exploring UKB RAP data and merging PA and protein data ###
### Lucy Goudswaard - 22nd July 2024 ###

## load packages
library(data.table)
library(psych)
library(dplyr)
library(tidyr)
library(ggplot2)
library(reshape2)
library(moosefun)
library(metaboprep)

## connect to project and setwd to script folder

# read in parameter file (specified on command line)
source("parameter_file/parameters_for_r.R")

## read in Olink data. Instance 0 is the only file with 3k proteins.
olink0 <- fread(paste0(input_dir, "olink_instance_0.csv")) ## 53016 participants
olink2 <- read.csv(file = paste0(input_dir, "olink_instance_2.csv"), sep = ",", na.strings = " ") ## 1172 participants
olink3 <- fread(paste0(input_dir, "olink_instance_3.csv"))

## read in physical activity file, 90012 – overall acceleration average, 40047 – sedentary overall average, 40049 – moderate to vigorous overall average,
PA_data_prep <- fread(paste0(input_dir, "MET_PA_covars_all_participant_technical_vars_participant.csv"))

## Drop variables no longer needed (MET minutes is duplicated in another dataframe)
PA_data_prep <- PA_data_prep[,-c("Summed MET minutes per week for all activity | Instance 0")]

## Read in additional MET instances
MET <- fread(paste0(input_dir, "MET_instances_with_ID.csv"))

## Merge
PA_data_prep2 <- merge(PA_data_prep, MET, by = "Participant ID", all = T)

## Merge in genetic PCs and linker
gen_PCs <-  fread(paste0(input_dir, "data.pca1-10.plink.txt"), header = F)
colnames(gen_PCs) <- c("ieu", "ieu2", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10")
linker <-  fread(paste0(input_dir, "16391_from_Becky/linker.csv"))
PC_linker <- merge(linker, gen_PCs, by = "ieu")

## Check withdrawals are gone
withdrawals <- fread(paste0(input_dir, "/16391_from_Becky/withdrawals/w16391_20241217.csv"))
colnames(withdrawals) <- "withdrawn"

## Merge PCs to PA_data
PA_data <- merge(PA_data_prep2, PC_linker, by.x = "Participant ID", by.y = "app", all.x = T)

## No withdrawals need to be made - these aren'r in the data
withdraw <- which(PA_data$`Participant ID` %in% withdrawals$withdrawn)
PA_data <- PA_data[-withdraw,]

## Rename PA columns
names(PA_data)[names(PA_data) == 'Overall acceleration average'] <- 'overall_acceleration_average'
names(PA_data)[names(PA_data) == 'Sedentary - Overall average | Instance 0'] <- 'sedentary_overall_average'
names(PA_data)[names(PA_data) == 'Moderate-Vigorous - Overall average | Instance 0'] <- 'moderate_to_vigorous_overall_average'
names(PA_data)[names(PA_data) == 'Summed MET minutes per week for all activity | Instance 0'] <- 'MET_mins_week_overall'
names(PA_data)[names(PA_data) == 'Summed MET minutes per week for all activity | Instance 1'] <- 'MET_mins_week_overall_1'
names(PA_data)[names(PA_data) == 'Summed MET minutes per week for all activity | Instance 2'] <- 'MET_mins_week_overall_2'
names(PA_data)[names(PA_data) == 'Summed MET minutes per week for all activity | Instance 3'] <- 'MET_mins_week_overall_3'

## Change the overall acceleration average to NA if the calibration and wear time were poor
PA_data$overall_acceleration_average <- ifelse(
  is.na(PA_data$`Data quality, good calibration`) | 
    is.na(PA_data$`Data quality, good wear time`) | 
    PA_data$`Data quality, good calibration` == "No" | 
    PA_data$`Data quality, good wear time` == "No",
  NA,
  PA_data$overall_acceleration_average
)

PA_data$sedentary_overall_average <- ifelse(
  is.na(PA_data$`Data quality, good calibration`) | 
    is.na(PA_data$`Data quality, good wear time`) | 
    PA_data$`Data quality, good calibration` == "No" | 
    PA_data$`Data quality, good wear time` == "No",
  NA,
  PA_data$sedentary_overall_average
)

PA_data$moderate_to_vigorous_overall_average <- ifelse(
  is.na(PA_data$`Data quality, good calibration`) | 
    is.na(PA_data$`Data quality, good wear time`) | 
    PA_data$`Data quality, good calibration` == "No" | 
    PA_data$`Data quality, good wear time` == "No",
  NA,
  PA_data$moderate_to_vigorous_overall_average
)

## List of PA column names 
PA_names <- grep("overall", names(PA_data), value = TRUE)

## Histograms of all physical activity data with summaries
#
pdf(file = paste0(results_dir, "figures/PA_raw_histogram.pdf"), height = 6, width = 8)
for (i in 1:length(PA_names)) {
  column_name <- PA_names[i]
  raw_data <- unlist(PA_data[, ..column_name])
  
  # Calculate mean and standard deviation
  meanvar <- mean(raw_data, na.rm = TRUE)
  sdvar <- sd(raw_data, na.rm = TRUE)
  
  # Exclude outliers based on 10*SD from the mean
  lower_limit <- meanvar - 10 * sdvar
  upper_limit <- meanvar + 10 * sdvar
  filtered_data <- raw_data[raw_data >= lower_limit & raw_data <= upper_limit]
  
  # Recalculate statistics for filtered data
  DataDescribed <- psych::describe(filtered_data)
  meanvar <- DataDescribed$mean
  medianvar <- DataDescribed$median
  minvar <- DataDescribed$min
  maxvar <- DataDescribed$max
  kurtosisvar <- DataDescribed$kurtosis
  skewnessvar <- DataDescribed$skew
  N <- length(filtered_data)
  missingness <- (sum(is.na(raw_data)) / length(raw_data)) * 100
  
  cols <- PA_names[i]
  hist_data <- hist(filtered_data, 
                    col = "red", 
                    main = cols, 
                    xlab = "peak area", 
                    breaks = 50) # Frequency histogram with 50 bins
  
  # Get plotting limits
  x_limits <- par("usr")[1:2]  # x_min and x_max
  y_limits <- par("usr")[3:4]  # y_min and y_max
  
  # Adjust text placement
  text_x <- x_limits[1] + (x_limits[2] - x_limits[1]) * 0.75  # Move to the right
  text_y <- y_limits[1] + (y_limits[2] - y_limits[1]) * 0.9   # Place near the top
  
  text(x = text_x, 
       y = text_y, 
       cex = 0.6, 
       paste("N=", N, "\npercent missing=", 
             signif(missingness, 3), "\nmin=", 
             signif(minvar, 3), " \nmax=",
             signif(maxvar, 3), 
             "\nmean=", 
             signif(meanvar, 3), " \nmedian=", signif(medianvar, 3), 
             "\nkurt=", 
             signif(kurtosisvar, 3), " \nskew=", 
             signif(skewnessvar, 3), sep = ''), 
       pos = 4, xpd = NA)
}
# Close the pdf file
dev.off()


## which Olink IDs occur more than once across the three dataframes
all_eids <- bind_rows(
  olink0 %>% dplyr::select(eid),
  olink2 %>% dplyr::select(eid),
  olink3 %>% dplyr::select(eid)
)

# Count occurrences of each eid
eid_counts <- all_eids %>%
  group_by(eid) %>%
  summarise(count = n())

# Filter eids that occur exactly twice or three times
eids_twice <- eid_counts %>% filter(count == 2)
eids_thrice <- eid_counts %>% filter(count == 3)

# Return the results as tables
ids_repeated <- list(
  eids_twice = eids_twice,
  eids_thrice = eids_thrice
)

# Identify common eids across the three dataframes
common_eids <- Reduce(intersect, list(olink0$eid, olink2$eid, olink3$eid))

# Subset the dataframes to only include common participants
olink0_common <- olink0 %>% filter(eid %in% common_eids)
olink2_common <- olink2 %>% filter(eid %in% common_eids)
olink3_common <- olink3 %>% filter(eid %in% common_eids)

# Ensure columns match and exclude the 'eid' column for correlation analysis
common_columns <- intersect(intersect(names(olink0_common), names(olink2_common)), names(olink3_common))
common_columns <- setdiff(common_columns, "eid")

# Add suffixes to all columns except for 'eid'
olink0_common <- olink0_common %>% rename_with(~ paste0(.x, "_olink0"), -eid)
olink2_common <- olink2_common %>% rename_with(~ paste0(.x, "_olink2"), -eid)
olink3_common <- olink3_common %>% rename_with(~ paste0(.x, "_olink3"), -eid)

# Combine the data frames by 'eid' for correlation calculation
combined_data <- full_join(olink0_common, olink2_common, by = "eid")
combined_data <- full_join(combined_data, olink3_common, by = "eid")

# Extract protein names (base names before _olink suffix)
protein_columns <- names(combined_data)[-1]  # Exclude 'eid'
protein_names <- unique(sub("_olink[0-9]", "", protein_columns))

# Count occurrences of each protein across dataframes
protein_counts <- sapply(protein_names, function(protein) {
  col_olink0 <- paste0(protein, "_olink0")
  col_olink2 <- paste0(protein, "_olink2")
  col_olink3 <- paste0(protein, "_olink3")
  
  # Check if the protein exists in at least two dataframes
  sum(c(col_olink0 %in% names(combined_data),
        col_olink2 %in% names(combined_data),
        col_olink3 %in% names(combined_data)))
})

# Filter proteins that are present in at least two dataframes
proteins_in_two_or_more <- protein_names[protein_counts >= 2]

# Initialize a list to store correlation matrices
correlation_matrices <- list()

# Calculate the correlation matrix for each protein and store it
for (protein in proteins_in_two_or_more) {
  col_olink0 <- paste0(protein, "_olink0")
  col_olink2 <- paste0(protein, "_olink2")
  col_olink3 <- paste0(protein, "_olink3")
  
  data_to_correlate <- combined_data %>%
    dplyr::select(eid, all_of(col_olink0), all_of(col_olink2), all_of(col_olink3)) %>%
    na.omit()
  
  if (ncol(data_to_correlate) == 4) {
    cor_matrix <- cor(data_to_correlate[, -1])  # Exclude 'eid' for correlation calculation
    correlation_matrices[[protein]] <- cor_matrix
  } else {
    correlation_matrices[[protein]] <- matrix(NA, 3, 3, dimnames = list(c("olink0", "olink2", "olink3"), c("olink0", "olink2", "olink3")))
  }
}

# Plot the correlation matrices with text annotations
#pdf(file = paste0(results_dir, "figures/olink_instance_corr_plot.pdf"), height = 6, width = 6)
for (protein in proteins_in_two_or_more) {
  cor_matrix <- correlation_matrices[[protein]]
  
  # Melt the correlation matrix for ggplot2
  cor_matrix_melted <- melt(cor_matrix)
  
  # Create the heatmap plot with annotations
  p <- ggplot(cor_matrix_melted, aes(x = Var1, y = Var2, fill = value)) +
    geom_tile() +
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, limit = c(-1, 1)) +
    geom_text(aes(label = sprintf("%.2f", value)), color = "black", size = 4) +
    labs(title = paste("Correlation Heatmap for", protein),
         x = "Dataset",
         y = "Dataset",
         fill = "Correlation") +
    theme_minimal()
  
  # Print the plot
  print(p)
}
dev.off()

## Olink proteins - how many independent proteins are there
## Run feature tree independence on a sample of 5000
olink0_feature_data <- olink0[,-c("eid")]
set.seed(123)
olink0_feature_data_sample <- olink0_feature_data %>% sample_n(5000)
tree_and_independent_features_ukb = tree_and_independent_features(olink0_feature_data_sample, minimum_samplesize = 30, tree_cut_height = 0.5, feature_names_2_exclude = NA)
length(tree_and_independent_features_ukb$independent_features) ## 1885

## save PA_data
write.table(PA_data, file = paste0(data_intermediate_dir, "PA_protein_data_with_MET.csv"), col.names = T, row.names = F, sep = ",")

## merge
olink_PA <- merge(PA_data, olink0, by.x = "Participant ID", by.y = "eid", all = FALSE) ## 42533 participants remain

## save protein data for metaboprep
protein_names <- read.table(paste0(input_dir, "field_names_instance_0.txt"), header = T, sep = "\t")
protein_names <- protein_names[[1]]
olink_PA <- as.data.frame(olink_PA)
names(olink_PA)[names(olink_PA) == 'Participant ID'] <- 'eid'
olink_protein_subset <- olink_PA[,c("eid", protein_names)]
#write.table(olink_protein_subset, file = paste0(data_intermediate_dir, "Olink_protein_data_subset_with_MET.txt"), col.names = T, row.names = F, sep = "\t")

