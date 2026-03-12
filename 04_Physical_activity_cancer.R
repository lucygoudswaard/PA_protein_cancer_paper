#### Physical activity and CRC association ####
## Started 28th October 2024

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
library(gtsummary)
library(stringr)
library(EpiViz)
library(ggrepel)
library(glue)
library(MASS)
library(ggforce)
library(survival)
library(lubridate)

## clear environment
rm(list=ls())

## connect to project and setwd to script folder

# read in parameter file (specified on command line)
source("parameter_file/parameters_for_r.R")

# read in CRC  data
cancer <- fread(file = paste0(input_dir, "cancer_all_ukb.csv"), header = T, sep = ",")
cancer <- as.data.frame(cancer)

# Define CRC status
## CRC case at any timepoint
# colorectal cancer cases. ICD10 = p40006. IC9 = 40013.
# Define the CRC codes to search for
CRC_codes <- c("C18", "C18.0", 'C18.1', 'C18.2', 'C18.3', 'C18.4', 'C18.5', 'C18.6', 'C18.7', 'C18.8', 'C18.9', "C19", "C20")

# Create a pattern that combines all CRC codes for detection
CRC_pattern <- paste(CRC_codes, collapse = "|")  # Create a regex pattern separated by "|"

# Update the CRC data frame by adding a new column
cancer <- cancer %>%
  mutate(CRC_ICD10_case = if_else(
    if_any(starts_with("p40006"), ~ str_detect(.x, CRC_pattern)),  # Check if any column starts with p40006 contains any of the CRC codes
    "case",  # If true, assign "case"
    "control"  # If false, assign "control"
  ))

# Define the CRC ICD-9 codes to search for
CRC_ICD9 <- c("1530", "1531", "1532", "1533", "1534", "1535", "1536", "1537", "1538", "1539", "1540", "1541", "1542", "1543", "1548")

# Create a pattern that combines all CRC ICD-9 codes for detection
CRC_ICD9_pattern <- paste(CRC_ICD9, collapse = "|")  # Create a regex pattern separated by "|"

# Update the CRC data frame by adding a new column for ICD-9 codes
cancer <- cancer %>%
  mutate(across(starts_with("p40013"), as.character)) %>%
  mutate(CRC_ICD9_case = if_else(
    if_any(starts_with("p40013"), ~ replace_na(str_detect(.x, CRC_ICD9_pattern), FALSE)),
    "case",
    "control"
  ))

## Check these new case/control variables

table(cancer$CRC_ICD10_case) ## 9042 cases, 493182 controls
table(cancer$CRC_ICD9_case) ## 291 cases, 501933 controls

## CRC date of diagnosis for ICD10 codes
# Identify all p40006 (ICD10 code) and p40005 (dates) columns
p40006_cols <- grep("^p40006_", names(cancer), value = TRUE)
p40005_cols <- grep("^p40005_", names(cancer), value = TRUE)

# Extract instance identifiers for all p40006 columns
instances <- str_extract(p40006_cols, "i[0-9]+")

# Initialize 'ICD10_date' column with NA values in date format
cancer$CRC_ICD10_date <- as.IDate(NA)

# Loop through each instance (i0, i1, ..., i21)
for (instance in unique(instances)) {
  # Get the specific p40006 and p40005 column for the current instance
  p40006_col <- paste0("p40006_", instance)
  p40005_col <- paste0("p40005_", instance)
  
  # Ensure both columns exist before proceeding
  if (p40006_col %in% names(cancer) & p40005_col %in% names(cancer)) {
    # Update 'ICD10_date' for rows where CRC code is detected in the current p40006 column
    cancer <- cancer %>%
      mutate(CRC_ICD10_date = if_else(
        is.na(CRC_ICD10_date) & CRC_ICD10_case == "case" & str_detect(get(p40006_col), CRC_pattern),
        get(p40005_col),
        CRC_ICD10_date
      ))
  }
}

## Repeat for ICD9
p40013_cols <- grep("^p40013_", names(cancer), value = TRUE)
p40005_cols <- grep("^p40005_", names(cancer), value = TRUE)

# Extract instance identifiers for all p40006 columns
instances2 <- str_extract(p40013_cols, "i[0-9]+")

# Initialize 'ICD10_date' column with NA values in date format
cancer$CRC_ICD9_date <- as.IDate(NA)

# Loop through each instance (i0, i1, ..., i21)
for (instance in unique(instances2)) {
  # Get the specific p40013 and p40005 column for the current instance
  p40013_col <- paste0("p40013_", instance)
  p40005_col <- paste0("p40005_", instance)
  
  # Ensure both columns exist before proceeding
  if (p40013_col %in% names(cancer) & p40005_col %in% names(cancer)) {
    # Update 'ICD10_date' for rows where CRC code is detected in the current p40006 column
    cancer <- cancer %>%
      mutate(CRC_ICD9_date = if_else(
        is.na(CRC_ICD9_date) & CRC_ICD9_case == "case" & str_detect(get(p40013_col), CRC_ICD9_pattern),
        get(p40005_col),
        CRC_ICD9_date
      ))
  }
}

## Which participants are cases according to both systems?
CRC_ICD9_ICD10 <- which(cancer$CRC_ICD10_case == "case" & cancer$CRC_ICD9_case == "case")
rows_consistent_CRC_case <- cancer[CRC_ICD9_ICD10,]

## Participants with diagnosis across either ICD10 code
CRC_ICD9_or_ICD10 <- which(cancer$CRC_ICD10_case == "case" | cancer$CRC_ICD9_case == "case")
rows_CRC_case_either <- cancer[CRC_ICD9_or_ICD10,]

## Exclude participants if self reported CRC but no ICD10 or ICD9 code - p20001 is cancer code, self reported, verbal interview
selected_cols <- grep("^p20001_", names(cancer), value = TRUE)

# Create the column for self-reported information
cancer$case_self_report_info <- NA_character_

# Identify cases
cases <- cancer$CRC_ICD10_case == "case" | cancer$CRC_ICD9_case == "case"

# Efficiently combine strings for cases
combined_strings <- apply(cancer[c(selected_cols)], 1, function(x) {
  paste(na.omit(x), collapse = "; ")
})

# Assign combined strings only to the case_self_report_info for identified cases
cancer$case_self_report_info <- combined_strings

# Define CRC_terms
CRC_terms <- c("colorectal", "caecum", "colon", "ileocaecal valve", "Appendix", 
               "Ascending colon", "Hepatic flexure", "Transverse colon", 
               "Splenic flexure", "Descending", "Sigmoid", 
               "Overlapping lesion of colon", "rectosigmoid", "rectum", "rectal")

# Create a regex pattern to match any of the terms in CRC_terms
CRC_pattern <- paste(CRC_terms, collapse = "|")

# Update the remove column based on the new criteria
cancer <- cancer %>%
  mutate(CRC_true_control = case_when(
    !cases & !str_detect(case_self_report_info, CRC_pattern) ~ "true_control",
    !cases & str_detect(case_self_report_info, CRC_pattern) ~ "control_mismatch",
    cases ~ "case"  # Assign "case" if they are a case
  ))

# Check the distribution of remove values
table(cancer$CRC_true_control)  # Check how many are "true_control", "control_mismatch", and "case"

# Explore mismatches with ICD10
mismatch_rows <- which(cancer$true_control == "control_mismatch")
mismatches <- cancer[mismatch_rows, ]

## Change mismatches to NA
cancer$CRC_true_control <- ifelse(cancer$CRC_true_control %in% c("case", "true_control"), 
                              cancer$CRC_true_control, 
                           NA)
                           
## Derive CRC variable for diagnosis after 2012-08-01. Changed to 2015-12-31 for sensitivity analysis.
cutoff_date <- as.IDate("2012-08-01")

# Create the CRC_incidental column based on the conditions
cancer <- cancer %>%
  mutate(CRC_incidental = case_when(
    # Check if any date is present and after the cutoff
    (CRC_ICD9_date >= cutoff_date | CRC_ICD10_date >= cutoff_date) ~ "incidental",
    
    # Check if any date is present and before the cutoff
    (CRC_ICD9_date < cutoff_date | CRC_ICD10_date < cutoff_date) ~ "prevalent",
    
    # If both dates are missing, classify as "control"
    is.na(CRC_ICD9_date) & is.na(CRC_ICD10_date) ~ "control"
  ))
table(cancer$CRC_incidental)

# Change mismatch to NA
cancer$CRC_incidental[mismatch_rows] <- NA

## Change prevalent to NA
w <- which(cancer$CRC_incidental == "prevalent")
cancer$CRC_incidental[w] <- NA
table(cancer$CRC_incidental)

## Repeat the above for breast cancer
## Breast cancer codes
BC_codes <- c("C50", "C50.0", "C50.1", "C50.2", "C50.3", "C50.4", "C50.5", "C50.6", "C50.7", "C50.8", "C50.9")

# Create a pattern that combines all CRC codes for detection
BC_pattern <- paste(BC_codes, collapse = "|")  # Create a regex pattern separated by "|"

# Update the CRC data frame by adding a new column
cancer <- cancer %>%
  mutate(BC_ICD10_case = if_else(
    if_any(starts_with("p40006"), ~ str_detect(.x, BC_pattern)),  # Check if any column starts with p40006 contains any of the CRC codes
    "case",  # If true, assign "case"
    "control"  # If false, assign "control"
  ))

# Define the CRC ICD-9 codes to search for - do I include 175 for male breast cancer?
BC_ICD9 <- c("174", "174.0", "174.1", "174.2", "174.3", "174.4", "174.5", "174.6", "174.7", "174.8", "174.9")

# Create a pattern that combines all CRC ICD-9 codes for detection
BC_ICD9_pattern <- paste(BC_ICD9, collapse = "|")  # Create a regex pattern separated by "|"

# Update the CRC data frame by adding a new column for ICD-9 codes
cancer <- cancer %>%
  mutate(across(starts_with("p40013"), as.character)) %>%
  mutate(BC_ICD9_case = if_else(
    if_any(starts_with("p40013"), ~ replace_na(str_detect(.x, BC_ICD9_pattern), FALSE)),
    "case",
    "control"
  ))

## Check these new case/control variables
table(cancer$BC_ICD10_case) ## 18844 cases, 483380 controls
table(cancer$BC_ICD9_case) ## 1767 cases, 500457 controls

## BC date of diagnosis for ICD10 codes
# Identify all p40006 (ICD10 code) and p40005 (dates) columns
p40006_cols <- grep("^p40006_", names(cancer), value = TRUE)
p40005_cols <- grep("^p40005_", names(cancer), value = TRUE)

# Extract instance identifiers for all p40006 columns
instances <- str_extract(p40006_cols, "i[0-9]+")

# Initialize 'ICD10_date' column with NA values in date format
cancer$BC_ICD10_date <- as.IDate(NA)

# Loop through each instance (i0, i1, ..., i21)
for (instance in unique(instances)) {
  # Get the specific p40006 and p40005 column for the current instance
  p40006_col <- paste0("p40006_", instance)
  p40005_col <- paste0("p40005_", instance)
  
  # Ensure both columns exist before proceeding
  if (p40006_col %in% names(cancer) & p40005_col %in% names(cancer)) {
    # Update 'BC_ICD10_date' for rows where BC code is detected in the current p40006 column
    cancer <- cancer %>%
      mutate(BC_ICD10_date = if_else(
        is.na(BC_ICD10_date) & BC_ICD10_case == "case" & str_detect(get(p40006_col), BC_pattern),
        get(p40005_col),
        BC_ICD10_date
      ))
  }
}

## Repeat for ICD9
p40013_cols <- grep("^p40013_", names(cancer), value = TRUE)
p40005_cols <- grep("^p40005_", names(cancer), value = TRUE)

# Extract instance identifiers for all p40006 columns
instances2 <- str_extract(p40013_cols, "i[0-9]+")

# Initialize 'ICD10_date' column with NA values in date format
cancer$BC_ICD9_date <- as.IDate(NA)

# Loop through each instance (i0, i1, ..., i21)
for (instance in unique(instances2)) {
  # Get the specific p40013 and p40005 column for the current instance
  p40013_col <- paste0("p40013_", instance)
  p40005_col <- paste0("p40005_", instance)
  
  # Ensure both columns exist before proceeding
  if (p40013_col %in% names(cancer) & p40005_col %in% names(cancer)) {
    # Update 'ICD10_date' for rows where CRC code is detected in the current p40006 column
    cancer <- cancer %>%
      mutate(BC_ICD9_date = if_else(
        is.na(BC_ICD9_date) & BC_ICD9_case == "case" & str_detect(get(p40013_col), BC_ICD9_pattern),
        get(p40005_col),
        BC_ICD9_date
      ))
  }
}

## Which participants are cases according to both systems?
BC_ICD9_ICD10 <- which(cancer$BC_ICD10_case == "case" & cancer$BC_ICD9_case == "case")
rows_consistent_BC_case <- cancer[BC_ICD9_ICD10,]

## Participants with diagnosis across either ICD10 code
BC_ICD9_or_ICD10 <- which(cancer$BC_ICD10_case == "case" | cancer$BC_ICD9_case == "case")
rows_BC_case_either <- cancer[BC_ICD9_or_ICD10,]

## Exclude participants if self reported BC but no ICD10 or ICD9 code - p20001 is cancer code, self reported, verbal interview
selected_cols <- grep("^p20001_", names(cancer), value = TRUE)

# Create the column for self-reported information
cancer$case_self_report_info <- NA_character_

# Identify cases
cases <- cancer$BC_ICD10_case == "case" | cancer$BC_ICD9_case == "case"

# Efficiently combine strings for cases
combined_strings <- apply(cancer[c(selected_cols)], 1, function(x) {
  paste(na.omit(x), collapse = "; ")
})

# Assign combined strings only to the case_self_report_info for identified cases
cancer$case_self_report_info <- combined_strings

# Define BC_terms
BC_terms <- c("breast")

# Create a regex pattern to match any of the terms in CRC_terms
BC_pattern <- paste(BC_terms, collapse = "|")

# Update the remove column based on the new criteria
cancer <- cancer %>%
  mutate(BC_true_control = case_when(
    !cases & !str_detect(case_self_report_info, BC_pattern) ~ "true_control",
    !cases & str_detect(case_self_report_info, BC_pattern) ~ "control_mismatch",
    cases ~ "case"  # Assign "case" if they are a case
  ))

# Check the distribution of remove values
table(cancer$BC_true_control)  # Check how many are "true_control", "control_mismatch", and "case"

# Explore mismatches with ICD10
mismatch_rows <- which(cancer$BC_true_control == "control_mismatch")
mismatches <- cancer[mismatch_rows, ]

## Change mismatches to NA
cancer$BC_true_control <- ifelse(cancer$BC_true_control %in% c("case", "true_control"), 
                              cancer$BC_true_control, 
                              NA)

## Derive BC variable for diagnosis after date
cutoff_date <- as.IDate("2012-08-01")

# Create the CRC_incidental column based on the conditions
cancer <- cancer %>%
  mutate(BC_incidental = case_when(
    # Check if any date is present and after the cutoff
    (BC_ICD9_date >= cutoff_date | BC_ICD10_date >= cutoff_date) ~ "incidental",
    
    # Check if any date is present and before the cutoff
    (BC_ICD9_date < cutoff_date | BC_ICD10_date < cutoff_date) ~ "prevalent",
    
    # If both dates are missing, classify as "control"
    is.na(BC_ICD9_date) & is.na(BC_ICD10_date) ~ "control"
  ))
table(cancer$BC_incidental)

# Change mismatch to NA
cancer$BC_incidental[mismatch_rows] <- NA

## Change prevalent to NA
w <- which(cancer$BC_incidental == "prevalent")
cancer$BC_incidental[w] <- NA
table(cancer$BC_incidental)

## Repeat for endometrial cancer
## Endometrial cancer codes
EC_ICD10 <- c("C54", "C54.1", "C54.2", "C54.3", "C54.9", "C55")
EC_ICD9 <- c("179", "1799", "180", "182", "1820", "1821", "1828")

# Create patterns for regex detection
EC_ICD10_pattern <- paste(EC_ICD10, collapse = "|")
EC_ICD9_pattern <- paste(EC_ICD9, collapse = "|")

# Update the cancer data frame by adding a new column for ICD10 codes
cancer <- cancer %>%
  mutate(EC_ICD10_case = if_else(
    if_any(starts_with("p40006"), ~ str_detect(.x, EC_ICD10_pattern)),
    "case",
    "control"
  ))

# Update the cancer data frame by adding a new column for ICD9 codes
cancer <- cancer %>%
  mutate(across(starts_with("p40013"), as.character)) %>%
  mutate(EC_ICD9_case = if_else(
    if_any(starts_with("p40013"), ~ replace_na(str_detect(.x, EC_ICD9_pattern), FALSE)),
    "case",
    "control"
  ))

## Check these new case/control variables
table(cancer$EC_ICD10_case)
table(cancer$EC_ICD9_case)

## EC date of diagnosis for ICD10 codes
p40006_cols <- grep("^p40006_", names(cancer), value = TRUE)
p40005_cols <- grep("^p40005_", names(cancer), value = TRUE)
instances <- str_extract(p40006_cols, "i[0-9]+")

# Initialize 'EC_ICD10_date' column with NA values in date format
cancer$EC_ICD10_date <- as.IDate(NA)

for (instance in unique(instances)) {
  p40006_col <- paste0("p40006_", instance)
  p40005_col <- paste0("p40005_", instance)
  if (p40006_col %in% names(cancer) & p40005_col %in% names(cancer)) {
    cancer <- cancer %>%
      mutate(EC_ICD10_date = if_else(
        is.na(EC_ICD10_date) & EC_ICD10_case == "case" & str_detect(get(p40006_col), EC_ICD10_pattern),
        get(p40005_col),
        EC_ICD10_date
      ))
  }
}

## Repeat for ICD9
p40013_cols <- grep("^p40013_", names(cancer), value = TRUE)
instances2 <- str_extract(p40013_cols, "i[0-9]+")

cancer$EC_ICD9_date <- as.IDate(NA)

for (instance in unique(instances2)) {
  p40013_col <- paste0("p40013_", instance)
  p40005_col <- paste0("p40005_", instance)
  if (p40013_col %in% names(cancer) & p40005_col %in% names(cancer)) {
    cancer <- cancer %>%
      mutate(EC_ICD9_date = if_else(
        is.na(EC_ICD9_date) & EC_ICD9_case == "case" & str_detect(get(p40013_col), EC_ICD9_pattern),
        get(p40005_col),
        EC_ICD9_date
      ))
  }
}

## Participants with diagnosis across either ICD10 or ICD9 codes
EC_ICD9_or_ICD10 <- which(cancer$EC_ICD10_case == "case" | cancer$EC_ICD9_case == "case")
rows_EC_case_either <- cancer[EC_ICD9_or_ICD10, ]

## Exclude participants if self-reported EC but no ICD10 or ICD9 code
selected_cols <- grep("^p20001_", names(cancer), value = TRUE)

cancer$case_self_report_info <- NA_character_
cases <- cancer$EC_ICD10_case == "case" | cancer$EC_ICD9_case == "case"

combined_strings <- apply(cancer[c(selected_cols)], 1, function(x) {
  paste(na.omit(x), collapse = "; ")
})

cancer$case_self_report_info <- combined_strings

EC_terms <- c("endometrial")
EC_pattern <- paste(EC_terms, collapse = "|")

cancer <- cancer %>%
  mutate(EC_true_control = case_when(
    !cases & !str_detect(case_self_report_info, EC_pattern) ~ "true_control",
    !cases & str_detect(case_self_report_info, EC_pattern) ~ "control_mismatch",
    cases ~ "case"
  ))

table(cancer$EC_true_control)

mismatch_rows <- which(cancer$EC_true_control == "control_mismatch")
mismatches <- cancer[mismatch_rows, ]

cancer$EC_true_control <- ifelse(cancer$EC_true_control %in% c("case", "true_control"), 
                              cancer$EC_true_control, 
                              NA)

## Derive EC variable for diagnosis after a cutoff date
cutoff_date <- as.IDate("2012-08-01")

cancer <- cancer %>%
  mutate(EC_incidental = case_when(
    (EC_ICD9_date >= cutoff_date | EC_ICD10_date >= cutoff_date) ~ "incidental",
    (EC_ICD9_date < cutoff_date | EC_ICD10_date < cutoff_date) ~ "prevalent",
    is.na(EC_ICD9_date) & is.na(EC_ICD10_date) ~ "control"
  ))

table(cancer$EC_incidental)

cancer$EC_incidental[mismatch_rows] <- NA
w <- which(cancer$EC_incidental == "prevalent")
cancer$EC_incidental[w] <- NA
table(cancer$EC_incidental)

## Repeat the above for prostate cancer
## Prostate cancer codes
PC_ICD10 <- c("C61")
PC_ICD9 <- c("185")

# Create patterns for regex detection
PC_ICD10_pattern <- paste(PC_ICD10, collapse = "|")
PC_ICD9_pattern <- paste(PC_ICD9, collapse = "|")

# Update the cancer data frame by adding a new column for ICD10 codes
cancer <- cancer %>%
  mutate(PC_ICD10_case = if_else(
    if_any(starts_with("p40006"), ~ str_detect(.x, PC_ICD10_pattern)),
    "case",
    "control"
  ))

# Update the cancer data frame by adding a new column for ICD9 codes
cancer <- cancer %>%
  mutate(across(starts_with("p40013"), as.character)) %>%
  mutate(PC_ICD9_case = if_else(
    if_any(starts_with("p40013"), ~ replace_na(str_detect(.x, PC_ICD9_pattern), FALSE)),
    "case",
    "control"
  ))

## Check these new case/control variables
table(cancer$PC_ICD10_case)
table(cancer$PC_ICD9_case)

## PC date of diagnosis for ICD10 codes
p40006_cols <- grep("^p40006_", names(cancer), value = TRUE)
p40005_cols <- grep("^p40005_", names(cancer), value = TRUE)
instances <- str_extract(p40006_cols, "i[0-9]+")

# Initialize 'PC_ICD10_date' column with NA values in date format
cancer$PC_ICD10_date <- as.IDate(NA)

for (instance in unique(instances)) {
  p40006_col <- paste0("p40006_", instance)
  p40005_col <- paste0("p40005_", instance)
  if (p40006_col %in% names(cancer) & p40005_col %in% names(cancer)) {
    cancer <- cancer %>%
      mutate(PC_ICD10_date = if_else(
        is.na(PC_ICD10_date) & PC_ICD10_case == "case" & str_detect(get(p40006_col), PC_ICD10_pattern),
        get(p40005_col),
        PC_ICD10_date
      ))
  }
}

## Repeat for ICD9
p40013_cols <- grep("^p40013_", names(cancer), value = TRUE)
instances2 <- str_extract(p40013_cols, "i[0-9]+")

cancer$PC_ICD9_date <- as.IDate(NA)

for (instance in unique(instances2)) {
  p40013_col <- paste0("p40013_", instance)
  p40005_col <- paste0("p40005_", instance)
  if (p40013_col %in% names(cancer) & p40005_col %in% names(cancer)) {
    cancer <- cancer %>%
      mutate(PC_ICD9_date = if_else(
        is.na(PC_ICD9_date) & PC_ICD9_case == "case" & str_detect(get(p40013_col), PC_ICD9_pattern),
        get(p40005_col),
        PC_ICD9_date
      ))
  }
}

## Participants with diagnosis across either ICD10 or ICD9 codes
PC_ICD9_or_ICD10 <- which(cancer$PC_ICD10_case == "case" | cancer$PC_ICD9_case == "case")
rows_PC_case_either <- cancer[PC_ICD9_or_ICD10, ]

## Exclude participants if self-reported PC but no ICD10 or ICD9 code
selected_cols <- grep("^p20001_", names(cancer), value = TRUE)

cancer$case_self_report_info <- NA_character_
cases <- cancer$PC_ICD10_case == "case" | cancer$PC_ICD9_case == "case"

combined_strings <- apply(cancer[c(selected_cols)], 1, function(x) {
  paste(na.omit(x), collapse = "; ")
})

cancer$case_self_report_info <- combined_strings

PC_terms <- c("prostate")
PC_pattern <- paste(PC_terms, collapse = "|")

cancer <- cancer %>%
  mutate(PC_true_control = case_when(
    !cases & !str_detect(case_self_report_info, PC_pattern) ~ "true_control",
    !cases & str_detect(case_self_report_info, PC_pattern) ~ "control_mismatch",
    cases ~ "case"
  ))

table(cancer$PC_true_control)

mismatch_rows <- which(cancer$PC_true_control == "control_mismatch")
mismatches <- cancer[mismatch_rows, ]

cancer$PC_true_control <- ifelse(cancer$PC_true_control %in% c("case", "true_control"), 
                              cancer$PC_true_control, 
                              NA)

## Derive PC variable for diagnosis after a cutoff date
cutoff_date <- as.IDate("2012-08-01")

cancer <- cancer %>%
  mutate(PC_incidental = case_when(
    (PC_ICD9_date >= cutoff_date | PC_ICD10_date >= cutoff_date) ~ "incidental",
    (PC_ICD9_date < cutoff_date | PC_ICD10_date < cutoff_date) ~ "prevalent",
    is.na(PC_ICD9_date) & is.na(PC_ICD10_date) ~ "control"
  ))

table(cancer$PC_incidental)

cancer$PC_incidental[mismatch_rows] <- NA
w <- which(cancer$PC_incidental == "prevalent")
cancer$PC_incidental[w] <- NA
table(cancer$PC_incidental)

# read in PA data. copied the .txt file as a .csv to try to make reading in easier.
PA_covar_data <- fread(paste0(data_intermediate_dir, "PA_protein_data_with_MET.csv"), header = T, sep = ",")
PA_protein <- fread(file = paste0(data_intermediate_dir, "PA_rnt_protein_data_with_MET2.txt"), header = T, sep = "\t")
PA_protein <- as.data.frame(PA_protein)
protein_names <- read.table(paste0(input_dir, "field_names_instance_0.txt"), header = T, sep = "\t")
protein_names <- protein_names[[1]]
rnt_protein_names <- paste0("rnt_", protein_names)
protein <- PA_protein[, c("Participant ID", rnt_protein_names), drop = FALSE]

## Merge CRC and PA and protein dataset
cancer_PA <- merge(cancer, PA_covar_data, by.x = "eid", by.y = "Participant ID")

## Remove outliers more than 5 SD from PA vars again
PA_cols_QC <- c("overall_acceleration_average", "sedentary_overall_average", 
                "moderate_to_vigorous_overall_average")

# Change outliers to NA in PA data (more than 5 SD away from the mean)
means <- colMeans(cancer_PA[PA_cols_QC], na.rm = TRUE)
sds <- apply(cancer_PA[PA_cols_QC], 2, sd, na.rm = TRUE)
cancer_PA[PA_cols_QC] <- apply(cancer_PA[PA_cols_QC], 2, function(column) {
  outliers <- abs(column - mean(column, na.rm = TRUE)) > (5 * sd(column, na.rm = TRUE))
  column[outliers] <- NA
  return(column)
})

## Derive quantile var again
cancer_PA <- cancer_PA %>%
  mutate(moderate_to_vigorous_quantile = ntile(moderate_to_vigorous_overall_average, 4) )

## Rank normal transform the three original vars 
PA_vars <- c("overall_acceleration_average", "sedentary_overall_average", "moderate_to_vigorous_overall_average")
for (i in 1:length(PA_vars)) {
  cancer_PA[, PA_vars[i]] <- rntransform(cancer_PA[, PA_vars[i]])  
}

## Rename covars again
names(cancer_PA)[names(cancer_PA) == "Age at recruitment"] <- "age"
names(cancer_PA)[names(cancer_PA) == "Sex"] <- "sex"
names(cancer_PA)[names(cancer_PA) == "Alcohol drinker status | Instance 0"] <- "alcohol"
cancer_PA$alcohol <- ifelse(cancer_PA$alcohol == "Current", 2,
                                ifelse(cancer_PA$alcohol == "Previous", 1,
                                       ifelse(cancer_PA$alcohol == "Never", 0, NA)))
names(cancer_PA)[names(cancer_PA) == "Ever smoked | Instance 0"] <- "smoking"
cancer_PA$smoking <- ifelse(cancer_PA$smoking == "Yes", 1,
                                ifelse(cancer_PA$smoking == "No", 0, NA))
names(cancer_PA)[names(cancer_PA) == "Townsend deprivation index at recruitment"] <- "townsend"
names(cancer_PA)[names(cancer_PA) == "Body mass index (BMI) | Instance 0"] <- "bmi"
names(cancer_PA)[names(cancer_PA) == "UK Biobank assessment centre | Instance 0"] <- "centre"
names(cancer_PA)[names(cancer_PA) == "Fasting time | Instance 0"] <- "fasting"
names(cancer_PA)[names(cancer_PA) == "Qualifications | Instance 0"] <- "education_raw"

# Recode education
education_levels <- c(
  "College or University degree" = 5,
  "Other professional qualifications eg: nursing, teaching" = 4,
  "NVQ or HND or HNC or equivalent" = 3,
  "A levels/AS levels or equivalent" = 2,
  "O levels/GCSEs or equivalent" = 1,
  "CSEs or equivalent" = 1,
  "None of the above" = 0,
  "Prefer not to answer" = NA
)

cancer_PA$education <- sapply(cancer_PA$education_raw, function(x) {
  # Split the entry into individual qualifications
  qualifications <- unlist(strsplit(x, "\\|"))
  
  # Find the numeric levels for the qualifications
  qualification_levels <- as.numeric(education_levels[qualifications])
  
  # Check if there are valid levels
  if (all(is.na(qualification_levels))) {
    return(NA) # Return NA if no valid qualifications are found
  } else {
    return(max(qualification_levels, na.rm = TRUE)) # Return the highest valid level
  }
})
cancer_PA$education <- as.factor(cancer_PA$education)

## # List of cancers to analyze
cancer_types <- c("CRC", "EC", "BC", "PC")

# Loop through each cancer type (using 'i' as the loop variable)
for (i in cancer_types) {
  # Create ICD10 case columns for each cancer
  cancer_PA[[paste0(i, "_ICD10_case")]] <- ifelse(cancer_PA[[paste0(i, "_ICD10_case")]] == "case", 1,
                                                  ifelse(cancer_PA[[paste0(i, "_ICD10_case")]] == "control", 0, NA))
  
  # Create status columns based on true_control for each cancer
  cancer_PA[[paste0(i, "_status")]] <- ifelse(cancer_PA[[paste0(i, "_true_control")]] == "case", 1,
                                              ifelse(cancer_PA[[paste0(i, "_true_control")]] == "true_control", 0, NA))
  
  # Create incidental columns for each cancer
  cancer_PA[[paste0(i, "_incidental")]] <- ifelse(cancer_PA[[paste0(i, "_incidental")]] == "incidental", 1,
                                                  ifelse(cancer_PA[[paste0(i, "_incidental")]] == "control", 0, NA))
  
  # Convert status columns to factors (case, control)
  cancer_PA[[paste0(i, "_status")]] <- factor(cancer_PA[[paste0(i, "_status")]], levels = c(0, 1), labels = c("control", "case"))
  cancer_PA[[paste0(i, "_incidental")]] <- factor(cancer_PA[[paste0(i, "_incidental")]], levels = c(0, 1), labels = c("control", "case"))
  
}

## Make a new case/control variable for all 4 cancer types (so someone is only a control if they don't have that cancer,
## or any of the other 3)
cancer_PA$combined_cancer_status <- ifelse(cancer_PA$CRC_status == "case" | cancer_PA$EC_status == "case" |
                                      cancer_PA$BC_status == "case" | cancer_PA$PC_status == "case", "case", "control" )
cancer_PA$combined_cancer_incidental <- ifelse(cancer_PA$CRC_incidental == "case" | cancer_PA$EC_incidental == "case" |
                                                 cancer_PA$BC_incidental == "case" | cancer_PA$PC_incidental == "case", "case", "control" )
cancer_PA$sex <- as.factor(cancer_PA$sex)

cancer_PA$CRC_status_combined_control <- factor(
  ifelse(cancer_PA$CRC_status == "case", 1,
         ifelse(cancer_PA$CRC_status == "control" & cancer_PA$combined_cancer_status == "control", 0, NA)),
  levels = c(0, 1), labels = c("control", "case")
)

cancer_PA$CRC_incidental_combined_control <- factor(
  ifelse(cancer_PA$CRC_incidental == "case", 1,
         ifelse(cancer_PA$CRC_incidental == "control" & cancer_PA$combined_cancer_incidental == "control", 0, NA)),
  levels = c(0, 1), labels = c("control", "case")
)

cancer_PA$BC_status_combined_control <- factor(
  ifelse(cancer_PA$BC_status == "case", 1,
         ifelse(cancer_PA$BC_status == "control" & cancer_PA$combined_cancer_status == "control", 0, NA)),
  levels = c(0, 1), labels = c("control", "case")
)

cancer_PA$BC_incidental_combined_control <- factor(
  ifelse(cancer_PA$BC_incidental == "case", 1,
         ifelse(cancer_PA$BC_incidental == "control" & cancer_PA$combined_cancer_incidental == "control", 0, NA)),
  levels = c(0, 1), labels = c("control", "case")
)

cancer_PA$EC_status_combined_control <- factor(
  ifelse(cancer_PA$EC_status == "case", 1,
         ifelse(cancer_PA$EC_status == "control" & cancer_PA$combined_cancer_status == "control", 0, NA)),
  levels = c(0, 1), labels = c("control", "case")
)

cancer_PA$EC_incidental_combined_control <- factor(
  ifelse(cancer_PA$EC_incidental == "case", 1,
         ifelse(cancer_PA$EC_incidental == "control" & cancer_PA$combined_cancer_incidental == "control", 0, NA)),
  levels = c(0, 1), labels = c("control", "case")
)

cancer_PA$PC_status_combined_control <- factor(
  ifelse(cancer_PA$PC_status == "case", 1,
         ifelse(cancer_PA$PC_status == "control" & cancer_PA$combined_cancer_status == "control", 0, NA)),
  levels = c(0, 1), labels = c("control", "case")
)

cancer_PA$PC_incidental_combined_control <- factor(
  ifelse(cancer_PA$PC_incidental == "case", 1,
         ifelse(cancer_PA$PC_incidental == "control" & cancer_PA$combined_cancer_incidental == "control", 0, NA)),
  levels = c(0, 1), labels = c("control", "case")
)
         
## Create a complete case dataset
cancer_PA_cc <- cancer_PA %>%
  filter(
    # Ensure all covariates are not missing
    !is.na(age) & 
      !is.na(sex) & 
      !is.na(smoking) & 
      !is.na(alcohol) & 
      !is.na(centre) & 
      !is.na(fasting) & 
      !is.na(bmi) & 
      !is.na(education) &
      # Retain rows that have at least one physical activity variable
      (!is.na(overall_acceleration_average) |
         !is.na(moderate_to_vigorous_quantile) |
         !is.na(sedentary_overall_average))
  )

## Regression function - PA and cancer adjusted for age and sex
perform_cancer_analysis_age_sex <- function(cancer_type, data) {
  
  # Variables for logistic regression
  outcome_vars_age_sex <- c(paste0(cancer_type, "_status_combined_control"), paste0(cancer_type, "_incidental_combined_control"))
  predictors_age_sex <- c("overall_acceleration_average", "moderate_to_vigorous_quantile", "sedentary_overall_average")
  covariates_age_sex <- c("age", "sex")  # Only using age and sex as covariates
  
  # Store results
  results_list_age_sex <- list()
  
  for (outcome in outcome_vars_age_sex) {
    for (predictor in predictors_age_sex) {
      
      # Filter data for the specific cancer type
      filtered_data_age_sex <- data %>%
        filter(!is.na(.data[[outcome]]) & !is.na(.data[[predictor]]))
      
      # Apply sex-specific filtering
      if (cancer_type == "EC" || cancer_type == "BC") {
        filtered_data_age_sex <- filtered_data_age_sex %>% filter(sex == "Female")
      } else if (cancer_type == "PC") {
        filtered_data_age_sex <- filtered_data_age_sex %>% filter(sex == "Male")
      }
      
      # Filter for missing values in covariates (only age and sex now)
      filtered_data_age_sex <- filtered_data_age_sex %>%
        filter(!is.na(age) & !is.na(sex))
      
      # Check number of rows after filtering
      cat("Rows for", outcome, "and", predictor, "after filtering: ", nrow(filtered_data_age_sex), "\n")
      
      # Skip the model if there are no rows
      if (nrow(filtered_data_age_sex) == 0) {
        cat("Skipping", outcome, "and", predictor, "due to no valid data.\n")
        next
      }
      
      # Check if the covariates (age and sex) have more than one level
      covariate_levels <- sapply(covariates_age_sex, function(covariate) {
        length(unique(filtered_data_age_sex[[covariate]]))
      })
      
      # If "sex" has fewer than 2 levels, exclude it from the model
      if (covariate_levels["sex"] < 2) {
        covariates_to_use <- covariates_age_sex[covariates_age_sex != "sex"]
        cat("Skipping 'sex' in the model due to less than 2 levels.\n")
      } else {
        covariates_to_use <- covariates_age_sex
      }
      
      # Convert outcome to factor if it's not already
      filtered_data_age_sex[[outcome]] <- as.factor(filtered_data_age_sex[[outcome]])
      
      # Count the cases and controls in the filtered data
      case_control_counts_age_sex <- filtered_data_age_sex %>%
        group_by(status = .data[[outcome]]) %>%
        summarise(count = n(), .groups = 'drop')
      
      num_cases_age_sex <- case_control_counts_age_sex$count[case_control_counts_age_sex$status == "case"]
      num_controls_age_sex <- case_control_counts_age_sex$count[case_control_counts_age_sex$status == "control"]
      
      # Fit the logistic regression model with the appropriate covariates
      model_age_sex <- glm(as.formula(paste(outcome, "~", paste(c(predictor, covariates_to_use), collapse = "+"))), 
                           data = filtered_data_age_sex, family = binomial)
      
      # Collect the results
      results_age_sex <- broom::tidy(model_age_sex, conf.int = TRUE) %>%
        filter(term != "(Intercept)") %>%
        mutate(
          cancer_type = cancer_type,
          outcome = outcome,
          predictor = predictor,
          N = nrow(filtered_data_age_sex),  # Number of observations used in regression
          case_control_label = paste("Cases: ", ifelse(is.na(num_cases_age_sex), 0, num_cases_age_sex), 
                                     ", Controls: ", ifelse(is.na(num_controls_age_sex), 0, num_controls_age_sex))
        )
      
      results_list_age_sex[[length(results_list_age_sex) + 1]] <- results_age_sex
    }
  }
  
  # Combine all results for this cancer
  combined_results_age_sex <- do.call(rbind, results_list_age_sex)
  
  return(combined_results_age_sex)
}

# List of cancer types to analyze
cancer_types_age_sex <- c("CRC", "BC", "PC", "EC")

all_cancer_results_age_sex <- list()

# Loop over each cancer and analyze
for (i in cancer_types_age_sex) {
  cancer_results_age_sex <- perform_cancer_analysis_age_sex(i, cancer_PA_cc)
  all_cancer_results_age_sex[[i]] <- cancer_results_age_sex
}

# Combine results for all cancers
final_results_age_sex <- do.call(rbind, all_cancer_results_age_sex)

# Save combined results
write.table(final_results_age_sex, file = paste0(results_dir, "tables/all_cancer_PA_regression_results_age_sex.txt"), col.names = TRUE, row.names = FALSE, sep = "\t")

# Restrict to relevant rows
final_results_filtered_age_sex <- final_results_age_sex %>%
  filter(term %in% c("overall_acceleration_average", "moderate_to_vigorous_quantile", "sedentary_overall_average")) %>%
  filter(outcome %in% c("CRC_status_combined_control", "CRC_incidental_combined_control", 
                        "BC_status_combined_control", "BC_incidental_combined_control", 
                        "PC_status_combined_control", "PC_incidental_combined_control", 
                        "EC_status_combined_control", "EC_incidental_combined_control"))

# Save restricted results
write.table(final_results_filtered_age_sex, file = paste0(results_dir, "tables/all_cancer_PA_regression_results_filtered_age_sex.txt"), 
            col.names = TRUE, row.names = FALSE, sep = "\t")


## Regression function - PA and cancer adjusted (without BMI)
perform_cancer_analysis <- function(cancer_type, data) {
  
  # Variables for logistic regression
  outcome_vars <- c(paste0(cancer_type, "_status_combined_control"), paste0(cancer_type, "_incidental_combined_control"))
  predictors <- c("overall_acceleration_average", "moderate_to_vigorous_quantile", "sedentary_overall_average")
  covariates <- c("age", "smoking", "alcohol", "centre", "education")
  
  # Store results
  results_list <- list()
  
  for (outcome in outcome_vars) {
    for (predictor in predictors) {
      
      # Filter data for the specific cancer type
      filtered_data <- data %>%
        filter(!is.na(.data[[outcome]]) & !is.na(.data[[predictor]]))
      
      # Apply sex-specific filtering
      if (cancer_type == "EC" || cancer_type == "BC") {
        filtered_data <- filtered_data %>% filter(sex == "Female")
      } else if (cancer_type == "PC") {
        filtered_data <- filtered_data %>% filter(sex == "Male")
      }
      
      # Filter for missing values in covariates
      filtered_data <- filtered_data %>%
        filter(!is.na(age) & !is.na(smoking) & !is.na(alcohol) & 
                 !is.na(centre) & !is.na(fasting) & !is.na(education) & !is.na(bmi))
      
      # Check number of rows after filtering
      cat("Rows for", outcome, "and", predictor, "after filtering: ", nrow(filtered_data), "\n")
      
      # Skip the model if there are no rows
      if (nrow(filtered_data) == 0) {
        cat("Skipping", outcome, "and", predictor, "due to no valid data.\n")
        next
      }
      
      # Count the cases and controls in the filtered data
      case_control_counts <- filtered_data %>%
        group_by(status = .data[[outcome]]) %>%
        summarise(count = n(), .groups = 'drop')
      
      num_cases <- case_control_counts$count[case_control_counts$status == "case"]
      num_controls <- case_control_counts$count[case_control_counts$status == "control"]
      
      # Fit the logistic regression model
      model <- glm(as.formula(paste(outcome, "~", paste(c(predictor, covariates), collapse = "+"))), 
                   data = filtered_data, family = binomial)
      
      # Collect the results
      results <- broom::tidy(model, conf.int = TRUE) %>%
        filter(term != "(Intercept)") %>%
        mutate(
          cancer_type = cancer_type,
          outcome = outcome,
          predictor = predictor,
          N = nrow(filtered_data),  # Number of observations used in regression
          case_control_label = paste("Cases: ", ifelse(is.na(num_cases), 0, num_cases), 
                                     ", Controls: ", ifelse(is.na(num_controls), 0, num_controls))
        )
      
      results_list[[length(results_list) + 1]] <- results
    }
  }
  
  # Combine all results for this cancer
  combined_results <- do.call(rbind, results_list)
  
  return(combined_results)
}

all_cancer_results_fully_adjusted <- list()

# Loop over each cancer and analyze
for (i in cancer_types) {
  cancer_results_adjusted <- perform_cancer_analysis(i, cancer_PA_cc)
  all_cancer_results_fully_adjusted[[i]] <- cancer_results_adjusted
}

# Combine results for all cancers
final_results_adjusted <- do.call(rbind, all_cancer_results_fully_adjusted)

# Save combined results
write.table(final_results_adjusted, file = paste0(results_dir, "tables/all_cancer_PA_regression_results_adjusted.txt"), col.names = TRUE, row.names = FALSE, sep = "\t")

# Restrict to relevant rows
final_results_adjusted_filtered <- final_results_adjusted %>%
  filter(term %in% c("overall_acceleration_average", "moderate_to_vigorous_quantile", "sedentary_overall_average")) %>%
  filter(outcome %in% c("CRC_status_combined_control", "CRC_incidental_combined_control", 
                        "BC_status_combined_control", "BC_incidental_combined_control", 
                        "PC_status_combined_control", "PC_incidental_combined_control", 
                        "EC_status_combined_control", "EC_incidental_combined_control"))

# Save restricted results
write.table(final_results_adjusted_filtered, file = paste0(results_dir, "tables/all_cancer_PA_regression_results_adjusted_filtered.txt"), 
  col.names = TRUE, row.names = FALSE, sep = "\t")

## Regression function - PA and cancer fully adjusted
perform_cancer_analysis <- function(cancer_type, data) {
  
  # Variables for logistic regression
  outcome_vars <- c(paste0(cancer_type, "_status_combined_control"), paste0(cancer_type, "_incidental_combined_control"))
  predictors <- c("overall_acceleration_average", "moderate_to_vigorous_quantile", "sedentary_overall_average")
  covariates <- c("age", "smoking", "alcohol", "centre", "education", "bmi")
  
  # Store results
  results_list <- list()
  
  for (outcome in outcome_vars) {
    for (predictor in predictors) {
      
      # Filter data for the specific cancer type
      filtered_data <- data %>%
        filter(!is.na(.data[[outcome]]) & !is.na(.data[[predictor]]))
      
      # Apply sex-specific filtering
      if (cancer_type == "EC" || cancer_type == "BC") {
        filtered_data <- filtered_data %>% filter(sex == "Female")
      } else if (cancer_type == "PC") {
        filtered_data <- filtered_data %>% filter(sex == "Male")
      }
      
      # Filter for missing values in covariates
      filtered_data <- filtered_data %>%
        filter(!is.na(age) & !is.na(smoking) & !is.na(alcohol) & 
                 !is.na(centre) & !is.na(fasting) & !is.na(education) & !is.na(bmi))
      
      # Check number of rows after filtering
      cat("Rows for", outcome, "and", predictor, "after filtering: ", nrow(filtered_data), "\n")
      
      # Skip the model if there are no rows
      if (nrow(filtered_data) == 0) {
        cat("Skipping", outcome, "and", predictor, "due to no valid data.\n")
        next
      }
      
      # Count the cases and controls in the filtered data
      case_control_counts <- filtered_data %>%
        group_by(status = .data[[outcome]]) %>%
        summarise(count = n(), .groups = 'drop')
      
      num_cases <- case_control_counts$count[case_control_counts$status == "case"]
      num_controls <- case_control_counts$count[case_control_counts$status == "control"]
      
      # Fit the logistic regression model
      model <- glm(as.formula(paste(outcome, "~", paste(c(predictor, covariates), collapse = "+"))), 
                   data = filtered_data, family = binomial)
      
      # Collect the results
      results <- broom::tidy(model, conf.int = TRUE) %>%
        filter(term != "(Intercept)") %>%
        mutate(
          cancer_type = cancer_type,
          outcome = outcome,
          predictor = predictor,
          N = nrow(filtered_data),  # Number of observations used in regression
          case_control_label = paste("Cases: ", ifelse(is.na(num_cases), 0, num_cases), 
                                     ", Controls: ", ifelse(is.na(num_controls), 0, num_controls))
        )
      
      results_list[[length(results_list) + 1]] <- results
    }
  }
  
  # Combine all results for this cancer
  combined_results <- do.call(rbind, results_list)
  
  return(combined_results)
}

all_cancer_results <- list()

# Loop over each cancer and analyze
for (i in cancer_types) {
  cancer_results <- perform_cancer_analysis(i, cancer_PA_cc)
  all_cancer_results[[i]] <- cancer_results
}

# Combine results for all cancers
final_results <- do.call(rbind, all_cancer_results)

# Save combined results
write.table(final_results, file = paste0(results_dir, "tables/all_cancer_PA_regression_results_fully_adjusted_bmi.txt"), col.names = TRUE, row.names = FALSE, sep = "\t")

# Restrict to relevant rows
final_results_filtered <- final_results %>%
  filter(term %in% c("overall_acceleration_average", "moderate_to_vigorous_quantile", "sedentary_overall_average")) %>%
  filter(outcome %in% c("CRC_status_combined_control", "CRC_incidental_combined_control", 
                        "BC_status_combined_control", "BC_incidental_combined_control", 
                        "PC_status_combined_control", "PC_incidental_combined_control", 
                        "EC_status_combined_control", "EC_incidental_combined_control"))

# Save restricted results
write.table(final_results_filtered, file = paste0(results_dir, "tables/all_cancer_PA_regression_results_filtered_fully_adjusted_bmi.txt"), 
            col.names = TRUE, row.names = FALSE, sep = "\t")

## Do the Cox PH model for PA and cancer
## Read in first date of wear variable
date_wear <- read.delim(
  file = paste0(input_dir, "16391_from_Becky/output_extracted_variables_16391_lucy_8.1.26.tsv"),
  stringsAsFactors = FALSE
)
date_wear <- date_wear[,c("projectID", "f.90003.0.0")]

cancer_PA_cc <- merge(cancer_PA_cc, date_wear, by.x = "eid", by.y = "projectID", all.x = T, all.y = F)
cancer_PA_cc$f.90003.0.0 <- as.Date(cancer_PA_cc$f.90003.0.0, format = "%Y-%m-%d %H:%M:%S")

## Function to get event date
get_event_date <- function(data, cancer) {
  icd10 <- paste0(cancer, "_ICD10_date")
  icd9  <- paste0(cancer, "_ICD9_date")
  
  as.Date(ifelse(
    !is.na(data[[icd10]]),
    data[[icd10]],
    data[[icd9]]
  ), origin = "1970-01-01")
}

## Cox PH with delayed entry (prevalent cases included)
perform_cancer_cox_analysis <- function(cancer_type, data, include_bmi = TRUE, ENTRY_DATE = as.Date("2012-08-01"), ADMIN_CENSOR_DATE = as.Date("2020-03-31")) {
  
  predictors <- c("overall_acceleration_average",
                  "moderate_to_vigorous_quantile",
                  "sedentary_overall_average")
  
  covariates <- c("age", "smoking", "alcohol", "centre", "education")
  if (include_bmi) covariates <- c(covariates, "bmi")
  
  # Get event date
  event_date <- get_event_date(data, cancer_type)
  
  results_list <- list()
  
  for (predictor in predictors) {
    
    df <- data %>%
      mutate(
        event_date = event_date,
        entry_date = ENTRY_DATE
      ) %>%
      
      # Sex-specific restriction
      {if (cancer_type %in% c("EC", "BC")) filter(., sex == "Female") else .} %>%
      {if (cancer_type == "PC") filter(., sex == "Male") else .} %>%
      
      # Keep complete cases
      filter(!is.na(.data[[predictor]]), complete.cases(across(all_of(covariates)))) %>%
      
      # Delayed entry: everyone starts at ENTRY_DATE
      mutate(
        t_entry = as.numeric(entry_date),
        t_exit  = as.numeric(ifelse(!is.na(event_date) & event_date > entry_date & event_date <= ADMIN_CENSOR_DATE,
                                    event_date,
                                    ADMIN_CENSOR_DATE)),
        event   = ifelse(!is.na(event_date) & event_date > entry_date & event_date <= ADMIN_CENSOR_DATE, 1, 0)
      )
    
    # Skip if too few participants/events
    if (nrow(df) < 50 || sum(df$event) < 10) next
    
    # Fit Cox PH model
    form <- reformulate(
      termlabels = c(predictor, covariates),
      response = "Surv(t_entry, t_exit, event)"
    )
    
    fit <- survival::coxph(form, data = df)
    
    # Store results
    res <- broom::tidy(fit, conf.int = TRUE) %>%
      filter(term == predictor) %>%
      mutate(
        cancer_type = cancer_type,
        predictor = predictor,
        N = nrow(df),
        events = sum(df$event)
      )
    
    results_list[[length(results_list) + 1]] <- res
  }
  
  bind_rows(results_list)
}

## -------------------------
## Run Cox PH for all cancer types
## -------------------------

## Run analysis for all cancer types (with BMI)
all_cancer_results_bmi <- lapply(cancer_types, function(ct) perform_cancer_cox_analysis(ct, cancer_PA_cc, include_bmi = TRUE))
final_results_bmi <- bind_rows(all_cancer_results_bmi)

write.table(
  final_results_bmi,
  file = paste0(results_dir, "tables/all_cancer_PA_cox_results_entry2012_fully_adjusted_bmi.txt"),
  col.names = TRUE,
  row.names = FALSE,
  sep = "\t"
)

## Run analysis for all cancer types (without BMI)
all_cancer_results_no_bmi <- lapply(cancer_types, function(ct) perform_cancer_cox_analysis(ct, cancer_PA_cc, include_bmi = FALSE))
final_results_no_bmi <- bind_rows(all_cancer_results_no_bmi)

write.table(
  final_results_no_bmi,
  file = paste0(results_dir, "tables/all_cancer_PA_cox_results_entry2012_no_bmi.txt"),
  col.names = TRUE,
  row.names = FALSE,
  sep = "\t"
)

## Change cut off date to 2015 (with BMI)
perform_cancer_cox_analysis_2015 <- function(cancer_type, data, include_bmi = TRUE, ENTRY_DATE = as.Date("2015-12-31"), ADMIN_CENSOR_DATE = as.Date("2020-03-31")) {
  
  predictors <- c("overall_acceleration_average",
                  "moderate_to_vigorous_quantile",
                  "sedentary_overall_average")
  
  covariates <- c("age", "smoking", "alcohol", "centre", "education")
  if (include_bmi) covariates <- c(covariates, "bmi")
  
  # Get event date
  event_date <- get_event_date(data, cancer_type)
  
  results_list <- list()
  
  for (predictor in predictors) {
    
    df <- data %>%
      mutate(
        event_date = event_date,
        entry_date = ENTRY_DATE
      ) %>%
      
      # Sex-specific restriction
      {if (cancer_type %in% c("EC", "BC")) filter(., sex == "Female") else .} %>%
      {if (cancer_type == "PC") filter(., sex == "Male") else .} %>%
      
      # Keep complete cases
      filter(!is.na(.data[[predictor]]), complete.cases(across(all_of(covariates)))) %>%
      
      # Delayed entry: everyone starts at ENTRY_DATE
      mutate(
        t_entry = as.numeric(entry_date),
        t_exit  = as.numeric(ifelse(!is.na(event_date) & event_date > entry_date & event_date <= ADMIN_CENSOR_DATE,
                                    event_date,
                                    ADMIN_CENSOR_DATE)),
        event   = ifelse(!is.na(event_date) & event_date > entry_date & event_date <= ADMIN_CENSOR_DATE, 1, 0)
      )
    
    # Skip if too few participants/events
    if (nrow(df) < 50 || sum(df$event) < 10) next
    
    # Fit Cox PH model
    form <- reformulate(
      termlabels = c(predictor, covariates),
      response = "Surv(t_entry, t_exit, event)"
    )
    
    fit <- survival::coxph(form, data = df)
    
    # Store results
    res <- broom::tidy(fit, conf.int = TRUE) %>%
      filter(term == predictor) %>%
      mutate(
        cancer_type = cancer_type,
        predictor = predictor,
        N = nrow(df),
        events = sum(df$event)
      )
    
    results_list[[length(results_list) + 1]] <- res
  }
  
  bind_rows(results_list)
}

all_cancer_results_bmi_2015 <- lapply(cancer_types, function(ct) perform_cancer_cox_analysis_2015(ct, cancer_PA_cc, include_bmi = TRUE))
final_results_bmi_2015 <- bind_rows(all_cancer_results_bmi_2015)

write.table(
  final_results_bmi_2015,
  file = paste0(results_dir, "tables/all_cancer_PA_cox_results_entry2015_bmi.txt"),
  col.names = TRUE,
  row.names = FALSE,
  sep = "\t"
)

## Change cut off date to 2015 (without BMI)
perform_cancer_cox_analysis_2015 <- function(cancer_type, data, include_bmi = FALSE, ENTRY_DATE = as.Date("2015-12-31"), ADMIN_CENSOR_DATE = as.Date("2020-03-31")) {
  
  predictors <- c("overall_acceleration_average",
                  "moderate_to_vigorous_quantile",
                  "sedentary_overall_average")
  
  covariates <- c("age", "smoking", "alcohol", "centre", "education")
  if (include_bmi) covariates <- c(covariates, "bmi")
  
  # Get event date
  event_date <- get_event_date(data, cancer_type)
  
  results_list <- list()
  
  for (predictor in predictors) {
    
    df <- data %>%
      mutate(
        event_date = event_date,
        entry_date = ENTRY_DATE
      ) %>%
      
      # Sex-specific restriction
      {if (cancer_type %in% c("EC", "BC")) filter(., sex == "Female") else .} %>%
      {if (cancer_type == "PC") filter(., sex == "Male") else .} %>%
      
      # Keep complete cases
      filter(!is.na(.data[[predictor]]), complete.cases(across(all_of(covariates)))) %>%
      
      # Delayed entry: everyone starts at ENTRY_DATE
      mutate(
        t_entry = as.numeric(entry_date),
        t_exit  = as.numeric(ifelse(!is.na(event_date) & event_date > entry_date & event_date <= ADMIN_CENSOR_DATE,
                                    event_date,
                                    ADMIN_CENSOR_DATE)),
        event   = ifelse(!is.na(event_date) & event_date > entry_date & event_date <= ADMIN_CENSOR_DATE, 1, 0)
      )
    
    # Skip if too few participants/events
    if (nrow(df) < 50 || sum(df$event) < 10) next
    
    # Fit Cox PH model
    form <- reformulate(
      termlabels = c(predictor, covariates),
      response = "Surv(t_entry, t_exit, event)"
    )
    
    fit <- survival::coxph(form, data = df)
    
    # Store results
    res <- broom::tidy(fit, conf.int = TRUE) %>%
      filter(term == predictor) %>%
      mutate(
        cancer_type = cancer_type,
        predictor = predictor,
        N = nrow(df),
        events = sum(df$event)
      )
    
    results_list[[length(results_list) + 1]] <- res
  }
  
  bind_rows(results_list)
}

all_cancer_results_no_bmi_2015 <- lapply(cancer_types, function(ct) perform_cancer_cox_analysis_2015(ct, cancer_PA_cc, include_bmi = FALSE))
final_results_no_bmi_2015 <- bind_rows(all_cancer_results_no_bmi_2015)

write.table(
  final_results_no_bmi_2015,
  file = paste0(results_dir, "tables/all_cancer_PA_cox_results_entry2015_no_bmi.txt"),
  col.names = TRUE,
  row.names = FALSE,
  sep = "\t"
)

## Redo PA-cancer Cox PH but use the wear date
perform_cancer_cox_analysis_cutoff <- function(cancer_type, data, include_bmi = TRUE, 
                                               cutoff_years = 0, 
                                               ADMIN_CENSOR_DATE = as.Date("2020-03-31")) {
  
  predictors <- c("overall_acceleration_average",
                  "moderate_to_vigorous_quantile",
                  "sedentary_overall_average")
  
  covariates <- c("age", "smoking", "alcohol", "centre", "education")
  if (include_bmi) covariates <- c(covariates, "bmi")
  
  # Get event date
  event_date <- get_event_date(data, cancer_type)
  
  results_list <- list()
  
  for (predictor in predictors) {
    
    df <- data %>%
      mutate(
        event_date = event_date,
        # Delayed entry: add cutoff_years to accelerometer date
        entry_date = f.90003.0.0 + years(cutoff_years)
      ) %>%
      
      # Sex-specific restriction
      {if (cancer_type %in% c("EC", "BC")) filter(., sex == "Female") else .} %>%
      {if (cancer_type == "PC") filter(., sex == "Male") else .} %>%
      
      # Keep complete cases and remove those with missing entry date
      filter(!is.na(.data[[predictor]]), complete.cases(across(all_of(covariates))),
             !is.na(entry_date)) %>%
      
      # Remove events that occur before the cutoff (t_entry)
      filter(is.na(event_date) | event_date >= entry_date) %>%
      
      # Compute time for Cox model
      mutate(
        t_entry = as.numeric(entry_date),
        t_exit  = as.numeric(ifelse(!is.na(event_date) & event_date <= ADMIN_CENSOR_DATE,
                                    event_date,
                                    ADMIN_CENSOR_DATE)),
        event   = ifelse(!is.na(event_date) & event_date <= ADMIN_CENSOR_DATE, 1, 0)
      )
    
    # Skip if too few participants/events
    if (nrow(df) < 50 || sum(df$event) < 10) next
    
    # Fit Cox PH model
    form <- reformulate(
      termlabels = c(predictor, covariates),
      response = "Surv(t_entry, t_exit, event)"
    )
    
    fit <- coxph(form, data = df)
    
    # Store results
    res <- broom::tidy(fit, conf.int = TRUE) %>%
      filter(term == predictor) %>%
      mutate(
        cancer_type = cancer_type,
        predictor = predictor,
        N = nrow(df),
        events = sum(df$event),
        cutoff_years = cutoff_years
      )
    
    results_list[[length(results_list) + 1]] <- res
  }
  
  bind_rows(results_list)
}


# 1-year cutoff
all_cancer_bmi_1yr <- lapply(cancer_types, function(ct) {
  perform_cancer_cox_analysis_cutoff(ct, cancer_PA_cc, include_bmi = TRUE, cutoff_years = 1)
})
final_bmi_1yr <- bind_rows(all_cancer_bmi_1yr)

# 2-year cutoff
all_cancer_bmi_2yr <- lapply(cancer_types, function(ct) {
  perform_cancer_cox_analysis_cutoff(ct, cancer_PA_cc, include_bmi = TRUE, cutoff_years = 2)
})
final_bmi_2yr <- bind_rows(all_cancer_bmi_2yr)

# 1-year cutoff
all_cancer_no_bmi_1yr <- lapply(cancer_types, function(ct) {
  perform_cancer_cox_analysis_cutoff(ct, cancer_PA_cc, include_bmi = FALSE, cutoff_years = 1)
})
final_no_bmi_1yr <- bind_rows(all_cancer_no_bmi_1yr)

# 2-year cutoff
all_cancer_no_bmi_2yr <- lapply(cancer_types, function(ct) {
  perform_cancer_cox_analysis_cutoff(ct, cancer_PA_cc, include_bmi = FALSE, cutoff_years = 2)
})
final_no_bmi_2yr <- bind_rows(all_cancer_no_bmi_2yr)

## Write all of these to file 
# Define output file paths
output_files <- list(
  final_bmi_1yr = paste0(results_dir, "tables/all_cancer_PA_cox_results_bmi_1yr_cutoff.txt"),
  final_bmi_2yr = paste0(results_dir, "tables/all_cancer_PA_cox_results_bmi_2yr_cutoff.txt"),
  final_no_bmi_1yr = paste0(results_dir, "tables/all_cancer_PA_cox_results_no_bmi_1yr_cutoff.txt"),
  final_no_bmi_2yr = paste0(results_dir, "tables/all_cancer_PA_cox_results_no_bmi_2yr_cutoff.txt")
)

# Create a named list of results
results_list <- list(
  final_bmi_1yr = final_bmi_1yr,
  final_bmi_2yr = final_bmi_2yr,
  final_no_bmi_1yr = final_no_bmi_1yr,
  final_no_bmi_2yr = final_no_bmi_2yr
)

# Write each to file
for (name in names(results_list)) {
  write.table(
    results_list[[name]],
    file = output_files[[name]],
    sep = "\t",
    col.names = TRUE,
    row.names = FALSE
  )
}

## Regression function - PA and cancer fully adjusted (with PCs)
perform_cancer_analysis <- function(cancer_type, data) {
  
  # Variables for logistic regression
  outcome_vars <- c(paste0(cancer_type, "_status_combined_control"), paste0(cancer_type, "_incidental_combined_control"))
  predictors <- c("overall_acceleration_average", "moderate_to_vigorous_quantile", "sedentary_overall_average")
  covariates <- c("age", "smoking", "alcohol", "centre", "education", "bmi", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7",
                  "PC8", "PC9", "PC10")
  
  # Store results
  results_list <- list()
  
  for (outcome in outcome_vars) {
    for (predictor in predictors) {
      
      # Filter data for the specific cancer type
      filtered_data <- data %>%
        filter(!is.na(.data[[outcome]]) & !is.na(.data[[predictor]]))
      
      # Apply sex-specific filtering
      if (cancer_type == "EC" || cancer_type == "BC") {
        filtered_data <- filtered_data %>% filter(sex == "Female")
      } else if (cancer_type == "PC") {
        filtered_data <- filtered_data %>% filter(sex == "Male")
      }
      
      # Filter for missing values in covariates
      filtered_data <- filtered_data %>%
        filter(!is.na(age) & !is.na(smoking) & !is.na(alcohol) & 
                 !is.na(centre) & !is.na(fasting) & !is.na(education) & !is.na(bmi))
      
      # Check number of rows after filtering
      cat("Rows for", outcome, "and", predictor, "after filtering: ", nrow(filtered_data), "\n")
      
      # Skip the model if there are no rows
      if (nrow(filtered_data) == 0) {
        cat("Skipping", outcome, "and", predictor, "due to no valid data.\n")
        next
      }
      
      # Count the cases and controls in the filtered data
      case_control_counts <- filtered_data %>%
        group_by(status = .data[[outcome]]) %>%
        summarise(count = n(), .groups = 'drop')
      
      num_cases <- case_control_counts$count[case_control_counts$status == "case"]
      num_controls <- case_control_counts$count[case_control_counts$status == "control"]
      
      # Fit the logistic regression model
      model <- glm(as.formula(paste(outcome, "~", paste(c(predictor, covariates), collapse = "+"))), 
                   data = filtered_data, family = binomial)
      
      # Collect the results
      results <- broom::tidy(model, conf.int = TRUE) %>%
        filter(term != "(Intercept)") %>%
        mutate(
          cancer_type = cancer_type,
          outcome = outcome,
          predictor = predictor,
          N = nrow(filtered_data),  # Number of observations used in regression
          case_control_label = paste("Cases: ", ifelse(is.na(num_cases), 0, num_cases), 
                                     ", Controls: ", ifelse(is.na(num_controls), 0, num_controls))
        )
      
      results_list[[length(results_list) + 1]] <- results
    }
  }
  
  # Combine all results for this cancer
  combined_results <- do.call(rbind, results_list)
  
  return(combined_results)
}

all_cancer_results_PCs <- list()

# Loop over each cancer and analyze
for (i in cancer_types) {
  cancer_results_PCs <- perform_cancer_analysis(i, cancer_PA_cc)
  all_cancer_results_PCs[[i]] <- cancer_results_PCs
}

# Combine results for all cancers
final_results_PCs <- do.call(rbind, all_cancer_results_PCs)

# Save combined results
write.table(final_results_PCs, file = paste0(results_dir, "tables/all_cancer_PA_regression_results_fully_adjusted_PCs_bmi.txt"), col.names = TRUE, row.names = FALSE, sep = "\t")

# Restrict to relevant rows
final_results_filtered_PCs <- final_results_PCs %>%
  filter(term %in% c("overall_acceleration_average", "moderate_to_vigorous_quantile", "sedentary_overall_average")) %>%
  filter(outcome %in% c("CRC_status_combined_control", "CRC_incidental_combined_control", 
                        "BC_status_combined_control", "BC_incidental_combined_control", 
                        "PC_status_combined_control", "PC_incidental_combined_control", 
                        "EC_status_combined_control", "EC_incidental_combined_control"))

# Save restricted results
write.table(final_results_filtered_PCs, file = paste0(results_dir, "tables/all_cancer_PA_regression_results_filtered_fully_adjusted_PCs_bmi.txt"), 
            col.names = TRUE, row.names = FALSE, sep = "\t")

# Define full names for each cancer type
# Mapping cancer abbreviations to full names
cancer_full_names <- c(
  "CRC" = "Colorectal Cancer",
  "EC" = "Endometrial Cancer",
  "BC" = "Breast Cancer",
  "PC" = "Prostate Cancer"
)

# List of physical activity variables
activity_vars <- c("overall_acceleration_average", "moderate_to_vigorous_quantile", "sedentary_overall_average")

# Map physical activity variable names to full descriptions
activity_full_names <- c(
  "overall_acceleration_average" = "overall activity",
  "moderate_to_vigorous_quantile" = "moderate to vigorous activity",
  "sedentary_overall_average" = "sedentary"
)

# Loop over each activity variable and plot
for (activity in activity_vars) {
  
  # Filter final results to include only the rows for the current activity
  filtered_activity_results <- final_results_filtered %>%
    filter(predictor == activity) %>%
    mutate(
      cancer = factor(cancer_full_names[cancer_type], levels = cancer_full_names),  # Add the full cancer name
      outcome_label = case_when(
        grepl("_status_combined_control", outcome) ~ "Prevalent",  # Rename _status to Prevalent
        grepl("_incidental_combined_control", outcome) ~ "Incident",
        TRUE ~ "Other"
      )
    )
  
  # Create the plot with dodged points
  position_dodge_val <- position_dodge(width = 0.5)
  
  plot <- ggplot(filtered_activity_results, aes(x = cancer, y = estimate, ymin = conf.low, ymax = conf.high, color = outcome_label)) +
    geom_pointrange(position = position_dodge_val) +  # Points are dodged
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +  # Red dashed line at y = 0
    coord_flip() +
    labs(
      x = "Cancer Type", 
      y = paste0("Log odds of cancer risk per unit higher physical activity measure (95% CI)"),
      title = activity_full_names[activity],
      color = "Outcome Type"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Save the plot to a PDF file with the activity variable name
  ggsave(filename = paste0(results_dir, "figures/forest_PA_", activity, "_cancer_inc_prev.pdf"), plot = plot, width = 7, height = 7)
}

## Plot the three models of physical activity and cancer
# Combine all three data frames into one
final_results_adjusted_filtered$model <- "Fully adjusted (excluding BMI)"
final_results_filtered_age_sex$model <- "Age and sex adjusted"
final_results_filtered$model <- "Fully adjusted (including BMI)"
final_results_filtered_PCs$model <- "Fully adjusted (including BMI and PCs)"

# Combine the three data frames
combined_results <- bind_rows(final_results_adjusted_filtered,
                              final_results_filtered_age_sex,
                              final_results_filtered,
                              final_results_filtered_PCs)

# Add a new column to distinguish between "_status" and "_incidental" outcomes
combined_results <- combined_results %>%
  mutate(
    outcome_type = ifelse(grepl("_status_combined_control$", outcome), "Prevalent", "Incident"),
    physical_activity_type = case_when(
      grepl("overall_acceleration_average", predictor) ~ "Overall acceleration",
      grepl("moderate_to_vigorous_quantile", predictor) ~ "Moderate to vigorous activity",
      grepl("sedentary_overall_average", predictor) ~ "Sedentary",
      TRUE ~ "Other"
    )
  )

combined_results <- combined_results %>%
  mutate(cancer_type = recode(cancer_type,
                              CRC = "Colorectal",
                              EC = "Endometrial",
                              BC = "Breast",
                              PC = "Prostate"))

# Order the outcomes (Prevalent first, then Incident)
combined_results$outcome_type <- factor(combined_results$outcome_type, 
                                        levels = c("Prevalent", "Incident"))

# Plot and save
pdf(paste0(results_dir, "figures/forest_PA_cancer_all_models.pdf"), height = 12, width = 12)
ggplot(combined_results, aes(x = estimate, y = cancer_type, color = model, shape = outcome_type)) +
  geom_point(position = position_dodge(width = 0.5), size = 3) +  # Separate points by cancer type and outcome type
  geom_errorbarh(aes(xmin = estimate - 1.96 * std.error, xmax = estimate + 1.96 * std.error),
                 position = position_dodge(width = 0.5), height = 0.2) +  # Add error bars for confidence intervals
  facet_wrap(~ physical_activity_type, scales = "free_y", ncol = 1) +  # One plot for each PA type
  geom_vline(xintercept = 0, linetype = "dotted", color = "black") +  # Add a dotted line at beta = 0
  theme_minimal() +
  labs(
    x = "Log odds of cancer risk per unit higher activity measure (±95% CI)",
    y = "Cancer site",
    title = "",
    color = "Model",
    shape = "Outcome Type",
    linetype = "Model Adjustment"
  ) +
  theme(
    legend.position = "top",
    strip.text = element_text(size = 12),  # Adjust facet labels font size
    axis.text.y = element_text(size = 10),  # Adjust y-axis labels font size
    legend.key = element_blank(),  # To avoid a box around legend keys
    axis.title.x = element_text(size = 12),  # Adjust x-axis title font size
    axis.title.y = element_text(size = 12)   # Adjust y-axis title font size
  )
dev.off()

## Proteins associated with all measures of PA
protein_PA_names <- fread(file = paste0(results_dir, "tables/overall_PA_protein_regression_fully_adjusted_bmi.txt"), header = T, sep = "\t")
threshold <- 0.05 / 1885
protein_PA_name_list <- protein_PA_names$protein[protein_PA_names$P_lm < threshold]
protein_PA_name_list_2 <- toupper(gsub("^rnt_", "", protein_PA_name_list))
write.table(protein_PA_name_list_2, 
            file = paste0(results_dir, "tables/proteins_for_MR.txt"), 
            col.names = FALSE, 
            row.names = FALSE, 
            sep = "\t", 
            quote = FALSE)

logistic_model <- function(wdata, dependent, independent, covariables) {
  ##############################
  ### 1. Define Model Data Frame
  ##############################
  if (is.na(covariables[1])) {
    model_variables = c(dependent, independent)
  } else {
    model_variables = c(dependent, covariables, independent)
  }
  
  mod_data = wdata[, c(model_variables)]
  
  ##############################
  ### 2. Perform Logistic Model 
  ##############################
  if (is.na(covariables[1])) {
    form = formula(paste0(dependent, " ~ ", independent))  
  } else {
    form = formula(paste0(dependent, " ~ ", independent, " + ", paste0(covariables, collapse = " + ")))  
  }
  
  glm_mod <- glm(form, data = mod_data, family = "binomial")
  
  #################
  ## 3. Summary stats
  #################
  protein_name = independent; names(protein_name) = "protein"  # Set protein name as independent variable
  s = summary(glm_mod)
  ## Sample size
  n = length(residuals(glm_mod)); names(n) = "n_glm"
  ## Dependent Variable effect estimates
  beta = s$coefficients[2,1]; names(beta) = "beta_glm"
  se = s$coefficients[2,2]; names(se) = "se_glm"
  pval = s$coefficients[2,4]; names(pval) = "P_glm"
  
  glm_results = c(protein_name, n, beta, se, pval)
  return(glm_results)
}

# Function to compute adjusted p-values after running logistic models
logistic_model_with_fdr <- function(wdata, dependent, proteins, covariables) {
  results <- lapply(proteins, function(protein) {
    logistic_model(wdata, dependent, protein, covariables)
  })
  
  # Convert list to a data frame
  results_df <- do.call(rbind, results)
  results_df <- as.data.frame(results_df)
  
  # Convert numeric columns back to numeric types
  results_df$n_glm <- as.numeric(results_df$n_glm)
  results_df$beta_glm <- as.numeric(results_df$beta_glm)
  results_df$se_glm <- as.numeric(results_df$se_glm)
  results_df$P_glm <- as.numeric(results_df$P_glm)
  
  # Add Benjamini-Hochberg FDR-adjusted p-value
  results_df$P_FDR <- p.adjust(results_df$P_glm, method = "BH")
  
  return(results_df)
}

## Summary table for participants across all analyses
# Merge datasets without renaming again
cancer_PA_protein_table <- dplyr::full_join(cancer_PA_cc, PA_protein, by = c("eid" = "Participant ID"))

## Bring in the raw PA variables for summary
PA_raw <- read.csv(paste0(data_intermediate_dir, "raw_PA.csv"), header = T, sep = ",")
PA_raw <- as.data.frame(PA_raw)

## Remove the outliers that were removed before transformation
## Remove outliers more than 5 SD
PA_cols_QC <- c("overall_acceleration_average", "sedentary_overall_average", 
                "moderate_to_vigorous_overall_average")

# Change outliers to NA in PA data (more than 5 SD away from the mean)
means <- colMeans(PA_raw[PA_cols_QC], na.rm = TRUE)
sds <- apply(PA_raw[PA_cols_QC], 2, sd, na.rm = TRUE)
PA_raw[PA_cols_QC] <- apply(PA_raw[PA_cols_QC], 2, function(column) {
  outliers <- abs(column - mean(column, na.rm = TRUE)) > (5 * sd(column, na.rm = TRUE))
  column[outliers] <- NA
  return(column)
})
names(PA_raw) <- c("eid", "moderate_to_vigorous_overall_average_raw", "sedentary_overall_average", 'overall_acceleration_average')

## Merge these together
PA_raw_restricted <- PA_raw %>%
  dplyr::semi_join(
    cancer_PA_protein_table %>% dplyr::select(eid),
    by = "eid"
  )

PA_summary <- PA_raw_restricted %>%
  summarise(
    N = n(),
    
    mean_overall_acc = mean(overall_acceleration_average, na.rm = TRUE),
    sd_overall_acc   = sd(overall_acceleration_average, na.rm = TRUE),
    
    mean_sedentary = mean(sedentary_overall_average, na.rm = TRUE),
    sd_sedentary   = sd(sedentary_overall_average, na.rm = TRUE),
    
    mean_mvpa = mean(moderate_to_vigorous_overall_average_raw, na.rm = TRUE),
    sd_mvpa   = sd(moderate_to_vigorous_overall_average_raw, na.rm = TRUE)
  )

PA_summary_Ns <- PA_raw_restricted %>%
  summarise(
    N_overall_acc = sum(!is.na(overall_acceleration_average)),
    N_sedentary  = sum(!is.na(sedentary_overall_average)),
    N_mvpa       = sum(!is.na(moderate_to_vigorous_overall_average_raw))
  )

# Use coalesce() to combine overlapping variables correctly
cancer_PA_protein_table <- cancer_PA_protein_table %>%
  mutate(
    education = coalesce(
      as.integer(as.character(education.x)),
      education.y
    )
  )

cancer_PA_protein_table <- cancer_PA_protein_table %>%
  mutate(
    age = coalesce(age.x, age.y),
    sex = coalesce(sex.x, sex.y),
    townsend = coalesce(townsend.x, townsend.y),
    bmi = coalesce(bmi.x, bmi.y),
    smoking = coalesce(smoking.x, smoking.y),
    alcohol = coalesce(alcohol.x, alcohol.y),
    centre = coalesce(centre.x, centre.y),
    moderate_to_vigorous_overall_average = coalesce(moderate_to_vigorous_overall_average.x, moderate_to_vigorous_overall_average.y),
    sedentary_overall_average = coalesce(sedentary_overall_average.x, sedentary_overall_average.y),
    overall_acceleration_average = coalesce(overall_acceleration_average.x, overall_acceleration_average.y)
  )

cat_characteristics_all <- cancer_PA_protein_table[,c("age", "sex", "townsend", "education", "bmi", "smoking", "alcohol", "centre", "moderate_to_vigorous_overall_average",
                                              "sedentary_overall_average", "overall_acceleration_average") ]

# Generate the summary table with updated code
table1 <- tbl_summary(
  cat_characteristics_all,
  statistic = list(
    all_continuous() ~ "{mean} ({sd})",
    all_categorical() ~ "{n} / {N} ({p}%)"
  ),
  missing = "ifany"  # Show missing data only if present
) %>%
  add_n() %>%  # Adds a column showing total N, including missing values
  add_stat_label()  # Optionally adds variable labels to the table

# Convert to Word document
table1 %>%
  as_flex_table() %>%
  flextable::save_as_docx(path = paste0(results_dir, "tables/table1_UKB_characteristics_all.docx"))

## Association between PA-associated proteins and cancer
cancer_PA_protein <- merge(cancer_PA, protein, all.x = F, all.y = T, by.x = "eid", by.y = "Participant ID")

## List of dependent variables
# Initialize an empty list to store results from each iteration
all_results_list <- list()
dependent_variables <- c("PC_incidental_combined_control", "EC_incidental_combined_control", "BC_incidental_combined_control", "CRC_incidental_combined_control")

# Loop over the dependent variables (e.g., PC_incidental, EC_incidental, etc.)
for (dependent_var in dependent_variables) {
  
  # Ensure the dependent variable is a factor (e.g., "PC", "EC", "BC", etc.)
  cancer_PA_protein[[dependent_var]] <- as.factor(cancer_PA_protein[[dependent_var]])
  
  # Run logistic models and add FDR-adjusted p-values
  prot_results <- logistic_model_with_fdr(
    wdata = cancer_PA_protein,
    dependent = dependent_var,
    proteins = protein_PA_name_list,  # List of protein traits
    covariables = c("age", "sex", "smoking", "alcohol", "centre", "bmi", "education")
  )
  
  # Sort results by raw p-values (adjusted p-values)
  prot_results <- prot_results[order(prot_results$P_FDR), ]
  
  # Add a column to indicate cancer type and status type
  cancer_type <- gsub("_incidental_combined_control", "", dependent_var)  # Extract cancer type from the dependent variable name
  prot_results$cancer_type <- cancer_type  # Add cancer type to results
  prot_results$variable_type <- "Incidental"  # Set status type as "Incidental"
  
  # Add the results for this cancer type to the list
  all_results_list[[dependent_var]] <- prot_results
}

all_results <- dplyr::bind_rows(all_results_list)
all_results$CI_lower <- all_results$beta_glm - 1.96 * all_results$se_glm
all_results$CI_upper <- all_results$beta_glm + 1.96 * all_results$se_glm
all_results$protein <- gsub("^rnt_", "", all_results$protein)
output_file <- paste0(results_dir, "tables/protein_cancer_incidence_reg_fully_adjusted_bmi_logistic.txt")

# Save results to a file
write.table(all_results, 
            file = output_file, 
            col.names = TRUE, row.names = FALSE, sep = "\t")

# Plot the forest plot using ggplot2
pdf(paste0(results_dir, "figures/forest_protein_incidence_reg_fully_adjusted_bmi_logistic.pdf"), height = 28, width = 8)
ggplot(all_results, aes(x = beta_glm, y = protein, color = cancer_type, shape = variable_type)) +
  geom_point(size = 3, position = position_dodge(width = 1)) +  # Add dodging to avoid overlap
  geom_errorbarh(aes(xmin = CI_lower, xmax = CI_upper), height = 0.2, position = position_dodge(width = 1)) +  # Apply dodging to error bars
  scale_color_manual(
    values = c("PC" = "blue", "EC" = "purple", "BC" = "red", "CRC" = "orange"),  # Colors for each cancer type
    labels = c("BC" = "Breast Cancer", "CRC" = "Colorectal Cancer", "EC" = "Endometrial Cancer", "PC" = "Prostate Cancer")  # Custom labels
  ) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "black") +
  labs(
    x = "Log odds of cancer per SD unit higher protein (±95%)", 
    y = "Proteins", 
    title = "Association between protein and cancer risk (incident)",
    color = "Cancer Site"  # Title of the color legend
  ) +
  theme_minimal() +
  theme(legend.position = "bottom") +  # Adjust legend position
  guides(shape = "none")  # Remove the legend for variable_type
dev.off()

## All proteins and incident cancer
protein_PA_name_list2 <- rnt_protein_names
all_results_list2 <- list()
dependent_variables2 <- c("PC_incidental_combined_control", "EC_incidental_combined_control", "BC_incidental_combined_control", "CRC_incidental_combined_control")

# Loop over the dependent variables (e.g., PC_incidental, EC_incidental, etc.)
for (dependent_var in dependent_variables2) {
  
  # Ensure the dependent variable is a factor (e.g., "PC", "EC", "BC", etc.)
  cancer_PA_protein[[dependent_var]] <- as.factor(cancer_PA_protein[[dependent_var]])
  
  # Run logistic models and add FDR-adjusted p-values
  prot_results <- logistic_model_with_fdr(
    wdata = cancer_PA_protein,
    dependent = dependent_var,
    proteins = protein_PA_name_list2,  # List of protein traits
    covariables = c("age", "sex", "smoking", "alcohol", "centre", "education")
  )
  
  # Sort results by raw p-values (adjusted p-values)
  prot_results <- prot_results[order(prot_results$P_FDR), ]
  
  # Add a column to indicate cancer type and status type
  cancer_type <- gsub("_incidental_combined_control", "", dependent_var)  # Extract cancer type from the dependent variable name
  prot_results$cancer_type <- cancer_type  # Add cancer type to results
  prot_results$variable_type <- "Incidental"  # Set status type as "Incidental"
  
  # Add the results for this cancer type to the list
  all_results_list[[dependent_var]] <- prot_results
}

all_results_proteins <- dplyr::bind_rows(all_results_list)
all_results_proteins$CI_lower <- all_results_proteins$beta_glm - 1.96 * all_results_proteins$se_glm
all_results_proteins$CI_upper <- all_results_proteins$beta_glm + 1.96 * all_results_proteins$se_glm
all_results_proteins$protein <- gsub("^rnt_", "", all_results_proteins$protein)
output_file2 <- paste0(results_dir, "tables/all_proteins_cancer_incidence_reg_fully_adjusted_logistic.txt")

# Save results to a file
write.table(all_results_proteins, 
            file = output_file2, 
            col.names = TRUE, row.names = FALSE, sep = "\t")

## Repeat regression and plot for prevalent cancer variable
# List of the cancer_status dependent variables (status variables)
# List of dependent variables for cancer status
dependent_variables2 <- c("PC_status_combined_control", "EC_status_combined_control", "BC_status_combined_control", "CRC_status_combined_control")

# Initialize an empty list to store results
all_results_list2 <- list()

# Loop over the dependent variables (e.g., PC_status, EC_status, etc.)
for (dependent_var in dependent_variables2) {
  
  # Ensure the dependent variable is a factor (e.g., "PC", "EC", "BC", etc.)
  cancer_PA_protein[[dependent_var]] <- as.factor(cancer_PA_protein[[dependent_var]])
  
  # Run logistic models and add FDR-adjusted p-values
  prot_results <- logistic_model_with_fdr(
    wdata = cancer_PA_protein,
    dependent = dependent_var,
    proteins = protein_PA_name_list,  # List of protein traits
    covariables = c("age", "sex", "smoking", "alcohol", "centre", "fasting", "bmi", "education")
  )
  
  # Sort results by raw p-values (adjusted p-values)
  prot_results <- prot_results[order(prot_results$P_FDR), ]
  
  # Add a column to indicate cancer type and status type
  cancer_type <- gsub("_status_combined_control", "", dependent_var)  # Extract cancer type from the dependent variable name
  prot_results$cancer_type <- cancer_type  # Add cancer type to results
  prot_results$variable_type <- "Status"  # Set status type as "Status"
  
  # Add the results for this cancer type to the list
  all_results_list2[[dependent_var]] <- prot_results
}

# Combine all the results into one data frame
all_results2 <- dplyr::bind_rows(all_results_list2)

# Calculate Confidence Intervals for the results
all_results2$CI_lower <- all_results2$beta_glm - 1.96 * all_results2$se_glm
all_results2$CI_upper <- all_results2$beta_glm + 1.96 * all_results2$se_glm

# Clean up the protein names
all_results2$protein <- gsub("^rnt_", "", all_results2$protein)

# Define the output file name
output_file2 <- paste0(results_dir, "tables/protein_cancer_status_reg_fully_adjusted_bmi_logistic.txt")

# Save results to a file
write.table(all_results2, 
            file = output_file2, 
            col.names = TRUE, row.names = FALSE, sep = "\t")

# Plot the forest plot using ggplot2
pdf(paste0(results_dir, "figures/forest_protein_prev_reg_fully_adjusted_bmi_logistic.pdf"), height = 28, width = 8)
ggplot(all_results2, aes(x = beta_glm, y = protein, color = cancer_type, shape = variable_type)) +
  geom_point(size = 3, position = position_dodge(width = 1)) +  # Add dodging to avoid overlap
  geom_errorbarh(aes(xmin = CI_lower, xmax = CI_upper), height = 0.2, position = position_dodge(width = 1)) +  # Apply dodging to error bars
  scale_color_manual(
    values = c("PC" = "blue", "EC" = "purple", "BC" = "red", "CRC" = "orange"),  # Colors for each cancer type
    labels = c("BC" = "Breast Cancer", "CRC" = "Colorectal Cancer", "EC" = "Endometrial Cancer", "PC" = "Prostate Cancer")  # Custom labels
  ) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "black") +
  labs(
    x = "Log odds of cancer per SD unit higher protein (±95%)", 
    y = "Proteins", 
    title = "Association between protein and cancer risk (prevalent)",
    color = "Cancer Site"  # Title of the color legend
  ) +
  theme_minimal() +
  theme(legend.position = "bottom") +  # Adjust legend position
  guides(shape = "none")  # Remove the legend for variable_type
dev.off()

# Add a new column to indicate whether the result is "Incidental" or "Status" 
all_results$variable_type <- "Incident"
all_results2$variable_type <- "Prevalent"

# Combine the results for both incidental and status variables
combined_results <- dplyr::bind_rows(all_results, all_results2)

# Now let's calculate the confidence intervals
combined_results$CI_lower <- combined_results$beta_glm - 1.96 * combined_results$se_glm
combined_results$CI_upper <- combined_results$beta_glm + 1.96 * combined_results$se_glm

# Clean up the protein names (removing the prefix 'rnt_')
combined_results$protein <- gsub("^rnt_", "", combined_results$protein)

# Plot the combined forest plot using ggplot2
# Number of proteins per page
proteins_per_page <- 20

# Calculate the number of pages needed
num_pages <- ceiling(n_distinct(combined_results$protein) / proteins_per_page)

# Save plots into a multi-page PDF
pdf(paste0(results_dir, "figures/forest_protein_cancer_combined_fully_adjusted_bmi_logistic.pdf"), height = 16, width = 8)

# Generate the paginated plots
for (page in 1:num_pages) {
  p <- ggplot(combined_results, aes(x = beta_glm, y = protein, color = cancer_type, shape = variable_type)) +
    geom_point(size = 3, position = position_dodge(width = 1)) +  # Add dodging to avoid overlap
    geom_errorbarh(aes(xmin = CI_lower, xmax = CI_upper), height = 0.2, position = position_dodge(width = 1)) +  # Apply dodging to error bars
    scale_color_manual(
      values = c("PC" = "blue", "EC" = "purple", "BC" = "red", "CRC" = "orange"),  # Colors for each cancer type
      labels = c("BC" = "Breast Cancer", "CRC" = "Colorectal Cancer", "EC" = "Endometrial Cancer", "PC" = "Prostate Cancer")  # Custom labels
    ) +
    scale_shape_manual(
      values = c("Prevalent" = 16, "Incident" = 17)  # Assign shapes: 16 = filled circle, 17 = triangle
    ) +
    geom_vline(xintercept = 0, linetype = "dotted", color = "black") +  # Dotted vertical line at x = 0
    labs(
      x = "Log odds of cancer risk (±95% CI) per SD higher protein", 
      y = "Proteins", 
      title = "",
      color = "Cancer Site",  # Title of the color legend
      shape = "Outcome Type"  # Title of the shape legend
    ) +
    theme_minimal() +
    theme(
      legend.position = "bottom",   # Place the legends at the bottom
      legend.box = "horizontal",    # Arrange the legends horizontally
      legend.box.margin = margin(10, 0, 0, 0),  # Add space around the legends
      legend.key.size = unit(1.5, "lines"),  # Adjust size of legend keys
      strip.text = element_blank(), # Remove the facet strip text (protein names above plots)
      strip.background = element_blank() # Remove strip background to clean up
    ) +
    guides(
      shape = guide_legend(override.aes = list(size = 4)),  # Adjust the size of the legend symbols for shapes
      color = guide_legend(ncol = 1)  # Make sure the cancer site legend is a single column
    ) +
    ggforce::facet_wrap_paginate(~protein, ncol = 1, nrow = proteins_per_page, page = page, scales = "free_y")
  
  print(p)
}

# Close the PDF device
dev.off()

## Do Cox PH for protein and cancer
get_event_date <- function(data, cancer) {
  icd10 <- paste0(cancer, "_ICD10_date")
  icd9  <- paste0(cancer, "_ICD9_date")
  
  as.Date(ifelse(
    !is.na(data[[icd10]]),
    data[[icd10]],
    data[[icd9]]
  ), origin = "1970-01-01")
}

# Administrative censor date
ADMIN_CENSOR_DATE <- as.Date("2020-03-31")

# Prepare protein dataset with entry age
cancer_PA_protein <- cancer_PA_protein %>%
  mutate(
    age_entry = age  # age at protein measurement
  )

# Cox PH model function with delayed entry cutoff
cox_model_cutoff <- function(data, cancer, protein, covariates, cutoff_years = 0) {
  
  # Compute event dates
  event_date <- get_event_date(data, cancer)
  
  df <- data %>%
    mutate(
      event_date = event_date,
      
      # Baseline entry date plus optional cutoff in years
      entry_date = as.Date(paste0(`Year of birth`, "-06-30")) +
        age_entry * 365.25 + years(cutoff_years)
    ) %>%
    
    # Exclude events before entry (delayed entry)
    filter(is.na(event_date) | event_date >= entry_date) %>%
    
    mutate(
      # Event indicator
      event = ifelse(!is.na(event_date) & event_date <= ADMIN_CENSOR_DATE, 1, 0),
      
      # Censoring date
      censor_date = ifelse(event == 1, event_date, ADMIN_CENSOR_DATE) %>%
        as.Date(origin = "1970-01-01"),
      
      # Age at exit (age timescale)
      age_exit = age_entry + as.numeric(censor_date - entry_date) / 365.25
    ) %>%
    
    # Keep only valid rows
    filter(
      age_exit > age_entry,
      !is.na(.data[[protein]]),
      complete.cases(across(all_of(covariates)))
    )
  
  # Skip if too few participants or events
  if (nrow(df) < 50 || sum(df$event) < 10) return(NULL)
  
  # Cox model formula
  form <- as.formula(
    paste0(
      "Surv(age_entry, age_exit, event) ~ ",
      protein, " + ",
      paste(covariates, collapse = " + ")
    )
  )
  
  fit <- coxph(form, data = df)
  
  # Return tidy results for the protein term
  broom::tidy(fit, conf.int = TRUE) %>%
    filter(term == protein) %>%
    mutate(
      protein = protein,
      cancer_type = cancer,
      N = nrow(df),
      events = sum(df$event),
      cutoff_years = cutoff_years
    )
}

# Covariates and cancers
covariates <- c("sex", "smoking", "alcohol", "centre", "bmi", "education")
cancer_types <- c("PC", "EC", "BC", "CRC")

# Run Cox models with a 2-year cutoff
cutoff_years <- 2
cox_results_cutoff <- list()

for (cancer in cancer_types) {
  for (prot in protein_PA_name_list) {
    res <- cox_model_cutoff(
      data = cancer_PA_protein,
      cancer = cancer,
      protein = prot,
      covariates = covariates,
      cutoff_years = cutoff_years
    )
    if (!is.null(res)) cox_results_cutoff[[length(cox_results_cutoff) + 1]] <- res
  }
}

cox_results_cutoff_df <- bind_rows(cox_results_cutoff)

## Add hazard ratios and 95% CI
cox_results_cutoff_df <- cox_results_cutoff_df %>%
  mutate(
    HR = exp(estimate),                 # hazard ratio
    HR_conf_low = exp(conf.low),        # lower 95% CI
    HR_conf_high = exp(conf.high)       # upper 95% CI
  )

write.table(
  cox_results_cutoff_df,
  file = paste0(results_dir,
                "tables/protein_cancer_incidence_cox_age_2years.txt"),
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE
)

## Redo PA-protein analysis in all controls
cancer_PA_protein_controls <- cancer_PA_protein %>%
  filter(combined_cancer_status == "control")

## Function to remove single level factors
remove_single_level_factors <- function(data) {
  factor_vars <- names(data)[sapply(data, is.factor)]  # Identify factor columns
  for (var in factor_vars) {
    if (length(unique(data[[var]])) < 2) {  # Check if only one level exists
      data[[var]] <- NULL  # Remove the factor column
    }
  }
  return(data)
}

### Run regression over all protein change cols
linear_model <- function(wdata, dependent, independent, covariables) {
  ##############################
  ### 1. Define Model Data Frame
  ##############################
  if (is.na(covariables[1])) {
    model_variables <- c(dependent, independent)
  } else {
    model_variables <- c(dependent, covariables, independent)
  }
  
  mod_data <- wdata[, model_variables, drop = FALSE]
  
  ##############################
  ### 2. Perform Linear Model 
  ##############################
  perform_lm <- function(data, dependent, independent, covariables) {
    # Remove single-level factors
    data <- remove_single_level_factors(data)
    
    # Ensure data is valid
    if (nrow(data) < 2 || length(unique(data[[dependent]])) < 2) {
      return(c(n_lm = NA, rsq_adj_lm = NA, beta_lm = NA, se_lm = NA, P_lm = NA))
    }
    
    # Construct formula
    if (is.na(covariables[1])) {
      form <- formula(paste0(dependent, " ~ ", independent))
    } else {
      form <- formula(paste0(dependent, " ~ ", independent, " + ", paste0(covariables, collapse = " + ")))
    }
    
    # Fit linear model
    lm_mod <- lm(form, data = data)
    
    # Extract summary stats
    s <- summary(lm_mod)
    n <- length(residuals(s)); names(n) <- "n_lm"
    rsq <- s$adj.r.squared; names(rsq) <- "rsq_adj_lm"
    beta <- s$coefficients[2, 1]; names(beta) <- "beta_lm"
    se <- s$coefficients[2, 2]; names(se) <- "se_lm"
    pval <- s$coefficients[2, 4]; names(pval) <- "P_lm"
    
    c(n, rsq, beta, se, pval)
  }
  
  # Perform Analysis with Updated Handling
  combined_results <- perform_lm(mod_data, dependent, independent, covariables)
  
  # Male Analysis
  male_data <- subset(mod_data, sex == "Male")
  male_data <- droplevels(male_data)
  male_data <- male_data[complete.cases(male_data), ]
  male_results <- if (nrow(male_data) > 0) {
    perform_lm(male_data, dependent, independent, covariables[!(covariables %in% "sex")])
  } else {
    c(n_lm = NA, rsq_adj_lm = NA, beta_lm = NA, se_lm = NA, P_lm = NA)
  }
  
  # Female Analysis
  female_data <- subset(mod_data, sex == "Female")
  female_data <- droplevels(female_data)
  female_data <- female_data[complete.cases(female_data), ]
  female_results <- if (nrow(female_data) > 0) {
    perform_lm(female_data, dependent, independent, covariables[!(covariables %in% "sex")])
  } else {
    c(n_lm = NA, rsq_adj_lm = NA, beta_lm = NA, se_lm = NA, P_lm = NA)
  }
  
  ##############################
  ### 4. Combine Results
  ##############################
  protein_name <- dependent; names(protein_name) <- "protein"
  
  results <- c(
    protein_name,
    combined_results,
    setNames(male_results, paste0("male_", names(male_results))),
    setNames(female_results, paste0("female_", names(female_results)))
  )
  
  return(results)
}

## Overall activity and protein in controls, fully adjusted
controls_overall_PA_protein3 = t( sapply(rnt_protein_names, function(trait){
  linear_model(wdata = cancer_PA_protein_controls, 
               dependent = trait, 
               independent = "overall_acceleration_average",
               covariables =  c("age", "sex", "smoking", "alcohol", "centre", "fasting", "education", "bmi"))
}) )
controls_overall_PA_protein3 <- as.data.frame(controls_overall_PA_protein3)
controls_overall_PA_protein3$P_lm <- as.numeric(controls_overall_PA_protein3$P_lm)
controls_overall_PA_protein3 <- controls_overall_PA_protein3[order(controls_overall_PA_protein3$P_lm), ]
write.table(controls_overall_PA_protein3, file = paste0(results_dir, "tables/controls_overall_PA_protein_regression_fully_adjusted_bmi.txt"), col.names = T, row.names = F, sep = "\t")

## Sedentary and protein in controls, fully adjusted 
controls_sedentary_PA_protein3 = t( sapply(rnt_protein_names, function(trait){
  linear_model(wdata = cancer_PA_protein_controls, 
               dependent = trait, 
               independent = "sedentary_overall_average",
               covariables =  c("age", "sex", "smoking", "alcohol", "centre", "fasting", "education", "bmi"))
}) )
controls_sedentary_PA_protein3 <- as.data.frame(controls_sedentary_PA_protein3)
controls_sedentary_PA_protein3$P_lm <- as.numeric(controls_sedentary_PA_protein3$P_lm)
controls_sedentary_PA_protein3 <- controls_sedentary_PA_protein3[order(controls_sedentary_PA_protein3$P_lm), ]
write.table(controls_sedentary_PA_protein3, file = paste0(results_dir, "tables/controls_sedentary_PA_protein_regression_full_adjusted_bmi.txt"), col.names = T, row.names = F, sep = "\t")

## Moderate_to_vigorous quantile and proteins in controls, fully adjusted
controls_mod_vig_PA_quant_numeric_protein3 = t( sapply(rnt_protein_names, function(trait){
  linear_model(wdata = cancer_PA_protein_controls, 
               dependent = trait, 
               independent = "moderate_to_vigorous_quantile",  # Changed to moderate_to_vigorous_quantile
               covariables =  c("age", "sex", "smoking", "alcohol", "centre", "fasting", "education", "bmi"))
}) )
controls_mod_vig_PA_quant_numeric_protein3 <- as.data.frame(controls_mod_vig_PA_quant_numeric_protein3)
controls_mod_vig_PA_quant_numeric_protein3$P_lm <- as.numeric(controls_mod_vig_PA_quant_numeric_protein3$P_lm)
controls_mod_vig_PA_quant_numeric_protein3 <- controls_mod_vig_PA_quant_numeric_protein3[order(controls_mod_vig_PA_quant_numeric_protein3$P_lm), ]
write.table(controls_mod_vig_PA_quant_numeric_protein3, file = paste0(results_dir, "tables/controls_mod_vig_PA_quantile_numeric_protein_regression_fully_adjusted_bmi.txt"), col.names = T, row.names = F, sep = "\t")

