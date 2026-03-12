### Physical activity and protein analysis ###
### Main analysis ###
### Lucy Goudswaard - 23rd July 2024 ###

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
library(glue)
library(MASS)
library(corrplot)

## clear environment
rm(list=ls())

## connect to project and setwd to script folder

# read in parameter file (specified on command line)
source("parameter_file/parameters_for_r.R")

## Read in physical activity and covar data, as well as olink dataframe
PA_data <- fread(paste0(data_intermediate_dir, "PA_protein_data_with_MET.csv"), header = T, sep = ",")
PA_data <- as.data.frame(PA_data)
olink <- fread(paste0(input_dir, "olink_instance_0.csv"))
olink <- as.data.frame(olink)

## Remove withdrawals from olink df (already removed from PA and covar data)
withdrawals <- fread(paste0(input_dir, "/16391_from_Becky/withdrawals/w16391_20241217.csv"))
colnames(withdrawals) <- "withdrawn"
olink_withdraw <- which(olink$eid %in% withdrawals$withdrawn)
olink <- olink[-olink_withdraw,]

## Read in protein names
protein_names <- read.table(paste0(input_dir, "field_names_instance_0.txt"), header = T, sep = "\t")
protein_names <- protein_names[[1]]

## Merge PA, covars and olink
#data <- merge(PA_data, olink, by.x = "Participant ID", by.y = "eid", all = TRUE) 
#data <- as.data.frame(data)

## protein name file - check this can be used to match.
## is there a specific ukb file?
olink_names <- read_xlsx(paste0(input_dir, "olink-explore-3072-assay-list-2023-06-08.xlsx"))

## transform olink data
rnt_data <- as.data.frame(matrix(nrow = nrow(olink), ncol = length(protein_names)+1))
for (i in 1:length(protein_names)) {
colnames(rnt_data)[1] <- "eid"
  colnames(rnt_data)[i+1] <- paste0("rnt_", protein_names[i])
  rnt_data[,i+1] <- rntransform(olink[, protein_names[i]])  
}
rnt_protein_names <- colnames(rnt_data)[2:ncol(rnt_data)]
rnt_data$eid <- olink$eid
rnt_data$eid <- as.integer(rnt_data$eid)

## bind this back to PA data, for those who just have protein levels
data_extended <- merge(PA_data, rnt_data, by.x = "Participant ID", by.y = "eid", all = F)
data_extended <- as.data.frame(data_extended)

## PA variables normality
hist(data_extended$overall_acceleration_average)
sampled_data_overall <- sample(data_extended$overall_acceleration_average, 5000)
shapiro.test(sampled_data_overall)
hist(data_extended$sedentary_overall_average)
sampled_data_sedentary <- sample(data_extended$sedentary_overall_average, 5000)
shapiro.test(sampled_data_sedentary)
hist(data_extended$moderate_to_vigorous_overall_average)
sampled_data_mod_vig <- sample(data_extended$moderate_to_vigorous_overall_average, 5000)
shapiro.test(sampled_data_mod_vig)

## Remove outliers more than 5 SD
PA_cols_QC <- c("overall_acceleration_average", "sedentary_overall_average", 
                      "moderate_to_vigorous_overall_average", "MET_mins_week_overall",
                "MET_mins_week_overall_1", "MET_mins_week_overall_2", "MET_mins_week_overall_3")

# Change outliers to NA in PA data (more than 5 SD away from the mean)
means <- colMeans(data_extended[PA_cols_QC], na.rm = TRUE)
sds <- apply(data_extended[PA_cols_QC], 2, sd, na.rm = TRUE)
data_extended[PA_cols_QC] <- apply(data_extended[PA_cols_QC], 2, function(column) {
  outliers <- abs(column - mean(column, na.rm = TRUE)) > (5 * sd(column, na.rm = TRUE))
  column[outliers] <- NA
  return(column)
})

## Derive MET in hours/day (currently it is MET minutes per week)
data_extended$MET_h_day <- (data_extended$MET_mins_week_overall / 7 ) / 60

## Check outliers have been removed - moderate to vigorous and MET have a negative skew - check transformation is appropriate
pdf(file = paste0(results_dir, "figures/hist_PA_outliers_excluded.pdf"))
par(mfrow = c(2, 2))
columns_to_plot <- c("overall_acceleration_average", "sedentary_overall_average", 
                     "moderate_to_vigorous_overall_average", "MET_h_day")
for (col in columns_to_plot) {
  hist(data_extended[[col]], 
       main = col, 
       xlab = col, 
       col = "lightblue", 
       border = "black")
}
dev.off()

## How many of Met_h_day and moderate_to_vigorous_overall are 0?
zero <- sum(data_extended$moderate_to_vigorous_overall_average == 0, na.rm=T) ## 577
zero_2 <- sum(data_extended$MET_h_day == 0, na.rm=T) ## 939

## Recode covariable data (smoking and alcohol)
names(data_extended)[names(data_extended) == "Age at recruitment"] <- "age"
names(data_extended)[names(data_extended) == "Sex"] <- "sex"
names(data_extended)[names(data_extended) == "Alcohol drinker status | Instance 0"] <- "alcohol"
data_extended$alcohol <- ifelse(data_extended$alcohol == "Current", 2,
                                ifelse(data_extended$alcohol == "Previous", 1,
                                       ifelse(data_extended$alcohol == "Never", 0, NA)))
names(data_extended)[names(data_extended) == "Ever smoked | Instance 0"] <- "smoking"
data_extended$smoking <- ifelse(data_extended$smoking == "Yes", 1,
                                ifelse(data_extended$smoking == "No", 0, NA))
names(data_extended)[names(data_extended) == "Townsend deprivation index at recruitment"] <- "townsend"
names(data_extended)[names(data_extended) == "Body mass index (BMI) | Instance 0"] <- "bmi"
names(data_extended)[names(data_extended) == "UK Biobank assessment centre | Instance 0"] <- "centre"
names(data_extended)[names(data_extended) == "Fasting time | Instance 0"] <- "fasting"
names(data_extended)[names(data_extended) == "Qualifications | Instance 0"] <- "education_raw"

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

data_extended$education <- sapply(data_extended$education_raw, function(x) {
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
data_extended$education <- as.factor(data_extended$education)

## Generate table of participant characteristics - just for those with accelerometry and protein
### Categorical table
# Create cat_characteristics by filtering rows where any of the specified columns is not NA
cat_characteristics <- data_extended %>%
  filter(
    !is.na(overall_acceleration_average) |
      !is.na(sedentary_overall_average) |
      !is.na(moderate_to_vigorous_overall_average)
  )
cat_characteristics <- cat_characteristics[,c("age", "sex", "townsend", "education", "bmi", "smoking", "alcohol", "centre", "moderate_to_vigorous_overall_average",
                                              "sedentary_overall_average", "overall_acceleration_average") ]
cat_characteristics <- cat_characteristics %>%
  mutate(
    age = as.numeric(age),
    bmi = as.numeric(bmi),
    sex = as.factor(sex),
    smoking = as.factor(smoking),
    alcohol = as.factor(alcohol),
    centre = as.factor(centre),
    education = as.factor(education)
  )

cat_characteristics <- cat_characteristics %>%
  mutate(
    smoking = factor(smoking, levels = c("0", "1")),
    alcohol = factor(alcohol, levels = c("0", "1", "2"))
  )

# Generate the summary table with updated code
table1 <- tbl_summary(
  cat_characteristics,
  statistic = list(
    all_continuous() ~ "{mean} ({sd})",
    all_categorical() ~ "{n} / {N} ({p}%)"
  ),
  missing = "no"
)

table1 %>%
  as_flex_table() %>%
  flextable::save_as_docx(path = paste0(results_dir, "tables/table1_UKB_characteristics.docx"))

## Correlation matrix across PA variables - those with protein data first
# Compute the correlation matrix
columns_of_interest <- c("overall_acceleration_average", 
                         "sedentary_overall_average", 
                         "moderate_to_vigorous_overall_average", 
                         "MET_mins_week_overall",
                         'MET_mins_week_overall_1',
                         'MET_mins_week_overall_2',
                         'MET_mins_week_overall_3')

data_subset <- PA_data[, columns_of_interest]
correlation_matrix <- cor(data_subset, use = "complete.obs", method = "pearson")

# Function to calculate the number of complete cases for each pair
calculate_n <- function(x) {
  outer(1:ncol(x), 1:ncol(x), Vectorize(function(i, j) {
    sum(complete.cases(x[, c(i, j)]))
  }))
}

# Calculate N for each pair of variables
n_matrix <- calculate_n(data_subset)

# Assign row and column names
rownames(n_matrix) <- colnames(data_subset)
colnames(n_matrix) <- colnames(data_subset)

# Display the N matrix
print(n_matrix)
write.table(n_matrix, file = paste0(results_dir, "tables/corr_matrix_PA_N.txt"), row.names = T, col.names = T, sep = "\t")

# Create the correlation plot
pdf(file = paste0(results_dir, "figures/corr_matrix_PA.pdf"), height =10, width = 10)
corrplot(correlation_matrix, method = "color", 
         type = "upper",  # Show only upper triangle
         addCoef.col = "black", # Add correlation coefficients
         tl.col = "black", tl.srt = 45, # Rotate labels
         col = colorRampPalette(c("blue", "white", "red"))(200))
dev.off()

## Categorise overall and MET variables as these were skewed so lm may not be most appropriate
data_extended <- data_extended %>%
  mutate(moderate_to_vigorous_quantile = ntile(moderate_to_vigorous_overall_average, 4) )
data_extended <- data_extended %>%
  mutate(MET_quantile = ntile(MET_h_day, 4) )

## Transform the PA data
PA_vars <- c("overall_acceleration_average", "sedentary_overall_average", "moderate_to_vigorous_overall_average", "MET_mins_week_overall")
for (i in 1:length(PA_vars)) {
  data_extended[, PA_vars[i]] <- rntransform(data_extended[, PA_vars[i]])  
}

## function to perform regression for each PA measure and each protein outcome
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


## Write this to file
write.table(data_extended, file = paste0(data_intermediate_dir, "PA_rnt_protein_data_with_MET2.txt"), row.names = F, col.names = T, sep = "\t")

## Restrict accelerometry dataset to complete case dataset
accelerometry_dataset <- data_extended %>%
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
         !is.na(moderate_to_vigorous_overall_average) |
         !is.na(moderate_to_vigorous_quantile) |
         !is.na(sedentary_overall_average))
  )
accelerometry_dataset$sex <- as.factor(accelerometry_dataset$sex)

## Regression to estimate association between overall physical activity and proteins
## Model 1
overall_PA_protein = t( sapply(rnt_protein_names, function(trait){
  linear_model(wdata = accelerometry_dataset, 
               dependent = trait, 
               independent = "overall_acceleration_average",
               covariables =  c("age", "sex"))
}) )
overall_PA_protein <- as.data.frame(overall_PA_protein)
overall_PA_protein$P_lm <- as.numeric(overall_PA_protein$P_lm)
overall_PA_protein <- overall_PA_protein[order(overall_PA_protein$P_lm), ]
write.table(overall_PA_protein, file = paste0(results_dir, "tables/overall_PA_protein_regression_age_sex.txt"), col.names = T, row.names = F, sep = "\t")

## Model 2
overall_PA_protein2 = t( sapply(rnt_protein_names, function(trait){
  linear_model(wdata = accelerometry_dataset, 
               dependent = trait, 
               independent = "overall_acceleration_average",
               covariables =  c("age", "sex", "smoking", "alcohol", "centre", "fasting", "education" ))
}) )
overall_PA_protein2 <- as.data.frame(overall_PA_protein2)
overall_PA_protein2$P_lm <- as.numeric(overall_PA_protein2$P_lm)
overall_PA_protein2 <- overall_PA_protein2[order(overall_PA_protein2$P_lm), ]
write.table(overall_PA_protein2, file = paste0(results_dir, "tables/overall_PA_protein_regression_fully_adjusted.txt"), col.names = T, row.names = F, sep = "\t")

## Model 3
overall_PA_protein3 = t( sapply(rnt_protein_names, function(trait){
  linear_model(wdata = accelerometry_dataset, 
               dependent = trait, 
               independent = "overall_acceleration_average",
               covariables =  c("age", "sex", "smoking", "alcohol", "centre", "fasting", "education", "bmi"))
}) )
overall_PA_protein3 <- as.data.frame(overall_PA_protein3)
overall_PA_protein3$P_lm <- as.numeric(overall_PA_protein3$P_lm)
overall_PA_protein3 <- overall_PA_protein3[order(overall_PA_protein3$P_lm), ]
write.table(overall_PA_protein3, file = paste0(results_dir, "tables/overall_PA_protein_regression_fully_adjusted_bmi.txt"), col.names = T, row.names = F, sep = "\t")

## Regression for mod-vigorous activity and proteins
## Model 1
mod_vig_PA_protein = t( sapply(rnt_protein_names, function(trait){
linear_model(wdata = accelerometry_dataset, 
             dependent = trait, 
             independent = "moderate_to_vigorous_overall_average",
             covariables =  c("age", "sex"))
}) )
mod_vig_PA_protein <- as.data.frame(mod_vig_PA_protein)
mod_vig_PA_protein$P_lm <- as.numeric(mod_vig_PA_protein$P_lm)
mod_vig_PA_protein <- mod_vig_PA_protein[order(mod_vig_PA_protein$P_lm), ]
write.table(mod_vig_PA_protein, file = paste0(results_dir, "tables/mod_vig_PA_protein_regression_age_sex.txt"), col.names = T, row.names = F, sep = "\t")

## Model 2
mod_vig_PA_protein2 = t( sapply(rnt_protein_names, function(trait){
  linear_model(wdata = accelerometry_dataset, 
               dependent = trait, 
               independent = "moderate_to_vigorous_overall_average",
               covariables =  c("age", "sex", "smoking", "alcohol", "centre", "fasting", "education"))
}) )
mod_vig_PA_protein2 <- as.data.frame(mod_vig_PA_protein2)
mod_vig_PA_protein2$P_lm <- as.numeric(mod_vig_PA_protein2$P_lm)
mod_vig_PA_protein2 <- mod_vig_PA_protein2[order(mod_vig_PA_protein2$P_lm), ]
write.table(mod_vig_PA_protein2, file = paste0(results_dir, "tables/mod_vig_PA_protein_regression_fully_adjusted.txt"), col.names = T, row.names = F, sep = "\t")

## Model 3
mod_vig_PA_protein3 = t( sapply(rnt_protein_names, function(trait){
  linear_model(wdata = accelerometry_dataset, 
               dependent = trait, 
               independent = "moderate_to_vigorous_overall_average",
               covariables =  c("age", "sex", "smoking", "alcohol", "centre", "fasting", "education", "bmi"))
}) )
mod_vig_PA_protein3 <- as.data.frame(mod_vig_PA_protein3)
mod_vig_PA_protein3$P_lm <- as.numeric(mod_vig_PA_protein3$P_lm)
mod_vig_PA_protein3 <- mod_vig_PA_protein3[order(mod_vig_PA_protein3$P_lm), ]
write.table(mod_vig_PA_protein3, file = paste0(results_dir, "tables/mod_vig_PA_protein_regression_fully_adjusted_bmi.txt"), col.names = T, row.names = F, sep = "\t")

## Regression sedentary average and proteins
## Model 1
sedentary_PA_protein = t( sapply(rnt_protein_names, function(trait){
  linear_model(wdata = accelerometry_dataset, 
               dependent = trait, 
               independent = "sedentary_overall_average",
               covariables =  c("age", "sex"))
}) )
sedentary_PA_protein <- as.data.frame(sedentary_PA_protein)
sedentary_PA_protein$P_lm <- as.numeric(sedentary_PA_protein$P_lm)
sedentary_PA_protein <- sedentary_PA_protein[order(sedentary_PA_protein$P_lm), ]
write.table(sedentary_PA_protein, file = paste0(results_dir, "tables/sedentary_PA_protein_regression.txt"), col.names = T, row.names = F, sep = "\t")

## Model 2
sedentary_PA_protein2 = t( sapply(rnt_protein_names, function(trait){
  linear_model(wdata = accelerometry_dataset, 
               dependent = trait, 
               independent = "sedentary_overall_average",
               covariables =  c("age", "sex", "smoking", "alcohol", "centre", "fasting", "education"))
}) )
sedentary_PA_protein2 <- as.data.frame(sedentary_PA_protein2)
sedentary_PA_protein2$P_lm <- as.numeric(sedentary_PA_protein2$P_lm)
sedentary_PA_protein2 <- sedentary_PA_protein2[order(sedentary_PA_protein2$P_lm), ]
write.table(sedentary_PA_protein2, file = paste0(results_dir, "tables/sedentary_PA_protein_regression_full_adjusted.txt"), col.names = T, row.names = F, sep = "\t")

## Model 3
sedentary_PA_protein3 = t( sapply(rnt_protein_names, function(trait){
  linear_model(wdata = accelerometry_dataset, 
               dependent = trait, 
               independent = "sedentary_overall_average",
               covariables =  c("age", "sex", "smoking", "alcohol", "centre", "fasting", "education", "bmi"))
}) )
sedentary_PA_protein3 <- as.data.frame(sedentary_PA_protein3)
sedentary_PA_protein3$P_lm <- as.numeric(sedentary_PA_protein3$P_lm)
sedentary_PA_protein3 <- sedentary_PA_protein3[order(sedentary_PA_protein3$P_lm), ]
write.table(sedentary_PA_protein3, file = paste0(results_dir, "tables/sedentary_PA_protein_regression_full_adjusted_bmi.txt"), col.names = T, row.names = F, sep = "\t")

## Create a new linear model for quantiles
# Adapted function for reporting estimates for quantiles 2, 3, and 4
linear_model_quantile <- function(wdata, dependent, independent, covariables) {
  ##############################
  ### 1. Define Model Data Frame
  ##############################
  if(is.na(covariables[1])){
    model_variables = c(dependent, independent)
  } else {
    model_variables = c(dependent, covariables, independent)
  }
  
  mod_data = wdata[, c(model_variables)]
  
  ##############################
  ### 2. Perform Linear Model 
  ##############################
  if(is.na(covariables[1])){
    form = formula(paste0(dependent, " ~ ", independent))  
  } else {
    form = formula(paste0(dependent, " ~ ", independent, " + ", paste0(covariables, collapse = " + ")))  
  }
  
  lm_mod <- lm(form, data = mod_data)
  
  ##############################
  ### 3. Extract summary stats
  ##############################
  protein_name = dependent; names(protein_name) = "protein"
  s = summary(lm_mod)
  
  ## Sample size
  n = length(residuals(s)); names(n) = "n_lm"
  
  ## Adjusted R-squared
  rsq = s$adj.r.squared; names(rsq) = "rsq_adj_lm"
  
  ## Dependent Variable effect estimates (for all quantiles)
  # Extract estimates for each quantile
  beta2 = s$coefficients["moderate_to_vigorous_quantile2", 1] # Quantile 2 estimate
  beta3 = s$coefficients["moderate_to_vigorous_quantile3", 1] # Quantile 3 estimate
  beta4 = s$coefficients["moderate_to_vigorous_quantile4", 1] # Quantile 4 estimate
  
  # Extract standard errors for each quantile
  se2 = s$coefficients["moderate_to_vigorous_quantile2", 2]
  se3 = s$coefficients["moderate_to_vigorous_quantile3", 2]
  se4 = s$coefficients["moderate_to_vigorous_quantile4", 2]
  
  # Extract p-values for each quantile
  pval2 = s$coefficients["moderate_to_vigorous_quantile2", 4]
  pval3 = s$coefficients["moderate_to_vigorous_quantile3", 4]
  pval4 = s$coefficients["moderate_to_vigorous_quantile4", 4]
  
  #######################
  ## Return results
  #######################
  lm_results = c(protein_name, n, rsq, 
                 beta2, se2, pval2, 
                 beta3, se3, pval3, 
                 beta4, se4, pval4)
  
  # Create a named data frame for the results
  lm_results_df <- data.frame(
    protein = protein_name,
    n_lm = n,
    rsq_adj_lm = rsq,
    beta2 = beta2, se2 = se2, pval2 = pval2,
    beta3 = beta3, se3 = se3, pval3 = pval3,
    beta4 = beta4, se4 = se4, pval4 = pval4,
    stringsAsFactors = FALSE
  )
  
  return(lm_results_df)
}


# Perform ordinal regression (polr) for each trait in rnt_protein_names with moderate to vigorous quantile var
# Model 1
accelerometry_dataset$moderate_to_vigorous_quantile <- as.factor(accelerometry_dataset$moderate_to_vigorous_quantile)
mod_vig_PA_quant_protein = t( sapply(rnt_protein_names, function(trait){
  linear_model_quantile(wdata = accelerometry_dataset, 
               dependent = trait, 
               independent = "moderate_to_vigorous_quantile",
               covariables =  c("age", "sex"))
}) )
mod_vig_PA_quant_protein <- as.data.frame(mod_vig_PA_quant_protein)
mod_vig_PA_quant_protein[] <- lapply(mod_vig_PA_quant_protein, function(x) {
  if (is.list(x)) {
    return(sapply(x, toString))  # Convert list elements to strings
  } else {
    return(x)  # Leave other columns unchanged
  }
})
mod_vig_PA_quant_protein$pval2 <- as.numeric(mod_vig_PA_quant_protein$pval2)
mod_vig_PA_quant_protein <- mod_vig_PA_quant_protein[order(mod_vig_PA_quant_protein$pval2), ]
write.table(mod_vig_PA_quant_protein, file = paste0(results_dir, "tables/mod_vig_PA_quant_protein_regression_age_sex.txt"), col.names = T, row.names = F, sep = "\t")

# Model 2
mod_vig_PA_quant_protein2 = t( sapply(rnt_protein_names, function(trait){
  linear_model_quantile(wdata = accelerometry_dataset, 
               dependent = trait, 
               independent = "moderate_to_vigorous_quantile",
               covariables =  c("age", "sex", "smoking", "alcohol", "centre", "fasting", "education"))
}) )
mod_vig_PA_quant_protein2 <- as.data.frame(mod_vig_PA_quant_protein2)
mod_vig_PA_quant_protein2$pval2 <- as.numeric(mod_vig_PA_quant_protein2$pval2)
mod_vig_PA_quant_protein2 <- mod_vig_PA_quant_protein2[order(mod_vig_PA_quant_protein2$pval2), ]
mod_vig_PA_quant_protein2[] <- lapply(mod_vig_PA_quant_protein2, function(x) {
  if (is.list(x)) {
    return(sapply(x, toString))  # Convert list elements to strings
  } else {
    return(x)  # Leave other columns unchanged
  }
})
write.table(mod_vig_PA_quant_protein2, file = paste0(results_dir, "tables/mod_vig_PA_quant_protein_regression_fully_adjusted.txt"), col.names = T, row.names = F, sep = "\t")

# Model 3
mod_vig_PA_quant_protein3 = t( sapply(rnt_protein_names, function(trait){
  linear_model_quantile(wdata = accelerometry_dataset, 
               dependent = trait, 
               independent = "moderate_to_vigorous_quantile",
               covariables =  c("age", "sex", "smoking", "alcohol", "centre", "fasting", "education", "bmi"))
}) )
mod_vig_PA_quant_protein3 <- as.data.frame(mod_vig_PA_quant_protein3)
mod_vig_PA_quant_protein3$pval2 <- as.numeric(mod_vig_PA_quant_protein3$pval2)
mod_vig_PA_quant_protein3 <- mod_vig_PA_quant_protein3[order(mod_vig_PA_quant_protein3$pval2), ]
mod_vig_PA_quant_protein3[] <- lapply(mod_vig_PA_quant_protein3, function(x) {
  if (is.list(x)) {
    return(sapply(x, toString))  # Convert list elements to strings
  } else {
    return(x)  # Leave other columns unchanged
  }
})
write.table(mod_vig_PA_quant_protein3, file = paste0(results_dir, "tables/mod_vig_PA_quant_protein_regression_fully_adjusted_bmi.txt"), col.names = T, row.names = F, sep = "\t")

## Regression for mod-vigorous activity and proteins
## Rerun as above but have PA quantiles as numeric
accelerometry_dataset$moderate_to_vigorous_quantile <- as.numeric(accelerometry_dataset$moderate_to_vigorous_quantile)

# Model 1: Linear regression with moderate_to_vigorous_quantile and age, sex as covariates
mod_vig_PA_quant_numeric_protein = t( sapply(rnt_protein_names, function(trait){
  linear_model(wdata = accelerometry_dataset, 
               dependent = trait, 
               independent = "moderate_to_vigorous_quantile",  # Changed to moderate_to_vigorous_quantile
               covariables =  c("age", "sex"))
}) )
mod_vig_PA_quant_numeric_protein <- as.data.frame(mod_vig_PA_quant_numeric_protein)
mod_vig_PA_quant_numeric_protein$P_lm <- as.numeric(mod_vig_PA_quant_numeric_protein$P_lm)
mod_vig_PA_quant_numeric_protein <- mod_vig_PA_quant_numeric_protein[order(mod_vig_PA_quant_numeric_protein$P_lm), ]
write.table(mod_vig_PA_quant_numeric_protein, file = paste0(results_dir, "tables/mod_vig_PA_quantile_numeric_protein_regression_age_sex.txt"), col.names = T, row.names = F, sep = "\t")

# Model 2: Linear regression with moderate_to_vigorous_quantile and age, sex, smoking, alcohol, centre, fasting as covariates
mod_vig_PA_quant_numeric_protein2 = t( sapply(rnt_protein_names, function(trait){
  linear_model(wdata = accelerometry_dataset, 
               dependent = trait, 
               independent = "moderate_to_vigorous_quantile",  # Changed to moderate_to_vigorous_quantile
               covariables =  c("age", "sex", "smoking", "alcohol", "centre", "fasting", "education"))
}) )
mod_vig_PA_quant_numeric_protein2 <- as.data.frame(mod_vig_PA_quant_numeric_protein2)
mod_vig_PA_quant_numeric_protein2$P_lm <- as.numeric(mod_vig_PA_quant_numeric_protein2$P_lm)
mod_vig_PA_quant_numeric_protein2 <- mod_vig_PA_quant_numeric_protein2[order(mod_vig_PA_quant_numeric_protein2$P_lm), ]
write.table(mod_vig_PA_quant_numeric_protein2, file = paste0(results_dir, "tables/mod_vig_PA_quantile_numeric_protein_regression_fully_adjusted.txt"), col.names = T, row.names = F, sep = "\t")

# Model 3: Linear regression with moderate_to_vigorous_quantile and age, sex, smoking, alcohol, centre, fasting, bmi as covariates
mod_vig_PA_quant_numeric_protein3 = t( sapply(rnt_protein_names, function(trait){
  linear_model(wdata = accelerometry_dataset, 
               dependent = trait, 
               independent = "moderate_to_vigorous_quantile",  # Changed to moderate_to_vigorous_quantile
               covariables =  c("age", "sex", "smoking", "alcohol", "centre", "fasting", "education", "bmi"))
}) )
mod_vig_PA_quant_numeric_protein3 <- as.data.frame(mod_vig_PA_quant_numeric_protein3)
mod_vig_PA_quant_numeric_protein3$P_lm <- as.numeric(mod_vig_PA_quant_numeric_protein3$P_lm)
mod_vig_PA_quant_numeric_protein3 <- mod_vig_PA_quant_numeric_protein3[order(mod_vig_PA_quant_numeric_protein3$P_lm), ]
write.table(mod_vig_PA_quant_numeric_protein3, file = paste0(results_dir, "tables/mod_vig_PA_quantile_numeric_protein_regression_fully_adjusted_bmi.txt"), col.names = T, row.names = F, sep = "\t")

