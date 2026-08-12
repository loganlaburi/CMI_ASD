# load package
library(sjPlot)
library(dplyr)
library(broom)
library(mice)

# Read data (it only has T-scores for the clinical cut-offs) 
asd_data <- read.csv('final_data.csv')
outcome <- read.csv('outcome.csv')
outcome[outcome == "."] <- NA

# only extract the one that is ASD 
data <- merge(asd_data, outcome, by = "Identifiers", all.x = TRUE)
colSums(is.na(data))

write.csv(data, 'test_run.csv')

## Run the imputation
dat_imp <- data
vars_outcomes <- c(
  "SDS.SDS_DA_T","SDS.SDS_DIMS_T","SDS.SDS_DOES_T",
  "SDS.SDS_SBD_T", "SDS.SDS_SHY_T", "SDS.SDS_SWDT_T", "SDS.SDS_Total_T",
  "SCARED_P.SCARED_P_GD" ,"SCARED_P.SCARED_P_PN","SCARED_P.SCARED_P_SC",
  "SCARED_P.SCARED_P_SH", "SCARED_P.SCARED_P_SP","SCARED_P.SCARED_P_Total",
  "CBCL.CBCL_Int_T","CBCL.CBCL_Ext_T","CBCL.CBCL_Total_T"       
)

vars_covars <- c("Identifiers", "hhi", "puberty_stage", "age", "bmi_percentile", "sleep_hr",
                 "BMI", "sex", "ethnicity", "race", "comorbidity_count")

vars_all    <- c(vars_outcomes, vars_covars)

# Subset the data
dat_imp <- dat_imp[, vars_all, drop = FALSE]

dat_imp <- dat_imp %>%
  dplyr::mutate(
    across(-c("Identifiers"), ~as.numeric(as.character(.x)))
  )

# Convert categorical variables to factors
dat_imp$sex <- factor(dat_imp$sex)
dat_imp$ethnicity <- factor(dat_imp$ethnicity)
dat_imp$race <- factor(dat_imp$race)
dat_imp$puberty_stage <- factor(dat_imp$puberty_stage)
dat_imp$hhi <- factor(dat_imp$hhi)

# Ensure numeric variables are numeric
dat_imp$age <- as.numeric(dat_imp$age)
dat_imp$bmi_percentile <- as.numeric(dat_imp$bmi_percentile)
dat_imp$BMI <- as.numeric(dat_imp$BMI)
dat_imp$sleep_hr <- as.numeric(dat_imp$sleep_hr)
dat_imp$comorbidity_count <- as.numeric(dat_imp$comorbidity_count)

# create the method 
ini  <- mice::mice(dat_imp, maxit = 0)
meth <- ini$method
pred <- ini$predictorMatrix
meth[vars_outcomes] <- "pmm"

# Methods for covariates
meth["age"]      <- "pmm"
meth["sex"]      <- "logreg"    
meth["bmi_percentile"]      <- "pmm"
meth["hhi"]   <- if (is.ordered(dat_imp$income)) "polr" else "polyreg"
meth["puberty_stage"]  <- if (is.ordered(dat_imp$puberty_stage)) "polr" else "polyreg"
meth["ethnicity"]  <- if (is.ordered(dat_imp$ethnicity)) "polr" else "polyreg"
meth["race"]  <- if (is.ordered(dat_imp$race)) "polr" else "polyreg"
meth["BMI"]      <- "pmm"
meth["comorbidity_count"] <- "pmm"
meth["sleep_hr"] <- "pmm"

# Do NOT impute ID
meth["Identifiers"] <- ""

# Exclude subject_id as predictor & target
pred[, "Identifiers"] <- 0
pred["Identifiers", ] <- 0

# Run MICE 
imp <- mice::mice(
  dat_imp,
  m = 30,
  method = meth,
  predictorMatrix = pred,
  maxit = 50,
  seed = 20260108
)

# Extract completed dataset
data_comp1 <- mice::complete(imp, action = "all")
saveRDS(data_comp1, file = "cutoff.rds")

# ---- Read the data
data <- readRDS("cutoff.rds")
set.seed(123)

# Check if there is missing value 
sapply(data, function(d) sum(is.na(as.data.frame(d))))

# Calculate the descripitive statistics
mutate_vars <- function(df) {
  df <- as.data.frame(df)
  
  # categorize the age group
  df <- df %>%
    mutate(
      age_group = case_when(
        age >=  9 & age < 11 ~ "pre", # pre adolescent
        age >= 11 & age < 15 ~ "early", # early adolescent
        age >= 15 & age < 20 ~ "late", # middle-to-late adolescents,
        TRUE ~ NA_character_
      )
    )
  
  # BMI category using BMI percentile 
  df$bmi_category <- with(df, ifelse(
    is.na(bmi_percentile), NA,
    ifelse(bmi_percentile < 5, "Underweight",
           ifelse(bmi_percentile < 85, "Normal Weight",
                  ifelse(bmi_percentile < 95, "Overweight", "Obesity")
           )
    )
  ))
  
  # Reorder into this order
  df$bmi_category <- factor(
    df$bmi_category,
    levels = c("Underweight", "Normal Weight", "Overweight", "Obesity"),
    ordered = TRUE
  )
  
  df
}

df_imp <- lapply(data, mutate_vars)

# Create the clinical cut off classifications
class_sdsc <- function(x) {
  case_when(
    is.na(x)        ~ NA_character_,
    x >= 70         ~ "Clinical",
    x >= 50         ~ "Borderline",
    TRUE            ~ "Non-clinical"
  )
}

class_cbcl <- function(x) {
  case_when(
    is.na(x)        ~ NA_character_,
    x >= 64         ~ "Clinical",
    x >= 60         ~ "Borderline",
    TRUE            ~ "Non-clinical"
  )
}

class_scared <- function(x, cutoff) {
  case_when(
    is.na(x)        ~ NA_character_,
    x >= cutoff    ~ "Clinical",
    TRUE            ~ "Non-clinical"
  )
}

df_imp <- lapply(df_imp, function(d) {
  d <- as.data.frame(d)
  
  ## SDSC
  sdsc_vars <- c("SDS.SDS_DA_T","SDS.SDS_DIMS_T","SDS.SDS_DOES_T",
                 "SDS.SDS_SBD_T","SDS.SDS_SHY_T","SDS.SDS_SWDT_T", "SDS.SDS_Total_T")
  for (v in sdsc_vars[sdsc_vars %in% names(d)]) {
    d[[paste0(v, "_class")]] <- class_sdsc(d[[v]])
  }
  
  ## CBCL
  cbcl_vars <- c("CBCL.CBCL_Int_T","CBCL.CBCL_Ext_T","CBCL.CBCL_Total_T")
  for (v in cbcl_vars[cbcl_vars %in% names(d)]) {
    d[[paste0(v, "_class")]] <- class_cbcl(d[[v]])
  }
  
  ## SCARED
  if ("SCARED_P.SCARED_P_Total" %in% names(d))
    d$SCARED_P.SCARED_P_Total <- class_scared(d$SCARED_P.SCARED_P_Total, 25)
  if ("SCARED_P.SCARED_P_GD" %in% names(d))
    d$SCARED_P.SCARED_P_GD <- class_scared(d$SCARED_P.SCARED_P_GD, 9)
  if ("SCARED_P.SCARED_P_PN" %in% names(d))
    d$SCARED_P.SCARED_P_PN <- class_scared(d$SCARED_P.SCARED_P_PN, 7)
  if ("SCARED_P.SCARED_P_SC" %in% names(d))
    d$SCARED_P.SCARED_P_SC <- class_scared(d$SCARED_P.SCARED_P_SC, 8)
  if ("SCARED_P.SCARED_P_SH" %in% names(d))
    d$SCARED_P.SCARED_P_SH <- class_scared(d$SCARED_P.SCARED_P_SH, 3)
  if ("SCARED_P.SCARED_P_SP" %in% names(d))
    d$SCARED_P.SCARED_P_SP <- class_scared(d$SCARED_P.SCARED_P_SP, 5)
  
  d
})


get_prop <- function(var) {
  
  # Compute proportions within each imputation
  results <- bind_rows(lapply(seq_along(df_imp), function(i) {
    df_imp[[i]] %>%
      filter(!is.na(age_group),
             !is.na(.data[[var]])) %>%
      count(age_group, Category = .data[[var]]) %>%
      group_by(age_group) %>%
      mutate(
        total_n = sum(n),
        prop = n / total_n,
        .imp = i
      )
  }))
  
  # Pool across imputations
  results <- results %>%
    group_by(age_group, Category) %>%
    summarise(
      prop = mean(prop, na.rm = TRUE),
      N = round(mean(total_n, na.rm = TRUE)),
      .groups = "drop"
    )
  
  # Compute counts that SUM EXACTLY to N
  results <- results %>%
    group_by(age_group) %>%
    mutate(
      n_raw = prop * N,
      n = floor(n_raw)   
    ) %>%
    mutate(
      remainder = N - sum(n)   
    ) %>%
    mutate(
      rank = rank(-n_raw, ties.method = "first"),
      n = n + ifelse(rank <= remainder, 1, 0)  
    ) %>%
    ungroup() %>%
    select(age_group, Category, prop, N, n)
  
  return(results)
}

# Extract the result
a <- data.frame(get_prop("SCARED_P.SCARED_P_Total"))
