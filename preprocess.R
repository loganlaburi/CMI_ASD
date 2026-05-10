
# load package
library(sjPlot)
library(dplyr)
library(broom)
library(mice) 

# Read data
asd_data <- read.csv('data_asd.csv')
diag <- read.csv('diagnosis.csv')
sleep <- read.csv('sleep_hr.csv')
diag <- diag[, c("Identifiers", "Has_Autism")]
cbcl <- read.csv('cbcl.csv')

# Attach the "has_autism" info based on dplyr
data <- asd_data %>%
  left_join(diag, by="Identifiers")

raw_data <- subset(data, data$Has_Autism==1)

# Combine all files together
merged_data <- raw_data %>%
  left_join(cbcl, by="Identifiers")

df <- merged_data[, c("Identifiers", "bmi_percentile", "FSQ.FSQ_04",
               "PPS.PPS_F_Score", "PPS.PPS_M_Score", 
               "Physical.Age", "Physical.BMI", "Physical.Sex",
               "PreInt_Demos_Fam.Child_Ethnicity",
               "PreInt_Demos_Fam.Child_Race",
               "SDS.SDS_DA_T", "SDS.SDS_DIMS_T",
               "SDS.SDS_DOES_T", "SDS.SDS_SBD_T",
               "SDS.SDS_SHY_T", "SDS.SDS_SWDT_T", "SDS.SDS_Total_T",
               "SCARED_P.SCARED_P_GD", "SCARED_P.SCARED_P_PN",
               "SCARED_P.SCARED_P_SC", "SCARED_P.SCARED_P_SH",
               "SCARED_P.SCARED_P_SP", "SCARED_P.SCARED_P_Total",
               "CBCL.CBCL_Int_T", "CBCL.CBCL_Ext_T", 
               "CBCL.CBCL_Total_T")]

# Change the column names
colnames(df) <- c("Identifiers", "BMI_perc", "hhi",
                  "pps_f_total", "pps_m_total", "age", "BMI", "sex",
                  "ethnicity", "race", 
                  "SDS_DA_T", "SDS_DIMS_T", "SDS_DOES_T",
                  "SDS_SBD_T", "SDS_SHY_T", "SDS_SWDT_T", "SDS_Total_T",
                  "SCARED_P_GD", "SCARED_P_PN", "SCARED_P_SC", 
                  "SCARED_P_SH","SCARED_P_SP", "SCARED_P_Total",
                  "CBCL_Int_T", "CBCL_Ext_T", "CBCL_Total_T")

# Put the puberty score into the categorical character, and convert 
# back to the numeric 
df <- df %>%
  mutate(puberty_stage = case_when(
    sex == "M" & pps_m_total == 3 ~ 1,                      # Prepubertal
    sex == "M" & pps_m_total %in% c(4, 5) ~ 2,              # Beginning Pubertal
    sex == "M" & pps_m_total %in% c(6, 7, 8) ~ 3,           # Mid-Pubertal
    sex == "M" & pps_m_total %in% c(9, 10, 11) ~ 4,         # Advanced Pubertal
    sex == "M" & pps_m_total == 12 ~ 5,
    sex == "F" & pps_f_total == 1 ~ 1,                      # Prepubertal
    sex == "F" & pps_f_total %in% 2 ~ 2,              # Beginning Pubertal
    sex == "F" & pps_f_total %in% 3 ~ 3,           # Mid-Pubertal
    sex == "F" & pps_f_total %in% 4 ~ 4,         # Advanced Pubertal
    sex == "F" & pps_f_total == 5 ~ 5 # Postpubertal
  ))

# convert everything except sex to numeric
df <- df %>%
  mutate(across(.cols = -c(sex, Identifiers), ~as.numeric(as.character(.x))))

# Add sleep information
dfs <- left_join(df, sleep, by = "Identifiers")

# sex coded M = 1 ; F = 2
dfs <- dfs %>%
  mutate(sex = recode(sex, "M" = 1, "F" = 2))

# Replace Unknown for NA 
dfs$hhi <- replace(dfs$hhi, dfs$hhi==12, NA)
dfs$sex <- replace(dfs$sex, dfs$sex==".", NA)
dfs$race <- replace(dfs$race, dfs$race==10 | dfs$race == 11, NA)
dfs$ethnicity <- replace(df$ethnicity, dfs$ethnicity==2 | dfs$ethnicity == 3, NA)

# Save the raw data file
write.csv(dfs, "raw_data.csv")

## IMPUTATION 
# make the copy of raw data
data <- dfs

dat_imp <- data

colSums(is.na(dat_imp)) # check the number of missingness for each variable

# Make sure what is needed to factorize is factorize
dat_imp$sex     <- factor(dat_imp$sex)              
dat_imp$hhi  <- if (!is.ordered(dat_imp$hhi)) ordered(dat_imp$hhi) else dat_imp$income
dat_imp$puberty_stage <- if (!is.ordered(dat_imp$puberty_stage)) ordered(dat_imp$puberty_stage) else dat_imp$puberty_stage
dat_imp$sleep_hr <- if (!is.ordered(dat_imp$sleep_hr)) ordered(dat_imp$sleep_hr) else dat_imp$sleep_hr
dat_imp$race <- if (!is.ordered(dat_imp$race)) ordered(dat_imp$race) else dat_imp$race
dat_imp$ethnicity <- if (!is.ordered(dat_imp$ethnicity)) ordered(dat_imp$ethnicity) else dat_imp$ethnicity

vars_outcomes <- c(
  "SDS_DA_T","SDS_DIMS_T","SDS_DOES_T","SDS_SBD_T", "SDS_SHY_T", "SDS_SWDT_T",
  "SCARED_P_GD" ,"SCARED_P_PN","SCARED_P_SC", "SCARED_P_SH", "SCARED_P_SP","SCARED_P_Total",
  "CBCL_Int_T","CBCL_Ext_T","CBCL_Total_T"       
)

vars_covars <- c("Identifiers", "BMI_perc", "hhi", "puberty_stage", "sleep_hr", "age",
                 "BMI", "sex", "ethnicity", "race", "comorbidity_count")

vars_all    <- c(vars_outcomes, vars_covars)

# Only extract the variable that will be imputed
dat_imp <- dat_imp[, vars_all, drop = FALSE]

ini  <- mice::mice(dat_imp, maxit = 0)
meth <- ini$method
pred <- ini$predictorMatrix

meth[vars_outcomes] <- "pmm"

meth["age"]      <- "pmm"
meth["sex"]      <- "logreg"    # binary
meth["hhi"]   <- if (is.ordered(dat_imp$income)) "polr" else "polyreg"
meth["sleep_hr"]   <- if (is.ordered(dat_imp$sleep_hr)) "polr" else "polyreg"
meth["puberty_stage"]  <- if (is.ordered(dat_imp$puberty_stage)) "polr" else "polyreg"
meth["ethnicity"]  <- if (is.ordered(dat_imp$ethnicity)) "polr" else "polyreg"
meth["race"]  <- if (is.ordered(dat_imp$race)) "polr" else "polyreg"
meth["BMI"]      <- "pmm"
meth["BMI_perc"] <- "pmm"  
meth["comorbidity_count"] <- "pmm"

# exclude ID for imputation (because we cannot impute the ID) 
meth["Identifiers"] <- ""
pred[, "Identifiers"] <- 0
pred["Identifiers", ] <- 0

# Run mice (make sure the 'mice package' is loaded
imp <- mice::mice(
  dat_imp,
  m = 30,
  method = meth,
  predictorMatrix = pred,
  maxit = 50,
  seed = 20260108
)

# Extract the complete set and save it
data_comp1 <- mice::complete(imp, action = "all")
saveRDS(data_comp1, file = "data_after_mice.rds")
