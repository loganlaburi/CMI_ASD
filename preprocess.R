
# load package
library(sjPlot)
library(dplyr)
library(broom)
library(mice)

# Read data
asd_data <- read.csv('asd_only.csv')

# Extract only the followings
df <- asd_data[, c("Identifiers", "bmi_percentile", "FSQ.FSQ_04", "SDS.SDS_01",
                      "PPS.PPS_F_Score", "PPS.PPS_M_Score", 
                      "Physical.Age", "Physical.BMI", "Physical.Sex",
                      "PreInt_Demos_Fam.Child_Ethnicity",
                      "PreInt_Demos_Fam.Child_Race",
                      "SDS.SDS_DA_Raw", "SDS.SDS_DIMS_Raw",
                      "SDS.SDS_DOES_Raw", "SDS.SDS_SBD_Raw",
                      "SDS.SDS_SHY_Raw", "SDS.SDS_SWTD_Raw", "SDS.SDS_Total_Raw",
                      "SCARED_P.SCARED_P_GD", "SCARED_P.SCARED_P_PN",
                      "SCARED_P.SCARED_P_SC", "SCARED_P.SCARED_P_SH",
                      "SCARED_P.SCARED_P_SP", "SCARED_P.SCARED_P_Total",
                      "CBCL.CBCL_Int", "CBCL.CBCL_Ext", 
                      "CBCL.CBCL_Total", "comorbidity_count")]

# Change the column names
colnames(df) <- c("Identifiers", "bmi_percentile", "hhi", "sleep_hr",
                  "pps_f_total", "pps_m_total", "age", "BMI", "sex",
                  "ethnicity", "race", 
                  "SDS_DA", "SDS_DIMS", "SDS_DOES",
                  "SDS_SBD", "SDS_SHY", "SDS_SWDT", "SDS_Total",
                  "SCARED_P_GD", "SCARED_P_PN", "SCARED_P_SC", 
                  "SCARED_P_SH","SCARED_P_SP", "SCARED_P_Total",
                  "CBCL_Int", "CBCL_Ext", "CBCL_Total", "comorbidity_count")

df <- df %>%
  mutate(
    pps_m_total = suppressWarnings(as.numeric(pps_m_total)),
    pps_f_total = suppressWarnings(as.numeric(pps_f_total))
  )

# Put the puberty score into the categorical character, and convert 
# back to the numeric 
df <- df %>%
  mutate(
    puberty_stage = case_when(
      sex == 1 & pps_m_total == 3 ~ 1,
      sex == 1 & pps_m_total %in% c(4, 5) ~ 2,
      sex == 1 & pps_m_total %in% c(6, 7, 8) ~ 3,
      sex == 1 & pps_m_total %in% c(9, 10, 11) ~ 4,
      sex == 1 & pps_m_total == 12 ~ 5,
      
      sex == 2 & pps_f_total == 1 ~ 1,
      sex == 2 & pps_f_total %in% c(2) ~ 2,
      sex == 2 & pps_f_total %in% c(3) ~ 3,
      sex == 2 & pps_f_total %in% c(4) ~ 4,
      sex == 2 & pps_f_total == 5 ~ 5,
      
      TRUE ~ NA_real_
    )
  )

# convert everything except sex to numeric
df <- df %>%
  dplyr::mutate(
    across(-c("Identifiers"), ~as.numeric(as.character(.x)))
  )

# Remove the columns
df <- df %>%
  dplyr::select(-c(pps_f_total, pps_m_total))

# Replace Unknown for NA 
df[df == "."] <- NA
df$hhi <- replace(df$hhi, df$hhi==12, NA)
df$sex <- replace(df$sex, df$sex==".", NA)
df$race <- replace(df$race, df$race==10 | df$race == 11, NA)
df$ethnicity <- replace(df$ethnicity, df$ethnicity==2 | df$ethnicity == 3, NA)

write.csv(df, "final_data.csv")

# Imputation
data <- df
dat_imp <- data

colSums(is.na(dat_imp)) 

dat_imp$sex     <- factor(dat_imp$sex)              
dat_imp$hhi  <- if (!is.ordered(dat_imp$hhi)) ordered(dat_imp$hhi) else dat_imp$income
dat_imp$puberty_stage <- if (!is.ordered(dat_imp$puberty_stage)) ordered(dat_imp$puberty_stage) else dat_imp$puberty_stage
dat_imp$race <- if (!is.ordered(dat_imp$race)) ordered(dat_imp$race) else dat_imp$race
dat_imp$ethnicity <- if (!is.ordered(dat_imp$ethnicity)) ordered(dat_imp$ethnicity) else dat_imp$ethnicity

# Define imputation variable
vars_outcomes <- c(
  "SDS_DA","SDS_DIMS","SDS_DOES","SDS_SBD", "SDS_SHY", "SDS_SWDT", "SDS_Total",
  "SCARED_P_GD" ,"SCARED_P_PN","SCARED_P_SC", "SCARED_P_SH", "SCARED_P_SP","SCARED_P_Total",
  "CBCL_Int","CBCL_Ext","CBCL_Total"       
)

vars_covars <- c("Identifiers", "hhi", "puberty_stage", "age", "bmi_percentile", "sleep_hr",
                 "BMI", "sex", "ethnicity", "race", "comorbidity_count")

vars_all    <- c(vars_outcomes, vars_covars)

# Subset to working frame
dat_imp <- dat_imp[, vars_all, drop = FALSE]

ini  <- mice::mice(dat_imp, maxit = 0)
meth <- ini$method
pred <- ini$predictorMatrix

# Methods for outcomes 
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

# Extract the completed set
data_comp1 <- mice::complete(imp, action = "all")
saveRDS(data_comp1, file = "data_after_mice.rds")
