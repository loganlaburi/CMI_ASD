# load package
library(sjPlot)
library(dplyr)
library(broom)

# Read data
data <- read.csv('data_asd.csv') # raw data
diag <- read.csv('diagnosis.csv') # Consensus diagnosis
sleep <- read.csv('sleep_hr.csv') # Sleep data
diag <- diag[, c("Identifiers", "Has_Autism")]

# Attach the "has_autism" info based on dplyr
data <- data %>%
  left_join(diag, by="Identifiers")

data <- subset(data, data$Has_Autism==1)

# Add uncorrected CBCL score in it
cbcl <- read.csv('cbcl.csv')

# Combine uncorrected CBCL score on the raw data set 
merged_data <- merge(data, cbcl, by = "Identifiers")

df <- merged_data[, c("Identifiers", "bmi_percentile", "FSQ.FSQ_04",
               "PPS.PPS_F_Score", "PPS.PPS_M_Score", 
               "Physical.Age", "Physical.BMI", "Physical.Sex",
               "PreInt_Demos_Fam.Child_Ethnicity",
               "PreInt_Demos_Fam.Child_Race",
               "SDS.SDS_DA_T", "SDS.SDS_DIMS_T",
               "SDS.SDS_DOES_T", "SDS.SDS_SBD_T",
               "SDS.SDS_SHY_T", "SDS.SDS_SWDT_T",
               "SCARED_P.SCARED_P_GD", "SCARED_P.SCARED_P_PN",
               "SCARED_P.SCARED_P_SC", "SCARED_P.SCARED_P_SH",
               "SCARED_P.SCARED_P_SP", "SCARED_P.SCARED_P_Total",
               "CBCL.CBCL_Int", "CBCL.CBCL_Ext", 
               "CBCL.CBCL_Total")]

# Change the column names
colnames(df) <- c("Identifiers", "BMI_perc", "hhi",
                  "pps_f_total", "pps_m_total", "age", "BMI", "sex",
                  "ethnicity", "race", 
                  "SDS_DA_T", "SDS_DIMS_T", "SDS_DOES_T",
                  "SDS_SBD_T", "SDS_SHY_T", "SDS_SWDT_T",
                  "SCARED_P_GD", "SCARED_P_PN", "SCARED_P_SC", 
                  "SCARED_P_SH","SCARED_P_SP", "SCARED_P_Total",
                  "CBCL_Int", "CBCL_Ext", "CBCL_Total")

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

write.csv(dfs, "final_data.csv")
