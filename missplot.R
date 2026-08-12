# Load packages
library(dplyr)
library(naniar)

# Read the raw data
data <- read.csv( "final_data.csv")

# Extract the age
df <- subset(data, data$age >=9 & data$age <=19.99)
df$sex[df$sex == ""] <- NA  
df$sex <- ifelse(df$sex == "M", 1, ifelse(df$sex == "F", 2, NA))

# select columns w/ missing data
missing_data <- df[, c(4, 28, 7, 6, 8, 11:26)]

# convert everything to numeric 
missing_data[] <- lapply(missing_data, function(x) {
  x[x == ""] <- NA
  as.numeric(x)
})

# rename columns for plot
colnames(missing_data) <- c(
  "Household Income", "Pubertal Stage", "BMI", "Age","Sex", 
  "SDSC Disorders of Arousal Score", 
  "SDSC Disorders of Initiating and Maintaining Sleep Score", 
  "SDSC Excessive Somnolence Score",
  "SDSC Sleep Breathing Disorder Score", 
  "SDSC Sleep Hyperhydrosis Score", 
  "SDSC Sleep Wake-Transition Disorders Score", 
  "SDSC Total Score",
  "SCARED Generalized Anxiety Score", "SCARED Somatic/Panic Score", 
  "SCARED Social Anxiety Score", 
  "SCARED School Avoidance Score",
  "SCARED Separation Anxiety Score", 
  "SCARED Total Score",
  "CBCL Internalizing Score", 
  "CBCL Externalizing Score", 
  "CBCL Total Score"
)

# create missing plot
gg_miss_upset(
  missing_data,
  nsets = ncol(missing_data),   
  order.by = "freq"            
)

