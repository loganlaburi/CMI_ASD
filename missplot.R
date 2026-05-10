# Load packages
library(dplyr)
library(naniar)

# Read the raw data file
data <- read.csv( "raw_data.csv")

# Extract the age of interest
df <- subset(data, data$age >=9 & data$age <=19.99)

df$sex[df$sex == ""] <- NA  

# numerical code 
df$sex <- ifelse(df$sex == "M", 1, ifelse(df$sex == "F", 2, NA))

#select columns w/ missing data
missing_data <- df[, c(3:4, 28, 7, 9, 12:27)]

#convert everything to numeric 
missing_data[] <- lapply(missing_data, function(x) {
  x[x == ""] <- NA
  as.numeric(x)
})

#rename columns for plot
colnames(missing_data) <- c(
  "BMI Percentile", "Household Income", "Pubertal Stage", 
  "Age","Sex",
  "SDSC Disorders of Arousal T-score", 
  "SDSC Disorders of Initiating and Maintaining Sleep T-score", 
  "SDSC Excessive Somnolence T-score",
  "SDSC Sleep Breathing Disorder T-score", 
  "SDSC Sleep Hyperhydrosis T-score", 
  "SDSC Sleep Wake-Transition Disorders T-score", 
  "SDSC Total T-score",
  "SCARED Generalized Anxiety Score", "SCARED Somatic/Panic Score", 
  "SCARED Social Anxiety Score", 
  "SCARED School Avoidance Score",
  "SCARED Separation Anxiety Score", 
  "SCARED Total Score",
  "CBCL Internalizing T-Score", 
  "CBCL Externalizing T-score", 
  "CBCL Total Score"
)

#create plot
gg_miss_upset(
  missing_data,
  nsets = ncol(missing_data),   # ensure all columns are included 
  order.by = "freq"             # order by frequency 
)
