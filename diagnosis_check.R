
# Load package
library(sjPlot)
library(dplyr)
library(broom)
library(emmeans)
library(car)
library(ggplot2)

# Read data
data <- read.csv("raw_data.csv")
df <- subset(data, as.numeric(data$Physical.Age) >=9 & 
                 as.numeric(data$Physical.Age) <= 19.99)

dx_cols <- grep("^Diagnosis_ClinicianConsensus\\.DX_(0[1-9]|10)$",
                names(df), value = TRUE)

norm_col <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x[x == ""] <- NA
  x
}

dx_df <- as.data.frame(lapply(df[dx_cols], norm_col), stringsAsFactors = FALSE)

# Unique diagnoses across all cells (excluding ASD)
all_dx <- unlist(dx_df, use.names = FALSE)
unique_dx_excl_asd <- sort(unique(all_dx[!is.na(all_dx) & all_dx != "Autism Spectrum Disorder"]))

ASD_LABEL  <- "Autism Spectrum Disorder"
asd_aliases <- c(ASD_LABEL, "ASD", "Autistic Disorder")  

# Detect diagnosis columns 
dx_pattern <- "^Diagnosis_ClinicianConsensus\\.DX_(0[1-9]|10)$"

dx_cols <- grep(dx_pattern, names(df), value = TRUE)
stopifnot(length(dx_cols) > 0)

dx_df <- as.data.frame(lapply(df[dx_cols], norm_col), stringsAsFactors = FALSE)

# Map each original label 
recode_map <- c(
  # ADHD 
  "ADHD-Combined Type"                                    = "ADHD",
  "ADHD-Hyperactive/Impulsive Type"                       = "ADHD",
  "ADHD-Inattentive Type"                                 = "ADHD",
  "Other Specified Attention-Deficit/Hyperactivity Disorder" = "ADHD",
  "Unspecified Attention-Deficit/Hyperactivity Disorder"  = "ADHD",
  
  # Anxiety 
  "Generalized Anxiety Disorder"                          = "Anxiety Disorder",
  "Other Specified Anxiety Disorder"                      = "Anxiety Disorder",
  "Unspecified Anxiety Disorder"                          = "Anxiety Disorder",
  "Separation Anxiety"                                    = "Anxiety Disorder",
  "Social Anxiety (Social Phobia)"                        = "Anxiety Disorder",
  "Panic Disorder"                                        = "Anxiety Disorder",
  "Agoraphobia"                                           = "Anxiety Disorder",
  "Specific Phobia"                                       = "Anxiety Disorder",
  
  # Depressive disorders
  "Other Specified Depressive Disorder"                   = "Depressive Disorder",
  "Major Depressive Disorder"                             = "Depressive Disorder",
  "Persistent Depressive Disorder (Dysthymia)"            = "Depressive Disorder",
  
  # Intellectual ability
  "Intellectual Disability-Mild"                          = "Intellectual Disability",
  "Intellectual Disability-Moderate"                      = "Intellectual Disability",
  "Borderline Intellectual Functioning"                   = "Borderline Intellectual Functioning", # keep separate (change if desired)
  
  # Learning disorders
  "Specific Learning Disorder with Impairment in Mathematics"        = "Specific Learning Disorder",
  "Specific Learning Disorder with Impairment in Reading"            = "Specific Learning Disorder",
  "Specific Learning Disorder with Impairment in Written Expression" = "Specific Learning Disorder",
  
  # Tic disorders
  "Provisional Tic Disorder"                               = "Tic Disorder",
  "Persistent (Chronic) Motor or Vocal Tic Disorder"       = "Tic Disorder",
  "Other Specified Tic Disorder"                           = "Tic Disorder",
  "Tourettes Disorder"                                     = "Tic Disorder",
  
  # OCD & related
  "Obsessive-Compulsive Disorder"                          = "Obsessive-Compulsive Disorder",
  "Trichotillomania (Hair-Pulling Disorder)"               = "OC-Related Disorder",
  "Excoriation (Skin-Picking) Disorder"                    = "OC-Related Disorder",
  
  # Trauma- & stressor-related / relational
  "Posttraumatic Stress Disorder"                          = "Trauma- and Stressor-Related Disorder",
  "Other Specified Trauma- and Stressor-Related Disorder"  = "Trauma- and Stressor-Related Disorder",
  "Adjustment Disorders"                                   = "Trauma- and Stressor-Related Disorder",
  "Parent-Child Relational Problem"                        = "Relational Problem",
  
  # Disruptive / impulse-control
  "Oppositional Defiant Disorder"                          = "Disruptive/Impulse-Control Disorder",
  "Intermittent Explosive Disorder"                        = "Disruptive/Impulse-Control Disorder",
  
  # Feeding & eating
  "Avoidant/Restrictive Food Intake Disorder"              = "Feeding and Eating Disorder",
  "Binge-Eating Disorder"                                  = "Feeding and Eating Disorder",
  "Pica in Children"                                       = "Feeding and Eating Disorder",
  
  # Sleep
  "Insomnia Disorder"                                      = "Sleep-Wake Disorder",
  
  # Elimination
  "Encopresis"                                             = "Elimination Disorder",
  "Enuresis"                                               = "Elimination Disorder",
  "Other Specified Elimination Disorder with Fecal Symptoms" = "Elimination Disorder",
  
  # Neurodevelopmental—language/speech/motor
  "Language Disorder"                                      = "Language Disorder",
  "Speech Sound Disorder"                                  = "Speech Sound Disorder",
  "Developmental Coordination Disorder"                    = "Developmental Coordination Disorder",
  
  # Psychotic spectrum
  "Schizophrenia"                                          = "Schizophrenia Spectrum and Other Psychotic Disorder",
  "Other Specified Schizophrenia Spectrum and Other Psychotic Disorder" = "Schizophrenia Spectrum and Other Psychotic Disorder",
  
  # Substance
  "Cannabis Use Disorder"                                  = "Substance Use Disorder"
)

# Recode
recode_vec <- function(x, map) {
  x0 <- as.character(x)
  xlow <- tolower(x0)
  keylow <- tolower(names(map))
  m <- match(xlow, keylow)
  out <- ifelse(is.na(m), x0, unname(map[m]))
  out
}

dx_rec <- as.data.frame(lapply(dx_df, recode_vec, map = recode_map),
                        stringsAsFactors = FALSE)

# Check the ASD cohort
has_asd <- apply(dx_rec, 1, function(r) any(tolower(r) %in% tolower(asd_aliases), na.rm = TRUE))
asd_n <- sum(has_asd)
cat("ASD cohort size:", asd_n, "of", nrow(dx_rec), "\n")

# Unique comorbidities excluding autism
all_labels_asd <- unlist(dx_rec[has_asd, , drop = FALSE], use.names = FALSE)
unique_comorbidities <- sort(unique(all_labels_asd[
  !is.na(all_labels_asd) & !(tolower(all_labels_asd) %in% tolower(asd_aliases))
]))
cat("\nUnique comorbidities (collapsed labels, excluding ASD):\n")
print(unique_comorbidities)

# Autism and other comorbidities
## Make sure ASD appears first in the combination string
make_combo_incl_asd <- function(row) {
  vals <- unique(na.omit(row))
  is_asd <- tolower(vals) %in% tolower(asd_aliases)
  vals[is_asd] <- ASD_LABEL
  vals <- unique(vals)
  
  if (any(vals == ASD_LABEL)) {
    others <- sort(vals[vals != ASD_LABEL])
    paste(c(ASD_LABEL, others), collapse = " + ")
  } else {
    paste(sort(vals), collapse = " + ")
  }
}

# Rebuild combos and table
combos <- apply(dx_rec[has_asd, , drop = FALSE], 1, make_combo_incl_asd)
combo_counts <- sort(table(combos), decreasing = TRUE)
combo_tbl <- data.frame(
  combo    = names(combo_counts),
  n_people = as.integer(combo_counts),
  stringsAsFactors = FALSE
)


write.csv(combo_tbl, "comorbid.csv")

# Recode table and map
recode_vec <- function(x, map) {
  x0 <- as.character(x)
  xlow <- tolower(x0)
  keylow <- tolower(names(map))
  m <- match(xlow, keylow)
  ifelse(is.na(m), x0, unname(map[m]))  
}

# Apply across columns of the df
dx_rec <- as.data.frame(lapply(dx_df, recode_vec, map = recode_map), stringsAsFactors = FALSE)

orig_long <- unlist(dx_df,  use.names = FALSE)
cat_long  <- unlist(dx_rec, use.names = FALSE)

keep <- !is.na(orig_long)
orig_long <- orig_long[keep]
cat_long  <- cat_long[keep]

diagnosis_map_tbl <- as.data.frame(table(original = orig_long, category = cat_long),
                                   stringsAsFactors = FALSE)
# Remove zeros 
diagnosis_map_tbl <- diagnosis_map_tbl[diagnosis_map_tbl$Freq > 0, ]
diagnosis_map_tbl <- diagnosis_map_tbl[order(diagnosis_map_tbl$category,
                                             -diagnosis_map_tbl$Freq,
                                             diagnosis_map_tbl$original), ]
diagnosis_map_tbl


## Comorbidities
asd_aliases <- c("Autism Spectrum Disorder", "ASD", "Autistic Disorder", "F84.0")

df$comorbidity_count <- apply(dx_rec, 1, function(row) {
  sum(!is.na(row) & !(tolower(row) %in% tolower(asd_aliases)))
})

df_asd_only <- df[has_asd, ]

write.csv(df_asd_only, "asd_only.csv")
