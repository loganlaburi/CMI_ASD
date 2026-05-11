
# load package
library(sjPlot)
library(dplyr)
library(broom)
library(emmeans)
library(car)
library(ggplot2)
library(tidyr)
library(purrr)  
library(miceadds)
library(mice)
library(emmeans)
library(stringr)
library(lmtest)
library(sandwich)

# Read data
data <- readRDS("data_after_mice.rds")

# Set seed
set.seed(123)

# Check if there is missing value 
sapply(data, function(d) sum(is.na(as.data.frame(d))))


## Descripitive statistics
mutate_vars <- function(df) {
  df <- as.data.frame(df)
  
  # categorize the age group
  df <- df %>%
    mutate(
      age_group = case_when(
        age >=  9 & age <= 10.99 ~ "pre", # pre adolescent
        age >= 11 & age <= 14.99 ~ "early", # early adolescent
        age >= 15 & age <= 19.99 ~ "late", # middle-to-late adolescents,
        TRUE ~ NA_character_
      )
    )
  
  # BMI category using BMI percentile 
  df$bmi_category <- with(df, ifelse(
    is.na(BMI_perc), NA,
    ifelse(BMI_perc < 5, "Underweight",
           ifelse(BMI_perc < 85, "Normal Weight",
                  ifelse(BMI_perc < 95, "Overweight", "Obesity")
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

var_int <- c("SDS_DA_T", "SDS_DIMS_T", "SDS_DOES_T",
             "SDS_SBD_T", "SDS_SHY_T", "SDS_SWDT_T", "SDS_Total_T",
             "SCARED_P_GD","SCARED_P_PN", "SCARED_P_SC",
             "SCARED_P_SH", "SCARED_P_SP", "SCARED_P_Total",
             "CBCL_Int_T", "CBCL_Ext_T", "CBCL_Total_T",
             "BMI_perc", "age", "BMI")  

summaries_by_imp <- lapply(seq_along(df_imp), function(i) {
  d <- df_imp[[i]]
  
  d %>%
    group_by(age_group) %>%
    summarise(
      across(
        all_of(var_int[var_int %in% names(d) &
                              sapply(d[var_int], is.numeric)]),
        list(mean = ~mean(.x, na.rm = TRUE),
             sd   = ~sd(.x, na.rm = TRUE)),
        .names = "{.col}__{.fn}"
      ),
      .groups = "drop"
    ) %>%
    mutate(.imp = i)
})

summary_all <- bind_rows(summaries_by_imp)

summary_pooled <- summary_all %>%
  pivot_longer(
    cols = -c(age_group, .imp),
    names_to = c("Variable", "Stat"),
    names_sep = "__",
    values_to = "Value"
  ) %>%
  group_by(age_group, Variable, Stat) %>%
  summarise(Value = mean(Value, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(
    names_from = Stat,
    values_from = Value
  ) %>%
  mutate(Mean_SD = sprintf("%.2f ± %.2f", mean, sd)) %>%
  select(age_group, Variable, Mean_SD, mean, sd)

# Calculating the descripitive stat in the categorical/factorize data
cat_var <- "sleep_hr" 

# For each imputation, change everything to numeric and get mean +/- sd
summaries_by_imp <- lapply(seq_along(df_imp), function(i) {
  d <- as.data.frame(df_imp[[i]])
  if (!cat_var %in% names(d)) return(NULL)
  
  x <- d[[cat_var]]
  if (is.factor(x)) x <- as.character(x)
  if (is.character(x)) x <- suppressWarnings(as.numeric(x))
  d[[cat_var]] <- x
  
  d %>%
    filter(!is.na(age_group)) %>%
    group_by(age_group) %>%
    summarise(
      mean = mean(.data[[cat_var]], na.rm = TRUE),
      sd   = sd(.data[[cat_var]],   na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(.imp = i, Variable = cat_var)
})

# Combine and average across imputations
summary_pooled <- bind_rows(summaries_by_imp) %>%
  group_by(age_group, Variable) %>%
  summarise(
    mean = mean(mean, na.rm = TRUE),
    sd   = mean(sd,   na.rm = TRUE),   # descriptive average of SDs
    .groups = "drop"
  ) %>%
  mutate(Mean_SD = sprintf("%.2f ± %.2f", mean, sd)) %>%
  select(age_group, Variable, Mean_SD, mean, sd)


cat_vars <- c("sex","race","ethnicity","bmi_category","puberty_stage",
              "hhi","meets_guideline","comorbidity_count")

levels_union <- lapply(cat_vars, function(v) {
  lev <- unique(unlist(lapply(df_imp, function(d) {
    if (v %in% names(d)) {
      x <- d[[v]]
      if (is.factor(x)) levels(x) else unique(as.character(x))
    } else NULL
  })))
  lev[!is.na(lev)]
})
names(levels_union) <- cat_vars

df_imp <- lapply(df_imp, function(d) {
  d <- as.data.frame(d)
  for (v in cat_vars) {
    if (v %in% names(d)) {
      if (!is.factor(d[[v]])) d[[v]] <- as.factor(d[[v]])
      d[[v]] <- factor(d[[v]], levels = levels_union[[v]])
    }
  }
  d
})

# average age group across denominator
group_sizes_all <- bind_rows(lapply(seq_along(df_imp), function(i) {
  as.data.frame(df_imp[[i]]) %>%
    filter(!is.na(age_group)) %>%
    count(age_group, name = "N_group") %>%
    mutate(.imp = i)
}))
group_sizes_avg <- group_sizes_all %>%
  group_by(age_group) %>%
  summarise(group_N_avg = mean(N_group), .groups = "drop")

# Average raw counts per age_group × Variable × Category across imputations
counts_by_imp <- lapply(seq_along(df_imp), function(i) {
  d <- as.data.frame(df_imp[[i]])
  d <- d[!is.na(d$age_group), , drop = FALSE]
  
  bind_rows(lapply(cat_vars[cat_vars %in% names(d)], function(v) {
    grid <- expand.grid(
      age_group = unique(d$age_group),
      Category  = levels(d[[v]]),
      KEEP.OUT.ATTRS = FALSE,
      stringsAsFactors = FALSE
    )
    cnt <- d %>%
      count(age_group, Category = .data[[v]], name = "N") %>%
      as_tibble()
    
    grid %>%
      left_join(cnt, by = c("age_group","Category")) %>%
      mutate(N = ifelse(is.na(N), 0L, N),
             Variable = v,
             .imp = i)
  }))
})
counts_long <- bind_rows(counts_by_imp)

counts_avg <- counts_long %>%
  group_by(age_group, Variable, Category) %>%
  summarise(N_mean = mean(N), .groups = "drop")

# Rescale averaged counts so they sum to group totals (and round exactly)
counts_adjusted <- counts_avg %>%
  left_join(group_sizes_avg, by = "age_group") %>%
  group_by(age_group, Variable) %>%
  mutate(
    total_mean = sum(N_mean, na.rm = TRUE),
    # proportional rescale within each Variable so categories sum to group_N_avg
    N_rescaled = ifelse(total_mean > 0, N_mean / total_mean * group_N_avg, 0),
    # integer rounding while preserving the exact group total
    N_floor = floor(N_rescaled),
    remainder = N_rescaled - N_floor,
    shortfall = round(unique(group_N_avg)) - sum(N_floor),
    rank_rem = rank(-remainder, ties.method = "first"),
    add_one = ifelse(rank_rem <= pmax(shortfall, 0), 1L, 0L),
    N = N_floor + add_one
  ) %>%
  ungroup() %>%
  select(age_group, Variable, Category, N) %>%
  arrange(Variable, age_group, Category)

## Correlation analysis
to_numeric_safe <- function(x) {
  if (is.numeric(x)) return(x)
  if (is.factor(x)) x <- as.character(x)
  if (is.character(x)) {
    y <- suppressWarnings(as.numeric(x))
    if (all(is.na(y))) return(x)
    return(y)
  }
  suppressWarnings(as.numeric(x))
}

# Choose the column that we want for the correlation
pick_cols <- c("Identifiers", names(df_imp[[1]])[c(17:21, 23, 1:15, 26:28)])

# Extract only existing variables, per imputation
dat_list <- lapply(df_imp, function(d) d[, intersect(pick_cols, names(d)), drop = FALSE])

pos_vars <- 2:25
dat_list_num <- lapply(dat_list, function(d) {
  d2 <- d
  
  # Recode sex to 1/2 consistently across imputations
  if ("sex" %in% names(d2)) {
    sx <- tolower(trimws(as.character(d2$sex)))
    d2$sex <- ifelse(sx %in% c("m"), 1,
                     ifelse(sx %in% c("f"), 2, NA_real_))
  }
  
  # convert everything to numeric except ID and sex
  cols <- intersect(pos_vars, seq_along(d2))
  for (j in cols) {
    nm <- names(d2)[j]
    if (nm %in% c("Identifiers", "sex")) next
    d2[[j]] <- to_numeric_safe(d2[[j]])
  }
  
  d2
  
})

# Keep only variables that are numeric in all imputations
common_names <- names(dat_list_num[[1]])
numeric_in_all <- common_names[
  sapply(common_names, function(v) all(sapply(dat_list_num, function(d) is.numeric(d[[v]]))))
]
vars_final <- setdiff(numeric_in_all, c("Identifiers"))

vars_final <- intersect(vars_final, names(dat_list_num[[1]])[pos_vars])

stopifnot(length(vars_final) >= 2)  

# Run pooled correlation
imp_data_num <- miceadds::datlist2mids(dat_list_num)
cor_out <- miceadds::micombine.cor(mi.res = imp_data_num, variables = vars_final, method = "pearson")

# save correlation result into csv
write.csv(cor_out, "correlation.csv")


## Run the linear regression
mi_lm_summary <- function(mids_obj, formula) {
  if (!inherits(formula, "formula")) formula <- stats::as.formula(formula)
  
  # Complete all imputations to a list of data frames
  imps <- mice::complete(mids_obj, action = "all")   # list length = m
  
  # Fit lm() in each completed dataset 
  fits <- lapply(imps, function(d) stats::lm(formula, data = d))
  
  # Convert to a 'mira' object and pool (Rubin's rules)
  mira_obj <- mice::as.mira(fits)
  pooled   <- mice::pool(mira_obj)
  
  # Coefficient table 
  coef_table <- as.data.frame(summary(pooled))
  
  # Adjusted R^2 
  adj_r2 <- mean(sapply(fits, function(m) summary(m)$adj.r.squared), na.rm = TRUE)
  
  # Pooled F, df1, df2, p 
  fstats <- do.call(rbind, lapply(fits, function(m) {
    fs <- summary(m)$fstatistic
    if (is.null(fs) || length(fs) < 3) return(c(F = NA, df1 = NA, df2 = NA))
    c(F = unname(fs[1]), df1 = unname(fs[2]), df2 = unname(fs[3]))
  }))
  valid <- stats::complete.cases(fstats)
  if (any(valid)) {
    F_pooled   <- mean(fstats[valid, "F"])
    df1_pooled <- mean(fstats[valid, "df1"])
    df2_pooled <- mean(fstats[valid, "df2"])
    p_pooled   <- stats::pf(F_pooled, df1_pooled, df2_pooled, lower.tail = FALSE)
  } else {
    F_pooled <- df1_pooled <- df2_pooled <- p_pooled <- NA_real_
  }
  
  list(
    coefficients = coef_table,  
    adj_r2       = adj_r2,      
    F            = F_pooled,     
    df1          = df1_pooled,   
    df2          = df2_pooled,
    p_value_F    = p_pooled,
    fstats_raw   = fstats       
  )
}

# Run the linear regression function 
res1 <- mi_lm_summary(imp_data_num, 
                      CBCL_Total_T ~ SDS_DA_T + SDS_DIMS_T + SDS_DOES_T + SDS_SBD_T + 
                        SDS_SHY_T + SDS_SWDT_T + BMI_perc + hhi + puberty_stage +
                        age + sex + comorbidity_count)

# Model summary
model_df <- data.frame(
  adj_r2   = res1$adj_r2,
  F        = res1$F,
  df1      = res1$df1,
  df2      = res1$df2,
  p_value  = res1$p_value_F
)

## MANOVA + post hoc test 
# Complete each imputation and align factor levels for group
complete_align <- function(mids_obj, group = "age_group", levels = c("pre","early","late")) {
  imps <- mice::complete(mids_obj, action = "all")
  lapply(imps, function(d) {
    stopifnot(group %in% names(d))
    d[[group]] <- factor(d[[group]], levels = levels)
    d
  })
}

# Manual Rubin pooling for a scalar 
rubin_pool_scalar <- function(est, se, df_complete = NULL) {
  est <- as.numeric(est); se <- as.numeric(se)
  ok  <- is.finite(est) & is.finite(se)
  est <- est[ok]; se <- se[ok]
  m   <- length(est)
  if (m < 1L) {
    return(data.frame(estimate = NA_real_, se = NA_real_, df = NA_real_,
                      t_ratio = NA_real_, p_value = NA_real_))
  }
  Qbar <- mean(est)
  Ubar <- mean(se^2)
  B    <- if (m > 1L) stats::var(est) else 0
  Tvar <- Ubar + (1 + 1/m) * B
  se_p <- sqrt(Tvar)
  
  if (B <= .Machine$double.eps) {
    df <- if (!is.null(df_complete)) mean(df_complete, na.rm = TRUE) else 1e9
  } else {
    r  <- ((1 + 1/m) * B) / Ubar
    df <- (m - 1) * (1 + 1/r)^2
  }
  
  tval <- Qbar / se_p
  pval <- 2 * stats::pt(abs(tval), df = df, lower.tail = FALSE)
  data.frame(estimate = Qbar, se = se_p, df = df, t_ratio = tval, p_value = pval)
}

# Normalize emmeans contrast rows 
normalize_contrast_rows <- function(pairs_df) {
  out <- pairs_df %>%
    tidyr::separate(contrast, into = c("g1","g2"), sep = " - ", remove = FALSE) %>%
    mutate(
      target = case_when(
        (g1 == "early" & g2 == "late") | (g1 == "late"  & g2 == "early") ~ "early-late",
        (g1 == "early" & g2 == "pre")  | (g1 == "pre"   & g2 == "early") ~ "early-pre",
        (g1 == "late"  & g2 == "pre")  | (g1 == "pre"   & g2 == "late")  ~ "late-pre",
        TRUE ~ NA_character_
      ),
      flip = case_when(                      # need to flip sign when reversed
        target == "early-late" & g1 == "late"  & g2 == "early" ~ TRUE,
        target == "early-pre"  & g1 == "pre"   & g2 == "early" ~ TRUE,
        target == "late-pre"   & g1 == "pre"   & g2 == "late"  ~ TRUE,
        TRUE ~ FALSE
      ),
      estimate = ifelse(flip, -estimate, estimate),
      `t.ratio` = ifelse(flip, -`t.ratio`, `t.ratio`),
      group1 = case_when(
        target == "early-late" ~ "early",
        target == "early-pre"  ~ "early",
        target == "late-pre"   ~ "late",
        TRUE ~ NA_character_
      ),
      group2 = case_when(
        target == "early-late" ~ "late",
        target == "early-pre"  ~ "pre",
        target == "late-pre"   ~ "pre",
        TRUE ~ NA_character_
      )
    ) %>%
    filter(!is.na(target)) %>%
    transmute(contrast = target, group1, group2,
              estimate, SE, df = df, t.ratio = `t.ratio`, p.value = `p.value`)
  out
}

# Tukey post-hoc pooled per DV
mi_tukey_pooled <- function(mids_obj, dv, group = "age_group",
                            levels = c("pre","early","late"),
                            adjust = "tukey",
                            robust = TRUE) {
  
  imps <- complete_align(mids_obj, group, levels)
  
  per_imp <- lapply(imps, function(d) {
    if (!all(c(dv, group) %in% names(d))) return(NULL)
    
    d <- d[is.finite(d[[dv]]) & !is.na(d[[group]]), , drop = FALSE]
    if (!nrow(d)) return(NULL)
    
    # Linear model
    m <- stats::lm(as.formula(paste(dv, "~", group)), data = d)
    
    # Robust covariance (HC3) to accomodate heterogenity 
    vc <- if (robust) {
      sandwich::vcovHC(m, type = "HC3")
    } else {
      stats::vcov(m)
    }
    
    # Estimated marginal means with robust SEs
    em <- emmeans::emmeans(
      m,
      as.formula(paste("~", group)),
      vcov. = vc
    )
    
    pr <- emmeans::contrast(em, method = "pairwise", adjust = adjust)
    pairs_df <- as.data.frame(pr)
    
    # Normalize contrasts
    norm_df <- normalize_contrast_rows(pairs_df)
    
    # Group means and SDs
    means_df <- d %>%
      dplyr::group_by(.group = .data[[group]]) %>%
      dplyr::summarise(
        mean = mean(.data[[dv]], na.rm = TRUE),
        sd   = sd(.data[[dv]],   na.rm = TRUE),
        .groups = "drop"
      )
    
    list(pairs = norm_df, means = means_df)
  })
  
  per_imp <- Filter(Negate(is.null), per_imp)
  if (!length(per_imp)) return(NULL)
  
  # Pool contrasts using Rubin's rules
  pairs_all <- dplyr::bind_rows(lapply(seq_along(per_imp), function(i) {
    x <- per_imp[[i]]$pairs
    x$.imp <- i
    x
  }))
  
  pooled_pairs <- dplyr::bind_rows(lapply(
    split(pairs_all, pairs_all$contrast),
    function(cdat) {
      pooled <- rubin_pool_scalar(
        est = cdat$estimate,
        se  = cdat$SE,
        df_complete = cdat$df
      )
      cbind(
        contrast = cdat$contrast[1],
        group1   = cdat$group1[1],
        group2   = cdat$group2[1],
        pooled
      )
    }
  ))
  
  # Pool means and SDs
  means_all <- dplyr::bind_rows(lapply(per_imp, function(x) x$means))
  
  means_pooled <- means_all %>%
    dplyr::group_by(.group) %>%
    dplyr::summarise(
      mean = mean(mean, na.rm = TRUE),
      sd   = mean(sd,   na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::rename(level = .group)
  
  final <- pooled_pairs %>%
    dplyr::left_join(means_pooled, by = c("group1" = "level")) %>%
    dplyr::rename(mean1 = mean, sd1 = sd) %>%
    dplyr::left_join(means_pooled, by = c("group2" = "level")) %>%
    dplyr::rename(mean2 = mean, sd2 = sd) %>%
    dplyr::mutate(DV = dv) %>%
    dplyr::select(
      DV, contrast, group1, group2,
      df, t_ratio, p_value,
      mean1, sd1, mean2, sd2
    ) %>%
    dplyr::arrange(
      factor(contrast, levels = c("early-late","early-pre","late-pre"))
    )
  
  final
}


## MANOVA (Pillai) across imputations
mi_manova_pillai <- function(mids_obj, dvs, group = "age_group",
                             levels = c("pre","early","late"),
                             stouffer_weights = c("equal","sqrt_df")) {
  stouffer_weights <- match.arg(stouffer_weights)
  imps <- complete_align(mids_obj, group, levels)
  
  rows <- lapply(imps, function(d) {
    if (!all(c(dvs, group) %in% names(d))) return(NULL)
    f <- as.formula(paste0("cbind(", paste(dvs, collapse = ","), ") ~ ", group))
    fit <- try(manova(f, data = d), silent = TRUE)
    if (inherits(fit, "try-error")) return(NULL)
    S <- summary(fit, test = "Pillai")$stats
    rn <- rownames(S)
    idx <- which(rn == group); if (length(idx) != 1) idx <- 1L
    
    cn <- colnames(S)
    pick <- function(pat) cn[grepl(pat, cn, ignore.case = TRUE)][1]
    
    data.frame(
      Pillai = as.numeric(S[idx, pick("^pillai")]),
      approxF= as.numeric(S[idx, pick("approx")]),
      numdf  = as.numeric(S[idx, pick("num")]),
      dendf  = as.numeric(S[idx, pick("den")]),
      p      = as.numeric(S[idx, pick("pr\\(>f\\)|p")])
    )
  })
  rows <- Filter(Negate(is.null), rows)
  if (!length(rows)) {
    return(data.frame(
      set = paste(dvs, collapse = "+"),
      Pillai_mean = NA, approxF_mean = NA, numdf_mean = NA, dendf_mean = NA,
      p_fisher = NA, p_stouffer = NA, p_median = NA, m_used = 0L
    ))
  }
  df_all <- dplyr::bind_rows(rows)

  eps <- .Machine$double.xmin
  df_all$p <- pmin(pmax(df_all$p, eps), 1 - eps)
  df_all <- df_all[is.finite(df_all$p), , drop = FALSE]
  m_used <- nrow(df_all)
  if (m_used == 0) {
    return(data.frame(
      set = paste(dvs, collapse = "+"),
      Pillai_mean = NA, approxF_mean = NA, numdf_mean = NA, dendf_mean = NA,
      p_fisher = NA, p_stouffer = NA, p_median = NA, m_used = 0L
    ))
  }
  
  # Fisher’s method
  X2 <- -2 * sum(log(df_all$p))
  p_fisher <- 1 - stats::pchisq(X2, df = 2 * m_used)
  
  # Stouffer 
  w <- if (stouffer_weights == "sqrt_df") sqrt(pmax(df_all$dendf, 1)) else rep(1, m_used)
  # Two-sided z from p
  z <- stats::qnorm(1 - df_all$p / 2)
  Zc <- sum(w * z) / sqrt(sum(w^2))
  p_stouffer <- 2 * (1 - stats::pnorm(abs(Zc)))
  
  data.frame(
    set          = paste(dvs, collapse = "+"),
    Pillai_mean  = mean(df_all$Pillai,  na.rm = TRUE),
    approxF_mean = mean(df_all$approxF, na.rm = TRUE),
    numdf_mean   = mean(df_all$numdf,   na.rm = TRUE),
    dendf_mean   = mean(df_all$dendf,   na.rm = TRUE),
    p_fisher     = p_fisher,
    p_stouffer   = p_stouffer,
    p_median     = stats::median(df_all$p),
    m_used       = m_used
  )
}

## RUN MANOVA and post hoc function 
mids_obj   <- imp_data_num              # make sure MIDs object
group_var  <- "age_group"
age_levels <- c("pre","early","late")

SCARED <- c("SCARED_P_GD","SCARED_P_PN","SCARED_P_SC","SCARED_P_SH","SCARED_P_SP")
SDS    <- c("SDS_DA_T","SDS_DIMS_T","SDS_DOES_T","SDS_SBD_T","SDS_SHY_T","SDS_SWDT_T")
CBCL   <- c("CBCL_Int_T","CBCL_Ext_T")

manova_scared <- mi_manova_pillai(mids_obj, SCARED, group = group_var,
                                  levels = age_levels, stouffer_weights = "sqrt_df")
manova_sds <- mi_manova_pillai(mids_obj, SDS, group = group_var,
                                  levels = age_levels, stouffer_weights = "sqrt_df")
manova_cbcl <- mi_manova_pillai(mids_obj, CBCL, group = group_var,
                               levels = age_levels, stouffer_weights = "sqrt_df")

tukey_scared <- dplyr::bind_rows(lapply(SCARED, function(dv) {
  mi_tukey_pooled(mids_obj, dv, group = group_var, levels = age_levels, adjust = "tukey")
}))

tukey_sds<- dplyr::bind_rows(lapply(SDS, function(dv) {
  mi_tukey_pooled(mids_obj, dv, group = group_var, levels = age_levels, adjust = "tukey")
}))

tukey_cbcl<- dplyr::bind_rows(lapply(CBCL, function(dv) {
  mi_tukey_pooled(mids_obj, dv, group = group_var, levels = age_levels, adjust = "tukey")
}))

## Chi-square test 
chis_by_imp <- lapply(df_imp, function(d) {
  d <- as.data.frame(d)
  
  if (!("age_group" %in% names(d))) stop("age_group not found in dataset")
  if (!("sex"       %in% names(d))) stop("sex not found in dataset")
  
  d$age_group <- as.factor(d$age_group)
  d$sex       <- as.factor(d$sex)
  
  # Build the 2-way table and run Pearson chi-square
  tab <- table(d$age_group, d$sex)
  
  if (any(dim(tab) == 0)) stop("Empty table: check factor levels or data.")
  
  cs <- suppressWarnings(chisq.test(tab, correct = FALSE))
  
  list(stat = unname(cs$statistic), df = unname(cs$parameter))
})


## Extract vectors and sanity-check
chisq_stats <- vapply(chis_by_imp, `[[`, numeric(1), "stat")
dfs         <- vapply(chis_by_imp, `[[`, numeric(1), "df")

# Make sure they are plain vectors with no dimensions
chisq_stats <- as.numeric(chisq_stats)
dfs         <- as.numeric(dfs)

# Drop any NA/Inf 
ok <- is.finite(chisq_stats) & is.finite(dfs)
chisq_stats <- chisq_stats[ok]
dfs         <- dfs[ok]

## Run chi-square on the imputed data
mi_res <- miceadds::micombine.chisquare(chisq_stats, dfs)

# Pull chi-square result
F_pool <- unname(mi_res[grep("^D\\d+$",  names(mi_res))][1])
p_pool <- unname(mi_res[grep("^p\\d+$",  names(mi_res))][1])
df_num <- unname(mi_res[grep("^df\\d+$", names(mi_res))][1])
df_den <- unname(mi_res[grep("^df2\\d+$", names(mi_res))][1])


## ANOVA 
fit <- with(
  imp_data_num,
  lm(SDS_Total_T ~ age_group)
)

stats <- lapply(fit$analyses, function(m) {
  
  a <- car::Anova(
    m,
    vcov. = vcovHC(m, type = "HC3")
  )
  
  data.frame(
    F   = a[1, "F"],
    df1 = a[1, "Df"],
    df2 = a["Residuals", "Df"]
  )
})

stats <- do.call(rbind, stats)

F_bar  <- mean(stats$F)
df1bar <- round(mean(stats$df1))
df2bar <- round(mean(stats$df2))

p_val <- pf(F_bar, df1bar, df2bar, lower.tail = FALSE)

c(F = F_bar, df1 = df1bar, df2 = df2bar, p = p_val)

# Classification of the clinical-range
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
  sdsc_vars <- c("SDS_DA_T","SDS_DIMS_T","SDS_DOES_T",
                 "SDS_SBD_T","SDS_SHY_T","SDS_SWDT_T", "SDS_Total_T")
  for (v in sdsc_vars[sdsc_vars %in% names(d)]) {
    d[[paste0(v, "_class")]] <- class_sdsc(d[[v]])
  }
  
  ## CBCL
  cbcl_vars <- c("CBCL_Int_T","CBCL_Ext_T","CBCL_Total_T")
  for (v in cbcl_vars[cbcl_vars %in% names(d)]) {
    d[[paste0(v, "_class")]] <- class_cbcl(d[[v]])
  }
  
  ## SCARED
  if ("SCARED_P_Total" %in% names(d))
    d$SCARED_P_Total_class <- class_scared(d$SCARED_P_Total, 25)
  if ("SCARED_P_GD" %in% names(d))
    d$SCARED_P_GD_class <- class_scared(d$SCARED_P_GD, 9)
  if ("SCARED_P_PN" %in% names(d))
    d$SCARED_P_PN_class <- class_scared(d$SCARED_P_PN, 7)
  if ("SCARED_P_SC" %in% names(d))
    d$SCARED_P_SC_class <- class_scared(d$SCARED_P_SC, 8)
  if ("SCARED_P_SH" %in% names(d))
    d$SCARED_P_SH_class <- class_scared(d$SCARED_P_SH, 3)
  if ("SCARED_P_SP" %in% names(d))
    d$SCARED_P_SP_class <- class_scared(d$SCARED_P_SP, 5)
  
  d
})


get_prop <- function(var) { 
  # compute proportions within each imputation
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
  
  # pool across imputations
  results <- results %>%
    group_by(age_group, Category) %>%
    summarise(
      prop = mean(prop, na.rm = TRUE),
      N = round(mean(total_n, na.rm = TRUE)),
      .groups = "drop"
    )
  
  # compute counts 
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

 # Run the function
get_prop("SDS_Total_T_class")
