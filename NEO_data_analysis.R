# =============================================================================
# NEO reduced-rank analysis
#   - predictors: pre-meal Nightingale metabolites
#   - outcomes  : CVA, CVD, DIA, STR, DTH, MI, TIA
# =============================================================================

library(haven)
library(tidyverse)
library(dplyr)
library(survival)
library(Rcpp)
library(RcppEigen)
library(glmnet)
library(riskRegression)
library(prodlim)
library(survC1)

sourceCpp("group_lasso_gd/grplasso_survrrr.cpp")
source("group_lasso_gd/grplasso_survrrr.R")
sourceCpp("lasso_gd/lasso_survrrr.cpp")
source("lasso_gd/lasso_survrrr.R")
sourceCpp("ridge_gd/ridge_survrrr.cpp")
source("ridge_gd/ridge_survrrr.R")
sourceCpp("mrcox/aligned.cpp")
source("mrcox/mrcox.R")
source("lasso_cc/survRRR_functions.R")
source("simdata_function_and_performance_functions.R")


# Load NEO and select post-meal metabolites

data <- read_dta("NEO_analyse_PPXX_metabolytes_morbidity_RodriquezGirondo_2025-02-18.dta")
data <- data %>% mutate(id = row_number())

metabolites <- data[, 64:ncol(data)]
pre_meal <- dplyr::select(metabolites, matches("_1"))   # pre-meal sample
post_meal <- dplyr::select(metabolites, matches("_3"))   # post-meal sample

## STATA labels (human-readable description)
metabolites_labels = lapply(metabolites, function(x) attributes(x)$label) %>% as.character()
metabolites_labels_unique = metabolites_labels %>%
  map_chr(~ sub("^\\S+\\s+", "", .)) %>%
  unique()

cleaned_labels = sapply(metabolites_labels_unique, function(text) {
  result = sub("\\ - .*$", "", text)
  result = gsub("\\S*\\.\\S*", "", result)
  result = trimws(result)
  return(result)
}) %>% unique()

cleaned_labels <- cleaned_labels[!cleaned_labels %in% c("NULL", NA, "")]

# data preprocessing
cat(sprintf("Raw : N = %d, p = %d metabolite columns\n",
            nrow(post_meal), ncol(post_meal)))
# Raw : N = 6671, p = 229 metabolite columns

## age >= 50
keep_age <- !is.na(data$leeftijd) & data$leeftijd >= 50
post_meal <- post_meal[keep_age, ]
data     <- data[keep_age, ]

## subjects with > 1 % missing across the metabolite block
p0      <- ncol(post_meal)
row_pct <- rowSums(is.na(post_meal)) / p0
keep    <- row_pct <= 0.01
post_meal <- post_meal[keep, ];  data <- data[keep, ]

## metabolites with > 1 % missing across the surviving subjects
col_pct  <- colSums(is.na(post_meal)) / nrow(post_meal)
post_meal <- post_meal[, col_pct <= 0.01]
cleaned_labels <- cleaned_labels[col_pct <= 0.01]

## subjects with any remaining NA
keep    <- complete.cases(post_meal)
post_meal <- post_meal[keep, ]; data <- data[keep, ]

cat(sprintf("Cleaned: N = %d, p = %d metabolites\n",
            nrow(post_meal), ncol(post_meal)))
# Cleaned: N = 3674, p = 224 metabolites

# Build the outcomes 
# outcome dataframe for CVA
outcome_CVA <- data %>%
  dplyr::select(id, visitdd, CVA2_Date_1, CVA2_Date_2, CVA2_Date_3, einddatum2) %>%
  mutate(
    # Get the earliest valid date greater than visitdd
    valid_dates = pmin(CVA2_Date_1, CVA2_Date_2, CVA2_Date_3, na.rm = TRUE),
    # Outcome is the earliest non NA date if it's greater than visitdd, otherwise 0
    outcome_CVA = ifelse(!is.na(valid_dates) & valid_dates > visitdd, 1, 0),
    # Store the valid date in outcome_time (if greater than visitdd, else NA)
    outcome_CVA_time = ifelse(!is.na(valid_dates) & valid_dates > visitdd, valid_dates, NA),
    # Censored is 1 if (all dates are NA) OR (the only valid date is <= visitdd)
    censored_CVA = ifelse(is.na(CVA2_Date_1) & is.na(CVA2_Date_2) & is.na(CVA2_Date_3), 
                          1, 
                          ifelse(!is.na(valid_dates) & valid_dates <= visitdd, 1, 0))
  ) %>% 
  dplyr::select(-valid_dates) 

# outcome dataframe for CVD
outcome_CVD <- data %>%
  dplyr::select(id, visitdd, CVD2_inc, CVD2_date_first, einddatum2) %>%
  mutate(
    outcome_CVD = ifelse(CVD2_inc == 1, 1, 0),
    outcome_CVD_time = ifelse(CVD2_inc == 1, CVD2_date_first, 0),
    censored_CVD = ifelse(CVD2_inc == 0, 1, 0)
  )

# outcome dataframe for diabetes
outcome_diabetes <- data %>%
  dplyr::select(id, visitdd, diab_prev, diabetes2, diabetes2_date, einddatum2) %>%
  mutate(outcome_DIA_date = ifelse(diab_prev == 0 & diabetes2 == 1, diabetes2_date, 0),
         outcome_DIA = ifelse(outcome_DIA_date > 0, 1, 0),
         censored_DIA = ifelse(outcome_DIA == 1, 0, 1)
  )


# Outcome dataframe for Stroke
outcome_stroke <- data %>%
  dplyr::select(id, visitdd, Stroke2_date_1, Stroke2_date_2, Stroke2_date_3, Stroke2_inc, einddatum2) %>%
  mutate(
    valid_dates = pmin(Stroke2_date_1, Stroke2_date_2, Stroke2_date_3, na.rm = TRUE),
    outcome_STR = ifelse(Stroke2_inc == 1 | (!is.na(valid_dates) & valid_dates > visitdd), 1, 0),
    outcome_STR_time = ifelse(!is.na(valid_dates) & valid_dates > visitdd, valid_dates, NA),
    censored_STR = ifelse(is.na(Stroke2_date_1) & is.na(Stroke2_date_2) & is.na(Stroke2_date_3), 
                          1, 
                          ifelse(!is.na(valid_dates) & valid_dates <= visitdd, 1, 0))
  ) %>%
  dplyr::select(-c(valid_dates))

# outcome dataframe for death
outcome_death <- data %>%
  dplyr::select(id, visitdd, einddatum2, eind2) %>% #4 = death
  mutate(outcome_DTH = ifelse(eind2 == 4, 1, 0),
         outcome_DTH_date = ifelse(outcome_DTH == 1, einddatum2, 0),
         censored_DTH = ifelse(outcome_DTH == 1, 0, 1))

# outcome dataframe for myocardial infarction (MI)
outcome_MI <- data %>%
  dplyr::select(id, visitdd, MI2_Date_1, einddatum2) %>%
  mutate(
    outcome_MI      = ifelse(!is.na(MI2_Date_1) & MI2_Date_1 > visitdd, 1, 0),
    outcome_MI_time = ifelse(!is.na(MI2_Date_1) & MI2_Date_1 > visitdd, MI2_Date_1, NA),
    censored_MI     = ifelse(outcome_MI == 1, 0, 1)
  )

# outcome dataframe for transient ischemic attack (TIA)
outcome_TIA <- data %>%
  dplyr::select(id, visitdd, TIA2_date_1, einddatum2) %>%
  mutate(
    outcome_TIA      = ifelse(!is.na(TIA2_date_1) & TIA2_date_1 > visitdd, 1, 0),
    outcome_TIA_time = ifelse(!is.na(TIA2_date_1) & TIA2_date_1 > visitdd, TIA2_date_1, NA),
    censored_TIA     = ifelse(outcome_TIA == 1, 0, 1)
  )


# ------ Build survival data for the outcomes ------ 

age_at_visit <- data$leeftijd   # age in years at study visit

build_outcome_age <- function(d, age, event_flag, event_date_num,
                              extra_exclude = rep(FALSE, nrow(d))) {
  ev_flag <- as.numeric(event_flag)
  ev_date <- as.numeric(event_date_num)
  end_date_d <- as.numeric(d$einddatum2)
  visit_d <- as.numeric(d$visitdd)
  
  flag_known <- !is.na(ev_flag)
  has <- flag_known & ev_flag == 1 & !is.na(ev_date) & ev_date > 0
  
  end_raw <- ifelse(has, ev_date, end_date_d)
  end <- as.Date(end_raw, origin = "1970-01-01")
  visit <- as.Date(visit_d, origin = "1970-01-01")
  dur <- as.numeric(difftime(end, visit, units = "days")) / 365.25
  
  prev <- has & dur < 0                  
  valid <- flag_known & !is.na(dur) & dur > 0 & !extra_exclude & !prev
  
  list(entry  = age,
       time   = age + dur,  
       status = as.integer(has & !prev),
       valid  = valid)
}

## Build all 7 outcomes
all_outcomes <- list(
  CVA = build_outcome_age(data, age_at_visit,
                          event_flag     = outcome_CVA$outcome_CVA,
                          event_date_num = outcome_CVA$outcome_CVA_time),
  CVD = build_outcome_age(data, age_at_visit,
                          event_flag     = outcome_CVD$outcome_CVD,
                          event_date_num = outcome_CVD$outcome_CVD_time),
  DIA = build_outcome_age(data, age_at_visit,
                          event_flag     = outcome_diabetes$outcome_DIA,
                          event_date_num = outcome_diabetes$outcome_DIA_date,
                          extra_exclude  = data$diab_prev == 1),
  STR = build_outcome_age(data, age_at_visit,
                          event_flag     = outcome_stroke$outcome_STR,
                          event_date_num = outcome_stroke$outcome_STR_time),
  DTH = build_outcome_age(data, age_at_visit,
                          event_flag     = outcome_death$outcome_DTH,
                          event_date_num = outcome_death$outcome_DTH_date),
  MI  = build_outcome_age(data, age_at_visit,
                          event_flag     = outcome_MI$outcome_MI,
                          event_date_num = outcome_MI$outcome_MI_time),
  TIA = build_outcome_age(data, age_at_visit,
                          event_flag     = outcome_TIA$outcome_TIA,
                          event_date_num = outcome_TIA$outcome_TIA_time)
)

# pick outcomes
selected_outcomes <- c("CVA", "DIA", "MI", "DTH", "CVD", "TIA")

stopifnot(all(selected_outcomes %in% names(all_outcomes)))

keep <- Reduce(`&`,
               lapply(selected_outcomes, function(o) all_outcomes[[o]]$valid))

cat(sprintf("Selected outcomes (%d): %s\n",
            length(selected_outcomes),
            paste(selected_outcomes, collapse = ", ")))
cat(sprintf("Subjects valid for ALL selected outcomes: N = %d\n", sum(keep)))
# Subjects valid for ALL selected outcomes: N = 3308

## Final design matrix
X <- as.matrix(post_meal[keep, ])
X <- scale(X, center = TRUE, scale = TRUE)
storage.mode(X) <- "double"

## Per-outcome lists
entry_list  <- setNames(
  lapply(selected_outcomes, function(o) all_outcomes[[o]]$entry [keep]),
  selected_outcomes)
y_list      <- setNames(
  lapply(selected_outcomes, function(o) all_outcomes[[o]]$time  [keep]),
  selected_outcomes)
status_list <- setNames(
  lapply(selected_outcomes, function(o) all_outcomes[[o]]$status[keep]),
  selected_outcomes)

## Follow-up scale (years from visit)
dur_list <- setNames(
  lapply(seq_along(y_list), function(k) y_list[[k]] - entry_list[[k]]),
  selected_outcomes)

cat("Event counts:\n")
for (o in selected_outcomes)
  cat(sprintf("  %-3s : %d\n", o, sum(status_list[[o]])))
cat(sprintf("Total N: %d\n", nrow(X)))

#   CVA : 85
#   DIA : 179
#   MI  : 62
#   DTH : 37
#   CVD : 139
#   TIA : 48


# ------ Group labels via Nightingale biomarker categories ------
classify_nightingale <- function(label) {
  s <- tolower(label %||% "")
  s <- gsub("[-_]", " ", s)
  
  ## --- (1) lipoprotein subclasses ---
  if (grepl("very large vldl|xxl vldl",  s)) return("VVLDL")
  if (grepl("very small vldl|xs vldl",   s)) return("VSVLDL")
  if (grepl("large vldl|xl vldl",        s)) return("LVLDL")
  if (grepl("medium vldl|m vldl",        s)) return("MVLDL")
  if (grepl("small vldl|s vldl",         s)) return("SVLDL")
  if (grepl("\\bidl\\b",                 s)) return("IDL")
  if (grepl("large ldl",                 s)) return("LLDL")
  if (grepl("medium ldl",                s)) return("MLDL")
  if (grepl("small ldl",                 s)) return("SLDL")
  if (grepl("very large hdl|xl hdl",     s)) return("VLHDL")
  if (grepl("large hdl",                 s)) return("LHDL")
  if (grepl("medium hdl",                s)) return("MHDL")
  if (grepl("small hdl",                 s)) return("SHDL")
  
  ## --- (2) apolipoproteins ---
  if (grepl("apolipoprotein|apo b|apo a", s)) return("Apolipoprotein")
  
  ## --- (3) fatty acid RATIOS ---
  if (grepl("ratio.*fatty|fatty.*ratio|fa %|percentage.*fatty|degree of unsat",
            s)) return("FA_ratio")
  
  ## --- (4) fatty acids ---
  if (grepl("fatty acid|omega|pufa|mufa|sfa|linoleic|docosahexaenoic|saturat|unsaturat",
            s)) return("Fatty_acid")
  
  ## --- (5) amino acids ---
  if (grepl("isoleucine|leucine|valine|branched", s))     return("BCAA")
  if (grepl("tyrosine|phenylalanine|aromatic",    s))     return("AAA")
  if (grepl("alanine|glutamine|glycine|histidine", s))    return("Other_AA")
  
  ## --- (6) glycolysis & related small molecules ---
  if (grepl("glucose|lactate|pyruvate|citrate|glycerol", s)) return("Glycolysis")
  
  ## --- (7) ketones ---
  if (grepl("hydroxybutyrate|acetoacetate|acetate|ketone", s)) return("Ketone")
  
  ## --- (8) inflammation ---
  if (grepl("glycoprotein|glyca", s)) return("Inflammation")
  
  ## --- (9) fluid balance / other small molecules ---
  if (grepl("creatinine|albumin", s)) return("Fluid_balance")
  
  ## --- (10) whole-serum lipid composition ---
  if (grepl("triglyceride",                                  s)) return("Triglyceride_total")
  if (grepl("phospholipid|phosphatidylcholine|sphingomyelin|choline|phosphoglyceride",
            s)) return("Phospholipid_total")
  if (grepl("cholesteryl ester",                             s)) return("CholEster_total")
  if (grepl("free cholesterol",                              s)) return("FreeChol_total")
  if (grepl("cholesterol",                                   s)) return("Cholesterol_total")
  if (grepl("total lipid|particle",                          s)) return("TotalLipid_total")
  
  return("Other")
}

metab_groups <- vapply(cleaned_labels, classify_nightingale, character(1))
group_labels <- as.integer(factor(metab_groups))      # 1..n_groups, integer
n_groups     <- max(group_labels)

cat(sprintf("\nNumber of groups: %d  (%d predictors)\n", n_groups, ncol(X)))
# Number of groups: 28  (224 predictors)
print(table(metab_groups))
#metab_groups
# AAA                  Apolipoprotein       BCAA              Cholesterol_total   FA_ratio           Fatty_acid 
# 2                    3                    3                 10                  8                  9
# Fluid_balance        FreeChol_total       Glycolysis        IDL                 Inflammation       Ketone 
# 2                    2                    3                 12                  1                  2 
# LHDL                 LLDL                 LVLDL             MHDL                MLDL               MVLDL 
# 12                   12                   19                12                  12                 12 
# Other_AA             Phospholipid_total   SHDL              SLDL                SVLDL              TotalLipid_total 
# 3                    4                    12                12                  12                 3 
# Triglyceride_total   VLHDL                VSVLDL            VVLDL 
# 6                    12                   12                12 

## any "Other" left?
if ("Other" %in% metab_groups) {
  cat("\nUnclassified metabolites (group = 'Other'):\n")
  print(cleaned_labels[metab_groups == "Other"])
}


#  ------ Predictive performance: 10-fold CV ------

# configuration
set.seed(7)
N <- nrow(X)
n_outcomes <- length(y_list)
n_folds    <- 10

methods_to_run <- c("rrr_grplasso", "rrr_lasso", "rrr_ridge", "pen", "mrcox")
rank_vector    <- c(1, 2, 3)

lambda_grids <- list(
  rrr_grplasso = seq(0.01, 0.25, length.out = 11),
  rrr_lasso    = seq(0.01, 0.25, length.out = 11),
  rrr_ridge    = seq(0.2,  1.0, length.out = 11),
  mrcox        = seq(0.01, 0.3, length.out = 11),
  pen          = seq(0.001, 0.05, length.out = 11)
)

fold_id <- sample(rep(seq_len(n_folds), length.out = N))
needs_followup <- function(method) identical(method, "mrcox")



# Uno C-index (IPCW); no entry support
compute_uno <- function(lp, time, status) {
  ev <- time[status == 1]
  if (length(ev) < 2) return(NA_real_)
  tau <- quantile(ev, 0.8, na.rm = TRUE)
  dat <- data.frame(time = time, event = status, risk = lp)
  tryCatch(
    survC1::Est.Cval(dat, tau = tau, nofit = TRUE)$Dhat,
    error = function(e) NA_real_)
}

# Harrell C-index;  Supports left truncation
compute_harrell <- function(lp, time, status, entry = NULL) {
  if (sum(status) < 2) return(NA_real_)
  df <- data.frame(time = time, status = status, lp = lp)
  fit <- tryCatch({
    if (!is.null(entry)) {
      df$entry <- entry
      survival::concordance(survival::Surv(entry, time, status) ~ lp,
                            data = df, reverse = TRUE)
    } else {
      survival::concordance(survival::Surv(time, status) ~ lp,
                            data = df, reverse = TRUE)
    }
  }, error = function(e) NULL)
  if (is.null(fit)) return(NA_real_)
  as.numeric(fit$concordance)
}

# Orient each column of B_hat so that "higher lp = higher risk" on TRAINING
orient_B_by_training <- function(B_hat, X_tr, y_tr, s_tr, e_tr = NULL) {
  K <- ncol(B_hat)
  eta_tr <- X_tr %*% B_hat
  for (k in seq_len(K)) {
    lp_k <- as.numeric(eta_tr[, k])
    if (stats::sd(lp_k) < 1e-12) next                # constant lp: skip
    df <- data.frame(time = y_tr[[k]], status = s_tr[[k]], lp = lp_k)
    fit_chk <- tryCatch({
      if (!is.null(e_tr)) {
        df$entry <- e_tr[[k]]
        survival::coxph(survival::Surv(entry, time, status) ~ lp,
                        data = df, method = "breslow",
                        control = survival::coxph.control(iter.max = 20))
      } else {
        survival::coxph(survival::Surv(time, status) ~ lp,
                        data = df, method = "breslow",
                        control = survival::coxph.control(iter.max = 20))
      }
    }, error = function(e) NULL)
    if (!is.null(fit_chk) && length(coef(fit_chk)) >= 1 &&
        is.finite(coef(fit_chk)[1]) && coef(fit_chk)[1] < 0) {
      B_hat[, k] <- -B_hat[, k]
    }
  }
  B_hat
}


# rank + lambda selection on a single dataset
select_R_and_lambda <- function(X, y_list, status_list, entry_list,
                                rank_vector, lambda_vector,
                                method, group_labels = NULL,
                                n_folds_inner = 5,
                                dlong_train = NULL,
                                data_wide_train = NULL,
                                early_stop_tol = 1e-3) {
  best <- list(R = NA, lambda = NA, cve = Inf, per_rank = list())
  pred_names <- paste0("x", seq_len(ncol(X)))
  k_out <- length(y_list)
  
  for (R in rank_vector) {
    cv <- tryCatch({
      if (method == "pen") {
        if (is.null(dlong_train) || is.null(data_wide_train))
          stop("pen needs dlong_train and data_wide_train")
        select_lambda_pen_vvh(
          dlong         = dlong_train,
          data_wide     = data_wide_train,
          R             = R,
          lambda_vector = lambda_vector,
          pred_names    = pred_names,
          k             = k_out,
          entry_list    = entry_list,
          n_folds       = n_folds_inner)
      } else {
        select_lambda_vvh(
          X = X, y_list = y_list, status_list = status_list,
          entry_list = entry_list, R = R,
          lambda_vector = lambda_vector,
          group_labels = group_labels,
          method = method, n_folds = n_folds_inner)
      }
    }, error = function(e) {
      cat(sprintf("  [%s R=%d] inner CV failed: %s\n", method, R, e$message))
      NULL
    })
    
    if (is.null(cv)) next
    
    best$per_rank[[as.character(R)]] <- cv
    cve_min <- min(cv$CVE_mean, na.rm = TRUE)
    if (!is.finite(cve_min)) next
    
    if (cve_min < best$cve) {
      ## new best
      best$R      <- R
      best$lambda <- cv$best_lambda
      best$cve    <- cve_min
      cat(sprintf("    R = %d : CVE = %.4f  [new best]\n", R, cve_min))
    } else if (cve_min > best$cve + early_stop_tol) {
      ## strictly worse beyond noise tolerance -> stop searching higher R
      cat(sprintf("    R = %d : CVE = %.4f  (worse than best R = %d, %.4f); stopping rank search.\n",
                  R, cve_min, best$R, best$cve))
      break
    } else {
      ## within tolerance -- record but keep going (could still drop later)
      cat(sprintf("    R = %d : CVE = %.4f  (within tol of best)\n", R, cve_min))
    }
  }
  best
}


# (only if pen is used) build dlong / wide
neo_wide <- dlong_neo <- NULL
if ("pen" %in% methods_to_run) {
  pred_names <- paste0("x", seq_len(ncol(X)))
  neo_wide <- data.frame(id = seq_len(N))
  for (k in seq_len(n_outcomes)) {
    neo_wide[[paste0("t", k)]] <- y_list[[k]]
    neo_wide[[paste0("d", k)]] <- status_list[[k]]
    neo_wide[[paste0("e", k)]] <- entry_list[[k]]
  }
  for (j in seq_len(ncol(X)))
    neo_wide[[pred_names[j]]] <- X[, j]
  dlong_neo <- build_dlong(neo_wide, k = n_outcomes, entry_list = entry_list)
}

### ------ JOINT ------
# Step A: pick (R, lambda) on full data per method 
selected_params <- list()
cat("\n========== Step A: inner CV for (R, lambda) ==========\n")
for (method in methods_to_run) {
  cat(sprintf("\n--- %s ---\n", method))
  if (needs_followup(method)) {
    y_for_fit <- dur_list; e_for_fit <- NULL
  } else {
    y_for_fit <- y_list;   e_for_fit <- entry_list
  }
  sel <- select_R_and_lambda(
    X = X, y_list = y_for_fit, status_list = status_list,
    entry_list = e_for_fit,
    rank_vector   = rank_vector,
    lambda_vector = lambda_grids[[method]],
    method = method, group_labels = group_labels,
    n_folds_inner = 5,
    dlong_train     = dlong_neo,
    data_wide_train = neo_wide)
  selected_params[[method]] <- sel
  cat(sprintf("  -> selected R = %s, lambda = %.4g  (CVE = %.4f)\n",
              sel$R, sel$lambda, sel$cve))
}


# Step B: 10-fold loop, store OOF lp 
lp_pooled <- array(NA_real_,
                   dim = c(N, length(methods_to_run), n_outcomes),
                   dimnames = list(NULL, methods_to_run, names(y_list)))

cat("\n========== Step B: 10-fold OOF predictions ==========\n")
for (fold in seq_len(n_folds)) {
  cat(sprintf("\n--- Outer fold %d/%d ---\n", fold, n_folds))
  tr <- which(fold_id != fold)
  te <- which(fold_id == fold)
  
  X_tr <- X[tr, , drop = FALSE]; X_te <- X[te, , drop = FALSE]
  y_tr   <- lapply(y_list,      `[`, tr)
  s_tr   <- lapply(status_list, `[`, tr)
  e_tr   <- lapply(entry_list,  `[`, tr)
  dur_tr <- lapply(dur_list,    `[`, tr)
  
  dlong_tr <- NULL
  if ("pen" %in% methods_to_run) {
    wide_tr  <- neo_wide[tr, , drop = FALSE]
    dlong_tr <- build_dlong(wide_tr, k = n_outcomes, entry_list = e_tr)
  }
  
  for (method in methods_to_run) {
    sel <- selected_params[[method]]
    if (is.na(sel$R) || is.na(sel$lambda)) {
      cat(sprintf("  [%s] no selected (R, lambda); skipping.\n", method))
      next
    }
    
    if (needs_followup(method)) {
      y_fit <- dur_tr; e_fit <- NULL
    } else {
      y_fit <- y_tr;   e_fit <- e_tr
    }
    
    fit <- tryCatch(
      fit_one_model(
        X = X_tr, y_list = y_fit, status_list = s_tr,
        entry_list = e_fit,
        dlong  = if (method == "pen") dlong_tr else NULL,
        method = method, r = sel$R, lambda = sel$lambda,
        group_labels = group_labels),
      error = function(e) {
        cat(sprintf("  [%s] fit failed: %s\n", method, e$message)); NULL
      })
    
    if (is.null(fit) || is.null(fit$B_hat) || !all(is.finite(fit$B_hat))) next

    fit$B_hat <- orient_B_by_training(
      fit$B_hat, X_tr,
      y_tr = if (needs_followup(method)) dur_tr else y_tr,
      s_tr = s_tr,
      e_tr = if (needs_followup(method)) NULL    else e_tr)

    eta_te <- X_te %*% fit$B_hat
    for (k in seq_len(n_outcomes)) lp_pooled[te, method, k] <- eta_te[, k]
    cat(sprintf("  [%s] R=%s lam=%.3g  done.\n",
                method, sel$R, sel$lambda))
  }
}

# Uno + Harrell C-index per (method, outcome)
cindex_uno <- matrix(NA_real_, length(methods_to_run), n_outcomes,
                     dimnames = list(methods_to_run, names(y_list)))
cindex_har <- matrix(NA_real_, length(methods_to_run), n_outcomes,
                     dimnames = list(methods_to_run, names(y_list)))

for (method in methods_to_run) {
  for (k in seq_len(n_outcomes)) {
    lp <- lp_pooled[, method, k]
    ok <- !is.na(lp)
    if (sum(ok) < 2) next
    cindex_uno[method, k] <- compute_uno(
      lp[ok], dur_list[[k]][ok], status_list[[k]][ok])
    if (needs_followup(method)) {
      cindex_har[method, k] <- compute_harrell(
        lp[ok], dur_list[[k]][ok], status_list[[k]][ok], entry = NULL)
    } else {
      cindex_har[method, k] <- compute_harrell(
        lp[ok], y_list[[k]][ok], status_list[[k]][ok],
        entry = entry_list[[k]][ok])
    }
  }
}

### ------ UNIVARIATE ------
uni_methods <- c("glmnet_lasso", "glmnet_ridge")
alpha_map   <- c(glmnet_lasso = 1, glmnet_ridge = 0)   

# Step A_uni: cv.glmnet on full data per (outcome, method)
selected_uni <- vector("list", n_outcomes)
names(selected_uni) <- names(y_list)

cat("\n========== Step A_uni: per-outcome cv.glmnet lambda selection ==========\n")
set.seed(7)
for (k in seq_len(n_outcomes)) {
  out_name <- names(y_list)[k]
  cat(sprintf("\n--- Outcome: %s (events = %d) ---\n",
              out_name, sum(status_list[[k]])))
  selected_uni[[k]] <- list()
  
  y_surv <- survival::Surv(entry_list[[k]], y_list[[k]], status_list[[k]])
  
  for (m in uni_methods) {
    cv_fit <- tryCatch(
      glmnet::cv.glmnet(
        x            = X,
        y            = y_surv,
        family       = "cox",
        alpha        = alpha_map[[m]],   
        nfolds       = 5,                
        type.measure = "deviance",      
        standardize  = FALSE
      ),
      error = function(e) {
        cat(sprintf("  [%s] cv.glmnet failed: %s\n", m, e$message))
        NULL
      })
    
    if (is.null(cv_fit)) {
      selected_uni[[k]][[m]] <- list(lambda = NA_real_, cve = NA_real_)
      next
    }
    
    selected_uni[[k]][[m]] <- list(
      lambda = cv_fit$lambda.min,
      cve    = min(cv_fit$cvm, na.rm = TRUE)
    )
    cat(sprintf("  %-13s lambda.min = %.4g  (cv deviance = %.4f)\n",
                m, cv_fit$lambda.min, min(cv_fit$cvm, na.rm = TRUE)))
  }
}

# Step B_uni: 10-fold OOF lp 
lp_uni_pooled <- array(NA_real_,
                       dim = c(N, length(uni_methods), n_outcomes),
                       dimnames = list(NULL, uni_methods, names(y_list)))

cat("\n========== Step B_uni: 10-fold univariate glmnet OOF ==========\n")
for (fold in seq_len(n_folds)) {
  cat(sprintf("\n--- Outer fold %d/%d ---\n", fold, n_folds))
  tr <- which(fold_id != fold); te <- which(fold_id == fold)
  X_tr <- X[tr, , drop = FALSE]; X_te <- X[te, , drop = FALSE]
  
  for (k in seq_len(n_outcomes)) {
    y_surv_tr <- survival::Surv(entry_list[[k]][tr],
                                y_list[[k]][tr],
                                status_list[[k]][tr])
    
    for (m in uni_methods) {
      sel <- selected_uni[[k]][[m]]
      if (is.na(sel$lambda)) next
      
      fit <- tryCatch(
        glmnet::glmnet(
          x           = X_tr,
          y           = y_surv_tr,
          family      = "cox",
          alpha       = alpha_map[[m]],
          lambda      = sel$lambda,
          standardize = FALSE
        ),
        error = function(e) {
          cat(sprintf("  [%s/%s] glmnet refit failed: %s\n",
                      names(y_list)[k], m, e$message)); NULL
        })
      
      if (is.null(fit)) next
      
      ## type="link" returns the linear predictor X_te %*% beta_hat
      lp_uni_pooled[te, m, k] <- as.numeric(
        predict(fit, newx = X_te, type = "link"))
    }
  }
}

# Step C_uni: Uno + Harrell C-index per (method, outcome)
cindex_uno_uni <- matrix(NA_real_, length(uni_methods), n_outcomes,
                         dimnames = list(uni_methods, names(y_list)))
cindex_har_uni <- matrix(NA_real_, length(uni_methods), n_outcomes,
                         dimnames = list(uni_methods, names(y_list)))

for (m in uni_methods) {
  for (k in seq_len(n_outcomes)) {
    lp <- lp_uni_pooled[, m, k]
    ok <- !is.na(lp)
    if (sum(ok) < 2) next
    cindex_uno_uni[m, k] <- compute_uno(
      lp[ok], dur_list[[k]][ok], status_list[[k]][ok])
    cindex_har_uni[m, k] <- compute_harrell(
      lp[ok], y_list[[k]][ok], status_list[[k]][ok],
      entry = entry_list[[k]][ok])
  }
}

# ------ Summary ------
cat("\n=== JOINT (multi-outcome) Pooled OOF Uno C-index (follow-up scale) ===\n")
print(round(cindex_uno, 3))

cat("\n=== JOINT Pooled OOF Harrell C-index (age scale; follow-up for mrcox) ===\n")
print(round(cindex_har, 3))

cat("\n=== Selected (R, lambda) per method (JOINT, Step A on full data) ===\n")
sel_tbl <- do.call(rbind, lapply(methods_to_run, function(m) {
  s <- selected_params[[m]]
  data.frame(method = m, R = s$R, lambda = signif(s$lambda, 4),
             CVE = signif(s$cve, 4))
}))
print(sel_tbl, row.names = FALSE)

cat("\n=== UNIVARIATE glmnet Pooled OOF Uno C-index (follow-up scale) ===\n")
print(round(cindex_uno_uni, 3))

cat("\n=== UNIVARIATE glmnet Pooled OOF Harrell C-index (age scale) ===\n")
print(round(cindex_har_uni, 3))

cat("\n=== Selected lambda per (outcome, method) (UNIVARIATE glmnet) ===\n")
sel_uni_tbl <- do.call(rbind, lapply(seq_len(n_outcomes), function(k) {
  do.call(rbind, lapply(uni_methods, function(m) {
    s <- selected_uni[[k]][[m]]
    data.frame(outcome = names(y_list)[k], method = m,
               lambda  = signif(s$lambda, 4),
               cv_dev  = signif(s$cve, 4))
  }))
}))
print(sel_uni_tbl, row.names = FALSE)

#Comparison
cat("\n=== Joint methods MINUS glmnet-lasso baseline (Uno) ===\n")
print(round(sweep(cindex_uno, 2, cindex_uno_uni["glmnet_lasso", ], "-"), 3))

cat("\n=== Joint methods MINUS glmnet-ridge baseline (Uno) ===\n")
print(round(sweep(cindex_uno, 2, cindex_uno_uni["glmnet_ridge", ], "-"), 3))

cat("\n=== Joint methods MINUS glmnet-lasso baseline (Harrell) ===\n")
print(round(sweep(cindex_har, 2, cindex_har_uni["glmnet_lasso", ], "-"), 3))

cat("\n=== Joint methods MINUS glmnet-ridge baseline (Harrell) ===\n")
print(round(sweep(cindex_har, 2, cindex_har_uni["glmnet_ridge", ], "-"), 3))


