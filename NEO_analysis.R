# ==== NEO reduced-rank survival analysis =====================================
#   - predictors: post-meal Nightingale metabolites
#   - outcomes  : CVA, CVD, DIA, STR, DTH, MI, TIA

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


# ============ Load NEO data ====================

data <- read_dta(
  "NEO_analyse_PPXX_metabolytes_morbidity_RodriquezGirondo_2025-02-18.dta"
)

data <- data %>%
  mutate(id = row_number())


# ================= Extract metabolite matrix ===============================
metabolites <- data[, 64:ncol(data)]

post_meal <- metabolites %>%
  dplyr::select(matches("_3"))


orig_colnames <- colnames(post_meal) # retain original colnames


# ================= Extract STATA labels ======================

metabolites_labels <- lapply(post_meal, function(x) {
  attributes(x)$label
}) %>% as.character()


# =========== Clean labels ==============
cleaned_labels <- sapply(metabolites_labels, function(text) {
  
  if (is.na(text) || text == "NULL") return(NA)
  result <- sub("^\\S+\\s+", "", text)
  result <- sub("\\ - .*$", "", result)
  result <- gsub("\\S*\\.\\S*", "", result)
  result <- gsub("\\b\\d+\\s*:\\s*\\d+\\b", "", result)
  result <- gsub(";", "", result)
  result <- gsub(",", "", result)
  result <- gsub("\\s+", " ", result)
  trimws(result)
  
})

cleaned_labels[cleaned_labels == ""] <- NA


# Initial dimensions

cat(sprintf(
  "Raw data: N = %d, p = %d metabolites\n",
  nrow(post_meal),
  ncol(post_meal)
))
# Raw data: N = 6671, p = 229 metabolites

# ================= Data cleaning ====================

# 1. age >= 50

keep_age <- !is.na(data$leeftijd) & data$leeftijd >= 50
data <- data[keep_age, ]
post_meal <- post_meal[keep_age, ]


# 2. delete subjects with missing rate of metabolite > 1% 

row_missing_pct <- rowSums(is.na(post_meal)) / ncol(post_meal)
keep_subject <- row_missing_pct <= 0.01
data <- data[keep_subject, ]
post_meal <- post_meal[keep_subject, ]


# 3. delete predictors with missing rate of metabolite > 1% 

col_missing_pct <- colSums(is.na(post_meal)) / nrow(post_meal)
keep_metabolite <- col_missing_pct <= 0.01
post_meal <- post_meal[, keep_metabolite]
cleaned_labels <- cleaned_labels[keep_metabolite]



# 4. delete remaining NA

keep_complete <- complete.cases(post_meal)
data <- data[keep_complete, ]
post_meal <- post_meal[keep_complete, ]


cat(sprintf(
  "After missing-value cleaning: N = %d, p = %d\n",
  nrow(post_meal),
  ncol(post_meal)
))
# After missing-value cleaning: N = 3674, p = 224

# ========================== Group labels via Nightingale biomarker categories =======================================
grp_def <- list(
  cholesterol = c(seq(11, 95, by = 7), 173:178),
  free_cholesterol = c(seq(13, 97, by = 7), 180),
  total_lipid = seq(9, 93, by = 7),
  concentration = seq(8, 92, by = 7),
  triglycerides = c(seq(14, 98, by = 7), 182:184),
  phospholipids = c(seq(10, 94, by = 7), 187),
  cholesterol_ester = c(seq(12, 96, by = 7), 179),
  other_lipid = c(190, 191),
  apolipoprotein = c(192, 193),
  fatty_acid = seq(195, 205)[-2][-5],
  amino_acid = 217:224
)


# original colnames -> group mapping
grp_colnames <- lapply(grp_def, function(idx) {
  orig_colnames[idx]
})

grouped_cols <- unique(unlist(grp_colnames))


# only keep metabolites in the group

keep_grp <- colnames(post_meal) %in% grouped_cols

post_meal <- post_meal[, keep_grp]
cleaned_labels <- cleaned_labels[keep_grp]


cat(sprintf(
  "After group filtering: p = %d metabolites\n",
  ncol(post_meal)
))
# After group filtering: p = 124 metabolites

# ==================== Construct group_labels ===================

group_labels <- rep(NA_integer_, ncol(post_meal))

for (g in seq_along(grp_def)) {
  grp_name <- names(grp_def)[g]
  idx <- which(
    colnames(post_meal) %in% grp_colnames[[grp_name]]
  )
  group_labels[idx] <- g
}


# ======================= check group membership ==================================

cat("\n========== Group membership ==========\n")
group_names <- rep(NA_character_, ncol(post_meal)) 
for (grp in names(grp_def)) { 
  idx <- which(colnames(post_meal) %in% grp_colnames[[grp]]) 
  group_names[idx] <- grp }

for (grp in names(grp_def)) {
  idx <- which(group_names == grp)
  cat(sprintf(
    "\n[%s]  n=%d\n",
    grp,
    length(idx)
  ))
  
  cat(
    paste(cleaned_labels[idx], collapse = " | "),
    "\n"
  )
}

cat(sprintf(
  "\nUngrouped metabolites: %d\n",
  sum(is.na(group_names))
))


# =============================== Build outcomes ===============================

# CVA

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


# CVD

outcome_CVD <- data %>%
  dplyr::select(id, visitdd, CVD2_inc, CVD2_date_first, einddatum2) %>%
  mutate(
    outcome_CVD = ifelse(CVD2_inc == 1, 1, 0),
    outcome_CVD_time = ifelse(CVD2_inc == 1, CVD2_date_first, 0),
    censored_CVD = ifelse(CVD2_inc == 0, 1, 0)
  )


# Diabetes

outcome_diabetes <- data %>%
  dplyr::select(id, visitdd, diab_prev, diabetes2, diabetes2_date, einddatum2) %>%
  mutate(outcome_DIA_date = ifelse(diab_prev == 0 & diabetes2 == 1, diabetes2_date, 0),
         outcome_DIA = ifelse(outcome_DIA_date > 0, 1, 0),
         censored_DIA = ifelse(outcome_DIA == 1, 0, 1)
  )



# Stroke

outcome_stroke <- data %>%
  dplyr::select(id, visitdd, Stroke2_date_1, Stroke2_date_2, Stroke2_date_3, Stroke2_inc, einddatum2) %>%
  mutate(
    valid_dates = pmin(Stroke2_date_1, Stroke2_date_2, Stroke2_date_3, na.rm = TRUE),
    outcome_STR = ifelse(!is.na(valid_dates) & valid_dates > visitdd, 1, 0),
    outcome_STR_time = ifelse(!is.na(valid_dates) & valid_dates > visitdd, valid_dates, NA),
    censored_STR = ifelse(is.na(Stroke2_date_1) & is.na(Stroke2_date_2) & is.na(Stroke2_date_3), 
                          1, 
                          ifelse(!is.na(valid_dates) & valid_dates <= visitdd, 1, 0))
  ) %>%
  dplyr::select(-c(valid_dates))


# Death

outcome_death <- data %>%
  dplyr::select(id, visitdd, einddatum2, eind2) %>% #4 = death
  mutate(outcome_DTH = ifelse(eind2 == 4, 1, 0),
         outcome_DTH_date = ifelse(outcome_DTH == 1, einddatum2, 0),
         censored_DTH = ifelse(outcome_DTH == 1, 0, 1))


# MI

outcome_MI <- data %>%
  dplyr::select(id, visitdd, MI2_Date_1, einddatum2) %>%
  mutate(
    outcome_MI      = ifelse(!is.na(MI2_Date_1) & MI2_Date_1 > visitdd, 1, 0),
    outcome_MI_time = ifelse(!is.na(MI2_Date_1) & MI2_Date_1 > visitdd, MI2_Date_1, NA),
    censored_MI     = ifelse(outcome_MI == 1, 0, 1))


# TIA

outcome_TIA <- data %>%
  dplyr::select(id, visitdd, TIA2_date_1, einddatum2) %>%
  mutate(
    outcome_TIA      = ifelse(!is.na(TIA2_date_1) & TIA2_date_1 > visitdd, 1, 0),
    outcome_TIA_time = ifelse(!is.na(TIA2_date_1) & TIA2_date_1 > visitdd, TIA2_date_1, NA),
    censored_TIA     = ifelse(outcome_TIA == 1, 0, 1))


# ======================== Build survival outcomes ===============================

age_at_visit <- data$leeftijd


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
       time  = age + dur,  
       status = as.integer(has & !prev),
       valid = valid)
}


all_outcomes <- list(
  
  CVA = build_outcome_age(
    data,
    age_at_visit,
    outcome_CVA$outcome_CVA,
    outcome_CVA$outcome_CVA_time
  ),
  
  CVD = build_outcome_age(
    data,
    age_at_visit,
    outcome_CVD$outcome_CVD,
    outcome_CVD$outcome_CVD_time
  ),
  
  DIA = build_outcome_age(
    data,
    age_at_visit,
    outcome_diabetes$outcome_DIA,
    outcome_diabetes$outcome_DIA_date,
    extra_exclude = data$diab_prev == 1
  ),
  
  STR = build_outcome_age(
    data,
    age_at_visit,
    outcome_stroke$outcome_STR,
    outcome_stroke$outcome_STR_time
  ),
  
  DTH = build_outcome_age(
    data,
    age_at_visit,
    outcome_death$outcome_DTH,
    outcome_death$outcome_DTH_date
  ),
  
  MI = build_outcome_age(
    data,
    age_at_visit,
    outcome_MI$outcome_MI,
    outcome_MI$outcome_MI_time
  ),
  
  TIA = build_outcome_age(
    data,
    age_at_visit,
    outcome_TIA$outcome_TIA,
    outcome_TIA$outcome_TIA_time
  )
)


# =============== Select outcomes ===========================

selected_outcomes <- c(
  "CVA",
  "DIA",
  "MI",
  "DTH",
  "CVD",
  "TIA",
  "STR"
)

stopifnot(
  all(selected_outcomes %in% names(all_outcomes))
)


# Keep subjects valid for ALL selected outcomes

keep <- Reduce(
  `&`,
  lapply(selected_outcomes, function(o) {
    all_outcomes[[o]]$valid
  })
)

cat(sprintf(
  "Subjects valid for all outcomes: N = %d\n",
  sum(keep)
))
# Subjects valid for all outcomes: N = 3308

# ================== Final design matrix==========================

X <- as.matrix(post_meal[keep, ])


# Standardize metabolites

X <- scale(X, center = TRUE, scale = TRUE)

storage.mode(X) <- "double"


cat(sprintf(
  "Final design matrix: N = %d, p = %d\n",
  nrow(X),
  ncol(X)
))
# Final design matrix: N = 3308, p = 124

# ================== Build survival lists ===================

entry_list <- setNames(
  lapply(selected_outcomes, function(o) {
    all_outcomes[[o]]$entry[keep]
  }),
  selected_outcomes
)

y_list <- setNames(
  lapply(selected_outcomes, function(o) {
    all_outcomes[[o]]$time[keep]
  }),
  selected_outcomes
)

status_list <- setNames(
  lapply(selected_outcomes, function(o) {
    all_outcomes[[o]]$status[keep]
  }),
  selected_outcomes
)


# Follow-up duration

dur_list <- setNames(
  lapply(seq_along(y_list), function(k) {
    y_list[[k]] - entry_list[[k]]
  }),
  selected_outcomes
)


# ================= Event counts =============

cat("\n========== Event counts ==========\n")

for (o in selected_outcomes) {
  
  cat(sprintf(
    "  %-3s : %d\n",
    o,
    sum(status_list[[o]])
  ))
}

#CVA : 85
#DIA : 179
#MI  : 62
#DTH : 37
#CVD : 139
#TIA : 48
#STR : 38

#########################################################

# =================== configuration ====================
set.seed(7)
N <- nrow(X)
n_outcomes <- length(y_list)
n_folds <- 10

methods_to_run <- c("rrr_grplasso", "rrr_lasso", "rrr_ridge", "pen", "mrcox")
rank_vector    <- c(1, 2, 3)

lambda_grids <- list(
  rrr_grplasso = seq(0.01, 0.25, length.out = 11),
  rrr_lasso    = seq(0.01, 0.25, length.out = 11),
  rrr_ridge    = seq(0.2,  0.8, length.out = 11),
  mrcox        = seq(0.1, 0.4, length.out = 11),
  pen          = seq(0.0005, 0.005, length.out = 11)
)

fold_id <- sample(rep(seq_len(n_folds), length.out = N))
needs_followup <- function(method) identical(method, "mrcox")

# ============ C-index functions ====================

compute_uno <- function(lp, time, status) {
  ev <- time[status == 1]
  if (length(ev) < 2) return(NA_real_)
  tau <- quantile(ev, 0.8, na.rm = TRUE)
  dat <- data.frame(time = time, event = status, risk = lp)
  tryCatch(
    survC1::Est.Cval(dat, tau = tau, nofit = TRUE)$Dhat,
    error = function(e) NA_real_)
}


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


# ================ Orient coefficient matrix ====================

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

# ================ Build pen format ======================
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

# ================= Step A: parameter tuning (rank & lambda) =====================

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

cat("\n========== Step A: parameter tuning ==========\n")

selected_params <- list()
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


# selected parameters
cat("\n=== Selected (R, lambda) per method (JOINT, Step A on full data) ===\n")
sel_tbl <- do.call(rbind, lapply(methods_to_run, function(m) {
  s <- selected_params[[m]]
  data.frame(method = m, R = s$R, lambda = signif(s$lambda, 4),
             CVE = signif(s$cve, 4))
}))
print(sel_tbl, row.names = FALSE)


# ================ Step B: OOF prediction ====================

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

# ==================== Univariate glmnet ========================

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


# ================ C-Index Summary =======================
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


cat("\n=== JOINT (multi-outcome) Pooled OOF Uno C-index (follow-up scale) ===\n")
print(round(cindex_uno, 3))

cat("\n=== JOINT Pooled OOF Harrell C-index (age scale; follow-up for mrcox) ===\n")
print(round(cindex_har, 3))

cat("\n=== UNIVARIATE glmnet Pooled OOF Uno C-index (follow-up scale) ===\n")
print(round(cindex_uno_uni, 3))

cat("\n=== UNIVARIATE glmnet Pooled OOF Harrell C-index (age scale) ===\n")
print(round(cindex_har_uni, 3))

# ================= Integrated Brier Score =======================
pooled_ibs_per_method <- function(lp_arr, methods_lab) {
  mat <- matrix(NA_real_, length(methods_lab), n_outcomes,
                dimnames = list(methods_lab, names(y_list)))
  for (i in seq_along(methods_lab)) {
    for (k in seq_len(n_outcomes)) {
      lp_oof <- lp_arr[, i, k]
      ok     <- !is.na(lp_oof)
      if (sum(ok) < 2) next
      dur_k <- dur_list[[k]][ok]
      s_k   <- status_list[[k]][ok]
      ev_k  <- dur_k[s_k == 1]
      if (length(ev_k) < 2) next
      tau_ibs <- as.numeric(quantile(ev_k, 0.8, na.rm = TRUE))
      mat[i, k] <- ibs_oos(
        lp_train = lp_oof[ok], y_train = dur_k, s_train = s_k,
        lp_test  = lp_oof[ok], y_test  = dur_k, s_test  = s_k,
        tau      = tau_ibs)
    }
  }
  mat
}

ibs_pooled     <- pooled_ibs_per_method(lp_pooled,     methods_to_run)
ibs_pooled_uni <- pooled_ibs_per_method(lp_uni_pooled, uni_methods)

cat("\n=== JOINT Pooled OOF IBS on [0, tau] (single value per cell) ===\n")
print(round(ibs_pooled, 4))

cat("\n=== UNIVARIATE glmnet Pooled OOF IBS on [0, tau] ===\n")
print(round(ibs_pooled_uni, 4))

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


# ================= Final model fitting =========================
selected_params[["pen"]]$lambda <- 0.0018
store_estimated_B <- setNames(vector("list", length(methods_to_run)),
                              methods_to_run)

for (method in methods_to_run) {
  
  # Look up selected (R, lambda) from sel_tbl by method name (safe)
  R_sel  <- selected_params[[method]]$R
  lam_sel <- selected_params[[method]]$lambda
  
  if (is.na(R_sel) || is.na(lam_sel)) {
    message(sprintf("[%s] NA R or lambda, skipping.", method)); next
  }
  
  # mrcox uses follow-up duration (no left truncation)
  if (needs_followup(method)) {
    y_fit     <- dur_list
    e_fit     <- NULL
  } else {
    y_fit     <- y_list
    e_fit     <- entry_list
  }
  
  fit1 <- tryCatch(
    fit_one_model(
      X            = X,
      y_list       = y_fit,
      status_list  = status_list,
      entry_list   = e_fit,                          
      dlong        = if (method == "pen") dlong_neo  else NULL,
      method       = method,
      r            = R_sel,
      lambda       = lam_sel,
      group_labels = group_labels
    ),
    error = function(e) {
      message(sprintf("[%s] fit_one_model failed: %s", method, e$message))
      NULL
    }
  )
  
  if (!is.null(fit1) && !is.null(fit1$B_hat)) {
    store_estimated_B[[method]] <- fit1$B_hat
    cat(sprintf("[%s] B_hat dim: %d x %d\n",
                method, nrow(fit1$B_hat), ncol(fit1$B_hat)))
  } else {
    message(sprintf("[%s] returned NULL B_hat", method))
  }
}



# Univariate glmnet 
store_estimated_B_uni <- setNames(vector("list", length(uni_methods)),
                                  uni_methods)

for (m in uni_methods) {
  B_uni <- matrix(0, nrow = ncol(X), ncol = n_outcomes,
                  dimnames = list(NULL, names(y_list)))
  for (k in seq_len(n_outcomes)) {
    sel <- selected_uni[[k]][[m]]
    if (is.na(sel$lambda)) next
    y_surv <- survival::Surv(entry_list[[k]], y_list[[k]], status_list[[k]])
    fit_uni <- tryCatch(
      glmnet::glmnet(
        x = X, y = y_surv, family = "cox",
        alpha       = alpha_map[[m]],
        lambda      = sel$lambda,
        standardize = FALSE),
      error = function(e) { message(sprintf("[%s/%s] glmnet failed: %s",
                                            m, names(y_list)[k], e$message)); NULL })
    if (!is.null(fit_uni))
      B_uni[, k] <- as.numeric(coef(fit_uni))
  }
  store_estimated_B_uni[[m]] <- B_uni
}


all_B <- c(store_estimated_B, store_estimated_B_uni)
# all_B is now a 7-element named list, each element is a p × K matrix


aggregate_by_group <- function(B_mat, group_vec) {
  grp_names <- sort(unique(group_vec))
  do.call(rbind, lapply(grp_names, function(g) {
    rows <- which(group_vec == g)
    colMeans(abs(B_mat[rows, , drop = FALSE]), na.rm = TRUE)
  })) |> `rownames<-`(grp_names)
}


# ========= Heatmap: estimated Beta per method ==============================
library(ggplot2)
library(reshape2)

grp_order <- names(grp_def)

pred_meta <- data.frame(
  Predictor = seq_len(ncol(X)),
  group = factor(group_names, levels = grp_order),
  stringsAsFactors = FALSE
)

# sort predictors by group
pred_meta <- pred_meta[
  order(pred_meta$group, pred_meta$Predictor),
]

pred_levels <- pred_meta$Predictor

group_sizes <- table(pred_meta$group)
group_breaks <- cumsum(group_sizes)
# remove last boundary
group_breaks <- group_breaks[-length(group_breaks)]

# group center positions (for labels)
group_centers <- cumsum(group_sizes) - group_sizes / 2

outcome_levels <- names(y_list)

# Heatmap

for (nm in names(all_B)) {
  
  B_mat <- all_B[[nm]]
  
  if (is.null(B_mat)) {
    message(sprintf("[%s] NULL B_hat", nm))
    next
  }
  
  colnames(B_mat) <- outcome_levels
  
  # matrix -> long format
  df <- reshape2::melt(
    B_mat,
    varnames = c("Predictor", "Outcome"),
    value.name = "Beta"
  )
  
  # ordering
  df$Predictor <- factor(
    df$Predictor,
    levels = pred_levels
  )
  
  df$Outcome <- factor(
    df$Outcome,
    levels = outcome_levels
  )
  
  # symmetric color limits
  lim <- max(abs(df$Beta), na.rm = TRUE)
  
  if (!is.finite(lim) || lim == 0) {
    lim <- 1
  }
  
  # ================= Plot =================
  
  p_hm <- ggplot(
    df,
    aes(x = Outcome,
        y = Predictor,
        fill = Beta)
  ) +
    
    geom_tile() +
    
    # dashed separators between groups
    geom_hline(
      yintercept = group_breaks + 0.5,
      linetype = "dashed",
      colour = "black",
      linewidth = 0.4
    ) +
    
    scale_fill_gradient2(
      low = "#2166ac",
      mid = "white",
      high = "#b2182b",
      midpoint = 0,
      limits = c(-lim, lim),
      name = expression(hat(beta))
    ) +
    
    # show group labels only
    scale_y_discrete(
      breaks = as.character(pred_levels[round(group_centers)]),
      labels = names(group_sizes)
    ) +
    
    labs(
      title = paste0("Estimated Beta: ", nm),
      x = "Outcome",
      y = "Metabolite groups"
    ) +
    
    theme_minimal(base_size = 10) +
    
    theme(
      axis.text.y = element_text(size = 8),
      axis.ticks.y = element_blank(),
      
      panel.grid = element_blank(),
      
      plot.title = element_text(
        face = "bold",
        size = 12
      )
    )
  
  print(p_hm)
}

