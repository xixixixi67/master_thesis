library(ggplot2)
library(survival)
library(pec)
library(prodlim)
library(riskRegression)


## Breslow cumulative hazard from training data; returns event times + steps
breslow_baseline <- function(y_train, s_train, lp_train) {
  ord    <- order(y_train)
  y_s    <- y_train[ord]
  s_s    <- s_train[ord]
  exp_lp <- exp(as.numeric(lp_train[ord]))
  rev_cs <- rev(cumsum(rev(exp_lp)))
  ev_mask <- s_s == 1
  list(time = y_s[ev_mask], H0_step = cumsum(1 / rev_cs[ev_mask]))
}

## H_0(t) — right-continuous step lookup
H0_at <- function(t, baseline) {
  sapply(t, function(tt) {
    idx <- sum(baseline$time <= tt)
    if (idx == 0) 0 else baseline$H0_step[idx]
  })
}

## KM censoring distribution from test data: G(t) = P(C > t)
G_at <- function(t, y_test, s_test) {
  cens_fit <- survival::survfit(survival::Surv(y_test, 1 - s_test) ~ 1)
  if (length(cens_fit$time) == 0) return(rep(1, length(t)))
  sapply(t, function(tt) {
    idx <- sum(cens_fit$time <= tt - 1e-12)
    if (idx == 0) 1 else cens_fit$surv[idx]
  })
}

## OOS IPCW Brier at a single time tau_single
brier_oos <- function(lp_train, y_train, s_train,
                      lp_test,  y_test,  s_test,
                      tau_single = 0.5) {
  tryCatch({
    bl      <- breslow_baseline(y_train, s_train, lp_train)
    H0_tau  <- H0_at(tau_single, bl)
    S_pred  <- exp(-H0_tau * exp(lp_test))
    
    G_y     <- pmax(G_at(y_test,    y_test, s_test), 1e-10)
    G_tau   <- pmax(G_at(tau_single, y_test, s_test), 1e-10)
    
    event_before <- (y_test <= tau_single) & (s_test == 1)
    alive_after  <- y_test  >  tau_single
    
    bs <- ifelse(event_before, (0 - S_pred)^2 / G_y,
                 ifelse(alive_after,  (1 - S_pred)^2 / G_tau, 0))
    mean(bs)
  }, error = function(e) {
    message(sprintf("[brier_oos] %s", e$message))
    NA_real_
  })
}



## OOS Integrated Brier Score on [0, tau]
ibs_oos <- function(lp_train, y_train, s_train,
                    lp_test,  y_test,  s_test,
                    tau) {
  tryCatch({
    eval_times <- sort(unique(y_test[s_test == 1 & y_test <= tau & y_test > 0]))
    if (length(eval_times) < 2) return(NA_real_)
    
    bl         <- breslow_baseline(y_train, s_train, lp_train)
    H0_eval    <- H0_at(eval_times, bl)
    S_mat      <- exp(-outer(exp(lp_test), H0_eval, "*"))   # N_te x n_eval
    
    G_y        <- pmax(G_at(y_test,     y_test, s_test), 1e-10)
    G_eval     <- pmax(G_at(eval_times, y_test, s_test), 1e-10)
    
    bs_vec <- numeric(length(eval_times))
    for (j in seq_along(eval_times)) {
      t_j          <- eval_times[j]
      S_t          <- S_mat[, j]
      event_before <- (y_test <= t_j) & (s_test == 1)
      alive_after  <- y_test  >  t_j
      bs_vec[j] <- mean(
        ifelse(event_before, (0 - S_t)^2 / G_y,
               ifelse(alive_after,  (1 - S_t)^2 / G_eval[j], 0))
      )
    }
    ## Trapezoid from 0 to tau, normalized by tau.
    t_grid <- c(0, eval_times)
    b_grid <- c(0, bs_vec)
    sum(diff(t_grid) * (b_grid[-1] + b_grid[-length(b_grid)]) / 2) / tau
  }, error = function(e) {
    message(sprintf("[ibs_oos] %s", e$message))
    NA_real_
  })
}


## Function to calculate C-index
cindex_simple <- function(lp, time, event, tau) {
  dat <- data.frame(time = time, event = event, risk = lp)
  tryCatch(
    survC1::Est.Cval(dat, tau = tau, nofit = TRUE)$Dhat,
    error = function(e) NA_real_
  )
}


## Function to evaluate ONE method on ONE training dataset
evaluate_one_model_simple <- function(sim_data_train,
                                      X_test, y_test, s_test,
                                      method, lambda, r,
                                      group_labels        = NULL,
                                      tau_quantile        = 0.8,
                                      tau_quantile_single = 0.5) {
  
  K <- length(y_test)
  
  X_train <- sim_data_train$X_list[[1]]
  y_train <- sim_data_train$y_list_all[[1]]
  s_train <- sim_data_train$status_list_all[[1]]
  
  ## pen-specific setup
  dlong_arg           <- NULL
  pen_gamma_start_arg <- NULL
  if (method == "pen") {
    if (is.null(sim_data_train$dlong_list))
      stop("method='pen' needs sim_data_train$dlong_list (regenerate with store_dlong = TRUE)")
    dlong_arg <- sim_data_train$dlong_list[[1]]
    pen_gamma_start_arg <- make_pen_gamma_init(
      X = X_train, y_list = y_train, status_list = s_train,
      R = r, k = K
    )
  }
  
  fit <- fit_one_model(
    X               = X_train,
    y_list          = y_train,
    status_list     = s_train,
    dlong           = dlong_arg,
    method          = method,
    r               = r,
    lambda          = lambda,
    group_labels    = group_labels,
    pen_gamma_start = pen_gamma_start_arg
  )
  
  if (is.null(fit) || is.null(fit$B_hat) || !all(is.finite(fit$B_hat))) {
    return(data.frame(outcome = seq_len(K),
                      cindex  = NA_real_,
                      bs      = NA_real_,
                      ibs     = NA_real_,
                      method  = method))
  }
  
  B_hat     <- fit$B_hat
  eta_train <- X_train %*% B_hat
  eta_test  <- X_test  %*% B_hat
  
  results <- data.frame(outcome = seq_len(K),
                        cindex  = NA_real_,
                        bs      = NA_real_,
                        ibs     = NA_real_)
  
  for (k in seq_len(K)) {
    lp_tr_k    <- eta_train[, k]
    lp_te_k    <- eta_test[,  k]
    time_te_k  <- y_test[[k]]
    event_te_k <- s_test[[k]]
    y_tr_k     <- y_train[[k]]
    s_tr_k     <- s_train[[k]]
    
    ev  <- time_te_k[event_te_k == 1]
    tau <- if (length(ev) > 0) quantile(ev, tau_quantile,        na.rm = TRUE)
    else                median(time_te_k,                  na.rm = TRUE)
    tau_single <- if (length(ev) > 0) quantile(ev, tau_quantile_single, na.rm = TRUE)
    else                median(time_te_k,                  na.rm = TRUE)
    
    results$cindex[k] <- cindex_simple(lp_te_k, time_te_k, event_te_k, tau)
    results$bs[k]     <- brier_oos(lp_tr_k, y_tr_k, s_tr_k,
                                   lp_te_k, time_te_k, event_te_k, tau_single)
    results$ibs[k]    <- ibs_oos  (lp_tr_k, y_tr_k, s_tr_k,
                                   lp_te_k, time_te_k, event_te_k, tau)
  }
  
  results$method <- method
  results
}


## Function to compare a list of methods on ONE training dataset
compare_methods_simple <- function(sim_data_train,
                                   X_test, y_test, s_test,
                                   method_configs,
                                   group_labels = NULL,
                                   tau_quantile = 0.8) {
  
  res_list <- lapply(method_configs, function(cfg) {
    evaluate_one_model_simple(
      sim_data_train = sim_data_train,
      X_test         = X_test,
      y_test         = y_test,
      s_test         = s_test,
      method         = cfg$method,
      lambda         = cfg$lambda,
      r              = cfg$r,
      group_labels   = group_labels,
      tau_quantile   = tau_quantile
    )
  })
  do.call(rbind, res_list)
}


## Function to compare methods over multiple training replicates
compare_methods_multi_sim <- function(sim_data_train, sim_data_test,
                                      method_configs,
                                      group_labels = NULL,
                                      tau_quantile = 0.8,
                                      n_sims       = NULL,
                                      verbose      = TRUE) {
  
  N_avail <- length(sim_data_train$X_list)
  if (is.null(n_sims)) n_sims <- N_avail
  n_sims <- min(n_sims, N_avail)
  
  X_test <- sim_data_test$X_list[[1]]
  y_test <- sim_data_test$y_list_all[[1]]
  s_test <- sim_data_test$status_list_all[[1]]
  
  res_list <- vector("list", n_sims)
  t0 <- proc.time()
  
  for (i in seq_len(n_sims)) {
    if (verbose) cat(sprintf("[sim %d/%d]\n", i, n_sims))
    
    sim_train_i <- list(
      X_list          = list(sim_data_train$X_list[[i]]),
      y_list_all      = list(sim_data_train$y_list_all[[i]]),
      status_list_all = list(sim_data_train$status_list_all[[i]]),
      dlong_list      = if (!is.null(sim_data_train$dlong_list))
        list(sim_data_train$dlong_list[[i]]) else NULL
    )
    
    res_i <- compare_methods_simple(
      sim_data_train = sim_train_i,
      X_test         = X_test,
      y_test         = y_test,
      s_test         = s_test,
      method_configs = method_configs,
      group_labels   = group_labels,
      tau_quantile   = tau_quantile
    )
    res_i$sim <- i
    res_list[[i]] <- res_i
  }
  
  if (verbose)
    cat(sprintf("Done. %.1f s elapsed.\n", (proc.time() - t0)[["elapsed"]]))
  
  do.call(rbind, res_list)
}


## Function to construct wide table "mean (sd)" per (method x outcome) 
make_metric_table_msd <- function(multi_res, metric = c("cindex", "bs", "ibs"),
                                  digits = 4) {
  metric <- match.arg(metric)
  fmt    <- paste0("%.", digits, "f (%.", digits, "f)")
  
  ## per (method, outcome) over sims
  agg <- aggregate(multi_res[[metric]],
                   by  = list(method  = multi_res$method,
                              outcome = multi_res$outcome),
                   FUN = function(x) c(m = mean(x, na.rm = TRUE),
                                       s = sd  (x, na.rm = TRUE),
                                       n = sum(!is.na(x))))
  agg$mean_val <- agg$x[, "m"]
  agg$sd_val   <- agg$x[, "s"]
  agg$cell     <- sprintf(fmt, agg$mean_val, agg$sd_val)
  
  wide <- reshape2::dcast(agg, method ~ outcome, value.var = "cell")
  colnames(wide)[-1] <- paste0("outcome_", colnames(wide)[-1])
  
  per_sim <- aggregate(multi_res[[metric]],
                       by  = list(sim = multi_res$sim, method = multi_res$method),
                       FUN = mean, na.rm = TRUE)
  overall <- aggregate(per_sim$x,
                       by  = list(method = per_sim$method),
                       FUN = function(x) c(m = mean(x, na.rm = TRUE),
                                           s = sd  (x, na.rm = TRUE)))
  overall$cell <- sprintf(fmt, overall$x[, "m"], overall$x[, "s"])
  wide$mean    <- overall$cell[match(wide$method, overall$method)]
  
  wide
}



## Function to generate a held-out test set whose true Beta matches the training scenario
make_test_set <- function(search_data,
                          n_test            = 5000,
                          k_outcomes        = 4,
                          target_event_rate = 0.6,
                          target_sd         = 1,
                          n_probe           = 300) {
  p_val <- nrow(search_data$Beta_matrix)
  r_val <- ncol(search_data$A_matrix)
  
  simulations_data(
    n                 = n_test,
    p                 = p_val,
    r                 = r_val,
    k                 = k_outcomes,
    n_simulations     = 1,
    case              = 2,
    target_sd         = target_sd,
    target_event_rate = target_event_rate,
    n_groups          = search_data$n_groups,
    n_groups_k        = 2,
    zero_blocks       = NULL,   
    store_dlong       = TRUE,
    n_probe           = n_probe,
    beta_obj          = search_data
  )
}


## Function to run the full prediction comparison for one scenario
run_prediction_for_scenario <- function(scenario_label,
                                        search_data, sim_data,
                                        auto_grplasso, auto_lasso,
                                        auto_ridge,    auto_mrcox,
                                        auto_pen     = NULL,
                                        group_labels = NULL,
                                        n_test       = 5000,
                                        n_sims       = 10,
                                        tau_quantile = 0.8) {
  cat(sprintf("\n=== Scenario: %s ===\n", scenario_label))
  
  ## (1) test set with the same fixed Beta
  sim_test <- make_test_set(search_data, n_test = n_test)
  
  ## (2) method_configs from auto-selected (lambda, rank)
  configs <- list(
    list(method = "rrr_grplasso", lambda = auto_grplasso$best_lambda,
         r      = auto_grplasso$best_rank),
    list(method = "rrr_lasso",    lambda = auto_lasso$best_lambda,
         r      = auto_lasso$best_rank),
    list(method = "rrr_ridge",    lambda = auto_ridge$best_lambda,
         r      = auto_ridge$best_rank),
    list(method = "mrcox",        lambda = auto_mrcox$best_lambda,
         r      = auto_mrcox$best_rank)
  )
  if (!is.null(auto_pen)) {
    configs <- c(configs,
                 list(list(method = "pen",
                           lambda = auto_pen$best_lambda,
                           r      = auto_pen$best_rank)))
  }
  
  ## (3) C-index, Brier @ tau_single, Integrated Brier on [0, tau]
  multi_res <- compare_methods_multi_sim(
    sim_data_train = sim_data,
    sim_data_test  = sim_test,
    method_configs = configs,
    group_labels   = group_labels,
    tau_quantile   = tau_quantile,
    n_sims         = n_sims
  )
  
  multi_res$scenario <- scenario_label
  multi_res
}



## ------ Run prediction for all 4 scenarios ------ 
N_TEST <- 5000      
N_SIMS <- 100       

set.seed(7)
res_S1 <- run_prediction_for_scenario(
  scenario_label = "S1: p=600, 10g",
  search_data    = search_data_1, sim_data = sim_data1,
  auto_grplasso  = auto_grplasso_1, auto_lasso = auto_lasso_1,
  auto_ridge     = auto_ridge_1,    auto_mrcox = auto_mrcox_1,
  auto_pen       = NULL,            
  group_labels   = group_labels_1,
  n_test = N_TEST, n_sims = N_SIMS
)

set.seed(7)
res_S2 <- run_prediction_for_scenario(
  scenario_label = "S2: p=600, 25g",
  search_data    = search_data_2, sim_data = sim_data2,
  auto_grplasso  = auto_grplasso_2, auto_lasso = auto_lasso_2,
  auto_ridge     = auto_ridge_2,    auto_mrcox = auto_mrcox_2,
  auto_pen       = NULL,            
  group_labels   = group_labels_2,
  n_test = N_TEST, n_sims = N_SIMS
)

set.seed(7)
res_S3 <- run_prediction_for_scenario(
  scenario_label = "S3: p=100, 10g",
  search_data    = search_data_3, sim_data = sim_data3,
  auto_grplasso  = auto_grplasso_3, auto_lasso = auto_lasso_3,
  auto_ridge     = auto_ridge_3,    auto_mrcox = auto_mrcox_3,
  auto_pen       = auto_pen_3,      
  group_labels   = group_labels_3,
  n_test = N_TEST, n_sims = N_SIMS
)

set.seed(7)
res_S4 <- run_prediction_for_scenario(
  scenario_label = "S4: p=100, 25g",
  search_data    = search_data_4, sim_data = sim_data4,
  auto_grplasso  = auto_grplasso_4, auto_lasso = auto_lasso_4,
  auto_ridge     = auto_ridge_4,    auto_mrcox = auto_mrcox_4,
  auto_pen       = auto_pen_4,      
  group_labels   = group_labels_4,
  n_test = N_TEST, n_sims = N_SIMS
)

all_res <- rbind(res_S1, res_S2, res_S3, res_S4)


## ------ (1) summary tables ------ 
make_scenario_method_table <- function(all_res,
                                       metric = c("cindex", "bs", "ibs"),
                                       digits = 4) {
  metric <- match.arg(metric)
  fmt    <- paste0("%.", digits, "f (%.", digits, "f)")
  
  per_sim <- aggregate(all_res[[metric]],
                       by  = list(sim      = all_res$sim,
                                  scenario = all_res$scenario,
                                  method   = all_res$method),
                       FUN = mean, na.rm = TRUE)
  
  agg <- aggregate(per_sim$x,
                   by  = list(scenario = per_sim$scenario,
                              method   = per_sim$method),
                   FUN = function(x) c(m = mean(x, na.rm = TRUE),
                                       s = sd  (x, na.rm = TRUE)))
  agg$cell <- sprintf(fmt, agg$x[, "m"], agg$x[, "s"])
  
  wide <- reshape2::dcast(agg, scenario ~ method, value.var = "cell")
  wide
}

cat("\n=== C-index | mean (sd) over sims, averaged across outcomes ===\n")
print(make_scenario_method_table(all_res, "cindex"), row.names = FALSE)

cat("\n=== Brier @ tau_single | mean (sd) ===\n")
print(make_scenario_method_table(all_res, "bs"),     row.names = FALSE)

cat("\n=== Integrated Brier on [0, tau] | mean (sd) ===\n")
print(make_scenario_method_table(all_res, "ibs"),    row.names = FALSE)



## ------ (2) Per-outcome detail tables ------ 


cat("\n--- Per-outcome detail tables ---\n")
for (s_lab in unique(all_res$scenario)) {
  cat(sprintf("\n>>> %s\n", s_lab))
  sub_res <- subset(all_res, scenario == s_lab)
  cat("\nC-index:\n")
  print(make_metric_table_msd(sub_res, "cindex"), row.names = FALSE)
  cat("\nBrier @ tau_single:\n")
  print(make_metric_table_msd(sub_res, "bs"),     row.names = FALSE)
  cat("\nIBS on [0, tau]:\n")
  print(make_metric_table_msd(sub_res, "ibs"),    row.names = FALSE)
}



## ------ validate manual Brier calculations against pec and riskRegression::Score on lung data ------

lung <- survival::lung
lung <- na.omit(lung)
lung$status <- as.integer(lung$status == 2)

set.seed(2002)
idx_train <- sample(seq_len(nrow(lung)), 150)
train <- lung[idx_train, ]
test  <- lung[-idx_train, ]

fit <- coxph(Surv(time, status) ~ age + sex + ph.ecog,
             data = train, x = TRUE, y = TRUE, method = "breslow")

lp_train <- predict(fit, newdata = train, type = "lp", reference = "zero")
lp_test  <- predict(fit, newdata = test,  type = "lp", reference = "zero")

tau <- 300

# manual brier_oos
manual_bs <- brier_oos(
  lp_train = lp_train, y_train = train$time, s_train = train$status,
  lp_test  = lp_test,  y_test  = test$time,  s_test  = test$status,
  tau_single = tau
)

# pec::pec
pec_res_bs <- pec::pec(
  object      = list(Cox = fit),
  formula     = Surv(time, status) ~ 1,
  data        = test,
  times       = tau,
  cens.model  = "marginal",
  splitMethod = "none",
  exact       = FALSE,
  reference   = FALSE,
  verbose     = FALSE
)

pec_bs <- as.numeric(pec_res_bs$AppErr$Cox[
  which.min(abs(pec_res_bs$time - tau))
])

# riskRegression::Score
S_test <- pec::predictSurvProb(fit, newdata = test, times = tau)
sc_bs  <- riskRegression::Score(
  object     = list(Cox = 1 - as.matrix(S_test)),
  formula    = Surv(time, status) ~ 1,
  data       = test, times = tau, metrics = "brier",
  null.model = FALSE, summary = NULL, contrasts = FALSE,
  cens.model = "km", se.fit = FALSE
)
sc_bs_df <- as.data.frame(sc_bs$Brier$score)
rr_bs <- as.numeric(sc_bs_df$Brier[sc_bs_df$model == "Cox"][1])

cat(sprintf("\n=== Brier @ tau = %g  (lung; n_train = %d, n_test = %d) ===\n",
            tau, nrow(train), nrow(test)))
cat(sprintf("  manual brier_oos      : %.6f\n", manual_bs))
cat(sprintf("  pec::pec              : %.6f\n", pec_bs))
cat(sprintf("  riskRegression::Score : %.6f\n", rr_bs))
cat(sprintf("  max abs gap           : %.2e\n",
            max(abs(c(manual_bs - pec_bs, manual_bs - rr_bs, pec_bs - rr_bs)))))
# 1.06e-04

# integrated Brier on [0, tau]
manual_ibs <- ibs_oos(
  lp_train = lp_train, y_train = train$time, s_train = train$status,
  lp_test  = lp_test,  y_test  = test$time,  s_test  = test$status,
  tau = tau
)

eval_times <- sort(unique(
  test$time[test$status == 1 & test$time <= tau & test$time > 0]
))

pec_res_ibs <- pec::pec(
  object      = list(Cox = fit),
  formula     = Surv(time, status) ~ 1,
  data        = test, times = eval_times,
  cens.model  = "marginal", splitMethod = "none",
  exact = FALSE, reference = FALSE, verbose = FALSE
)
pec_ibs <- as.numeric(pec::crps(pec_res_ibs, times = tau))

cat(sprintf("\n=== Integrated Brier on [0, %g] ===\n", tau))
cat(sprintf("  manual ibs_oos        : %.6f   (trapezoid 0 -> tau, /tau)\n",
            manual_ibs))
cat(sprintf("  pec::crps             : %.6f   (trapezoid 0 -> tau, /tau)\n",
            pec_ibs))
cat(sprintf("  abs gap               : %.2e\n",
            abs(manual_ibs - pec_ibs)))
# 1.88e-02
