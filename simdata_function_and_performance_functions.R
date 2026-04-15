library(MASS)
library(survival)
library(ggplot2)
library(gridExtra)
library(reshape2)
library(grid)
library(glmnet)
library(dplyr)
library(scales)
library(RColorBrewer)

## Function to check the rank of the simulated beta matrix
rank_check <- function(Beta_m) {
  sv    <- svd(Beta_m)$d
  ratio <- sv[2] / sv[1]
  cat(sprintf("Singular values: %s\nRatio s2/s1: %.4f\n",
              paste(round(sv, 4), collapse = ", "), ratio))
  if (ratio < 0.01) cat("Matrix is close to rank 1.\n") else
    cat("Matrix is NOT close to rank 1.\n")
}

## Function to generate covariate matrix, alpha, Gamma and beta
add_zeros <- function(matrix, percentage = 0.10) {
  total_components <- length(matrix)
  zeros_count <- floor(total_components * percentage)
  
  indexes <- 1:total_components
  
  zero_indexes <- sample(indexes, zeros_count)
  
  matrix[zero_indexes] <- 0
  
  return(matrix)
}

## Functions to simulated beta matrix
simulate_Beta_matrix <- function(p, r, k, case = 1, scale = 1,
                                 n_groups = NULL, zero_blocks = NULL,
                                 n_groups_k = NULL,
                                 rho_w = 0.7,   rho_b = 0.1,
                                 rho_w_k = 0.7, rho_b_k = 0.1,
                                 mu_l = 0.2,    mu_u = 1.0,
                                 use_kms_beta  = TRUE,
                                 use_kms_gamma = TRUE) {
  
  make_kms_corr <- function(dim_total, n_g, pg, rho_w, rho_b) {
    Sigma <- matrix(rho_b, dim_total, dim_total)
    for (g in seq_len(n_g)) {
      idx <- ((g - 1) * pg + 1):(g * pg)
      Sigma[idx, idx] <- toeplitz(
        pmax(rho_w ^ ((0:(pg - 1)) / (pg - 1)), rep(rho_b, pg))
      )
    }
    D_inv <- diag(1 / sqrt(diag(Sigma)))
    D_inv %*% Sigma %*% D_inv
  }
  
  sample_copula_col <- function(R_corr, mu_vec = NULL) {
    d   <- nrow(R_corr)
    z   <- MASS::mvrnorm(1, mu = rep(0, d), Sigma = R_corr)
    u   <- pnorm(z)
    if (!is.null(mu_vec))
      u <- pmin(pmax(u + mu_vec - 0.5, 0), 1)
    u
  }
  
  if (case == 1) {
    
    A_matrix     <- matrix(runif(p * r, 0, 1), p, r)
    Gamma_matrix <- matrix(runif(k * r, 0, 1), k, r)
    
  } else if (case == 2) {
    
    if (is.null(n_groups) || p %% n_groups != 0)
      stop(sprintf("p (%d) must be divisible by n_groups (%d)", p, n_groups))
    p_g <- p / n_groups

    zero_group_ids <- if (!is.null(zero_blocks)) sapply(zero_blocks, `[`, 1) else integer(0)
    active_groups  <- setdiff(seq_len(n_groups), zero_group_ids)
    
    if (use_kms_beta) {
      R_beta <- make_kms_corr(p, n_groups, p_g, rho_w, rho_b)

      mus     <- numeric(n_groups)
      mus[active_groups] <- runif(length(active_groups), mu_l, mu_u)
      mu_full <- rep(mus, each = p_g)  
      
      A_matrix <- matrix(NA, p, r)
      for (col in seq_len(r)) {
        a_col <- sample_copula_col(R_beta, mu_vec = mu_full)
        zero_idx <- unlist(lapply(zero_group_ids,
                                  function(g) ((g-1)*p_g + 1):(g*p_g)))
        if (length(zero_idx)) a_col[zero_idx] <- 0
        A_matrix[, col] <- a_col
      }
    } else {
      A_matrix <- matrix(runif(p * r, 0, 1), p, r)
      for (blk in zero_blocks) {
        rows <- ((blk[1]-1)*p_g + 1):(blk[1]*p_g)
        A_matrix[rows, ] <- 0
      }
    }
    
    if (!is.null(n_groups_k) && use_kms_gamma) {
      if (k %% n_groups_k != 0)
        stop(sprintf("k (%d) must be divisible by n_groups_k (%d)", k, n_groups_k))
      k_g         <- k / n_groups_k
      R_gamma     <- make_kms_corr(k, n_groups_k, k_g, rho_w_k, rho_b_k)
      
      Gamma_matrix <- matrix(NA, k, r)
      for (col in seq_len(r)) {
        Gamma_matrix[, col] <- sample_copula_col(R_gamma)     # ← Copula
      }
    } else {
      Gamma_matrix <- matrix(runif(k * r, 0, 1), k, r)
    }
    
  } else if (case == 3) {
    
    A_matrix     <- add_zeros(matrix(runif(p * r, 0, 1), p, r), 0.10)
    Gamma_matrix <- add_zeros(matrix(runif(k * r, 0, 1), k, r), 0.10)
    
  }
  
  Beta_matrix <- A_matrix %*% t(Gamma_matrix)
  
  return(list(
    Beta_matrix  = Beta_matrix * scale,
    A_matrix     = A_matrix,
    Gamma_matrix = Gamma_matrix,
    group_size   = if (case == 2) p / n_groups              else NULL,
    n_groups     = if (case == 2) n_groups                  else NULL,
    group_size_k = if (case == 2 && !is.null(n_groups_k)) k / n_groups_k else NULL,
    n_groups_k   = if (case == 2) n_groups_k                else NULL
  ))
}

## Function to generate group_labels
gen_group_labels <- function(p,n_groups) {
  group_size <- p/n_groups
  rep(1:n_groups, each=group_size)
}

## Function to simulate survival data

# X ~ N(0,1), T ~ Exp(exp(X %*% B)), C ~ Uniform(0,1)
# T_obs = min(T, C), delta = I(T <= C)
simulate_data_model = function(Beta_matrix, n, corr=0, group_corr = 0, group_labels = NULL, 
                               use_kms_X=TRUE, n_groups=NULL, rho_w=0.7, rho_b=0.1, Print_flug=FALSE){
  p <- nrow(Beta_matrix)   
  k <- ncol(Beta_matrix)
  
  # KMS Toeplitz
  if (use_kms_X && !is.null(n_groups)) {
    if (p %% n_groups != 0)
      stop(sprintf("p (%d) must be divisible by n_groups (%d)", p, n_groups))
    p_g   <- p / n_groups
    Sigma <- matrix(rho_b, p, p)  
    for (g in seq_len(n_groups)) {
      idx <- ((g - 1) * p_g + 1):(g * p_g)
      
      within <- toeplitz(
        pmax(rho_w^((0:(p_g - 1)) / (p_g - 1)), rep(rho_b, p_g))
      )
      Sigma[idx, idx] <- within
    }
    X <- MASS::mvrnorm(n, mu = rep(0, p), Sigma = Sigma)
    
    # only simple correlation within group
  } else if (group_corr > 0 && !is.null(group_labels)) {
    n_groups_old <- max(group_labels)
    Sigma        <- diag(p)
    for (g in seq_len(n_groups_old)) {
      idx                 <- which(group_labels == g)
      Sigma[idx, idx]     <- group_corr
      diag(Sigma)[idx]    <- 1
    }
    X <- mvtnorm::rmvnorm(n, mean = rep(0, p), sigma = Sigma)
    
  } else if (corr > 0) {
    sigma        <- matrix(corr, p, p)
    diag(sigma)  <- 1
    X <- mvtnorm::rmvnorm(n, mean = rep(0, p), sigma = sigma)
    
  } else {
    X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  }
  
  eta <- X%*%Beta_matrix
  exp_eta <- exp(eta)
  exp_eta[exp_eta<=1e-300] <- 1e-300
  
  #GENERATE T(times)
  # later: mean_age = 60 
  # Outcome times T from exp distribution
  T_sim <- matrix(NA, nrow = n, ncol = k)
  for (i in 1:k) {
    T_sim[, i] <- rexp(n, rate =  exp(eta)[,i])
  }
  
  #GENERATE censoring 
  #Random censoring time from Uniform distribution
  C <- runif(n, min = 0, max = 1)
  
  #DETERMINE the observed time and censoring indicator
  T_obs <- matrix(NA, nrow = n, ncol = k)
  delta <- matrix(NA, nrow = n, ncol = k)
  
  for (i in 1:k) {
    T_obs[,i] <- pmin(T_sim[,i],C)
    delta[,i] <- as.integer(T_sim[,i]<=C)
  }
  
  # Check the proportion of censored vs event occurrences
  prop_censored <- sum(delta == 0)/n/k
  prop_event <- sum(delta == 1)/n/k
  # Output proportions 
  if (Print_flug) {
    print(paste('Proportion of censoring',prop_censored))
    print(paste('Proportion of events',prop_event))
  }
  y_list <- lapply(1:k, function(i) T_obs[,i])
  status_list <- lapply(1:k, function(i) delta[,i])
  
  #CREATE data frame
  data <- data.frame(id = 1:n)
  for (i in 1:k) {
    data[[paste0("t", i)]] <- T_obs[, i]
    data[[paste0("d", i)]] <- delta[, i]
  }
  # Add predictors 
  for (i in 1:p) {
    data[[paste0("x", i)]] <- X[, i]
  }
  
  return (list(
    y_list = y_list,
    status_list = status_list,
    X = X,
    T_sim = T_sim,
    T_obs = T_obs,
    data = data
  ))
}


## initialization for "pen"
make_pen_gamma_init <- function(X, y_list, status_list, R, k,
                                rho = 0.01) {
  p  <- ncol(X)
  B0 <- matrix(0, p, k)
  
  for (kk in seq_len(k)) {
    fit_ridge <- tryCatch(
      glmnet::glmnet(X,
                     survival::Surv(y_list[[kk]], status_list[[kk]]),
                     family  = "cox",
                     alpha   = 0,
                     lambda  = rho),
      error = function(e) NULL
    )
    if (!is.null(fit_ridge))
      B0[, kk] <- as.vector(coef(fit_ridge))
  }
  
  sv      <- svd(B0, nu = 0, nv = min(R, ncol(B0)))
  Gamma0  <- sv$v[, 1:R, drop = FALSE]   # k x R
  
  return(t(Gamma0))   # R x K
}

## Function to preprocess data to be ready to enter the redrank model
build_dlong = function(data,k){
  library(mstate)
  #subjects
  n <- nrow(data)
  # Transition matrix(k outcomes)
  tmat <- trans.comprisk(k)
  # Column names for survival times and event indicators
  tnames <- paste("t", 1:k, sep="")
  dnames <- paste("d", 1:k, sep="")
  # Create dlong
  dlong <- data.frame(
    id = rep(data$id, each = k),
    Tstart = 0,
    Tstop = c(t(matrix(unlist(c(data[, tnames])), n, k))),
    status = c(t(matrix(unlist(c(data[, dnames])), n, k))),
    from = 1,
    to = rep(2:(k+1), n),
    trans = rep(1:k, n)
  )
  
  #add predicotrs
  predictor_names <- grep("^x", names(data), value = TRUE)
  # Add each predictor to the dlong data frame
  for (predictor in predictor_names) {
    dlong[[predictor]] <- rep(data[[predictor]], each = k)
  }
  
  # Calculate the time
  dlong$time <- dlong$Tstop - dlong$Tstart
  # Define the transition matrix and class
  class(dlong) <- c("msdata", "data.frame")
  attr(dlong, "trans") <- tmat
  return(dlong)
}

## Function to fit specified model
fit_one_model <- function(X, y_list, status_list, dlong = NULL,
                          method = c("rrr_grplasso", "rrr_lasso", "rrr_ridge", "pen", "mrcox"),
                          r, lambda, group_labels = NULL,
                          pen_gamma_start = NULL, ...) {
  method <- match.arg(method)
  
  if (method == "rrr_grplasso") {
    raw   <- solve_RR_GrpLasso(X=X, y_list=y_list,
                               status_list=status_list,
                               R=r, lambda_alpha=lambda,
                               group_labels=group_labels, ...)
    res   <- raw$result[[1]]
    Alpha <- res$alpha; Gamma <- res$Gamma
    B_hat <- Alpha %*% t(Gamma)
    
  } else if (method == "rrr_lasso") {
    raw   <- solve_RR_Lasso(X=X, y_list=y_list,
                            status_list=status_list,
                            R=r, lambda_alpha=lambda, ...)
    res   <- raw$result[[1]]
    Alpha <- res$alpha; Gamma <- res$Gamma
    B_hat <- Alpha %*% t(Gamma)
    
  } else if (method == "rrr_ridge") {
    raw   <- solve_RR_Ridge(X=X, y_list=y_list,
                            status_list=status_list,
                            R=r, lambda_alpha=lambda, ...)
    res   <- raw$result[[1]]
    Alpha <- res$alpha; Gamma <- res$Gamma
    B_hat <- Alpha %*% t(Gamma)
    
  } else if (method == "mrcox") {
    raw <- solve_aligned(X=X, y_list=y_list, status_list=status_list,
                         lambda_1=0, lambda_2=lambda)
    B_hat <- raw$result[[1]]
    Alpha <- NULL; Gamma <- NULL
    
  } else {
    if (is.null(dlong)) stop("method='pen' needs dlong")
    
    if (is.null(pen_gamma_start)) {
      k_outcomes      <- length(y_list)
      pen_gamma_start <- make_pen_gamma_init(X, y_list, status_list,
                                             R=r, k=k_outcomes)
    }
    
    pred_names  <- grep("^x", names(dlong), value=TRUE)
    formula_str <- paste("Surv(Tstop, status) ~",
                         paste(pred_names, collapse=" + "))
    
    raw <- tryCatch({
      suppressWarnings(
        invisible(capture.output( 
          raw_inner <- pen.survrrr(
            as.formula(formula_str),
            dat          = dlong,
            R            = r,
            Gamma.iter   = pen_gamma_start,
            lambda.alpha = lambda,
            lambda.gamma = lambda,
            eps = 1e-3, maxit = 1e6, thresh = 1e-5,
            standardize.opt = FALSE,
            alpha = 1, ...)
        ))
      )
      raw_inner
    },
    error = function(e) NULL)
    
    if (is.null(raw)) return(NULL) 
    
    Alpha <- raw$Alpha; Gamma <- raw$Gamma
    B_hat <- Alpha %*% Gamma
  }
  
  return(list(alpha=Alpha, Gamma=Gamma, B_hat=B_hat, raw=raw))
}

## Function to check performance of given model
performance_model=function(B_hat, B_true, Print_flug=FALSE) {
  
  k <- ncol(B_true)
  
  #correlation coeff
  cor_coeff = cor(as.vector(B_hat),as.vector(B_true))
  
  # Each coefficient separately
  # bias
  bias <- (B_hat - B_true)
  # MSE
  mse<- (B_hat - B_true)^2
  bias_each_coeff = mean(abs(bias))
  mse_each_coeff = mean(mse)
  
  if (Print_flug==TRUE){
    print(paste('Mean bias taking each coeff seperatly',bias_each_coeff))
    print(paste('Mean mse taking each coeff seperatly',mse_each_coeff))
    print(paste('correlation between estimated Beta and the true Beta',cor_coeff))
  }
  
  
  #For each outcome separately
  
  b_pred=numeric(k)
  mse_pred=numeric(k)
  
  for (i in 1:k){
    mse_pred[i] = mean((B_hat[,i]-B_true[,i])^2)
    b_pred[i] = mean(B_hat[,i]-B_true[,i])
    if (Print_flug==TRUE) cat(sprintf("outcome: %d | bias: %.4f | MSE: %.4f\n",
                                      i, b_pred[i], mse_pred[i]))
  }
  
  
  return(list(
    bias_out = bias_each_coeff,
    mse_out = mse_each_coeff,
    b_pred = b_pred,
    mse_pred = mse_pred,
    cor_coeff = cor_coeff,
    bias_for_all_coeffs = bias,
    mse_for_all_coeffs = mse
    
  ))
  
}


## Function to simulate n_simulation datasets of given rank r
simulations_data <- function(n, p, r, k, n_simulations,
                             case = 1, scale = 1, corr = 0, group_corr = 0,
                             n_groups = NULL, zero_blocks = NULL,
                             store_dlong = FALSE, beta_obj = NULL,
                             random_beta = FALSE) {
  
  # shared-beta path: generate once (or use supplied beta_obj)
  if (!random_beta) {
    if (is.null(beta_obj))
      beta_obj <- simulate_Beta_matrix(p, r, k, case = case, scale = scale,
                                       n_groups = n_groups,
                                       zero_blocks = zero_blocks)
    fixed_Beta <- beta_obj$Beta_matrix
  }
  
  group_labels_vec <- if (case == 2) gen_group_labels(p, n_groups) else NULL
  
  y_list_all <- status_list_all <- X_list <-
    T_sim_list <- T_obs_list <- data_wide_list <-
    Beta_matrix_list <-                     
    vector("list", n_simulations)
  dlong_list <- if (store_dlong) vector("list", n_simulations) else NULL
  
  for (i in seq_len(n_simulations)) {
    
    if (random_beta) {
      beta_i   <- simulate_Beta_matrix(p, r, k, case = case, scale = scale,
                                       n_groups = n_groups,
                                       zero_blocks = zero_blocks)
      Beta_i   <- beta_i$Beta_matrix
    } else {
      Beta_i   <- fixed_Beta
    }
    Beta_matrix_list[[i]] <- Beta_i
    
    sim <- simulate_data_model(Beta_i, n, corr = corr,
                               group_corr = group_corr,
                               group_labels = group_labels_vec)
    y_list_all[[i]]      <- sim$y_list
    status_list_all[[i]] <- sim$status_list
    X_list[[i]]          <- sim$X
    T_sim_list[[i]]      <- sim$T_sim
    T_obs_list[[i]]      <- sim$T_obs
    data_wide_list[[i]]  <- sim$data
    if (store_dlong) dlong_list[[i]] <- build_dlong(sim$data, k)
  }
  
  # for backward compatibility: Beta_matrix = first replicate's beta
  Beta_ref <- Beta_matrix_list[[1]]
  A_ref    <- if (!random_beta) beta_obj$A_matrix    else NULL
  G_ref    <- if (!random_beta) beta_obj$Gamma_matrix else NULL
  
  return(list(
    y_list_all        = y_list_all,
    status_list_all   = status_list_all,
    X_list            = X_list,
    Beta_matrix       = Beta_ref,  
    Beta_matrix_list  = Beta_matrix_list, 
    random_beta       = random_beta,
    A_matrix          = A_ref,
    Gamma_matrix      = G_ref,
    T_sim_list        = T_sim_list,
    T_obs_list        = T_obs_list,
    data_wide_list    = data_wide_list,
    dlong_list        = dlong_list,
    group_labels      = group_labels_vec,
    n_groups          = if (case == 2) n_groups       else NULL,
    group_size        = if (case == 2) p / n_groups   else NULL
  ))
}

## Function to compute CVE (V&VH, 1993) for one fold
compute_vvh_fold <- function(X_train, X_test,
                             y_train, y_test,
                             s_train, s_test,
                             B_hat) {
  K <- ncol(B_hat)
  N_test <- nrow(X_test)
  X_full <- rbind(X_train, X_test)
  p <- ncol(X_train)
  pred_names <- paste0("xv", seq_len(p))
  
  fml <- as.formula(
    paste("survival::Surv(y, s) ~", paste(pred_names, collapse = "+")))
  
  loglik_full  <- 0
  loglik_train <- 0
  
  for (k in seq_len(K)) {
    beta_k <- B_hat[, k]
    
    if (!all(is.finite(beta_k))) next 
    
    dat_full <- cbind(
      data.frame(y=c(y_train[[k]], y_test[[k]]),
                 s=c(s_train[[k]], s_test[[k]])),
      setNames(as.data.frame(X_full), pred_names)
    )
    
    fit_full <- tryCatch(
      survival::coxph(fml, data = dat_full, init = beta_k,
                      control = survival::coxph.control(iter.max = 0,
                                                        timefix  = FALSE)),
      error = function(e) NULL)
    
    dat_train <- cbind(
      data.frame(y=y_train[[k]], s=s_train[[k]]),
      setNames(as.data.frame(X_train), pred_names)
    )
    
    dat_train <- cbind(
      data.frame(y = y_train[[k]], s = s_train[[k]]),
      setNames(as.data.frame(X_train), pred_names))
    
    fit_train <- tryCatch(
      survival::coxph(fml, data = dat_train, init = beta_k,
                      control = survival::coxph.control(iter.max = 0,
                                                        timefix  = FALSE)),
      error = function(e) NULL)
    
    if (!is.null(fit_full)  && is.finite(fit_full$loglik[2]) &&
        !is.null(fit_train) && is.finite(fit_train$loglik[2])) {
      loglik_full  <- loglik_full  + fit_full$loglik[2]
      loglik_train <- loglik_train + fit_train$loglik[2]
    }
  }
  
  result <- -2 * (loglik_full - loglik_train) / N_test
  if (!is.finite(result)) return(NA) 
  return(result)
}

## Function to select the best lambda using V&VH criterion
select_lambda_vvh <- function(X, y_list, status_list, R,
                              lambda_vector,
                              group_labels = NULL,          
                              method  = c("rrr_grplasso", "rrr_lasso", "rrr_ridge", "mrcox"),
                              n_folds = 5, seed = 123, ...) {
  method  <- match.arg(method)
  set.seed(seed)
  N <- nrow(X)
  n_lam <- length(lambda_vector)
  fold_id <- sample(rep(seq_len(n_folds), length.out=N))
  CVE_mat <- matrix(NA, n_folds, n_lam)
  
  for (fold in seq_len(n_folds)) {
    idx_te <- which(fold_id == fold)
    idx_tr <- which(fold_id != fold)
    X_tr <- X[idx_tr, , drop = FALSE]
    X_te <- X[idx_te, , drop = FALSE]
    y_tr <- lapply(y_list,      `[`, idx_tr)
    y_te <- lapply(y_list,      `[`, idx_te)
    s_tr <- lapply(status_list, `[`, idx_tr)
    s_te <- lapply(status_list, `[`, idx_te)
    
    for (lam in seq_along(lambda_vector)) {
      fit <- tryCatch(
        fit_one_model(X = X_tr, y_list = y_tr, status_list = s_tr,
                      method = method, r = R,
                      lambda = lambda_vector[lam],
                      group_labels = group_labels, ...),
        error = function(e) NULL)
      
      if (is.null(fit) || !all(is.finite(fit$B_hat))) next
      
      cve <- compute_vvh_fold(X_tr, X_te, y_tr, y_te, s_tr, s_te, fit$B_hat)
      if (!is.na(cve)) CVE_mat[fold, lam] <- cve
    }
    cat(sprintf("  [%s] fold %d/%d done\n", method, fold, n_folds))
  }
  
  CVE_mean <- colMeans(CVE_mat, na.rm = TRUE)
  CVE_se <- apply(CVE_mat, 2, sd, na.rm = TRUE) / sqrt(n_folds)
  best_lambda <- lambda_vector[which.min(CVE_mean)]
  
  cat(sprintf("[%s VVH] best lambda: %g  (CVE = %.4f)\n",
              method, best_lambda, min(CVE_mean, na.rm = TRUE)))
  
  return(list(CVE_mean = CVE_mean, CVE_se = CVE_se, CVE_mat = CVE_mat,
              lambda = lambda_vector, best_lambda = best_lambda))
}

## Function to select the best lambda wrt. "pen" method
select_lambda_pen_vvh <- function(dlong, data_wide, R,
                                  lambda_vector, pred_names, k,
                                  n_folds = 5, seed = 123, ...) {
  set.seed(seed)
  subject_ids <- unique(dlong$id)
  N_subj <- length(subject_ids)
  fold_id <- setNames(
    sample(rep(seq_len(n_folds), length.out = N_subj)),
    as.character(subject_ids))
  
  n_lam <- length(lambda_vector)
  CVE_mat <- matrix(NA, n_folds, n_lam)
  
  for (fold in seq_len(n_folds)) {
    test_ids <- names(fold_id)[fold_id == fold]
    train_ids <- names(fold_id)[fold_id != fold]
    N_test <- length(test_ids)
    
    data_train <- data_wide[data_wide$id %in% train_ids, ]
    dlong_train_ms <- build_dlong(data_train, k)
    dlong_full <- dlong
    
    X_tr <- as.matrix(data_train[, pred_names])
    y_tr <- lapply(seq_len(k), function(i) data_train[[paste0("t", i)]])
    s_tr <- lapply(seq_len(k), function(i) data_train[[paste0("d", i)]])
    gamma_init <- make_pen_gamma_init(X_tr, y_tr, s_tr, R = R, k = k)
    
    formula_rr <- as.formula(
      paste("Surv(Tstop, status) ~", paste(pred_names, collapse = "+")))
    
    for (lam in seq_along(lambda_vector)) {
      raw_fit <- tryCatch({
        suppressWarnings(invisible(capture.output(
          raw_inner <- pen.survrrr(
            formula_rr, dat = dlong_train_ms, R = R,
            Gamma.iter   = gamma_init,
            lambda.alpha = lambda_vector[lam],
            lambda.gamma = lambda_vector[lam],
            eps = 1e-3, maxit = 1e6, thresh = 1e-5,
            standardize.opt = FALSE, alpha = 1, ...)
        )))
        raw_inner
      }, error = function(e) NULL)
      
      if (is.null(raw_fit)) next
      B_hat <- raw_fit$Alpha %*% raw_fit$Gamma
      # Use dlong_train_ms for both to ensure consistency
      loglik_full <- calculate.log.partial.lik(dlong_full,     B_hat, pred_names)
      loglik_train <- calculate.log.partial.lik(dlong_train_ms, B_hat, pred_names)
      CVE_mat[fold, lam] <- -2 * (loglik_full - loglik_train) / N_test
    }
    cat(sprintf("  [pen] fold %d/%d done\n", fold, n_folds))
  }
  
  CVE_mean <- colMeans(CVE_mat, na.rm = TRUE)
  CVE_se  <- apply(CVE_mat, 2, sd, na.rm = TRUE) / sqrt(n_folds)
  best_lambda <- lambda_vector[which.min(CVE_mean)]
  
  cat(sprintf("[pen VVH] best lambda: %g  (CVE = %.4f)\n",
              best_lambda, min(CVE_mean, na.rm = TRUE)))
  
  return(list(CVE_mean = CVE_mean, CVE_se = CVE_se, CVE_mat = CVE_mat,
              lambda = lambda_vector, best_lambda = best_lambda))
}

## Function to fit data n_simulations times and access overall performance 
simulations_fit_and_performance <- function(datas, p, k, n_simulations, r,
                                            lambda_alpha, group_labels = NULL,
                                            method, pen_gamma_start = NULL, ...) {
  if (method == "pen" && is.null(datas$dlong_list))
    stop("method='pen' needs dlong_list; use store_dlong=TRUE.")
  
  use_random_beta <- isTRUE(datas$random_beta)
  
  N_obs <- nrow(datas$X_list[[1]])
  bias1 <- mse1 <- cor_coeff_vec <- c()
  bias_list <- mse_list <- bias_for_all_coeffs <- mse_for_all_coeffs <-
    vector("list", n_simulations)
  simulation_times <- numeric(n_simulations)
  
  y_list_all <- datas$y_list_all
  status_list_all <- datas$status_list_all
  X_list <- datas$X_list
  dlong_list<- datas$dlong_list
  
  for (i in seq_len(n_simulations)) {
    t0 <- Sys.time()
    
    # pick the correct Beta_true for this replicate
    Beta_true_i <- if (use_random_beta) datas$Beta_matrix_list[[i]] else
      datas$Beta_matrix
    
    dlong <- if (method == "pen") dlong_list[[i]] else NULL
    
    if (method == "pen") {
      gamma_init <- if (!is.null(pen_gamma_start)) pen_gamma_start else
        if (p > 50 || N_obs > 500)
          make_pen_gamma_init(X_list[[i]], y_list_all[[i]],
                              status_list_all[[i]], R = r, k = k)
      else matrix(rnorm(r * k), r, k)
    } else {
      gamma_init <- NULL
    }
    
    fit <- tryCatch(
      fit_one_model(X = X_list[[i]], y_list = y_list_all[[i]],
                    status_list = status_list_all[[i]], dlong = dlong,
                    method = method, r = r, lambda = lambda_alpha,
                    group_labels = group_labels,
                    pen_gamma_start = gamma_init, ...),
      error = function(e) {
        cat(sprintf("error(sim %d): %s\n", i, e$message)); NULL })
    
    if (is.null(fit)) next
    
    res  <- performance_model(fit$B_hat, Beta_true_i)
    bias1 <- c(bias1, res$bias_out)
    mse1 <- c(mse1, res$mse_out)
    cor_coeff_vec <- c(cor_coeff_vec, res$cor_coeff)
    bias_list[[i]] <- res$b_pred
    mse_list[[i]] <- res$mse_pred
    bias_for_all_coeffs[[i]] <- res$bias_for_all_coeffs
    mse_for_all_coeffs[[i]] <- res$mse_for_all_coeffs
    simulation_times[i] <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  }
  
  cat(sprintf("[%s] %d sims done, total %.1f sec\n",
              method, n_simulations, sum(simulation_times)))
  
  valid_bias <- Filter(Negate(is.null), bias_list)
  valid_mse  <- Filter(Negate(is.null), mse_list)
  
  if (length(valid_bias) == 0) {
    warning("All simulations failed, returning NULL.")
    return(NULL)
  }
  
  bias_df <- do.call(rbind, valid_bias)
  mse_df  <- do.call(rbind, valid_mse)
  colnames(bias_df) <- colnames(mse_df) <- paste0("Outcome_", 1:k)
  
  return(list(bias_df   = bias_df, mse_df  = mse_df,
              bias1     = bias1,   mse1    = mse1,
              cor_coeff_vec     = cor_coeff_vec,
              simulation_times  = simulation_times,
              bias_for_all_coeffs_and_simulations = bias_for_all_coeffs,
              mse_for_all_coeffs_and_simulation   = mse_for_all_coeffs))
}

## Function to fit the models for multiple datasets and find the best lambda (cv)
simulations_fit_find_best_lambda <- function(datasets_list,
                                             lambda_vector,
                                             p, k, r,
                                             group_labels = NULL,
                                             method,
                                             n_folds_0 = 5,
                                             ...) {
  if (method == "pen" && is.null(datasets_list$dlong_list))
    stop("method='pen' needs dlong_list, use store_dlong=TRUE.")
  
  n_datasets      <- length(datasets_list$y_list_all)
  optimal_lambdas <- numeric(0)
  results_df      <- data.frame()
  dataset_times <- numeric(n_datasets)
  
  for (j in seq_len(n_datasets)) {
    cat(sprintf("\n=== Dataset %d/%d ===\n", j, n_datasets))
    t0 <- Sys.time()
    X           <- datasets_list$X_list[[j]]
    y_list      <- datasets_list$y_list_all[[j]]
    status_list <- datasets_list$status_list_all[[j]]
    
    if (method == "pen") {
      o_res <- select_lambda_pen_vvh(
        dlong         = datasets_list$dlong_list[[j]],
        data_wide     = datasets_list$data_wide_list[[j]],
        R             = r,
        lambda_vector = lambda_vector,
        pred_names    = paste0("x", 1:p),
        k             = k,
        n_folds       = n_folds_0, ...)
    } else {
      o_res <- select_lambda_vvh(
        X             = X,
        y_list        = y_list,
        status_list   = status_list,
        R             = r,
        lambda_vector = lambda_vector,
        group_labels  = group_labels,
        method        = method,
        n_folds       = n_folds_0, ...)
    }
    
    optimal_lambdas <- c(optimal_lambdas, o_res$best_lambda)
    results_df      <- rbind(results_df,
                             data.frame(dataset  = j,
                                        lambda   = lambda_vector,
                                        CVE_mean = o_res$CVE_mean,
                                        CVE_se   = o_res$CVE_se))
    dataset_times[j] <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    
    cat(sprintf("  -> dataset %d | best_lambda: %g | time: %.2f sec\n", 
                j, o_res$best_lambda, dataset_times[j]))
  }
  
  best_lambda <- lambda_vector[which.min(abs(lambda_vector - median(optimal_lambdas)))]
  cat(sprintf("\n[%s VVH] global best lambda: %g, total CV time %.1f sec\n", 
              method, best_lambda, sum(dataset_times)))
  return(list(results_df      = results_df,
              best_lambda     = best_lambda,
              optimal_lambdas = optimal_lambdas))
}

## Helper: subset a datas object to first n datasets
subset_datas <- function(datas, n) {
  n <- min(n, length(datas$y_list_all))
  datas$y_list_all       <- datas$y_list_all[1:n]
  datas$status_list_all  <- datas$status_list_all[1:n]
  datas$X_list           <- datas$X_list[1:n]
  datas$T_sim_list       <- datas$T_sim_list[1:n]
  datas$T_obs_list       <- datas$T_obs_list[1:n]
  datas$data_wide_list   <- datas$data_wide_list[1:n]
  datas$Beta_matrix_list <- datas$Beta_matrix_list[1:n]
  if (!is.null(datas$dlong_list))
    datas$dlong_list <- datas$dlong_list[1:n]
  return(datas)
}

## Auto rank selection with optional separate CV dataset
simulations_fit_auto_rank <- function(datas,
                                      datas_cv      = NULL,
                                      p, k,
                                      n_simulations,
                                      lambda_vector,
                                      group_labels  = NULL,
                                      method,
                                      rank_vector   = c(1, 2, 3),
                                      n_folds_cv    = 5,
                                      n_datasets_cv = 10,
                                      verbose       = TRUE,
                                      ...) {
  
  extract_best_cve <- function(cv_res) {
    if (!is.null(cv_res$CVE_mean) && !is.null(cv_res$lambda)) {
      idx <- which(cv_res$lambda == cv_res$best_lambda)
      return(mean(cv_res$CVE_mean[idx], na.rm = TRUE))
    }
    if (!is.null(cv_res$results_df)) {
      sub <- cv_res$results_df[cv_res$results_df$lambda == cv_res$best_lambda, ]
      return(mean(sub$CVE_mean, na.rm = TRUE))
    }
    return(Inf)
  }
  
  if (!is.null(datas_cv)) {
    cv_data <- datas_cv
    if (verbose)
      cat(sprintf("Lambda selection: using supplied datas_cv (%d datasets).\n",
                  length(cv_data$y_list_all)))
  } else {
    cv_data <- subset_datas(datas, n_datasets_cv)
    if (verbose)
      cat(sprintf("Lambda selection: using first %d of %d datasets.\n",
                  length(cv_data$y_list_all), length(datas$y_list_all)))
  }
  

  if (method == "mrcox") {
    if (verbose) cat("\n[mrcox] Full-rank method — skipping rank selection.\n")
    
    cv_res <- simulations_fit_find_best_lambda(
      datasets_list = cv_data,
      lambda_vector = lambda_vector,
      p = p, k = k, r = 1,    
      group_labels  = group_labels,
      method        = method,
      n_folds_0     = n_folds_cv, ...)
    
    lam <- cv_res$best_lambda
    cve <- extract_best_cve(cv_res)
    if (verbose)
      cat(sprintf("[mrcox] best lambda = %g  |  best CVE = %.4f\n", lam, cve))
    
    perf <- simulations_fit_and_performance(
      datas         = datas,
      p = p, k = k,
      n_simulations = n_simulations,
      r             = 1,    
      lambda_alpha  = lam,
      group_labels  = group_labels,
      method        = method, ...)
    
    if (verbose && !is.null(perf))
      cat(sprintf("[mrcox] mean MSE = %.6f\n", mean(perf$mse1, na.rm = TRUE)))
    
    return(list(
      best_rank   = NULL,  
      best_lambda = lam,
      best_cve    = cve,
      best_perf   = perf,
      all_results = list(),
      cv_result   = cv_res     
    ))
  }
  
  
  best_rank   <- NULL
  best_perf   <- NULL
  best_lambda <- NULL
  best_cve    <- Inf
  all_results <- list()
  
  for (r in rank_vector) {
    
    if (verbose)
      cat(sprintf("\n========== [%s] Trying rank %d ==========\n", method, r))
    
    if (method == "pen") {
      if (is.null(cv_data$dlong_list))
        stop("method='pen' needs dlong_list; use store_dlong=TRUE.")
      
      n_cv <- length(cv_data$y_list_all)
      optimal_lambdas_pen <- numeric(0)
      results_df_pen      <- data.frame()
      
      for (j in seq_len(n_cv)) {
        cat(sprintf("\n=== [pen] Dataset %d/%d ===\n", j, n_cv))
        t0 <- Sys.time()
        o_res_j <- select_lambda_pen_vvh(
          dlong         = cv_data$dlong_list[[j]],
          data_wide     = cv_data$data_wide_list[[j]],
          R             = r,
          lambda_vector = lambda_vector,
          pred_names    = paste0("x", seq_len(p)),
          k             = k,
          n_folds       = n_folds_cv, ...)
        optimal_lambdas_pen <- c(optimal_lambdas_pen, o_res_j$best_lambda)
        results_df_pen <- rbind(results_df_pen,
                                data.frame(dataset  = j,
                                           lambda   = lambda_vector,
                                           CVE_mean = o_res_j$CVE_mean,
                                           CVE_se   = o_res_j$CVE_se))
        cat(sprintf("  -> dataset %d | best_lambda: %g | time: %.2f sec\n",
                    j, o_res_j$best_lambda,
                    as.numeric(difftime(Sys.time(), t0, units = "secs"))))
      }
      
      best_lam_pen <- lambda_vector[which.min(
        abs(lambda_vector - median(optimal_lambdas_pen)))]
      cat(sprintf("\n[pen VVH] global best lambda: %g\n", best_lam_pen))
      
      cv_res <- list(
        results_df      = results_df_pen,
        best_lambda     = best_lam_pen,
        optimal_lambdas = optimal_lambdas_pen,
        lambda          = lambda_vector,
        CVE_mean        = tapply(results_df_pen$CVE_mean, results_df_pen$lambda,
                                 mean, na.rm = TRUE)[as.character(lambda_vector)],
        CVE_se          = tapply(results_df_pen$CVE_mean, results_df_pen$lambda, 
                                 function(x) sd(x, na.rm=TRUE)/sqrt(sum(!is.na(x)))
        )[as.character(lambda_vector)]
      )
      
    } else {
      cv_res <- simulations_fit_find_best_lambda(
        datasets_list = cv_data,
        lambda_vector = lambda_vector,
        p = p, k = k, r = r,
        group_labels  = group_labels,
        method        = method,
        n_folds_0     = n_folds_cv, ...)
    }
    
    lam_r <- cv_res$best_lambda
    cve_r <- extract_best_cve(cv_res)
    
    if (verbose)
      cat(sprintf("[%s rank %d] best lambda = %g  |  best CVE = %.4f\n",
                  method, r, lam_r, cve_r))
    
    perf_r <- simulations_fit_and_performance(
      datas         = datas,
      p = p, k = k,
      n_simulations = n_simulations,
      r             = r,
      lambda_alpha  = lam_r,
      group_labels  = group_labels,
      method        = method, ...)
    
    if (is.null(perf_r)) {
      cat(sprintf("[%s rank %d] all simulations failed, skipping.\n", method, r))
      next
    }
    
    mse_r <- mean(perf_r$mse1, na.rm = TRUE)
    if (verbose)
      cat(sprintf("[%s rank %d] mean MSE = %.6f  (for reference)\n", method, r, mse_r))
    
    all_results[[as.character(r)]] <- list(
      rank        = r,
      lambda      = lam_r,
      best_cve    = cve_r,
      mean_mse    = mse_r,
      performance = perf_r,
      cv_result   = cv_res
    )
    
    if (cve_r < best_cve) {
      best_cve    <- cve_r
      best_rank   <- r
      best_perf   <- perf_r
      best_lambda <- lam_r
      if (verbose)
        cat(sprintf("[%s rank %d] CVE improved, continuing.\n", method, r))
    } else {
      if (verbose)
        cat(sprintf("[%s rank %d] CVE %.4f >= best %.4f at rank %d — stopping.\n",
                    method, r, cve_r, best_cve, best_rank))
      break
    }
  }
  
  if (verbose)
    cat(sprintf(
      "\n[%s] Auto rank done: best rank = %d | lambda = %g | CVE = %.4f | MSE = %.6f\n",
      method, best_rank, best_lambda, best_cve,
      mean(best_perf$mse1, na.rm = TRUE)))
  
  return(list(
    best_rank   = best_rank,
    best_lambda = best_lambda,
    best_cve    = best_cve,
    best_perf   = best_perf,
    all_results = all_results
  ))
}

## Extract comparison df at a specific rank
comparison_at_rank <- function(auto_list, r, method_names) {
  rows <- lapply(seq_along(auto_list), function(i) {
    if (length(auto_list[[i]]$all_results) == 0)
      return(NULL)
    entry <- auto_list[[i]]$all_results[[as.character(r)]]
    if (is.null(entry))
      return(data.frame(method = method_names[i], rank = r,
                        mean_bias = NA, mean_mse = NA, mean_corr = NA))
    perf <- entry$performance
    data.frame(method    = method_names[i],
               rank      = r,
               mean_bias = mean(perf$bias1,         na.rm = TRUE),
               mean_mse  = mean(perf$mse1,          na.rm = TRUE),
               mean_corr = mean(perf$cor_coeff_vec, na.rm = TRUE))
  })
  do.call(rbind, Filter(Negate(is.null), rows))
}

## Extract comparison df at best (auto-selected) rank
comparison_best_rank <- function(auto_list, method_names) {
  rows <- lapply(seq_along(auto_list), function(i) {
    res  <- auto_list[[i]]
    perf <- res$best_perf
    if (is.null(perf))
      return(data.frame(method    = method_names[i], best_rank = NA,
                        mean_bias = NA, mean_mse  = NA, mean_corr = NA))
    
    rank_label <- if (is.null(res$best_rank)) "full" else res$best_rank
    
    data.frame(method    = method_names[i],
               best_rank = rank_label,
               mean_bias = mean(perf$bias1,         na.rm = TRUE),
               mean_mse  = mean(perf$mse1,          na.rm = TRUE),
               mean_corr = mean(perf$cor_coeff_vec, na.rm = TRUE))
  })
  do.call(rbind, rows)
}

print_rank_comparisons <- function(auto_list, method_names, scenario_name) {
  cat(sprintf("\n============================================================\n"))
  cat(sprintf("  %s\n", scenario_name))
  cat(sprintf("============================================================\n"))
  
  ranks_tried <- sort(unique(unlist(
    lapply(auto_list, function(x)
      if (length(x$all_results) > 0) as.integer(names(x$all_results)) else NULL
    ))))
  
  for (r in ranks_tried) {
    cat(sprintf("\n--- Rank %d ---\n", r))
    df <- comparison_at_rank(auto_list, r, method_names)
    if (!is.null(df) && nrow(df) > 0) print(df, row.names = FALSE)
  }
  
  cat("\n--- Best Rank (auto-selected per method) ---\n")
  print(comparison_best_rank(auto_list, method_names), row.names = FALSE)
}

## Function to plot simulated data
generate_data_plots <- function(dlong_list, n, r, p, k) {
  
  combined_data <- if (is.data.frame(dlong_list)) {
    dlong_list
  } else {
    do.call(rbind, dlong_list)
  }
  
  predictor_names <- grep("^x", names(combined_data), value = TRUE)
  predictors_long <- reshape2::melt(combined_data, id.vars = "id",
                                    measure.vars = predictor_names)
  
  p1 <- ggplot(predictors_long, aes(x = value, color = variable, fill = variable)) +
    geom_density(alpha = 0.05) +
    labs(title = "Distribution of Predictors", x = "Value", y = "Density") +
    theme_minimal(base_size = 13) +
    theme(legend.title = element_blank())
  
  p2 <- ggplot(combined_data, aes(x = Tstop)) +
    geom_density(fill = "steelblue", alpha = 0.4) +
    labs(title = "Distribution of Tstop", x = "Tstop", y = "Density") +
    theme_minimal(base_size = 13)
  
  top <- textGrob(
    paste0("n = ", n, "  |  p = ", p, "  |  k = ", k, " outcomes"),
    gp = gpar(fontsize = 16, fontface = "bold")
  )
  grid_plots <- grid.arrange(p1, p2, ncol = 2, top = top)
  
  status_pct <- prop.table(table(combined_data$status)) * 100
  cat(sprintf("Censored: %.1f%%   Event: %.1f%%\n",
              status_pct["0"], status_pct["1"]))
  invisible(grid_plots)
}


## Heatmaps
beta_heatmap = function(title,datas1, performance_results11,minb,maxb,minmax=TRUE) {
  
  bias_list = performance_results11$bias_for_all_coeffs_and_simulations
  Beta_true = datas1$Beta_matrix
  Beta_list_estimated = list()
  
  for(i in 1:length(datas1$y_list_all)){
    Beta_list_estimated[[i]] = bias_list[[i]] + Beta_true
  }
  
  sum_betas_estimated <- Reduce("+", Beta_list_estimated) #sum element wise all matrices in the list
  
  # Calculate the mean
  mean_beta <- sum_betas_estimated / length(Beta_list_estimated)
  
  colnames(mean_beta) <- paste("k=", 1:ncol(mean_beta), sep="")
  rownames(mean_beta) <- paste("X", 1:nrow(mean_beta), sep="")
  
  beta_df <- reshape2::melt(mean_beta)
  
  
  
  
  if(minmax==TRUE){
    breaks <- seq(minb, maxb, length.out = 101)
  }
  else{
    breaks <- seq(min(mean_beta), max(mean_beta), length.out = 101) 
  }
  
  norm_breaks <- (breaks - min(breaks)) / (max(breaks) - min(breaks))
  my_palette <- colorRampPalette(c("lightgrey", "orange", "black"))(100)
  # Generate the heatmap
  # heatmap.2(mean_beta, Rowv=F, Colv=F, scale = "none", col = my_palette,
  #           density.info = "none", trace = "none", dendrogram = "none",breaks=breaks)
  
  p <- ggplot(beta_df, aes(x = Var2, y = Var1, fill = value)) +
    geom_tile() +
    scale_fill_gradientn(colors = my_palette, 
                         values = norm_breaks, # specify normalized breaks here
                         limits = range(breaks)) + # Set limits to the min and max of your breaks
    theme_minimal() +
    labs(x = "", y = "", fill = "Value") +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, size = 15),                        
      axis.title.x = element_text(size = 15),                       
      axis.title.y = element_text(size = 15),                       
      legend.text = element_text(size = 20),                         
      legend.title = element_text(size = 20),                       
      plot.title = element_text(size = 30),           
      strip.text = element_text(size = 20)                           
    )+
    ggtitle(title)
  
  print(min(mean_beta))
  print(max(mean_beta))
  print(mean_beta)
  return(p)
}

compare_beta_heatmap <- function(title, datas, perf_results, minb, maxb) {
  
  # True Beta
  Beta_true <- datas$Beta_matrix  
  colnames(Beta_true) <- paste0("k=", 1:ncol(Beta_true))
  rownames(Beta_true) <- paste0("X",  1:nrow(Beta_true))
  
  breaks      <- seq(minb, maxb, length.out=101)
  norm_breaks <- (breaks - min(breaks)) / (max(breaks) - min(breaks))
  my_palette  <- colorRampPalette(c("lightgrey","orange","black"))(100)
  
  p_true <- ggplot(reshape2::melt(Beta_true), aes(x=Var2, y=Var1, fill=value)) +
    geom_tile() +
    scale_fill_gradientn(colors=my_palette, values=norm_breaks,
                         limits=range(breaks)) +
    theme_minimal() +
    labs(x="", y="", fill="Value", title="True Beta") +
    theme(axis.text.x  = element_text(angle=90, hjust=1, size=12),
          axis.text.y  = element_blank(), 
          axis.ticks.y = element_blank(),
          plot.title   = element_text(size=14, face="bold"),
          legend.text  = element_text(size=12),
          legend.title = element_text(size=12))
  
  # Estimated Beta 
  p_est <- beta_heatmap("Estimated Beta (mean)",
                        datas, perf_results,
                        minb=minb, maxb=maxb, minmax=TRUE) +
    theme(axis.text.y  = element_blank(),
          axis.ticks.y = element_blank(),
          plot.title   = element_text(size=14, face="bold"))
  
  combined_plot <- p_true + p_est + 
    plot_annotation(title = title, 
                    theme = theme(plot.title = element_text(size=15, face="bold", hjust=0.5)))
  
  return(combined_plot)
}

## Function to plot lambdas
lambda_plot <- function(o_result, title = "V&VH Lambda Selection",
                        ylim = NULL) {  
  
  if (is.data.frame(o_result)) {
    df          <- o_result
    best_lambda <- df$lambda[which.min(
      tapply(df$CVE_mean, df$lambda, mean, na.rm=TRUE))]
  } else {
    df          <- o_result$results_df
    best_lambda <- o_result$best_lambda
  }
  
  df$dataset <- as.factor(df$dataset)
  
  best_per_dataset <- df %>%
    group_by(dataset) %>%
    filter(CVE_mean == min(CVE_mean, na.rm=TRUE)) %>%
    ungroup()
  
  p <- ggplot(df, aes(x=lambda, y=CVE_mean,
                      color=dataset, group=dataset)) +
    geom_line(alpha=0.6) +
    geom_point(size=1, alpha=0.6) +
    geom_point(data=best_per_dataset,
               aes(x=lambda, y=CVE_mean),
               color="red", size=3, shape=16) +
    geom_vline(xintercept=best_lambda,
               linetype="dashed", color="black", linewidth=1) +
    labs(title   = title,
         x       = expression(lambda),
         y       = "Mean CVE (V&VH)",
         color   = "Dataset",
         caption = sprintf("Red: per-dataset min  |  Black dashed: global best λ = %g",
                           best_lambda)) +
    theme_minimal() +
    theme(legend.position="right",
          plot.caption=element_text(size=10))
  
  if (!is.null(ylim))
    p <- p + coord_cartesian(ylim=ylim)  
  
  return(p)
}

## Function to compare rank
select_lambda_by_rank <- function(X, y_list, status_list,
                                  rank_vector,      # e.g. c(1, 2, 3)
                                  lambda_vector,
                                  group_labels,
                                  method = "rrr_grplasso",
                                  n_folds = 5, seed = 123, ...) {
  set.seed(seed)
  N       <- nrow(X)
  fold_id <- sample(rep(seq_len(n_folds), length.out = N))
  results <- list()
  for (R in rank_vector) {
    cat(sprintf("  Rank %d ...\n", R))
    n_lam   <- length(lambda_vector)
    CVE_mat <- matrix(NA, n_folds, n_lam)
    for (fold in seq_len(n_folds)) {
      idx_te <- which(fold_id == fold)
      idx_tr <- which(fold_id != fold)
      X_tr <- X[idx_tr, , drop=FALSE]; X_te <- X[idx_te, , drop=FALSE]
      y_tr <- lapply(y_list,      `[`, idx_tr); y_te <- lapply(y_list,      `[`, idx_te)
      s_tr <- lapply(status_list, `[`, idx_tr); s_te <- lapply(status_list, `[`, idx_te)
      for (lam in seq_along(lambda_vector)) {
        fit <- tryCatch(
          fit_one_model(X=X_tr, y_list=y_tr, status_list=s_tr,
                        method=method, r=R,
                        lambda=lambda_vector[lam],
                        group_labels=group_labels, ...),
          error = function(e) NULL)
        if (is.null(fit) || !all(is.finite(fit$B_hat))) next
        cve <- compute_vvh_fold(X_tr, X_te, y_tr, y_te, s_tr, s_te, fit$B_hat)
        if (!is.na(cve)) CVE_mat[fold, lam] <- cve
      }
    }
    CVE_mean <- colMeans(CVE_mat, na.rm=TRUE)
    CVE_se   <- apply(CVE_mat, 2, sd, na.rm=TRUE) / sqrt(n_folds)
    results[[as.character(R)]] <- data.frame(
      R        = R,
      lambda   = lambda_vector,
      CVE_mean = CVE_mean,
      CVE_se   = CVE_se
    )
  }
  do.call(rbind, results)
}

select_lambda_by_rank_multi <- function(datas, rank_vector, lambda_vector,
                                        group_labels, method="rrr_grplasso",
                                        n_datasets = 10, 
                                        n_folds=5, seed=123, ...) {
  all_results <- list()
  
  for (j in seq_len(n_datasets)) {
    cat(sprintf("=== Dataset %d/%d ===\n", j, n_datasets))
    res_j <- select_lambda_by_rank(
      X           = datas$X_list[[j]],
      y_list      = datas$y_list_all[[j]],
      status_list = datas$status_list_all[[j]],
      rank_vector  = rank_vector,
      lambda_vector = lambda_vector,
      group_labels  = group_labels,
      method        = method,
      n_folds       = n_folds,
      seed          = seed + j,  
      ...
    )
    res_j$dataset <- j
    all_results[[j]] <- res_j
  }
  
  do.call(rbind, all_results) %>%
    group_by(R, lambda) %>%
    summarise(CVE_mean = mean(CVE_mean, na.rm=TRUE),
              CVE_se   = sd(CVE_mean,  na.rm=TRUE) / sqrt(n()),
              .groups  = "drop")
}


plot_rank_comparison <- function(rank_cve_df, title = "Group Lasso: CVE by rank and lambda") {
  rank_cve_df$R <- factor(rank_cve_df$R)
  best_per_rank <- rank_cve_df %>%
    group_by(R) %>%
    filter(CVE_mean == min(CVE_mean, na.rm=TRUE)) %>%
    ungroup()
  ggplot(rank_cve_df, aes(x=log(lambda), y=CVE_mean,
                          color=R, group=R)) +
    geom_line(linewidth=0.9) +
    geom_point(size=1.5) +
    geom_point(data=best_per_rank,
               aes(x=log(lambda), y=CVE_mean),
               shape=23, size=4, fill="white", stroke=1.5) +
    scale_color_manual(values=c("1"="tomato","2"="steelblue","3"="darkgreen"),
                       name="Fitted rank") +
    labs(title    = title,
         subtitle = "True rank = 2  |  Diamond: best lambda per rank",
         x        = expression(log(lambda)),
         y        = "Mean CVE (V&VH)") +
    theme_minimal(base_size=12) +
    theme(legend.position="right")
}

lambda_plot_by_rank <- function(auto_result, method_name,
                                title = NULL, ylim = NULL) {
  
  rank_colors <- c("Rank 1" = "tomato", "Rank 2" = "steelblue",
                   "Rank 3" = "darkgreen", "Rank 4" = "purple",
                   "Full"   = "black")
  
  # ── mrcox: all_results 为空，直接从 best_perf 对应的 cv_result 画 ──────
  if (length(auto_result$all_results) == 0) {
    cv_res <- auto_result$cv_result  # mrcox 返回时没有存，需补充（见下注）
    
    # 如果 cv_res 不存在就无法画图
    if (is.null(cv_res)) {
      cat("No CV results to plot for", method_name, "\n")
      return(invisible(NULL))
    }
    
    # 提取 df
    if (!is.null(cv_res$results_df)) {
      df_all <- cv_res$results_df %>%
        group_by(lambda) %>%
        summarise(CVE_mean = mean(CVE_mean, na.rm = TRUE),
                  CVE_se   = mean(CVE_se,   na.rm = TRUE),
                  .groups  = "drop") %>%
        mutate(rank = "Full") %>%
        as.data.frame()
    } else if (!is.null(cv_res$lambda)) {
      df_all <- data.frame(rank     = "Full",
                           lambda   = cv_res$lambda,
                           CVE_mean = cv_res$CVE_mean,
                           CVE_se   = cv_res$CVE_se)
    } else {
      cat("No plottable CV data found for", method_name, "\n")
      return(invisible(NULL))
    }
    
    df_all$rank      <- factor(df_all$rank, levels = "Full")
    best_per_rank    <- df_all %>% filter(CVE_mean == min(CVE_mean, na.rm = TRUE)) %>% slice(1)
    global_best_lam  <- auto_result$best_lambda
    
    p <- ggplot(df_all, aes(x = lambda, y = CVE_mean, color = rank, group = rank)) +
      geom_line(linewidth = 0.9) +
      geom_point(size = 1.5) +
      geom_point(data = best_per_rank, aes(x = lambda, y = CVE_mean),
                 shape = 23, size = 4, fill = "white", stroke = 1.5) +
      geom_vline(xintercept = global_best_lam,
                 linetype = "dashed", color = "black", linewidth = 1) +
      scale_color_manual(values = rank_colors, name = "Rank") +
      labs(title    = title %||% paste0("[", method_name, "] CVE by lambda"),
           subtitle = paste0("Diamond: best lambda  |  Black dashed: \u03bb = ",
                             round(global_best_lam, 4)),
           x = expression(lambda), y = "Mean CVE (V&VH)") +
      theme_minimal(base_size = 12) +
      theme(legend.position = "right")
    
    if (!is.null(ylim)) p <- p + coord_cartesian(ylim = ylim)
    return(p)
  }
  
  # ── 其他方法：正常 rank 循环结果 ─────────────────────────────────────
  ranks_tried <- names(auto_result$all_results)
  
  df_all <- do.call(rbind, lapply(ranks_tried, function(r) {
    cv_res <- auto_result$all_results[[r]]$cv_result
    
    if (!is.null(cv_res$lambda) && length(cv_res$CVE_mean) > 1) {
      return(data.frame(rank     = paste0("Rank ", r),
                        lambda   = cv_res$lambda,
                        CVE_mean = cv_res$CVE_mean,
                        CVE_se   = cv_res$CVE_se))
    }
    
    if (!is.null(cv_res$results_df)) {
      return(
        cv_res$results_df %>%
          group_by(lambda) %>%
          summarise(CVE_mean = mean(CVE_mean, na.rm = TRUE),
                    CVE_se   = mean(CVE_se,   na.rm = TRUE),
                    .groups  = "drop") %>%
          mutate(rank = paste0("Rank ", r)) %>%
          as.data.frame()
      )
    }
    
    warning(sprintf("rank %s: unrecognised cv_result structure, skipping.", r))
    return(NULL)
  }))
  
  if (is.null(df_all) || nrow(df_all) == 0) {
    cat("No plottable CV data found.\n"); return(invisible(NULL))
  }
  
  df_all$rank <- factor(df_all$rank,
                        levels = paste0("Rank ", sort(as.integer(ranks_tried))))
  
  best_per_rank   <- df_all %>%
    group_by(rank) %>%
    filter(CVE_mean == min(CVE_mean, na.rm = TRUE)) %>%
    slice(1) %>%
    ungroup()
  
  global_best_lam <- auto_result$best_lambda
  best_rank_label <- if (is.null(auto_result$best_rank)) "full" else auto_result$best_rank
  
  p <- ggplot(df_all, aes(x = lambda, y = CVE_mean, color = rank, group = rank)) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 1.5) +
    geom_point(data = best_per_rank, aes(x = lambda, y = CVE_mean),
               shape = 23, size = 4, fill = "white", stroke = 1.5) +
    geom_vline(xintercept = global_best_lam,
               linetype = "dashed", color = "black", linewidth = 1) +
    scale_color_manual(values = rank_colors, name = "Rank") +
    labs(title    = title %||% paste0("[", method_name, "] CVE by lambda and rank"),
         subtitle = paste0("Diamond: best lambda per rank  |  ",
                           "Black dashed: global best  (rank ",
                           best_rank_label, ", \u03bb = ",
                           round(global_best_lam, 4), ")"),
         x = expression(lambda), y = "Mean CVE (V&VH)") +
    theme_minimal(base_size = 12) +
    theme(legend.position = "right")
  
  if (!is.null(ylim)) p <- p + coord_cartesian(ylim = ylim)
  return(p)
}

`%||%` <- function(a, b) if (!is.null(a)) a else b

## plot beta
plot_B_true_boxplot <- function(datas,
                                outcomes     = "all",
                                group_labels = NULL,
                                title        = NULL) {
  
  Beta_list <- datas$Beta_matrix_list
  n_sims    <- length(Beta_list)
  p         <- nrow(Beta_list[[1]])
  K         <- ncol(Beta_list[[1]])
  
  if (identical(outcomes, "all")) outcomes <- seq_len(K)
  
  # group_labels: caller > datas > uniform group-1
  if (is.null(group_labels)) {
    group_labels <- if (!is.null(datas$group_labels))
      datas$group_labels
    else
      rep(1L, p)
  }
  
  grp_ids <- sort(unique(group_labels))
  n_grps  <- length(grp_ids)
  
  # ── 1. Build long data frame, aggregated by group ─────────────────────────
  df <- do.call(rbind, lapply(seq_len(n_sims), function(i) {
    B <- Beta_list[[i]]
    do.call(rbind, lapply(outcomes, function(k) {
      # For each group, collect all Beta values across predictors in that group
      data.frame(
        group     = factor(paste0("G", group_labels),
                           levels = paste0("G", grp_ids)),
        group_num = group_labels,
        B_true    = B[, k],
        outcome   = factor(paste0("Outcome ", k),
                           levels = paste0("Outcome ", outcomes)),
        sim       = i,
        stringsAsFactors = FALSE
      )
    }))
  }))
  
  # ── 2. Summary for group labels (mean β per group, for annotation) ────────
  grp_summary <- df %>%
    dplyr::group_by(group, outcome) %>%
    dplyr::summarise(
      mean_beta = mean(B_true, na.rm = TRUE),
      .groups = "drop"
    )
  
  # ── 3. Build plot ─────────────────────────────────────────────────────────
  is_random <- isTRUE(datas$random_beta)
  
  p_out <- ggplot(df, aes(x = group, y = B_true)) +
    
    # zero-reference
    geom_hline(yintercept = 0, linetype = "dashed",
               linewidth  = 0.6, color = "grey50") +
    
    # main layer (boxplot for random, point for fixed)
    { if (is_random) {
      geom_boxplot(outlier.size  = 0.8,
                   outlier.alpha = 0.5,
                   fill          = "white",
                   color         = "black",
                   alpha         = 0.8,
                   width         = 0.7)
    } else {
      geom_point(size = 2.5, color = "steelblue", alpha = 0.8)
    }
    } +
    
    # group means as points (optional, matches figure style)
    geom_point(data = grp_summary,
               aes(x = group, y = mean_beta),
               color = "red", size = 2.5, shape = 16) +
    
    # outcome faceting (if multiple outcomes)
    { if (length(outcomes) > 1L)
      facet_wrap(~outcome, ncol = 1, scales = "free_y")
    } +
    
    labs(
      title    = title %||% "Distribution of true β-values",
      subtitle = paste0(
        "g = ", n_grps,
        if (!is.null(datas$g_nz)) paste0(", g_nz = ", datas$g_nz),
        ", p = ", p,
        if (!is.null(datas$rho_w)) paste0(", ρ_w = ", datas$rho_w),
        if (!is.null(datas$rho_b)) paste0(", ρ_b = ", datas$rho_b),
        if (!is.null(datas$t_snr)) paste0(", t_snr = ", datas$t_snr)
      ),
      x = "Group Number",
      y = "True coefficient value"
    ) +
    
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x        = element_text(size = 10, face = "bold"),
      axis.title.x       = element_text(size = 11, margin = margin(t = 10)),
      axis.title.y       = element_text(size = 11, margin = margin(r = 10)),
      panel.grid.major.x = element_blank(),
      panel.grid.minor   = element_blank(),
      strip.text         = element_text(face = "bold", size = 10),
      plot.title         = element_text(face = "bold", size = 14, hjust = 0.5),
      plot.subtitle      = element_text(size = 9, color = "grey30", hjust = 0.5),
      plot.margin        = margin(t = 10, r = 10, b = 10, l = 10)
    )
  
  return(p_out)
}
