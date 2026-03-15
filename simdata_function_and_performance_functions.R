library(MASS)
library(survival)
library(ggplot2)
library(gridExtra)
library(reshape2)
library(grid)
library(glmnet)
library(dplyr)
library(scales)

rank_check <- function(Beta_m) {
  sv    <- svd(Beta_m)$d
  ratio <- sv[2] / sv[1]
  cat(sprintf("Singular values: %s\nRatio s2/s1: %.4f\n",
              paste(round(sv, 4), collapse = ", "), ratio))
  if (ratio < 0.01) cat("Matrix is close to rank 1.\n") else
    cat("Matrix is NOT close to rank 1.\n")
}
## Function to generate covariate matrix, A,$\Gamma$ and B

add_zeros <- function(matrix, percentage = 0.10) {
  total_components <- length(matrix)
  zeros_count <- floor(total_components * percentage)
  
  indexes <- 1:total_components
  
  zero_indexes <- sample(indexes, zeros_count)
  
  matrix[zero_indexes] <- 0
  
  return(matrix)
}


simulate_Beta_matrix = function(p,r,k,case=1,scale=1,n_groups=NULL,zero_blocks=NULL){
  if (case == 1){
    #GENERATE A, \Gamma and B matrix
    #A matrix
    A_matrix <- matrix(NA, nrow = p, ncol = r)
    for (i in 1:r) {
      A_matrix [, i] <- runif(p, 0,1)
    }
    #Gamma matrix
    Gamma_matrix <- matrix(NA, nrow = k, ncol = r)
    for (i in 1:r) {
      Gamma_matrix [, i] <- runif(k, 0,1)
    }
  }
  else if(case ==2){
    if (p %% n_groups != 0) {
      stop(sprintf("p (%d) must be divisible by n_groups (%d)", p, n_groups))
    }
    group_size <- p/n_groups
    
    A_matrix <- matrix(runif(p*r,0,1),p,r)
    if (!is.null(zero_blocks)) {
      for (blk in zero_blocks) {
        g <- blk[1]
        rank_r <- blk[2]
        if (g<1 || g>n_groups)
          stop(sprintf("zero_blocks: group index %d out of range [1, %d]", g, n_groups))
        if (rank_r<1 || rank_r>r)
          stop(sprintf("zero_blocks: rank index %d out of range [1, %d]", rank_r, r))
        rows <- ((g-1)*group_size+1):(g*group_size)
        A_matrix[rows, rank_r] <- 0
      }
    }
    Gamma_matrix <- matrix(runif(k*r,0,1), nrow = k, ncol = r)
  }
  
  else if(case ==3){
    #GENERATE A, \Gamma and B matrix
    #A matrix
    A_matrix <- matrix(NA, nrow = p, ncol = r)
    for (i in 1:r) {
      A_matrix [, i] <- runif(p, 0,1)
    }
    #Gamma matrix
    Gamma_matrix <- matrix(NA, nrow = k, ncol = r)
    for (i in 1:r) {
      Gamma_matrix [, i] <- runif(k, 0,1)
    }
    A_matrix <- add_zeros(A_matrix, percentage = 0.10)
    Gamma_matrix <- add_zeros(Gamma_matrix, percentage = 0.10)
  }
  
  #Coefficient matrix B (p x k)
  Beta_matrix = A_matrix %*% t(Gamma_matrix)
  return(list(
    Beta_matrix = Beta_matrix*scale,
    A_matrix = A_matrix,
    Gamma_matrix = Gamma_matrix,
    group_size = if (case==2) p/n_groups else NULL,
    n_groups = if (case==2) n_groups else NULL
  ))
}

## Function to generate group_labels
gen_group_labels <- function(p,n_groups) {
  group_size <- p/n_groups
  rep(1:n_groups, each=group_size)
}

## Function to simulate Dataset

# X ~ N(0,1), T ~ Exp(exp(X %*% B)), C ~ Uniform(0,1)
# T_obs = min(T, C), delta = I(T <= C)
simulate_data_model = function(Beta_matrix, n, corr=0, group_corr = 0, group_labels = NULL, Print_flug=FALSE){
  
  if (group_corr > 0 && !is.null(group_labels)) {
    # 组内相关，组间独立
    n_groups   <- max(group_labels)
    group_size <- p / n_groups
    Sigma      <- diag(p)
    
    for (g in seq_len(n_groups)) {
      idx <- which(group_labels == g)
      # 组内用复合对称结构：对角为1，组内交叉项为 group_corr
      Sigma[idx, idx] <- group_corr
      diag(Sigma)[idx] <- 1
    }
    X <- mvtnorm::rmvnorm(n, mean = rep(0, p), sigma = Sigma)
    
  } else if (corr > 0) {
    # 原来的全局相关
    sigma <- matrix(corr, p, p); diag(sigma) <- 1
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


## Extract coefficients from outputs

extract_coef <- function(model, lambda_idx = 1) {
  if (!inherits(model, "survRRR")) stop("model must be survRRR object")
  res <- model$result[[lambda_idx]]
  Alpha <- res$alpha   # p x R
  Gamma <- res$Gamma   # K x R
  return(Alpha %*% t(Gamma))  # p x K
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
  
  # SVD 分解，取前 R 个右奇异向量作为 Gamma 初始值
  sv      <- svd(B0, nu = 0, nv = min(R, ncol(B0)))
  Gamma0  <- sv$v[, 1:R, drop = FALSE]   # k x R
  
  # pen 方法的 Gamma 维度是 R x K，需转置
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
                          method = c("rrr_grplasso", "rrr_lasso", "pen"),
                          r, lambda, group_labels,
                          pen_gamma_start = NULL, ...) {
  method <- match.arg(method)
  
  if (method == "rrr_grplasso") {
    raw   <- solve_reduced_rank(X = X, y_list = y_list,
                                status_list = status_list,
                                R = r, lambda_alpha = lambda,
                                group_labels = group_labels, ...)
    res   <- raw$result[[1]]
    Alpha <- res$alpha; Gamma <- res$Gamma
    B_hat <- Alpha %*% t(Gamma)
    
  } else if (method == "rrr_lasso") {
    raw   <- solve_reduced_rank_lasso(X = X, y_list = y_list,
                                      status_list = status_list,
                                      R = r, lambda_alpha = lambda, ...)
    res   <- raw$result[[1]]
    Alpha <- res$alpha; Gamma <- res$Gamma
    B_hat <- Alpha %*% t(Gamma)
    
  } else {
    if (is.null(dlong)) stop("method='pen' needs dlong")
    
    # 如果没有提供初始值，用 SVD 初始化（比随机更稳定）
    if (is.null(pen_gamma_start)) {
      k_outcomes  <- length(y_list)
      pen_gamma_start <- make_pen_gamma_init(X, y_list, status_list,
                                             R = r, k = k_outcomes)
    }
    
    pred_names  <- grep("^x", names(dlong), value = TRUE)
    formula_str <- paste("Surv(Tstop, status) ~",
                         paste(pred_names, collapse = " + "))
    raw   <- pen.survrrr(as.formula(formula_str),
                         dat          = dlong,
                         R            = r,
                         Gamma.iter   = pen_gamma_start,
                         lambda.alpha = lambda,
                         lambda.gamma = lambda,
                         eps          = 1e-4,
                         maxit = 1e6,
                         standardize.opt = FALSE,
                         alpha        = 1, ...)
    Alpha <- raw$Alpha; Gamma <- raw$Gamma
    B_hat <- Alpha %*% Gamma
  }
  
  return(list(alpha = Alpha, Gamma = Gamma, B_hat = B_hat, raw = raw))
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

simulations_data = function(n, p, r, k, n_simulations,
                            case = 1, scale = 1, corr = 0,group_corr = 0,
                            n_groups = NULL, zero_blocks = NULL,
                            store_dlong = FALSE) {
  beta_obj <- simulate_Beta_matrix(p,r,k,case=case,scale=scale,n_groups=n_groups,zero_blocks=zero_blocks)
  Beta_matrix <- beta_obj$Beta_matrix
  group_labels_vec <- if (case == 2) gen_group_labels(p, n_groups) else NULL
  
  y_list_all <- vector("list", n_simulations)
  status_list_all <- vector("list", n_simulations)
  X_list <- vector("list", n_simulations)
  T_sim_list <- vector("list", n_simulations)
  T_obs_list <- vector("list", n_simulations)
  data_wide_list  <- vector("list", n_simulations)
  dlong_list <- if (store_dlong) vector("list", n_simulations) else NULL
  
  for (i in 1:n_simulations) {
    sim <- simulate_data_model(Beta_matrix,n,corr=corr,group_corr=group_corr,group_labels=group_labels_vec)
    y_list_all[[i]] <- sim$y_list
    status_list_all[[i]] <- sim$status_list
    X_list[[i]] <- sim$X
    T_sim_list[[i]] <- sim$T_sim
    T_obs_list[[i]] <- sim$T_obs
    data_wide_list[[i]] <- sim$data
    
    if (store_dlong) {
      dlong_list[[i]] <- build_dlong(sim$data, k)
    }
  }
  
  return(list(
    y_list_all = y_list_all,
    status_list_all = status_list_all,
    X_list = X_list,
    Beta_matrix = Beta_matrix,
    A_matrix = beta_obj$A_matrix,
    Gamma_matrix = beta_obj$Gamma_matrix,
    T_sim_list = T_sim_list,
    T_obs_list = T_obs_list,
    data_wide_list = data_wide_list,
    dlong_list = dlong_list,
    group_labels = if (case == 2) gen_group_labels(p, n_groups) else NULL,
    n_groups = beta_obj$n_groups,
    group_size = beta_obj$group_size
  ))
}

#################
compute_vvh_fold <- function(X_train, X_test,
                             y_train, y_test,
                             s_train, s_test,
                             B_hat) {
  K      <- ncol(B_hat)
  N_test <- nrow(X_test)
  X_full <- rbind(X_train, X_test)
  
  loglik_full  <- 0
  loglik_train <- 0
  
  for (k in seq_len(K)) {
    y_full <- c(y_train[[k]], y_test[[k]])
    s_full <- c(s_train[[k]], s_test[[k]])
    beta_k <- B_hat[, k]
    
    # ℓ(β̂^{-k})：全数据上的偏似然，β̂ 来自训练集
    dat_full <- cbind(
      data.frame(y = y_full, s = s_full),
      as.data.frame(X_full)
    )
    pred_names <- paste0("xv", seq_len(ncol(X_full)))
    colnames(dat_full)[3:ncol(dat_full)] <- pred_names
    
    formula_full <- paste(
      "survival::coxph(survival::Surv(y, s) ~",
      paste(pred_names, collapse = "+"),
      ", data = dat_full, init = beta_k,",
      "control = survival::coxph.control(iter.max = 0, timefix = FALSE))"
    )
    fit_full <- tryCatch(eval(parse(text = formula_full)), error = function(e) NULL)
    
    # ℓ^{-k}(β̂^{-k})：训练集上的偏似然，同一个 β̂
    dat_train <- cbind(
      data.frame(y = y_train[[k]], s = s_train[[k]]),
      as.data.frame(X_train)
    )
    colnames(dat_train)[3:ncol(dat_train)] <- pred_names
    
    formula_train <- paste(
      "survival::coxph(survival::Surv(y, s) ~",
      paste(pred_names, collapse = "+"),
      ", data = dat_train, init = beta_k,",
      "control = survival::coxph.control(iter.max = 0, timefix = FALSE))"
    )
    fit_train <- tryCatch(eval(parse(text = formula_train)), error = function(e) NULL)
    
    if (!is.null(fit_full) && !is.null(fit_train)) {
      loglik_full  <- loglik_full  + fit_full$loglik[2]
      loglik_train <- loglik_train + fit_train$loglik[2]
    }
  }
  
  # V&VH 公式：-2 * (ℓ_full - ℓ_train)，除以测试集大小归一化
  return(-2 * (loglik_full - loglik_train) / N_test)
}
## ============================================================
## select_lambda_vvh：V&VH 准则选 lambda（修正版）
## ============================================================
select_lambda_vvh <- function(X, y_list, status_list, R,
                              lambda_vector, group_labels,
                              method  = c("rrr_grplasso", "rrr_lasso"),
                              n_folds = 5, seed = 123, ...) {
  method  <- match.arg(method)
  set.seed(seed)
  N       <- nrow(X)
  n_lam   <- length(lambda_vector)
  fold_id <- sample(rep(seq_len(n_folds), length.out = N))
  CVE_mat <- matrix(NA, n_folds, n_lam)
  
  for (fold in seq_len(n_folds)) {
    idx_te  <- which(fold_id == fold)
    idx_tr  <- which(fold_id != fold)
    
    X_tr  <- X[idx_tr, , drop = FALSE]
    X_te  <- X[idx_te, , drop = FALSE]
    y_tr  <- lapply(y_list,      `[`, idx_tr)
    y_te  <- lapply(y_list,      `[`, idx_te)
    s_tr  <- lapply(status_list, `[`, idx_tr)
    s_te  <- lapply(status_list, `[`, idx_te)
    
    for (lam in seq_along(lambda_vector)) {
      fit <- tryCatch(
        fit_one_model(
          X = X_tr, y_list = y_tr, status_list = s_tr,
          method = method, r = R,
          lambda = lambda_vector[lam],
          group_labels = group_labels, ...
        ),
        error = function(e) NULL
      )
      if (is.null(fit)) next
      
      CVE_mat[fold, lam] <- compute_vvh_fold(
        X_tr, X_te, y_tr, y_te, s_tr, s_te, fit$B_hat
      )
    }
    cat(sprintf("  fold %d/%d done\n", fold, n_folds))
  }
  
  CVE_mean    <- colMeans(CVE_mat, na.rm = TRUE)
  CVE_se      <- apply(CVE_mat, 2, sd, na.rm = TRUE) / sqrt(n_folds)
  best_lambda <- lambda_vector[which.min(CVE_mean)]  # 越小越好
  
  cat(sprintf("[%s VVH] best lambda: %g  (CVE = %.4f)\n",
              method, best_lambda, min(CVE_mean, na.rm = TRUE)))
  
  return(list(
    CVE_mean    = CVE_mean,
    CVE_se      = CVE_se,
    CVE_mat     = CVE_mat,
    lambda      = lambda_vector,
    best_lambda = best_lambda
  ))
}
## ============================================================
## select_lambda_pen_vvh：pen 方法的 V&VH 版本
## 利用已有的 calculate.log.partial.lik
## ============================================================
select_lambda_pen_vvh <- function(dlong, data_wide, R,
                                  lambda_vector, pred_names, k,
                                  n_folds = 5, seed = 123, ...) {
  set.seed(seed)
  subject_ids <- unique(dlong$id)
  N_subj      <- length(subject_ids)
  fold_id     <- setNames(
    sample(rep(seq_len(n_folds), length.out = N_subj)),
    as.character(subject_ids)
  )
  
  n_lam   <- length(lambda_vector)
  CVE_mat <- matrix(NA, n_folds, n_lam)
  
  for (fold in seq_len(n_folds)) {
    test_ids  <- names(fold_id)[fold_id == fold]
    train_ids <- names(fold_id)[fold_id != fold]
    N_test    <- length(test_ids)
    
    dlong_train <- dlong[dlong$id %in% train_ids, ]
    dlong_full  <- dlong
    
    data_train     <- data_wide[data_wide$id %in% train_ids, ]
    dlong_train_ms <- build_dlong(data_train, k)
    
    # ── SVD 初始化（大规模必须，每折算一次，所有 lambda 共用）──
    X_tr   <- as.matrix(data_train[, pred_names])
    y_tr   <- lapply(seq_len(k), function(i) data_train[[paste0("t", i)]])
    s_tr   <- lapply(seq_len(k), function(i) data_train[[paste0("d", i)]])
    gamma_init <- make_pen_gamma_init(X_tr, y_tr, s_tr, R = R, k = k)
    
    formula_str <- paste("Surv(Tstop, status) ~",
                         paste(pred_names, collapse = "+"))
    
    for (lam in seq_along(lambda_vector)) {
      raw_fit <- tryCatch(
        pen.survrrr(
          as.formula(formula_str),
          dat          = dlong_train_ms,
          R            = R,
          Gamma.iter   = gamma_init,      # ← SVD 初始化
          lambda.alpha = lambda_vector[lam],
          lambda.gamma = lambda_vector[lam],
          eps = 1e-4, standardize.opt = FALSE, alpha = 1, ...
        ),
        error = function(e) NULL
      )
      if (is.null(raw_fit)) next
      
      B_hat        <- raw_fit$Alpha %*% raw_fit$Gamma
      loglik_full  <- calculate.log.partial.lik(dlong_full,  B_hat, pred_names)
      loglik_train <- calculate.log.partial.lik(dlong_train, B_hat, pred_names)
      CVE_mat[fold, lam] <- -2 * (loglik_full - loglik_train) / N_test
    }
    cat(sprintf("  pen VVH fold %d/%d done\n", fold, n_folds))
  }
  
  CVE_mean    <- colMeans(CVE_mat, na.rm = TRUE)
  CVE_se      <- apply(CVE_mat, 2, sd, na.rm = TRUE) / sqrt(n_folds)
  best_lambda <- lambda_vector[which.min(CVE_mean)]
  
  cat(sprintf("[pen VVH] best lambda: %g  (CVE = %.4f)\n",
              best_lambda, min(CVE_mean, na.rm = TRUE)))
  
  return(list(CVE_mean = CVE_mean, CVE_se = CVE_se, CVE_mat = CVE_mat,
              lambda = lambda_vector, best_lambda = best_lambda))
}
## Function to fit data n_simulations times and access overall performance 

simulations_fit_and_performance <- function(datas, p, k, n_simulations, r,
                                            lambda_alpha, group_labels, method,
                                            pen_gamma_start = NULL, ...) {
  if (method == "pen" && is.null(datas$dlong_list))
    stop("method='pen' needs dlong_list, use store_dlong=TRUE.")
  
  bias1 <- mse1 <- cor_coeff_vec <- c()
  bias_list <- mse_list <- bias_for_all_coeffs <- mse_for_all_coeffs <-
    vector("list", n_simulations)
  simulation_times <- numeric(n_simulations)
  
  Beta_matrix     <- datas$Beta_matrix
  y_list_all      <- datas$y_list_all
  status_list_all <- datas$status_list_all
  X_list          <- datas$X_list
  dlong_list      <- datas$dlong_list
  
  for (i in seq_len(n_simulations)) {
    t0    <- Sys.time()
    dlong <- if (method == "pen") dlong_list[[i]] else NULL
    
    # pen 方法：根据规模自动选择初始化策略
    if (method == "pen") {
      if (!is.null(pen_gamma_start)) {
        gamma_init <- pen_gamma_start          # 用户手动指定
      } else if (p > 50 || n > 500) {
        gamma_init <- make_pen_gamma_init(     # 大规模用 SVD
          X_list[[i]], y_list_all[[i]],
          status_list_all[[i]], R = r, k = k)
      } else {
        gamma_init <- matrix(rnorm(r * k), r, k)  # 小规模随机即可
      }
    } else {
      gamma_init <- NULL
    }
    
    fit <- tryCatch(
      fit_one_model(
        X = X_list[[i]], y_list = y_list_all[[i]],
        status_list = status_list_all[[i]], dlong = dlong,
        method = method, r = r, lambda = lambda_alpha,
        group_labels = group_labels,
        pen_gamma_start = gamma_init, ...),
      error = function(e) {
        cat(sprintf("error(sim %d): %s\n", i, e$message)); NULL })
    if (is.null(fit)) next
    
    res  <- performance_model(fit$B_hat, Beta_matrix)
    bias1          <- c(bias1,          res$bias_out)
    mse1           <- c(mse1,           res$mse_out)
    cor_coeff_vec  <- c(cor_coeff_vec,  res$cor_coeff)
    bias_list[[i]] <- res$b_pred
    mse_list[[i]]  <- res$mse_pred
    bias_for_all_coeffs[[i]] <- res$bias_for_all_coeffs
    mse_for_all_coeffs[[i]]  <- res$mse_for_all_coeffs
    simulation_times[i] <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  }
  
  cat(sprintf("[%s] %d sims done, total %.1f sec\n",
              method, n_simulations, sum(simulation_times)))
  
  bias_df <- do.call(rbind, bias_list)
  mse_df  <- do.call(rbind, mse_list)
  colnames(bias_df) <- colnames(mse_df) <- paste0("Outcome_", 1:k)
  
  return(list(bias_df = bias_df, mse_df = mse_df,
              bias1 = bias1, mse1 = mse1,
              cor_coeff_vec = cor_coeff_vec,
              simulation_times = simulation_times,
              bias_for_all_coeffs_and_simulations = bias_for_all_coeffs,
              mse_for_all_coeffs_and_simulation    = mse_for_all_coeffs))
}

## Function to fit the models for multiple datasets and find the best lambda (cv)

simulations_fit_find_best_lambda <- function(datasets_list,
                                             lambda_vector,
                                             p, k, r,
                                             group_labels,
                                             method,
                                             n_folds_0 = 5,
                                             ...) {
  if (method == "pen" && is.null(datasets_list$dlong_list))
    stop("method='pen' needs dlong_list, use store_dlong=TRUE.")
  
  n_datasets      <- length(datasets_list$y_list_all)
  optimal_lambdas <- numeric(0)
  results_df      <- data.frame()
  
  for (j in seq_len(n_datasets)) {
    cat(sprintf("\n=== Dataset %d/%d ===\n", j, n_datasets))
    
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
    cat(sprintf("  -> dataset %d | best_lambda: %g\n", j, o_res$best_lambda))
  }
  
  best_lambda <- lambda_vector[which.min(abs(lambda_vector - median(optimal_lambdas)))]
  cat(sprintf("\n[%s VVH] global best lambda: %g\n", method, best_lambda))
  
  return(list(results_df      = results_df,
              best_lambda     = best_lambda,
              optimal_lambdas = optimal_lambdas))
}
## Function for Calculating Monte Carlo Error for MSE and Bias

MC_errors_calc = function(mse_for_all_sims,bias_for_all_sims){
  
  MC_error_mse = sd(mse_for_all_sims)/sqrt(length(mse_for_all_sims))
  MC_error_bias = sd(bias_for_all_sims)/sqrt(length(bias_for_all_sims))
  
  return(list(
    MC_error_mse = MC_error_mse,
    MC_error_bias = MC_error_bias
  ))
}


## Function to generate histograms

generate_histograms <- function(performance_results_list, k,r,barslength=0.005) {
  
  
bias_df = performance_results_list$bias_df
mse_df = performance_results_list$mse_df
bias1 = performance_results_list$bias1
mse1 = performance_results_list$mse1
cor_coeff_vec = performance_results_list$cor_coeff_vec


  # Helper function to create individual histograms
  create_histogram <- function(data, title, xlab, fill_color,binw) {
    ggplot(data.frame(value = data), aes(x = value)) +
      geom_histogram(aes(y = ..density..), binwidth = binw, fill = fill_color, color = "black", alpha = 0.7) +
      geom_density(color = "blue", linetype = "dashed", linewidth = 1) +
      geom_vline(aes(xintercept = mean(data)), color = "red", linetype = "dotted", linewidth = 1) +
      labs(title = title, x = xlab, y = "Density") +
      theme_minimal()
  }
  
  # Plot 1: Histogram of Mean Bias
  p1 <- create_histogram(bias1, "", "Mean Absolute Bias", "pink",barslength)
  
  # Plot 2: Histogram of Mean MSE
  p2 <- create_histogram(mse1, "", "Mean MSE", "orange",barslength/10)
  
  # Plot 3: Histograms of Mean Bias for each outcome
  bias_plots <- lapply(1:k, function(i) {
    create_histogram(bias_df[,i], paste("Outcome", i), "Mean Bias", "pink",barslength)
  })
  
  # Plot 4: Histograms of Mean MSE for each outcome
  mse_plots <- lapply(1:k, function(i) {
    create_histogram(mse_df[,i], paste("Outcome", i), "Mean MSE", "orange",barslength/10)
  })
  
  # Plot 5: Histogram of Correlation coefficient between B est and B true
  corr_plot <- create_histogram(cor_coeff_vec, "", "Corelation Coefficient between B true and B", "green",barslength/10)
  
  # Arrange plots in a grid 
  top <- textGrob(
  paste('Fitted rank:',r),
  gp = gpar(fontsize = 20)  # Adjust the fontsize as needed
)
  
  
  pL1 = grid.arrange(corr_plot, ncol = 1, top = top)
                    #"Histogrm of Corelation Coefficient between B true and B estimate across 200 simulations")
  pL2 =grid.arrange(p1, p2, ncol = 2, top = top)
                     #'Histograms of Mean Bias(left) and Mean MSE(right) across 200 Simulations')
  pL3 =grid.arrange(grobs = bias_plots, ncol = k/2, top = top)
                     #"Histograms of Mean Bias for each outcome across 200 Simulations")
  pL4 =grid.arrange(grobs = mse_plots, ncol = k/2, top = top)
  #"Histograms of Mean MSE for each outcome across 200 Simulations")
  
  return(list(
    meanbias1 = mean(bias1),
    meanmse1 = mean(mse1),
    meancorr =mean(cor_coeff_vec),
    pL1=pL1,
    pL2=pL2,
    pL3=pL3,
    pL4=pL4
  ))
}


## Function to plot simulated data

generate_data_plots <- function(dlong_list,n,r,p,k,T_sim=NULL,T_obs=NULL) {

combined_data=dlong_list
T_sim_combined=T_sim
T_obs_combined=T_obs
# Combine all data sets into a single data frame
if (!is.data.frame(dlong_list)){
combined_data <- do.call(rbind, dlong_list)
print(summary(combined_data))
}
  
if (!is.null(T_sim)){
T_sim_combined=do.call(rbind, T_sim)
T_obs_combined=do.call(rbind, T_obs)
}

predictor_names <- grep("^x", names(combined_data), value = TRUE)
predictors_long <- reshape2::melt(combined_data, id.vars = 'id', measure.vars = predictor_names)

# Plot the distribution of predictors on the same plot
p1 = ggplot(predictors_long, aes(x = value, color = variable, fill = variable)) +
  geom_density(alpha = 0.05) +
  labs(title = 'Distribution of \n Predictors', x = 'Value', y = 'Density') +
  theme(legend.title = element_blank())+
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 15),  
    axis.text.y = element_text(size = 15),                         
    axis.title.x = element_text(size = 15),                       
    axis.title.y = element_text(size = 15),                       
    legend.text = element_text(size = 15),                         
    legend.title = element_text(size = 15),                       
    plot.title = element_text(size = 20),           
    strip.text = element_text(size = 15)                           
  )
theme_minimal() 

# Plot the distribution of Tstop
p2 = ggplot(combined_data, aes(x = Tstop)) +
  geom_density(fill = 'blue', alpha = 0.05) +
  labs(title = 'Distribution of Tstop', x = 'Tstop', y = 'Density') +
  theme(legend.title = element_blank())+
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 15),  
    axis.text.y = element_text(size = 15),                         
    axis.title.x = element_text(size = 15),                       
    axis.title.y = element_text(size = 15),                       
    legend.text = element_text(size = 15),                         
    legend.title = element_text(size = 15),                       
    plot.title = element_text(size = 20),           
    strip.text = element_text(size = 15)                           
  )
theme_minimal() 

# T_obs_long <- as.data.frame(T_obs_combined) %>%
#   pivot_longer(
#     cols = everything(), 
#     names_to = "variable",  
#     values_to = "value"     
#   )

# p3 = ggplot(T_obs_long, aes(x = value, color = variable)) +
#   geom_density() +
#   labs(title = 'T_obs Plot for Each Outcome',
#        x = 'T_obs',
#        y = 'Density') +
#   theme_minimal()
# 
# T_sim_long <- as.data.frame(T_sim_combined) %>%
#   pivot_longer(
#     cols = everything(),  
#     names_to = "variable",  
#     values_to = "value"     
#   )
# 
# p4 = ggplot(T_sim_long, aes(x = value, color = variable)) +
#   geom_density() +
#   labs(title = 'T_sim Plot for Each Outcome',
#        x = 'T_simulated',
#        y = 'Density') +
#   theme_minimal()

top = paste('Datasets with ',p,' covariates and ',k,' outcomes for ',n,' individuals')
top <- textGrob(
  paste('Datasets with ',p,' covariates, \n',k,' outcomes, ',n,' individuals'),
  gp = gpar(fontsize = 25)  
)
grid_plots = grid.arrange(p1, p2, ncol = 2, top = top)
 
status_counts <- table(combined_data$status)
status_percentages <- prop.table(status_counts) * 100
# Print the percentages
print(paste('Percentage of Censored data:', status_percentages[1],' Percentage of Uncensored data', status_percentages[2]))
return(grid_plots)
}


## Calculate the Linear Predictor error

linear_predictor_err <- function(datas1, performance_results11) {
    
    X_list = datas1$X_list
    
    bias_list = performance_results11$bias_for_all_coeffs_and_simulations
    Beta_true = datas1$Beta_matrix
    
    Beta_list_estimated = list()
    eta_list_estimated = list()
    eta_list_real = list()
    
    lp_errors=list()
    
    ave_lp_errors = list()
    
    cor_lp=c()
    
    for(i in 1:length(datas1$y_list_all)){
      Beta_list_estimated[[i]] = bias_list[[i]] + Beta_true
      eta_list_estimated[[i]] = X_list[[i]] %*% Beta_list_estimated[[i]]
      
      eta_list_real[[i]]=X_list[[i]]  %*% Beta_true
      
      cor_lp=c(cor_lp,cor(as.vector(eta_list_estimated[[i]]),as.vector(eta_list_real[[i]])))
      
      K=dim(eta_list_estimated[[i]])[2]
      
      lp_errors[[i]]=rowSums((eta_list_estimated[[i]] - eta_list_real[[i]])^2)/K
      
      ave_lp_errors[[i]] = mean(lp_errors[[i]])
      
    }
    
    return(list(
    ave_lp_errors = ave_lp_errors,
    lp_errors = lp_errors,
    cor_lp=cor_lp))
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
    axis.text.y = element_text(size = 15),                         
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





## Function to plot lambdas

lambda_plot <- function(o_result, title = "V&VH Lambda Selection") {
  df <- data.frame(
    lambda   = o_result$lambda,
    CVE_mean = o_result$CVE_mean,
    CVE_se   = o_result$CVE_se
  )
  ggplot(df, aes(x = log(lambda), y = CVE_mean)) +
    geom_line() + geom_point() +
    geom_errorbar(aes(ymin = CVE_mean - CVE_se,
                      ymax = CVE_mean + CVE_se), width = 0.05) +
    geom_vline(xintercept = log(o_result$best_lambda),
               linetype = "dashed", color = "red") +
    labs(title   = title,
         x       = expression(log(lambda)),
         y       = "Mean CVE (V&VH)",
         caption = "Red: min CVE") +
    theme_minimal()
}