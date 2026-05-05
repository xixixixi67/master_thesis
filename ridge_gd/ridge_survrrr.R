library(glmnet)
library(survival)

#' Fit Penalized Reduced Rank Multi-Outcome Cox Model
#'
#' Both alpha (p x R) and Gamma (K x R) are penalized with Ridge.
#' Gamma has NO orthogonality constraint.
#' A single lambda sequence is used for both: lambda_gamma = lambda_alpha.
#' @return Object of class "survRRR_Ridge"
solve_RR_Ridge <- function(X, y_list, entry_list=NULL, status_list, R,
                           lambda_alpha,
                           rho = 0.01,
                           alpha0 = NULL,
                           Gamma0 = NULL,
                           noise_sd = 0.1,
                           step_size = 0.01,
                           max_iter = 500,
                           tol = 1e-4,
                           verbose = FALSE)
{
  if (!is.matrix(X)) {
    stop("X must be a matrix")
  }
  if (!is.list(y_list) || !is.list(status_list)) {
    stop("y_list and status_list must be lists")
  }
  if (length(y_list) != length(status_list)) {
    stop("y_list and status_list must have the same length")
  }
  K <- length(y_list)
  N <- nrow(X)
  p <- ncol(X)
  
  # lambda_gamma is always equal to lambda_alpha (single shared sequence)
  lambda_gamma <- lambda_alpha
  
  if (is.null(alpha0) || is.null(Gamma0)) {
    B0 <- matrix(0, p, K)
    for (k in 1:K) {
      fit_ridge <- glmnet(X,
                          Surv(y_list[[k]], status_list[[k]]),
                          family = "cox",
                          alpha = 0,
                          lambda = rho)
      B0[, k] <- as.vector(coef(fit_ridge))
    }
    svd_B   <- svd(B0)
    U_R     <- svd_B$u[, 1:R, drop = FALSE]
    Sigma_R <- diag(svd_B$d[1:R], nrow = R, ncol = R)
    V_R     <- svd_B$v[, 1:R, drop = FALSE]
    
    alpha_base <- U_R %*% Sigma_R
    Gamma_base <- V_R
    
    scale_alpha <- mean(abs(alpha_base))
    scale_Gamma <- mean(abs(Gamma_base))
    
    alpha0_perturbed <- alpha_base + matrix(rnorm(p*R,0,noise_sd*scale_alpha), p, R)
    Gamma0_perturbed <- Gamma_base + matrix(rnorm(K*R,0,noise_sd*scale_Gamma), K, R)
    
    if (is.null(alpha0)) alpha0 <- alpha0_perturbed   # p x R
    if (is.null(Gamma0)) Gamma0 <- Gamma0_perturbed   # K x R
  }
  
  order_list <- vector("list", K)
  order_list_entry <- vector("list", K)
  status_mat <- matrix(0.0, nrow = N, ncol = K)
  rankmin_exit    <- matrix(0L,  nrow = N, ncol = K)
  rankmax_exit    <- matrix(0L,  nrow = N, ncol = K)
  rankmin_entry    <- matrix(0L,  nrow = N, ncol = K)
  exit_sorted_mat <- matrix(0.0, nrow = N, ncol = K)
  entry_sorted_mat <- matrix(0.0, nrow = N, ncol = K) 
  
  if (is.null(entry_list)) {
    if (verbose) message("entry_list not provided. Assuming all entries at time 0.")
    entry_list <- lapply(y_list, function(y) rep(0, length(y)))
  }
  
  for (k in 1:K) {
    y <- y_list[[k]]
    s <- status_list[[k]]
    e <- entry_list[[k]]
    o_exit <- order(y)
    y_sorted <- y[o_exit]
    s_sorted <- s[o_exit]
    e_sorted_by_exit <- e[o_exit]
    
    order_list[[k]] <- order(o_exit) - 1L
    
    n_events <- sum(s_sorted)
    if (n_events == 0) stop(paste("Outcome", k, "has no events"))
    status_mat[, k] <- s_sorted / n_events
    
    rankmin_exit[, k] <- rank(y_sorted, ties.method = "min") - 1L
    rankmax_exit[, k] <- rank(y_sorted, ties.method = "max") - 1L
    
    o_entry <- order(e_sorted_by_exit)
    order_list_entry[[k]] <- order(o_entry) - 1L
    e_double_sorted <- e_sorted_by_exit[o_entry]
    rankmin_entry[, k] <- rank(e_double_sorted, ties.method = "min") - 1L
  }
  
  if (ncol(alpha0) != R) stop("alpha0 must have R columns")
  if (nrow(alpha0) != p) stop("alpha0 must have p rows")
  if (ncol(Gamma0) != R) stop("Gamma0 must have R columns")
  if (nrow(Gamma0) != K) stop("Gamma0 must have K rows")
  
  result <- fit_RR_Ridge(
    X                = X,
    status           = status_mat,
    rankmin          = rankmin_exit,
    rankmax          = rankmax_exit,
    rankmin_entry    = rankmin_entry,
    entry_sorted_mat = entry_sorted_mat,  
    exit_sorted_mat  = exit_sorted_mat,
    order_list       = order_list,
    order_list_entry = order_list_entry,
    alpha0           = alpha0,
    Gamma0           = Gamma0,
    lambda_alpha_all = lambda_alpha,
    lambda_gamma_all = lambda_gamma,   # same as lambda_alpha
    step_size        = step_size,
    niter            = max_iter,
    tol              = tol,
    linesearch_beta  = 0.5,
    verbose          = verbose
  )
  
  result$call         <- match.call()
  result$lambda_alpha <- lambda_alpha
  result$lambda_gamma <- lambda_gamma
  result$dimensions   <- list(N = N, p = p, K = K, R = R)
  
  class(result) <- "survRRR_ridge"
  return(result)
}



#' Compute Residuals for Reduced Rank Cox Model
get_residual_RR_Ridge <- function(X, y_list, entry_list=NULL,status_list, alpha, Gamma)
{
  K <- length(y_list)
  N <- nrow(X)
  p <- ncol(X)
  R <- ncol(alpha)
  
  if (nrow(alpha) != p) stop("alpha must be p x R")
  if (nrow(Gamma) != K) stop("Gamma must be K x R")
  if (ncol(Gamma) != R) stop("Gamma must have R columns")
  
  order_list <- vector("list", K)
  order_list_entry <- vector("list", K)
  status_mat <- matrix(0.0, nrow = N, ncol = K)
  rankmin_exit    <- matrix(0L,  nrow = N, ncol = K)
  rankmax_exit    <- matrix(0L,  nrow = N, ncol = K)
  rankmin_entry    <- matrix(0L,  nrow = N, ncol = K)
  exit_sorted_mat <- matrix(0.0, nrow = N, ncol = K)
  entry_sorted_mat <- matrix(0.0, nrow = N, ncol = K) 
  
  if (is.null(entry_list)) {
    entry_list <- lapply(y_list, function(y) rep(0, length(y)))
  }
  
  for (k in 1:K) {
    y <- y_list[[k]]
    s <- status_list[[k]]
    e <- entry_list[[k]]
    o_exit <- order(y)
    y_sorted <- y[o_exit]
    s_sorted <- s[o_exit]
    e_sorted_by_exit <- e[o_exit]
    order_list[[k]] <- order(o_exit) - 1L
    n_events <- sum(s_sorted)
    if (n_events == 0) stop(paste("Outcome", k, "has no events"))
    status_mat[, k] <- s_sorted / n_events
    rankmin_exit[, k] <- rank(y_sorted, ties.method = "min") - 1L
    rankmax_exit[, k] <- rank(y_sorted, ties.method = "max") - 1L
    
    o_entry <- order(e_sorted_by_exit)
    order_list_entry[[k]] <- order(o_entry) - 1L
    e_double_sorted <- e_sorted_by_exit[o_entry]
    rankmin_entry[, k] <- rank(e_double_sorted, ties.method = "min") - 1L
  }
  
  compute_residual_RR_Ridge(X, status_mat, rankmin_exit, rankmax_exit, rankmin_entry, entry_sorted_mat, 
                            exit_sorted_mat, order_list, order_list_entry, alpha, Gamma)
}

#' Print method for survRRR_Ridge objects
print.survRRR_ridge <- function(x, ...) {
  cat("Reduced Rank Multi-Outcome Cox Model (Ridge on alpha and Gamma)\n")
  cat("=================================================================\n\n")
  cat("Dimensions:\n")
  cat("  Samples (N):", x$dimensions$N, "\n")
  cat("  Predictors (p):", x$dimensions$p, "\n")
  cat("  Outcomes (K):", x$dimensions$K, "\n")
  cat("  Rank (R):", x$dimensions$R, "\n\n")
  cat("Lambdas fitted:", length(x$result), "\n")
  cat("Range of lambda (alpha = Gamma):",
      round(min(x$lambda_alpha), 4), "to",
      round(max(x$lambda_alpha), 4), "\n")
  cat("Convergence:\n")
  cat("  Mean iterations:", round(mean(x$num_iterations), 1), "\n")
  cat("  Range:", min(x$num_iterations), "to", max(x$num_iterations), "\n")
  invisible(x)
}


#' Summary method for survRRR_ridge objects
summary.survRRR_ridge <- function(object, ...) {
  print.survRRR_ridge(object)
  cat("\nObjective Function Values:\n")
  print(summary(object$objective_values))
  invisible(object)
}
