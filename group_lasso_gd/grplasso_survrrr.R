library(glmnet)
library(survival)

solve_RR_GrpLasso <- function(X, y_list, status_list, R,
                               lambda_alpha,
                               group_labels,
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
  if (length(group_labels) != ncol(X)) {
    stop("group_labels must have length p (one entry per predictor)")
  }
  
  K <- length(y_list)
  N <- nrow(X)
  p <- ncol(X)
  
  # Ridge initialization: fit separate ridge Cox for each outcome
  if (is.null(alpha0) || is.null(Gamma0)) {
    B0 <- matrix(0, p, K)
    for (k in 1:K) {
      fit_ridge <- glmnet(X,
                          Surv(y_list[[k]], status_list[[k]]),
                          family = "cox",
                          alpha = 0,       # ridge
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
  
  # Preprocess survival data for C++ 
  order_list <- vector("list", K)
  status_mat <- matrix(0.0, nrow = N, ncol = K)
  rankmin    <- matrix(0L,  nrow = N, ncol = K)
  rankmax    <- matrix(0L,  nrow = N, ncol = K)
  
  for (k in 1:K) {
    y <- y_list[[k]]
    s <- status_list[[k]]
    o <- order(y)              
    y_sorted <- y[o]
    s_sorted <- s[o]
    
    # restore original order from sorted
    order_list[[k]] <- order(o) - 1L
    
    # normalize status by number of events
    n_events <- sum(s_sorted)
    if (n_events == 0) stop(paste("Outcome", k, "has no events"))
    status_mat[, k] <- s_sorted / n_events
    
    # tied events
    rankmin[, k] <- rank(y_sorted, ties.method = "min") - 1L
    rankmax[, k] <- rank(y_sorted, ties.method = "max") - 1L
  }
  
  if (ncol(alpha0) != R)  stop("alpha0 must have R columns")
  if (nrow(alpha0) != p)  stop("alpha0 must have p rows")
  if (ncol(Gamma0) != R)  stop("Gamma0 must have R columns")
  if (nrow(Gamma0) != K)  stop("Gamma0 must have K rows")
  
  # Call C++ optimizer
  result <- fit_RR_GrpLasso(
    X              = X,
    status         = status_mat,
    rankmin        = rankmin,
    rankmax        = rankmax,
    order_list     = order_list,
    alpha0         = alpha0,
    Gamma0         = Gamma0,
    lambda_alpha_all = lambda_alpha,
    group_labels   = as.integer(group_labels),
    step_size      = step_size,
    niter          = max_iter,
    tol            = tol,
    linesearch_beta = 0.5,
    verbose        = verbose
  )
  
  result$call        <- match.call()
  result$lambda_alpha <- lambda_alpha
  result$R           <- R
  result$dimensions  <- list(N = N, p = p, K = K, R = R)
  
  class(result) <- "survRRR" # name can be changed
  return(result)
}


# compute residuals
get_residual_RR_GrpLasso <- function(X, y_list, status_list, alpha, Gamma)
{
  K <- length(y_list)
  N <- nrow(X)
  p <- ncol(X)
  R <- ncol(alpha)
  
  if (nrow(alpha) != p)  stop("alpha must be p x R")
  if (nrow(Gamma) != K)  stop("Gamma must be K x R")
  if (ncol(Gamma) != R)  stop("Gamma must have R columns")
  
  order_list <- vector("list", K)
  status_mat <- matrix(0.0, nrow = N, ncol = K)
  rankmin    <- matrix(0L,  nrow = N, ncol = K)
  rankmax    <- matrix(0L,  nrow = N, ncol = K)
  
  for (k in 1:K) {
    y <- y_list[[k]]
    s <- status_list[[k]]
    o <- order(y)
    y_sorted <- y[o]
    s_sorted <- s[o]
    order_list[[k]] <- order(o) - 1L
    n_events <- sum(s_sorted)
    if (n_events == 0) stop(paste("Outcome", k, "has no events"))
    status_mat[, k] <- s_sorted / n_events
    rankmin[, k] <- rank(y_sorted, ties.method = "min") - 1L
    rankmax[, k] <- rank(y_sorted, ties.method = "max") - 1L
  }
  
  compute_residual_RR_GrpLasso(X, status_mat, rankmin, rankmax, order_list, alpha, Gamma)
}


# Print method for survRRR objects
print.survRRR_GrpLasso <- function(x, ...) {
  cat("Reduced Rank Multi-Outcome Cox Model\n")
  cat("=====================================\n\n")
  cat("Dimensions:\n")
  cat("  Samples (N):", x$dimensions$N, "\n")
  cat("  Predictors (p):", x$dimensions$p, "\n")
  cat("  Outcomes (K):", x$dimensions$K, "\n")
  cat("  Rank (R):", x$dimensions$R, "\n\n")
  cat("Lambdas fitted:", length(x$result), "\n")
  cat("Range of lambda_alpha:",
      round(min(x$lambda_alpha), 4), "to",
      round(max(x$lambda_alpha), 4), "\n")
  cat("Convergence:\n")
  cat("  Mean iterations:", round(mean(x$num_iterations), 1), "\n")
  cat("  Range:", min(x$num_iterations), "to", max(x$num_iterations), "\n")
  invisible(x)
}


# Summary method for survRRR objects
summary.survRRR_GrpLasso <- function(object, ...) {
  print(object)
  cat("\nObjective Function Values:\n")
  print(summary(object$objective_values))
  invisible(object)
}
