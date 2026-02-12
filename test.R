solve_reduced_rank <- function(X, y_list, status_list, R,
                               lambda_alpha, lambda_gamma,
                               pfac_alpha = NULL,
                               pfac_gamma = NULL,
                               alpha0 = NULL,
                               Gamma0 = NULL,
                               step_size = 0.01,
                               max_iter = 500,
                               tol = 1e-4,
                               verbose = TRUE)
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
  
  order_list <- list()
  status <- matrix(nrow = N, ncol = K)
  rankmin <- matrix(0L, nrow = N, ncol = K)
  rankmax <- matrix(0L, nrow = N, ncol = K)
  
  for (k in 1:K) {
    y <- y_list[[k]]
    o <- order(y)
    y <- y[o]

    order_list[[k]] <- order(o) - 1L
    
    status[, k] <- status_list[[k]][o] / sum(status_list[[k]])

    rankmin[, k] <- rank(y, ties.method = "min") - 1L
    rankmax[, k] <- rank(y, ties.method = "max") - 1L
  }
  
  # set default penalty weight
  if (is.null(pfac_alpha)) {
    pfac_alpha <- rep(1, p)
  }
  if (is.null(pfac_gamma)) {
    pfac_gamma <- rep(1, K)
  }
  
  if (is.null(alpha0)) {
    alpha0 <- matrix(rnorm(p * R, 0, 0.1), p, R)
  }
  if (is.null(Gamma0)) {
    Gamma0 <- matrix(rnorm(K * R, 0, 0.1), K, R)
  }
  
  if (ncol(alpha0) != R) {
    stop("alpha0 must have R columns")
  }
  if (nrow(alpha0) != p) {
    stop("alpha0 must have p rows")
  }
  if (ncol(Gamma0) != R) {
    stop("Gamma0 must have R columns")
  }
  if (nrow(Gamma0) != K) {
    stop("Gamma0 must have K rows")
  }
  
  # call C++ function
  result <- fit_reduced_rank(
    X = X,
    status = status,
    rankmin = rankmin,
    rankmax = rankmax,
    order_list = order_list,
    alpha0 = alpha0,
    Gamma0 = Gamma0,
    lambda_alpha_all = lambda_alpha,
    lambda_Gamma_all = lambda_gamma,
    pfac_alpha = pfac_alpha,
    pfac_Gamma = pfac_gamma,
    step_size = step_size,
    niter = max_iter,
    tol = tol,
    linesearch_beta = 0.5,
    verbose = verbose
  )
  
  result$call <- match.call()
  result$lambda_alpha <- lambda_alpha
  result$lambda_gamma <- lambda_gamma
  result$R <- R
  result$dimensions <- list(N = N, p = p, K = K, R = R)
  
  class(result) <- "survRRR"
  
  return(result)
}


# Compute Residuals for Reduced Rank Cox Model

get_residual_rr <- function(X, y_list, status_list, alpha, Gamma)
{
  K <- length(y_list)
  N <- nrow(X)
  p <- ncol(X)
  R <- ncol(alpha)
  
  if (nrow(alpha) != p) {
    stop("alpha must be p × R")
  }
  if (nrow(Gamma) != K) {
    stop("Gamma must be K × R")
  }
  if (ncol(Gamma) != R) {
    stop("Gamma must have R columns")
  }
  
  order_list <- list()
  status <- matrix(nrow = N, ncol = K)
  rankmin <- matrix(0L, nrow = N, ncol = K)
  rankmax <- matrix(0L, nrow = N, ncol = K)
  
  for (k in 1:K) {
    y <- y_list[[k]]
    o <- order(y)
    y <- y[o]
    order_list[[k]] <- order(o) - 1L
    status[, k] <- status_list[[k]][o] / sum(status_list[[k]])
    rankmin[, k] <- rank(y, ties.method = "min") - 1L
    rankmax[, k] <- rank(y, ties.method = "max") - 1L
  }
  
  # call C++ function
  compute_residual_rr(X, status, rankmin, rankmax, order_list, alpha, Gamma)
}



# Print Method for survRRR Objects
print.survRRR <- function(x, ...) {
  cat("Reduced Rank Multi-Outcome Cox Model\n")
  cat("=====================================\n\n")
  
  cat("Dimensions:\n")
  cat("  Samples (N):", x$dimensions$N, "\n")
  cat("  Predictors (p):", x$dimensions$p, "\n")
  cat("  Outcomes (K):", x$dimensions$K, "\n")
  cat("  Rank (R):", x$dimensions$R, "\n\n")
  
  cat("Lambda pairs fitted:", length(x$result), "\n")
  cat("Range of lambda_alpha:", round(min(x$lambda_alpha), 4), "to", round(max(x$lambda_alpha), 4), "\n")
  cat("Range of lambda_gamma:", round(min(x$lambda_gamma), 4), "to", round(max(x$lambda_gamma), 4), "\n\n")
  
  cat("Convergence:\n")
  cat("  Mean iterations:", round(mean(x$num_iterations), 1), "\n")
  cat("  Range:", min(x$num_iterations), "to", max(x$num_iterations), "\n")
  
  invisible(x)
}


# Summary Method for survRRR Objects

summary.survRRR <- function(object, ...) {
  print(object)
  
  cat("\n")
  cat("Objective Function Values:\n")
  print(summary(object$objective_values))
  
  invisible(object)
}