library(Rcpp)
library(RcppEigen)

sourceCpp("src/test.cpp")
source("test.R")


# simulate data
set.seed(20020607)
n <- 200
p <- 30
K <- 3
R_true <- 2

X <- matrix(rnorm(n * p), n, p)

alpha_true <- matrix(runif(p * R_true, 0, 1), p, R_true)
Gamma_true <- matrix(runif(K * R_true, 0, 1), K, R_true)
B_true <- alpha_true %*% t(Gamma_true)

eta <- X %*% B_true

y_list <- list()
status_list <- list()

for (k in 1:K) {
  lambda_0 <- 0.01
  U <- runif(n)
  lambda_i <- lambda_0 * exp(eta[, k])
  T_event <- -log(U) / lambda_i
  
  C <- runif(n, 0, quantile(T_event, 0.7))
  
  Y <- pmin(T_event, C)
  delta <- as.numeric(T_event <= C)
  
  y_list[[k]] <- Y
  status_list[[k]] <- delta
}

cat("Event rates:\n")
for (k in 1:K) {
  cat("  Event", k, ":", round(mean(status_list[[k]]), 3), "\n")
}


# single lambda pair
fit1 <- solve_reduced_rank(
  X = X,
  y_list = y_list,
  status_list = status_list,
  R = 2,
  lambda_alpha = 0.1,
  lambda_gamma = 0.1,
  step_size = 0.01,
  max_iter = 100,
  verbose = TRUE
)

fit1$result[[1]]
cat("  Objective value:", fit1$objective_values[1], "\n")
cat("  Iterations:", fit1$num_iterations[1], "\n\n")


