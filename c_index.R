library(survC1)
library(parallel)

C_index_datasets <- function(
    sim_data,                     
    model_configs,                 # method, r, lambda, group_labels
    tau_years = NULL,              
    n_cores = 10,                 
    verbose = TRUE
) {
  
  # extract data
  X_list <- sim_data$X_list
  y_list_all <- sim_data$y_list_all
  status_list_all <- sim_data$status_list_all
  
  n_datasets <- length(X_list)          
  K <- length(y_list_all[[1]])          
  
  if (verbose) {
    cat("\n============================================================\n")
    cat(sprintf("Datasets: %d | Outcomes: %d\n", n_datasets, K))
    cat("============================================================\n")
  }
  
  fit_and_cindex <- function(dataset_idx, method, r, lambda, group_labels) {

    X <- X_list[[dataset_idx]]
    y_list <- y_list_all[[dataset_idx]]
    status_list <- status_list_all[[dataset_idx]]
    
    fit <- tryCatch({
      fit_one_model(
        X = X,
        y_list = y_list,
        status_list = status_list,
        method = method,
        r = r,
        lambda = lambda,
        group_labels = group_labels,
        dlong = NULL,          
        verbose = FALSE
      )
    }, error = function(e) {
      if (verbose) cat(sprintf("  Dataset %d, method %s: fitting error: %s\n", 
                               dataset_idx, method, e$message))
      return(NULL)
    })
    
    if (is.null(fit) || is.null(fit$B_hat)) {
      return(rep(NA, K))
    }
    
    B_hat <- fit$B_hat
    eta <- X %*% B_hat 
    
    cindices <- numeric(K)
    for (k in 1:K) {
      if (is.null(tau_years)) {
        event_times <- y_list[[k]][status_list[[k]] == 1]
        if (length(event_times) > 0) {
          tau <- quantile(event_times, 0.8, na.rm = TRUE)
        } else {
          tau <- median(y_list[[k]], na.rm = TRUE)
        }
      } else {
        tau <- tau_years[k] * 365.25
      }
      
      mydata <- data.frame(
        time = y_list[[k]],
        event = status_list[[k]],
        risk = eta[, k]
      )
      
      cval <- tryCatch({
        Est.Cval(mydata, tau = tau, nofit = TRUE)$Dhat
      }, error = function(e) NA)
      
      cindices[k] <- cval
    }
    return(cindices)
  }
  

  all_results <- list()
  
  for (cfg in model_configs) {
    method <- cfg$method
    r <- cfg$r
    lambda <- cfg$lambda
    group_labels <- if (!is.null(cfg$group_labels)) cfg$group_labels else NULL
    
    if (verbose) {
      cat(sprintf("\n===== Method: %s, r=%s, lambda=%.6f =====\n", 
                  method, ifelse(is.null(r), "NULL", r), lambda))
    }
    
    cindex_matrix <- mclapply(1:n_datasets, function(i) {
      if (verbose && i %% 20 == 0) cat(sprintf("  Dataset %d/%d\n", i, n_datasets))
      fit_and_cindex(i, method, r, lambda, group_labels)
    }, mc.cores = n_cores)
    
    # n_datasets x K
    cindex_mat <- do.call(rbind, cindex_matrix)
    colnames(cindex_mat) <- paste0("Outcome", 1:K)
    rownames(cindex_mat) <- paste0("Dataset", 1:n_datasets)
    
    # average C-index per outcome
    avg_cindex_per_outcome <- colMeans(cindex_mat, na.rm = TRUE)
    overall_avg <- mean(avg_cindex_per_outcome)
    
    if (verbose) {
      cat("\nAverage C-index per outcome:\n")
      print(round(avg_cindex_per_outcome, 4))
      cat(sprintf("Overall average: %.4f\n", overall_avg))
    }
    
    all_results[[method]] <- list(
      method = method,
      parameters = list(r = r, lambda = lambda),
      cindex_matrix = cindex_mat,                # 100 x K
      avg_cindex_per_outcome = avg_cindex_per_outcome,
      overall_avg_cindex = overall_avg
    )
  }

  
  return(all_results)
}

