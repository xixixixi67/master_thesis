library(survC1)

# ?????? scale ?????? sd(eta) ???????????????????????????
# test set with 5000 samples
set.seed(7)
sim_data1_test <- simulations_data(5000, 2000, 2, 4,
                              n_simulations = 1,
                              n_groups      = n_groups,
                              zero_blocks   = zero_blocks,
                              store_dlong   = FALSE,
                              case          = 2,
                              group_corr    = 0.7,
                              beta_obj      = search_data_1,
                              scale         = 0.15)

# model1 <-> rrr_grplasso
X_list_train <- search_data_1$X_list
y_list_all_train <- search_data_1$y_list_all
status_list_all_train <- search_data_1$status_list_all
K <- 4
model1 <- fit_one_model(
  X = X_list_train[[1]],
  y_list = y_list_all_train[[1]],
  status_list = status_list_all_train[[1]],
  method = "mrcox",
  group_labels = group_labels_1,
  lambda = 0.215714 ,
  r = 1,

  dlong = NULL,          
  verbose = TRUE
)

X_test <- sim_data1_test$X_list[[1]]
Y_list_test <- sim_data1_test$y_list_all[[1]]
status_list_test <- sim_data1_test$status_list_all[[1]]

B_hat1 <- model1$B_hat
eta1 <- X_test %*% B_hat1

cindices1 <- numeric(K)
for (k in 1:K) {
  event_times <- Y_list_test[[k]][status_list_test[[k]] == 1]
  if (length(event_times) > 0) {
      tau <- quantile(event_times, 0.8, na.rm = TRUE)
  } else {
    tau <- median(Y_list_test[[k]], na.rm = TRUE)
  }
  
  mydata1 <- data.frame(
    time = Y_list_test[[k]],
    event = status_list_test[[k]],
    risk = eta1[, k]
  )
  
  cval1 <- Est.Cval(mydata1, tau = tau, nofit = TRUE)$Dhat
  
  cindices1[k] <- cval1
}
cindices1
