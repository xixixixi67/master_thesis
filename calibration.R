library(CalibrationCurves)
library(survival)
library(parallel)


compute_calibration_metrics <- function(predicted_probs, actual_status, 
                                        prob_intervals = list(c(0.95, 1), 
                                                              c(0.9, 0.95), 
                                                              c(0.85, 0.9), 
                                                              c(0.8, 0.85))) {
  
  results <- list()
  
  for (interval in prob_intervals) {
    lower <- interval[1]
    upper <- interval[2]
    
    # 筛选处于当前概率区间内的样本
    idx <- which(predicted_probs >= lower & predicted_probs <= upper)
    
    if (length(idx) == 0) {
      results[[paste0(lower, "-", upper)]] <- c(intercept = NA, slope = NA)
      next
    }
    
    # 准备数据并拟合逻辑回归模型
    df <- data.frame(
      status = actual_status[idx],
      logit_pred = log(predicted_probs[idx] / (1 - predicted_probs[idx]))
    )
    
    fit <- tryCatch({
      glm(status ~ 1 + offset(logit_pred), data = df, family = binomial)
    }, error = function(e) NULL)
    
    if (is.null(fit)) {
      results[[paste0(lower, "-", upper)]] <- c(intercept = NA, slope = NA)
      next
    }
    

    intercept <- coef(fit)[1]               # 校准截距 (α_c)
    slope <- 1                              # 校准斜率 (ζ) 在offset模型中固定为1
    # 注意：上述模型斜率固定为1，仅估计截距。若需同时估计截距和斜率，应使用：
    # fit2 <- glm(status ~ logit_pred, data = df, family = binomial)
    # slope <- coef(fit2)[2]                # 校准斜率 (ζ)
    # intercept <- coef(fit2)[1]            # 校准截距 (α)
    
    results[[paste0(lower, "-", upper)]] <- c(intercept = intercept, slope = slope)
  }
  
  return(do.call(rbind, results))
}

# 主函数：处理100个数据集
process_100_datasets <- function(data_list, time_point, prob_intervals) {
  n_datasets <- length(data_list)
  K <- length(prob_intervals)
  
  # 存储所有数据集的截距和斜率
  intercepts <- matrix(NA, nrow = n_datasets, ncol = K)
  slopes <- matrix(NA, nrow = n_datasets, ncol = K)
  
  for (i in 1:n_datasets) {
    data <- data_list[[i]]
    
    # 假设 data 包含以下元素：
    # - X: 协变量矩阵
    # - y: 生存时间
    # - status: 事件指示器
    # - surv_prob: 模型预测的生存概率矩阵（N x T）
    
    # 提取指定时间点的预测生存概率和实际状态
    pred_probs <- data$surv_prob[, which(time_points == time_point)]
    actual_status <- ifelse(data$y <= time_point & data$status == 1, 1, 0)
    
    # 调用校准函数
    metrics <- compute_calibration_metrics(pred_probs, actual_status, prob_intervals)
    
    intercepts[i, ] <- metrics[, "intercept"]
    slopes[i, ] <- metrics[, "slope"]
  }
  
  # 计算平均值和标准差
  avg_intercepts <- colMeans(intercepts, na.rm = TRUE)
  sd_intercepts <- apply(intercepts, 2, sd, na.rm = TRUE)
  avg_slopes <- colMeans(slopes, na.rm = TRUE)
  sd_slopes <- apply(slopes, 2, sd, na.rm = TRUE)
  
  # 整理结果
  results <- data.frame(
    Interval = sapply(prob_intervals, function(x) paste0("(", x[1], ",", x[2], ")")),
    Intercept_Mean = avg_intercepts,
    Intercept_SD = sd_intercepts,
    Slope_Mean = avg_slopes,
    Slope_SD = sd_slopes
  )
  
  return(results)
}