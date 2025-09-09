set.seed(123)

# Set working directory
current_working_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(current_working_dir)
getwd()

# Load libraries
library(RGCCA)
library(qgraph)
library(MASS) 

###########################################################################################
# Evaluation Metrics

evaluate_variable_selection <- function(W_true, W_estimated) {
  
  # Variable Selection Metrics
  W_true_bin <- ifelse(W_true != 0, 1, 0)
  W_est_bin <- ifelse(W_estimated != 0, 1, 0)
  
  TP <- sum(W_true_bin == 1 & W_est_bin == 1)
  FP <- sum(W_true_bin == 0 & W_est_bin == 1)
  FN <- sum(W_true_bin == 1 & W_est_bin == 0)
  TN <- sum(W_true_bin == 0 & W_est_bin == 0)
  
  precision <- TP / (TP + FP + 1e-8)
  recall <- TP / (TP + FN + 1e-8)
  f1_score <- 2 * (precision * recall) / (precision + recall + 1e-8)
  accuracy <- (TP + TN) / (TP + FP + FN + TN)
  
  return(list(
    precision = precision,
    recall = recall,
    f1 = f1_score,
    recovery = accuracy
  ))
}

reconstruction_metrics <- function(X, W, P_T) {
  
  # Reconstruction Metrics
  X_hat <- X %*% W %*% t(P_T)
  error_matrix <- X - X_hat
  mse <- mean(error_matrix^2)
  var_explained <- 1 - (sum(error_matrix^2) / sum((X - mean(X))^2))
  return(list(mse = mse, R2 = var_explained))
}

score_metrics <- function(est, true) {
  
  # Compute MAE and RMSE
  mae <- mean(abs(est - true))
  rmse <- sqrt(mean((est - true)^2))
  return(list(mae = mae, rmse = rmse))
}

compute_bias_variance_mse <- function(W_true, W_est) {
  
  # compute bias, variance and mse
  W_true_vec <- as.vector(W_true)
  W_est_vec <- as.vector(W_est)
  bias <- mean(W_est_vec - W_true_vec)
  variance <- var(W_est_vec - W_true_vec)
  mse <- mean((W_est_vec - W_true_vec)^2)
  
  return(list(bias = bias, variance = variance, mse = mse))
}


sparsity_level <- function(W) {
  
  # Compute sparsity in parameter
  total_elements <- length(W)
  zero_elements <- sum(W == 0)
  return(zero_elements / total_elements)
}

compute_vaf <- function(X, W, P) {
  # Variance Accounted For calculation
  X_hat <- X %*% W %*% t(P)
  sum_sq_error <- sum((X - X_hat)^2)
  total_variance <- sum(X^2)
  vaf <- 1 - (sum_sq_error / total_variance)
  return(vaf)
}


