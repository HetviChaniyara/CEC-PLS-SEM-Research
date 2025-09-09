set.seed(123)

# Set working directory
current_working_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(current_working_dir)
getwd()

# Load libraries
library(RGCCA)
library(qgraph)
library(MASS) 
library(dplyr)

# Load functions from SGCCA functions file
source("~/BEP-Project-CEC-PLS-SEM/Scripts/SGCCA_Functions.R")

# Load simulation info
load("~/BEP-Project-CEC-PLS-SEM/Scripts/Data-R-W-Sparse/Info_simulation.RData")
Info_matrix <- Info_simulation$design_matrix_replication
Ndatasets <- Info_simulation$n_data_sets

results_list <- list()

for (i in 1:Ndatasets) {

  # Load data file
  filename <- paste0("~/BEP-Project-CEC-PLS-SEM/Scripts/DATA-R-W-Sparse/Wsparse", i, ".RData")
  load(filename) 
  
  # Extract true data
  X <- out$X
  W_true <- out$W
  P_true <- out$P
  T_true <- out$Z 
  
  R <- 2 # number of latent variables, can be changed to n_components in info_matrix
  I <- nrow(X)
  J <- ncol(X)
  conditions <- Info_matrix[i, ]
  
  # Define list of data
  blocks <- list(bl1 = X)
  
  n_vars <- ncol(blocks[[1]])
  var_names <- colnames(blocks[[1]])
  
  # Fit model based on the defined sparisty in info_matrix
  if (as.numeric(conditions[3]) == 0) {
    
    # Regular RGCCA (no sparsity)
    fit_sgcca <- rgcca(
      blocks = blocks,
      method = "rgcca",
      superblock = FALSE,
      ncomp = R,
      scheme = "factorial",
      comp_orth = FALSE,
      verbose = FALSE
    )
    
  } else {
    
    # Sparse SGCCA
    fit_sgcca <- rgcca(
      blocks = blocks,
      method = "sgcca",
      sparsity = as.numeric(conditions[3])*sqrt(n_vars),
      superblock = FALSE,
      ncomp = R,
      scheme = "factorial",
      comp_orth = FALSE, # Not used in this study, could be investigated in future work if this is required
      verbose = FALSE
    )
  }
  
  # Extract weights
  estimatedW <- matrix(0, ncol = R, nrow = J)
  block_weights <- fit_sgcca$a$bl1

  # SGCCA sometimes drops variables that are fully zero so you need to reconstruct the weights matrix
  if (is.null(rownames(block_weights))) {
    if (is.vector(block_weights)) {
      estimatedW[, 1:length(block_weights)] <- block_weights
    } else if (is.matrix(block_weights)) {
      estimatedW[,] <- block_weights
    }
  } else {
    # Extract row indices from block_weights row names like "V1_10"
    row_indices <- as.integer(sub("V1_", "", rownames(block_weights)))

    # Number of columns in block_weights
    n_col_bw <- ncol(block_weights)

    # Assign block_weights into estimatedW by row index and column index
    for (l in seq_along(row_indices)) {
      row_idx <- row_indices[l]
      estimatedW[row_idx, 1:n_col_bw] <- block_weights[l, ]
    }
  }
  
  # Calculate Scores and loadings 
  T_scores <- X%*%estimatedW
  P_estimated<- t(solve(t(T_scores) %*% T_scores) %*% t(T_scores) %*% X)  # R x p

  # Reconstruction loss
  X_hat <- T_scores %*% t(P_estimated)
  recon_mse <- sum((X - X_hat)^2)
  
  # Store results
  best_loss <- recon_mse
  
  # Similarity metrics
  sim_weights_loadings <- score_metrics(estimatedW, P_estimated)
  sim_weights_true <- score_metrics(estimatedW, W_true)
  sim_scores_true <- score_metrics(T_scores, T_true)
  sim_p_est_p_true <- score_metrics(P_estimated, P_true)
  
  # Reconstruction metrics
  recon_metrics <- reconstruction_metrics(X, estimatedW, P_estimated)
  
  # Variable selection evaluation
  selection_eval <- evaluate_variable_selection(W_true, estimatedW)
  
  # Bias-variance-MSE
  bias_var_mse <- compute_bias_variance_mse(W_true, estimatedW)
  
  # Sparsity
  sparsity <- sparsity_level(estimatedW)
  # Reconstruction and selection metrics
  bias_var_mse_W <- compute_bias_variance_mse(W_true, estimatedW)
  bias_var_mse_P <- compute_bias_variance_mse(P_true, P_estimated)
  
  #VAF
  vaf = compute_vaf(X, estimatedW, P_estimated)
  
  # Correlation comp 1 and comp 2
  inter_comp_corr <- cor(T_scores)[1, 2]
  
  # Save results for this dataset
  results_list[[i]] <- data.frame(
    Dataset = i,
    Sample_Size = I,
    Num_Items = J,
    Num_Latent = R,
    Condition1 = conditions[1],
    Condition2 = conditions[2],
    Condition3 = conditions[3],
    Final_Loss = best_loss,
    VAFx = conditions$VAFx,
    
    # Similarity metrics
    P_vs_Ptrue_MAE = sim_p_est_p_true$mae,
    P_vs_Ptrue_RMSE = sim_p_est_p_true$rmse,
    
    W_vs_Loadings_MAE = sim_weights_loadings$mae,
    W_vs_Loadings_RMSE = sim_weights_loadings$rmse,
    
    W_vs_Wtrue_MAE = sim_weights_true$mae,
    W_vs_Wtrue_RMSE = sim_weights_true$rmse,
    
    Score_vs_True_MAE = sim_scores_true$mae,
    Score_vs_True_RMSE = sim_scores_true$rmse,
    
    # Reconstruction
    MSE_Recon = recon_metrics$mse,
    R2_Recon = recon_metrics$R2,
    
    # Selection metrics
    Precision = selection_eval$precision,
    Recall = selection_eval$recall,
    F1_Score = selection_eval$f1,
    Recovery = selection_eval$recovery,
    
    # Bias, variance, mse weights
    Bias_W = bias_var_mse_W$bias,
    Variance_W = bias_var_mse_W$variance,
    MSE_W = bias_var_mse_W$mse,
    Bias_P = bias_var_mse_P$bias,
    Variance_P = bias_var_mse_P$variance,
    MSE_P = bias_var_mse_P$mse,
    
    #Vaf and corr
    VAF = vaf,
    Correlation_Comp1_Comp2 = inter_comp_corr,
    
    Sparsity = sparsity
    
  )
  cat("Dataset", i, "completed with best loss =", best_loss, "\n")
}

# Combine all results
results_table <- do.call(rbind, results_list)

# Average over the 5 repetition of conditions
results_table <- results_table %>%
  mutate(
    n_variables = n_variables,
    s_size = s_size,
    n_components = Num_Latent,
    p_sparse = p_sparse,
    VAFx = VAFx
  )

# Group by the specified columns and compute mean of all others
summary_table <- results_table %>%
  select(-Dataset) %>%
  group_by(n_variables,Num_Latent,s_size, p_sparse, VAFx) %>%
  summarise(across(everything(), mean, na.rm = TRUE), .groups = "drop")

# View or save result
print(summary_table)
write.csv(summary_table, "SGGCA_W_Sparse.csv", row.names = FALSE)

