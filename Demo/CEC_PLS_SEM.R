# Cardinality and Equality Constrained PlS-SEM
# Hetvi Chaniyara
# Bachelor End Project

current_working_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(current_working_dir)
getwd()

# install.packages("MASS")
# install.packages("gtools")
# install.packages("dplyr")
library(MASS)
library(gtools)
library(dplyr)

# Load functions from functions file
source("~/BEP-Project-CEC-PLS-SEM/Scripts/CEC_PLS_SEM_Functions.R")

n_multistarts <- 3 # Change for more multistarts
base_seed <- 100

# Load the design Matrix
load("~/BEP-Project-CEC-PLS-SEM/Scripts/Data-R-P-Sparse-W-not-P/Info_simulation.RData")
Info_matrix <- Info_simulation$design_matrix_replication
Ndatasets <- Info_simulation$n_data_sets

# Initalize list to get results
results_list <- list()

# Repeat loop for every file in the data folder
for (i in 1:Ndatasets) {
  
  # load the data folder
  filename <- paste0("~/BEP-Project-CEC-PLS-SEM/Scripts/Data-R-P-Sparse-W-not-P/Psparse", i, ".RData") # Change data directory accordingly
  load(filename)
  
  # Extract true data
  X <- out$X
  W_true <- out$W
  P_true <- out$P
  T_true <- out$Z
  R <- dim(W_true)[2]
  conditions <- Info_matrix[i, ]
  phi <- as.numeric((1 - conditions[3]) * dim(X)[2])
  rho <- 0.5
  
  # Multistart implementation
  best_result <- NULL
  best_loss <- Inf
  
  for (m in 1:n_multistarts) {
    set.seed(base_seed + m) # For reproducibility
    result <- CEC_PLS_SEM(X, R, epsilon = 1e-5, phi, rho)
    
    if (result$Residual < best_loss) {
      best_loss <- result$Residual
      best_result <- result
    }
  }
  
  
  # Align estimated matrices
  aligned_weights <- align_components(best_result$weights, W_true)
  aligned_loadings <- align_components(t(best_result$loadings), P_true)
  aligned_T <- align_components(best_result$Scores, T_true)
  
  # Compute similarity metrics
  sim_weights_true <- score_metrics(aligned_weights, W_true)
  sim_p_est_p_true <- score_metrics(aligned_loadings, P_true)
  sim_scores_true <- score_metrics(aligned_T, T_true)
  sim_weights_loadings <- score_metrics(best_result$weights, t(best_result$loadings))
  
  
  # Reconstruction and selection metrics
  bias_var_mse_W <- compute_bias_variance_mse(W_true, aligned_weights)
  recon_metrics <- reconstruction_metrics(X, best_result$weights, best_result$loadings)
  selection_eval <- evaluate_variable_selection(W_true, best_result$weights)
  bias_var_mse_P <- compute_bias_variance_mse(P_true, aligned_loadings)
  
  #VAF
  vaf <- compute_vaf(X, best_result$weights, best_result$loadings)
  
  # Correlation comp 1 and comp 2
  inter_comp_corr <- cor(aligned_T)[1, 2]
  
  # Sparsity
  sparsity <- sparsity_level(best_result$weights)
  
  # Store all results in consistent format
  results_list[[i]] <- data.frame(
    Dataset = i,
    n_variables = conditions$n_variables,
    s_size = conditions$s_size,
    p_sparse = conditions$p_sparse,
    n_components = conditions$n_components,
    VAFx = conditions$VAFx,
    Final_Loss = best_result$Residual,
    Num_Iterations = best_result$n_iterations,
    
    # Similarity metrics
    P_vs_Ptrue_MAE = sim_p_est_p_true$mae,
    P_vs_Ptrue_RMSE = sim_p_est_p_true$rmse,
    P_vs_Ptrue_Corr = sim_p_est_p_true$correlation,
    
    W_vs_Loadings_MAE = sim_weights_loadings$mae,
    W_vs_Loadings_RMSE = sim_weights_loadings$rmse,
    W_vs_Loadings_Corr = sim_weights_loadings$correlation,
    
    W_vs_Wtrue_MAE = sim_weights_true$mae,
    W_vs_Wtrue_RMSE = sim_weights_true$rmse,
    W_vs_Wtrue_Corr = sim_weights_true$correlation,
    
    Score_vs_True_MAE = sim_scores_true$mae,
    Score_vs_True_RMSE = sim_scores_true$rmse,
    Score_vs_True_Corr = sim_scores_true$correlation,
    
    # Reconstruction
    MSE_Recon = recon_metrics$mse,
    R2_Recon = recon_metrics$R2,
    
    # Selection
    Precision = selection_eval$precision,
    Recall = selection_eval$recall,
    F1_Score = selection_eval$f1,
    Recovery = selection_eval$recovery,
    Sparsity = sparsity,
    
    # Bias-variance-MSE
    Bias_W = bias_var_mse_W$bias,
    Variance_W = bias_var_mse_W$variance,
    MSE_W = bias_var_mse_W$mse,
    Bias_P = bias_var_mse_P$bias,
    Variance_P = bias_var_mse_P$variance,
    MSE_P = bias_var_mse_P$mse,
    
    #Vaf and corr
    VAF = vaf,
    Correlation_Comp1_Comp2 = inter_comp_corr
  )
  cat("Dataset", i, "completed with best loss =", best_result$Residual, "\n")
}

# Collect all results in a dataframe
results_table <- do.call(rbind, results_list)

# Average of the 5 iteration for each condition
results_table <- results_table %>%
  mutate(
    n_variables = n_variables,
    s_size = s_size,
    n_components = n_components,
    p_sparse = p_sparse,
    VAFx = VAFx
  )

# Group by the specified columns and compute mean of all others
summary_table <- results_table %>%
  select(-Dataset) %>%
  group_by(n_variables,n_components,s_size, p_sparse, VAFx) %>%
  summarise(across(everything(), mean, na.rm = TRUE), .groups = "drop")

# Save and view table
print(summary_table)
write.csv(summary_table, "CECPLSSEM-Psparse-W-not-P.csv", row.names = FALSE)
