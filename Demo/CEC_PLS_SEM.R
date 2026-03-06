# CEC-PLS-SEM Research
# Updated Version March 2026 Hetvi Chaniyara
# Runs CEC-PLS-SEM for various data folders
# Incorporates the changes as proposed in Katrijn's December 2025 version

current_working_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(current_working_dir)
getwd()

source("../Scripts/CEC_PLS_SEM_Functions.R")
library(dplyr)

# list of all data folders we want to run CEC-PLS-SEM for
folders <- c("../Scripts/DATA-R-W-Sparse", "../Scripts/DATA-R-P-Sparse")
results_list <- list()

for (f in folders) {
  # load the design parameters for the data
  load(file.path(f, "Info_simulation.RData"))
  design <- Info_simulation$design_matrix_replication
  prefix <- ifelse(grepl("W-Sparse", f), "Wsparse", "Psparse")
  
  for (i in 1:36) {
    data_file <- file.path(f, paste0(prefix, i, ".RData"))
    if(!file.exists(data_file)) {
      cat("File missing:", data_file, "\n")
      next
    }
    load(data_file)
    
    # initialize phi and rho
    X <- out$X; R <- out$k; J <- ncol(X)
    phi <- round((1 - design$p_sparse[i]) * J * R)
    rho <- sum(X^2) / R
    
    # run with multistart
    best_res <- NULL; best_loss <- Inf
    for (m in 1:2) {
      set.seed(100 + m)
      res <- CEC_PLS_SEM(X, R, 1e-8, phi, rho, constrained=0, MaxIter=100)
      if (res$Residual < best_loss) { best_loss <- res$Residual; best_res <- res }
    }
    
    # metrics calculation
    W_aligned <- align_components(best_res$weights, out$W)
    P_aligned <- align_components(best_res$loadings, out$P)
    selection <- evaluate_variable_selection(out$W, W_aligned)
    bvm_W <- compute_bias_variance_mse(out$W, W_aligned)
    bvm_P <- compute_bias_variance_mse(out$P, P_aligned) # Added for P comparison
    
    results_list[[length(results_list)+1]] <- data.frame(
      Folder = f, 
      Dataset = i, 
      design[i,],
      Loss = best_res$Residual,
      VAF = compute_vaf(X, best_res$weights, best_res$loadings),
      Recovery_Rate = selection$recovery,
      MSE_W = bvm_W$mse,
      MSE_P = bvm_P$mse, 
      W_Corr = diag(cor(W_aligned, out$W)) %>% mean(),
      P_Corr = diag(cor(P_aligned, out$P)) %>% mean(), 
      Iterations = best_res$n_iterations
    )
    cat("Folder:", f, "Dataset:", i, "Complete\n")
  }
}

# summarise results
final_summary <- do.call(rbind, results_list) %>%
  group_by(Folder, n_variables, s_size, p_sparse, VAFx) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE), .groups = "drop")

write.csv(final_summary, "Simulation_Results_Consolidated_penalised.csv", row.names = FALSE)
