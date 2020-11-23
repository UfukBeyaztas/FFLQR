source("auxiliary_f_qr_drosophila.R")

load("zygotic.RData")
load("muscle.RData")
load("eye.RData")

embriyo = zygotic_genes[,1:31]
larva = zygotic_genes[,32:41]
pupa = zygotic_genes[,42:58]

# Number of simulations
nsim = 2
# Number of predictors
n_pred = 2
# Number of functions in training sample
n_train = 10
# Number of functions in test sample
n_test = 11
# Number of bootstrap simulation
B = 100
# Nominal level
alpha = 0.05

# Number of basis functions for predictor
nbasis_est_x = c(6,4)
# Number of basis functions for response
nbasis_est_y = 6


# MSPE
mspe_full_qr = numeric()
mspe_selected_qr = numeric()


# RMSPE
rmspe_full_qr = numeric()
rmspe_selected_qr = numeric()

#MAPE
mape_full_qr = numeric()
mape_selected_qr = numeric()

# For prediction interval
score_rq_full = matrix(NA, nrow = nsim, ncol = 2)
score_rq_selected = matrix(NA, nrow = nsim, ncol = 2)

for(sim in 1:nsim){
  train_index = sample(1: (n_train+n_test), n_train, replace = FALSE)
  Y_train= pupa[train_index,]
  Y_test= pupa[-(train_index),]
  
  X_train_list = list(embriyo[train_index,], larva[train_index,])
  X_test_list = list(embriyo[-(train_index),], larva[-(train_index),])
  
  qr_model = qr_fun(response = Y_train, predictor = X_train_list, predictor_test = X_test_list, nbf_response = nbasis_est_y, nbf_vec_predictors = nbasis_est_x,
                     tau = 0.50)

  # MSPE
  mspe_full_qr[sim] = mean((Y_test - qr_model$preds_full)^2)
  mspe_selected_qr[sim] = mean((Y_test - qr_model$preds_selected)^2)

  # RMSPE
  rmspe_full_qr[sim] = sqrt(mean(((Y_test - qr_model$preds_full)/Y_test)^2))
  rmspe_selected_qr[sim] = sqrt(mean(((Y_test - qr_model$preds_selected)/Y_test)^2))

  #MAPE
  mape_full_qr[sim] = mean(abs((Y_test - qr_model$preds_full)/Y_test))
  mape_selected_qr[sim] = mean(abs((Y_test - qr_model$preds_selected)/Y_test))

  boot_preds_qr_full = vector("list",)
  boot_preds_qr_selected = vector("list",)

  for(boot in 1:B){
    boot_index = sample(1:n_train, n_train, replace = TRUE)
    
    Y_boot_train = Y_train[boot_index,]
    X_boot_train = vector("list",)
    for(ib in 1:n_pred)
      X_boot_train[[ib]] = X_train_list[[ib]][boot_index,]
    
    boot_qr_model = qr_fun(response = Y_boot_train, predictor = X_boot_train, predictor_test = X_test_list, nbf_response = nbasis_est_y, nbf_vec_predictors = nbasis_est_x,
                            tau = 0.50)
  
    boot_preds_qr_full[[boot]] = boot_qr_model$preds_full
    boot_preds_qr_selected[[boot]] = boot_qr_model$preds_selected
  }
  
  
  lq_matrix_qr_full = matrix(NA, nrow = nrow(Y_test), ncol = ncol(Y_test))
  uq_matrix_qr_full = matrix(NA, nrow = nrow(Y_test), ncol = ncol(Y_test))

  lq_matrix_qr_selected = matrix(NA, nrow = nrow(Y_test), ncol = ncol(Y_test))
  uq_matrix_qr_selected = matrix(NA, nrow = nrow(Y_test), ncol = ncol(Y_test))

  for(iq in 1:n_test){
    
    iq_matrix_qr_full = matrix(NA, nrow = B, ncol = ncol(Y_test))
    iq_matrix_qr_selected = matrix(NA, nrow = B, ncol = ncol(Y_test))

    for(jq in 1:B){
      
      iq_matrix_qr_full[jq,] = boot_preds_qr_full[[jq]][iq,]
      iq_matrix_qr_selected[jq,] = boot_preds_qr_selected[[jq]][iq,]
    }
    
    lq_matrix_qr_full[iq,] = apply(iq_matrix_qr_full, 2, quantile, probs = alpha/2)
    lq_matrix_qr_selected[iq,] = apply(iq_matrix_qr_selected, 2, quantile, probs = alpha/2)

    uq_matrix_qr_full[iq,] = apply(iq_matrix_qr_full, 2, quantile, probs = (1-alpha/2))
    uq_matrix_qr_selected[iq,] = apply(iq_matrix_qr_selected, 2, quantile, probs = (1-alpha/2))
  }
  
  int_score_rq_full = matrix(NA, nrow = nrow(Y_test), ncol = 2)
  int_score_rq_selected = matrix(NA, nrow = nrow(Y_test), ncol = 2)

  for(is in 1:nrow(Y_test)){
    
    int_score_rq_full[is,] = interval_score(Y_test[is,], lq_matrix_qr_full[is,], uq_matrix_qr_full[is,], alpha)
    int_score_rq_selected[is,] = interval_score(Y_test[is,], lq_matrix_qr_selected[is,], uq_matrix_qr_selected[is,], alpha)
  }
  
  score_rq_full[sim,] = apply(int_score_rq_full, 2, mean)
  score_rq_selected[sim,] = apply(int_score_rq_selected, 2, mean)

}
