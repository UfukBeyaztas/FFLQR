source("auxiliary_f_qr_dti.R")
library(refund)

# Number of simulations
nsim = 2
# Number of functions in training sample
n_train = 33
# Number of functions in test sample
n_test = 33
# Number of bootstrap simulation
B = 100
# Nominal level
alpha = 0.05


# Number of basis functions for predictor
nbasis_est_x = 10
# Number of basis functions for response
nbasis_est_y = 15


# MSPE
mspe_qr = numeric()
# RMSPE
rmspe_qr = numeric()
# MAPE
mape_qr = numeric()
# For prediction interval
score_rq = matrix(NA, nrow = nsim, ncol = 2)

for(sim in 1:nsim){
  
  train_index = sample(1: (n_train+n_test), n_train, replace = FALSE)
  DTI.complete <- subset(DTI, complete.cases(DTI))
  DTI.complete = DTI.complete[DTI.complete$visit == 1,]
  
  Y_train= DTI.complete$cca[train_index,]
  Y_test= DTI.complete$cca[-(train_index),]
  
  X_train_list=list(DTI.complete$rcst[train_index,])
  X_test_list=list(DTI.complete$rcst[-(train_index),])
  
  qr_model = qr_fun(response = Y_train, predictor = X_train_list, predictor_test = X_test_list, nbf_response = nbasis_est_y, nbf_predictors = nbasis_est_x,
                    tau = 0.5)
  
  # MSPE
  mspe_qr[sim] = mean((Y_test - qr_model$preds)^2)

  # RMSPE
  rmspe_qr[sim] = sqrt(mean(((Y_test - qr_model$preds)/Y_test)^2))

  # MAPE
  mape_qr[sim] = mean(abs((Y_test - qr_model$preds)/Y_test))

  boot_preds_qr = vector("list",)

  for(boot in 1:B){
    boot_index = sample(1:n_train, n_train, replace = TRUE)
    
    Y_boot_train = Y_train[boot_index,]
    X_boot_train = list(X_train_list[[1]][boot_index,])
    
    boot_qr_model = qr_fun(response = Y_boot_train, predictor = X_boot_train, predictor_test = X_test_list, nbf_response = nbasis_est_y, nbf_predictors = nbasis_est_x,
                           tau = 0.5)
    
    boot_preds_qr[[boot]] = boot_qr_model$preds
  }
  
  
  lq_matrix_qr = matrix(NA, nrow = nrow(Y_test), ncol = ncol(Y_test))
  uq_matrix_qr = matrix(NA, nrow = nrow(Y_test), ncol = ncol(Y_test))

  for(iq in 1:n_test){
    
    iq_matrix_qr = matrix(NA, nrow = B, ncol = ncol(Y_test))

    for(jq in 1:B){
      
      iq_matrix_qr[jq,] = boot_preds_qr[[jq]][iq,]
    }
    
    lq_matrix_qr[iq,] = apply(iq_matrix_qr, 2, quantile, probs = alpha/2)
    uq_matrix_qr[iq,] = apply(iq_matrix_qr, 2, quantile, probs = (1-alpha/2))
  }
  
  int_score_rq = matrix(NA, nrow = nrow(Y_test), ncol = 2)

  for(is in 1:nrow(Y_test)){
    
    int_score_rq[is,] = interval_score(Y_test[is,], lq_matrix_qr[is,], uq_matrix_qr[is,], alpha)
  }
  
  score_rq[sim,] = apply(int_score_rq, 2, mean)
}

