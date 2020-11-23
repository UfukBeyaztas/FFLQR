##############Packages##############
library(fda)                       #
library(quantreg)                  #
library(matrixStats)               #
####################################

# Quantile regression function
qr_fun = function(response, predictor, predictor_test, nbf_response, nbf_vec_predictors, tau){
  # response: a matrix containing functional response variable
  # predictor: a list containing functional predictor variables (training sample)
  # predictor_test: a list containing functional predictor variables (test sample)
  # nbf_response: number of basis functions to approximate the functional response
  # nbf_predictors: a vector containing numbers of basis functions to approximate the functional predictors
  # tau: quantile level

  # Number of predictors
  np = length(predictor)
  
  # Discrete time points
  dtpy = seq(0, 1, length=ncol(response))
  dptx = list()
  for(i in 1:np)
    dptx[[i]] = seq(0, 1, length=ncol(predictor[[i]]))
  
  # B-spline for response
  B_spline_basis_y = create.bspline.basis(c(0,1), nbasis = nbf_response)
  B_spline_basis_fun_y = eval.basis(dtpy, B_spline_basis_y)
  
  # B-spline for predictors
  B_spline_basis_x = vector("list",)
  B_spline_basis_funs_x = vector("list",)
  
  for(i in 1:np){
    B_spline_basis_x[[i]] = create.bspline.basis(c(0,1), nbasis = nbf_vec_predictors[i])
    B_spline_basis_funs_x[[i]] = eval.basis(dptx[[i]], B_spline_basis_x[[i]])
  }
  
  # Inner products
  Inner_prod_y = inprod(B_spline_basis_y, B_spline_basis_y)
  Inner_prod_x = vector("list",)
  
  for(i in 1:np)
    Inner_prod_x[[i]] = inprod(B_spline_basis_x[[i]], B_spline_basis_x[[i]])

  # Weight arguments
  w_arg = matrix(dtpy, nrow = nrow(response), ncol = ncol(response), byrow=T)
  w_arg_x_train = list()
  for(i in 1:np)
    w_arg_x_train[[i]] = matrix(dptx[[i]], nrow = nrow(predictor[[i]]), ncol = ncol(predictor[[i]]), byrow=T)

  # Weight matrices of the response and predictors
  W_y = t(smooth.basis(argvals=t(w_arg), y=t(response), fdParobj=B_spline_basis_y)$fd$coefs)
  W_x = vector("list",)

  for(i in 1:np){
    W_x[[i]] = t(smooth.basis(argvals=t(w_arg_x_train[[i]]), y=t(predictor[[i]]), fdParobj=B_spline_basis_x[[i]])$fd$coefs)
  }

  
  ##################################################Variable selection##################################################
  all_effects = vector("list", length = np)
  for(all in 1:np){
    all_effects[[all]] = W_x[[all]] %*% Inner_prod_x[[all]]
  }

  # Define the importances of the predictors
  error_all = numeric()
  for(order in 1:np){
    # Regression matrix
    Reg_mat_x = all_effects[[order]]
    
    # Model estimation
    coef_q = rqs.fit(x = cbind(Reg_mat_x), y = W_y, tau = tau, tol = 0.001)
    
    # Estimation of the coefficient function
    coef_hat = t(B_spline_basis_fun_y %*% coef_q %*% t(B_spline_basis_funs_x[[order]]))
    
    # Estimation of the intercept function
    int_hat = apply(response, 2, quantile, probs = tau) - (apply(predictor[[order]], 2, quantile, probs = tau) %*% coef_hat) / length(dtpy)
    
    # Predicted functions
    pred_fun_q = (predictor[[order]] %*% coef_hat) / length(dtpy)
    for(i in 1:dim(pred_fun_q)[1])
      pred_fun_q[i,] = pred_fun_q[i,] + int_hat
    
    error_all[order] = mean((response - pred_fun_q)^2)
  }
  
  # Order of variables according to their MSEs
  imp_var = order(error_all)
  
  # Forward procedure
  # Starting model
  Reg_mat_forw_start = all_effects[[imp_var[1]]]
  forward_error = min(error_all)
  
  np_forw = c(imp_var[1], rep(NA, (length(imp_var)-1)))
  selected_main = c(which.min(error_all))
  
  for(forw1 in 2:np){
    error_for_forward_selection = rbind(subset(imp_var, !(imp_var %in% selected_main)), NA)
    for(forw in 1:ncol(error_for_forward_selection)){
      # Regression matrix
      Reg_mat_forw_model = cbind(Reg_mat_forw_start, all_effects[[error_for_forward_selection[1,forw]]])
      
      # Model estimation
      coef_q = rqs.fit(x = Reg_mat_forw_model, y = W_y, tau = tau, tol = 0.001)
      
      # Estimation of the coefficient function
      frxs = c(np_forw[!(np_forw %in% NA)], error_for_forward_selection[1,forw])
      fs = list()
      for(i in 1:length(frxs))
        fs[[i]] = 1:nbf_vec_predictors[frxs[i]]
      cnms = list()
      k = 0
      for(i in 1:length(np_forw)){
        cnms[[i]] = fs[[i]]+k
        k = cnms[[i]][length(cnms[[i]])]
      }
      
      lv = length(cnms)
      
      coef_hat = list()
      for(i in 1: lv)
        coef_hat[[i]] = t(B_spline_basis_fun_y %*% coef_q[,cnms[[i]]] %*% t(B_spline_basis_funs_x[[frxs[i]]]))
      
      # Estimation of the intercept function
      int_hat = apply(response, 2, quantile, probs = tau) - 
        Reduce("+", lapply(1:lv,  function(k){apply(predictor[[frxs[k]]], 2, quantile, probs = tau) %*% coef_hat[[k]] / length(dtpy)}))
      
      
      # Predicted functions
      pred_fun_q = Reduce("+", lapply(1:lv,  function(k){predictor[[frxs[k]]] %*% coef_hat[[k]] / length(dtpy)}))
      for(i in 1:dim(pred_fun_q)[1])
        pred_fun_q[i,] = pred_fun_q[i,] + int_hat
      
      error_for_forward_selection[2,forw] = mean((response - pred_fun_q)^2)
    }
    
    err_next = error_for_forward_selection[2,][which.min(error_for_forward_selection[2,])]
    var_next = error_for_forward_selection[1,][which.min(error_for_forward_selection[2,])]
    
    selected_main = c(selected_main,var_next)
    
    if(err_next < forward_error){
      Reg_mat_forw_start = cbind(Reg_mat_forw_start, all_effects[[var_next]])
      forward_error = err_next
      np_forw[forw1] = var_next
    }else if(err_next > forward_error){
      Reg_mat_forw_start = Reg_mat_forw_start
      forward_error = forward_error
      np_forw[forw1] = np_forw[forw1]
    }
  }
  
  np_forw = sort(subset(np_forw, !(np_forw %in% NA)))
  
  ##################################################Variable selection##################################################
  
  # Matrices for the regressions (selected variables)
  Reg_mat_q_selected = vector("list",)

  for(i in np_forw){
    Reg_mat_q_selected[[i]] = W_x[[i]] %*% Inner_prod_x[[i]]
  }
  
  Reg_mat_q_selected = do.call(cbind, Reg_mat_q_selected)
  
  # Matrices for the regressions (all variables)
  Reg_mat_q_all = vector("list",)

  for(i in 1:np){
    Reg_mat_q_all[[i]] = W_x[[i]] %*% Inner_prod_x[[i]]
  }
  
  Reg_mat_q_all = do.call(cbind, Reg_mat_q_all)

  # Model estimation
  # Selected model
  coef_q_selected = rqs.fit(x = Reg_mat_q_selected, y = W_y, tau = tau, tol = 0.001)
  # Full model
  coef_q_all = rqs.fit(x = Reg_mat_q_all, y = W_y, tau = tau, tol = 0.001)

  # Estimation of the coefficient function
  # Selected model
  fs = list()
  for(i in 1:length(np_forw))
    fs[[i]] = 1:nbf_vec_predictors[np_forw[i]]
  cnms = list()
  k = 0
  for(i in 1:length(np_forw)){
    cnms[[i]] = fs[[i]]+k
    k = cnms[[i]][length(cnms[[i]])]
  }
  
  lv = length(cnms)
  
  coef_hat_selected = list()
  for(i in 1: lv)
    coef_hat_selected[[i]] = t(B_spline_basis_fun_y %*% coef_q_selected[,cnms[[i]]] %*% t(B_spline_basis_funs_x[[np_forw[i]]]))
  
  # Estimation of the intercept function
  int_hat_selected = apply(response, 2, quantile, probs = tau) - 
    Reduce("+", lapply(1:lv,  function(k){apply(predictor[[np_forw[k]]], 2, quantile, probs = tau) %*% coef_hat_selected[[k]] / length(dtpy)}))
  
  
  # Fitted functions
  fit_fun_q_selected = Reduce("+", lapply(1:lv,  function(k){predictor[[np_forw[k]]] %*% coef_hat_selected[[k]] / length(dtpy)}))
  for(i in 1:dim(fit_fun_q_selected)[1])
    fit_fun_q_selected[i,] = fit_fun_q_selected[i,] + int_hat_selected
  
  # Full model
  fs = list()
  for(i in 1:length(predictor))
    fs[[i]] = 1:nbf_vec_predictors[i]
  cnms_full = list()
  k = 0
  for(i in 1:length(predictor)){
    cnms_full[[i]] = fs[[i]]+k
    k = cnms_full[[i]][length(cnms_full[[i]])]
  }
  
  lv_full = length(cnms_full)
  
  coef_hat_full = list()
  for(i in 1: lv_full)
    coef_hat_full[[i]] = t(B_spline_basis_fun_y %*% coef_q_all[,cnms_full[[i]]] %*% t(B_spline_basis_funs_x[[i]]))
  
  # Intercept function
  int_hat_full = apply(response, 2, quantile, probs = tau) - 
    Reduce("+", lapply(1:lv_full,  function(k){apply(predictor[[k]], 2, quantile, probs = tau) %*% coef_hat_full[[k]] / length(dtpy)}))
  
  # Fitted functions
  fit_fun_q_full = Reduce("+", lapply(1:lv_full,  function(k){predictor[[k]] %*% coef_hat_full[[k]] / length(dtpy)}))
  for(i in 1:dim(fit_fun_q_full)[1])
    fit_fun_q_full[i,] = fit_fun_q_full[i,] + int_hat_full
  
  # Predicted functions
  # Selected model
  pred_fun_q_selected = Reduce("+", lapply(1:lv,  function(k){predictor_test[[np_forw[k]]] %*% coef_hat_selected[[k]] / length(dtpy)}))
  for(i in 1:dim(pred_fun_q_selected)[1])
    pred_fun_q_selected[i,] = pred_fun_q_selected[i,] + int_hat_selected
  
  # Full model
  pred_fun_q_full = Reduce("+", lapply(1:lv_full,  function(k){predictor_test[[k]] %*% coef_hat_full[[k]] / length(dtpy)}))
  for(i in 1:dim(pred_fun_q_full)[1])
    pred_fun_q_full[i,] = pred_fun_q_full[i,] + int_hat_full

  
  return(list("fits_selected"=fit_fun_q_selected, "fits_full"=fit_fun_q_full,
              "preds_selected"=pred_fun_q_selected, "preds_full"=pred_fun_q_full))
}

interval_score <- function(holdout, lb, ub, alpha){
  lb_ind = ifelse(holdout < lb, 1, 0)
  ub_ind = ifelse(holdout > ub, 1, 0)
  score = (ub - lb) + 2/alpha * ((lb - holdout) * lb_ind + (holdout - ub) * ub_ind)
  cover = 1 - (length(which(lb_ind == 1)) + length(which(ub_ind == 1)))/length(holdout)
  cpd = abs(cover - (1 - alpha))
  return(c(mean(score), cpd))
}
