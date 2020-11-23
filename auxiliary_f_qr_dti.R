##############Packages##############
library(fda)                       #
library(quantreg)                  #
library(matrixStats)               #
####################################

# Quantile regression function
qr_fun = function(response, predictor, predictor_test, nbf_response, nbf_predictors, tau){
  # response: a matrix containing functional response variable
  # predictor: a list containing functional predictor variables (training sample)
  # predictor_test: a list containing functional predictor variables (test sample)
  # nbf_response: number of basis functions to approximate the functional response
  # nbf_predictors: number of basis functions to approximate the functional predictors
  # tau: quantile level


  # Discrete time points
  dtpy = seq(0, 1, length=ncol(response))
  dtpx = seq(0, 1, length=ncol(predictor[[1]]))
  
  # B-spline for response
  B_spline_basis_y = create.bspline.basis(c(0,1), nbasis = nbf_response)
  B_spline_basis_fun_y = eval.basis(dtpy, B_spline_basis_y)
  
  # B-spline for predictors
  B_spline_basis_x = create.bspline.basis(c(0,1), nbasis = nbf_predictors)
  B_spline_basis_funs_x = eval.basis(dtpx, B_spline_basis_x)
  
  # Inner products
  Inner_prod_y = inprod(B_spline_basis_y, B_spline_basis_y)
  Inner_prod_x = inprod(B_spline_basis_x, B_spline_basis_x)

  # Weight arguments
  w_arg = matrix(dtpy, nrow = nrow(response), ncol = ncol(response), byrow=T)
  w_arg_x_train = matrix(dtpx, nrow = nrow(predictor[[1]]), ncol = ncol(predictor[[1]]), byrow=T)

  # Weight matrices of the response and predictors
  W_y = t(smooth.basis(argvals=t(w_arg), y=t(response), fdParobj=B_spline_basis_y)$fd$coefs)
  W_x = t(smooth.basis(argvals=t(w_arg_x_train), y=t(predictor[[1]]), fdParobj=B_spline_basis_x)$fd$coefs)

  # Matrix for the regression
  Reg_mat_q = W_x %*% Inner_prod_x

  # Model estimation
  coef_q = rqs.fit(x = Reg_mat_q, y = W_y, tau = tau, tol = 0.001)

  # Estimation of the coefficient function
  coef_hat = t(B_spline_basis_fun_y %*% coef_q %*% t(B_spline_basis_funs_x))

  # Estimation of the intercept function
  int_hat = apply(response, 2, quantile, probs = tau) - apply(predictor[[1]], 2, quantile, probs = tau) %*% coef_hat / length(dtpy)

  # Fitted functions
  fit_fun_q = predictor[[1]] %*% coef_hat / length(dtpy)
  for(i in 1:dim(fit_fun_q)[1])
    fit_fun_q[i,] = fit_fun_q[i,] + int_hat
  
  # Predicted functions
  pred_fun_q = predictor_test[[1]] %*% coef_hat / length(dtpy)
  for(i in 1:dim(pred_fun_q)[1])
    pred_fun_q[i,] = pred_fun_q[i,] + int_hat

  
  return(list("fits" = fit_fun_q, "preds" = pred_fun_q))
}

interval_score <- function(holdout, lb, ub, alpha){
  lb_ind = ifelse(holdout < lb, 1, 0)
  ub_ind = ifelse(holdout > ub, 1, 0)
  score = (ub - lb) + 2/alpha * ((lb - holdout) * lb_ind + (holdout - ub) * ub_ind)
  cover = 1 - (length(which(lb_ind == 1)) + length(which(ub_ind == 1)))/length(holdout)
  cpd = abs(cover - (1 - alpha))
  return(c(mean(score), cpd))
}
