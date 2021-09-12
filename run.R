rm(list=ls())
source("dgp1.R")
source("auxiliary_functions.R")


# Number of discrete time points
t_j = 100
# Equally spaced time points in the interval [0, 1]
tp_seq = seq(0, 1, length = t_j)
# Number of predictors
n_pred = 5
# Number of functions in training sample
n_train = 200
# Number of functions in test sample
n_test = 300

# Lag
nlag = 4
# Nois level
sigma2.y = 1
# Tau level
ltau = 0.5
# Significance level
alpha = 0.05

true_index = c(2,4,5)


# Number of basis functions for predictor
K_X = rep(20, n_pred)
# Number of basis functions for response
K_Y = 20

# Rangeval values for the response and predictors
rangeval_Y = c(0,1)
rangeval_X = list(c(0,1),c(0,1),c(0,1),c(0,1),c(0,1))

set.seed(12345) # for reproduction
# Generate the data under Gaussian errors
sim_dat = DGP(n.curves = n_pred, ntest = n_test, ntrain = n_train, lengthX = t_j, lengthY = t_j, sigma2.e = sigma2.y, Lag = nlag)

Y = sim_dat$Y_train
Y_test = sim_dat$fun_Y_test
X = sim_dat$X_train
X_test = sim_dat$X_test

# FFLQR fit under full model
model_qr_full = fpca(fY = Y, fX = X, fX_test = X_test, fmodel = "full", tau = ltau, rangevalY = rangeval_Y, rangevalX = rangeval_X,
                     fnbasisY = K_Y, fnbasisX = K_X, fncomp_Y = 5, fncomp_X =5)

# FFLQR fit under true model
model_qr_true = fpca(fY = Y, fX = X, fX_test = X_test, fmodel = "true", tau = ltau, rangevalY = rangeval_Y, rangevalX = rangeval_X,
                     mindex = true_index, fnbasisY = K_Y, fnbasisX = K_X, fncomp_Y = 5, fncomp_X =5)
# FFLQR fit under selected model
model_qr_selected = fpca(fY = Y, fX = X, fX_test = X_test, fmodel = "selected", tau = ltau, rangevalY = rangeval_Y, rangevalX = rangeval_X,
                         fnbasisY = K_Y, fnbasisX = K_X, fncomp_Y = 5, fncomp_X =5)


# MSPE under true model
mean((Y_test - model_qr_full$pred)^2) # 0.04557994
# MSPE under true model
mean((Y_test - model_qr_true$pred)^2) # 0.03053818
# MSPE under true model
mean((Y_test - model_qr_selected$pred)^2) # 0.03906812

# Plot of the functional response in the test sample
ts.plot(t(Y_test))
# Plots of the predicted functional response in the test sample (selected model)
ts.plot(t(model_qr_selected$pred))


# To construct prediction intervals directly using the proposed method fit two different models at tau levels alpha/2 and (1-alpha/2)

model_qr_true_q1 = fpca(fY = Y, fX = X, fX_test = X_test, fmodel = "true", tau = alpha/2, rangevalY = rangeval_Y, rangevalX = rangeval_X,
                        mindex = true_index, fnbasisY = K_Y, fnbasisX = K_X, fncomp_Y = 5, fncomp_X =5)
model_qr_full_q1 = fpca(fY = Y, fX = X, fX_test = X_test, fmodel = "full", tau = alpha/2, rangevalY = rangeval_Y, rangevalX = rangeval_X,
                        fnbasisY = K_Y, fnbasisX = K_X, fncomp_Y = 5, fncomp_X =5)
model_qr_selected_q1 = fpca(fY = Y, fX = X, fX_test = X_test, fmodel = "selected", tau = alpha/2, rangevalY = rangeval_Y, rangevalX = rangeval_X,
                            fnbasisY = K_Y, fnbasisX = K_X, fncomp_Y = 5, fncomp_X =5)
model_qr_true_q2 = fpca(fY = Y, fX = X, fX_test = X_test, fmodel = "true", tau = (1-alpha/2), rangevalY = rangeval_Y, rangevalX = rangeval_X,
                        mindex = true_index, fnbasisY = K_Y, fnbasisX = K_X, fncomp_Y = 5, fncomp_X =5)
model_qr_full_q2 = fpca(fY = Y, fX = X, fX_test = X_test, fmodel = "full", tau = (1-alpha/2), rangevalY = rangeval_Y, rangevalX = rangeval_X,
                        fnbasisY = K_Y, fnbasisX = K_X, fncomp_Y = 5, fncomp_X =5)
model_qr_selected_q2 = fpca(fY = Y, fX = X, fX_test = X_test, fmodel = "selected", tau = (1-alpha/2), rangevalY = rangeval_Y, rangevalX = rangeval_X,
                            fnbasisY = K_Y, fnbasisX = K_X, fncomp_Y = 5, fncomp_X =5)


# Compute the score and CPD values
int_score_qrq_full = matrix(NA, nrow = nrow(Y_test), ncol = 2)
int_score_qrq_true = matrix(NA, nrow = nrow(Y_test), ncol = 2)
int_score_qrq_selected = matrix(NA, nrow = nrow(Y_test), ncol = 2)

for(is in 1:nrow(Y_test)){
  int_score_qrq_full[is,] = interval_score(Y_test[is,], model_qr_full_q1$pred[is,], model_qr_full_q2$pred[is,], alpha)
  int_score_qrq_true[is,] = interval_score(Y_test[is,], model_qr_true_q1$pred[is,], model_qr_true_q2$pred[is,], alpha)
  int_score_qrq_selected[is,] = interval_score(Y_test[is,], model_qr_selected_q1$pred[is,], model_qr_selected_q2$pred[is,], alpha)
}

# score and CPD values under full model
apply(int_score_qrq_full, 2, mean) # 2.57259747 0.06736667
# score and CPD values under full model
apply(int_score_qrq_true, 2, mean) # 2.998603 0.056200
# score and CPD values under full model
apply(int_score_qrq_selected, 2, mean) # 2.77259709 0.04993333



# To obtain the results when the errors follow skewed chi-square(1) distribution
source("dgp2.R")


# Generate the data under chi-square(1) errors
sim_dat = DGP(n.curves = n_pred, ntest = n_test, ntrain = n_train, lengthX = t_j, lengthY = t_j, sigma2.e = sigma2.y, Lag = nlag)

Y = sim_dat$Y_train
Y_test = sim_dat$fun_Y_test
X = sim_dat$X_train
X_test = sim_dat$X_test

# FFLQR fit under full model
model_qr_full = fpca(fY = Y, fX = X, fX_test = X_test, fmodel = "full", tau = ltau, rangevalY = rangeval_Y, rangevalX = rangeval_X,
                     fnbasisY = K_Y, fnbasisX = K_X, fncomp_Y = 5, fncomp_X =5)

# FFLQR fit under true model
model_qr_true = fpca(fY = Y, fX = X, fX_test = X_test, fmodel = "true", tau = ltau, rangevalY = rangeval_Y, rangevalX = rangeval_X,
                     mindex = true_index, fnbasisY = K_Y, fnbasisX = K_X, fncomp_Y = 5, fncomp_X =5)
# FFLQR fit under selected model
model_qr_selected = fpca(fY = Y, fX = X, fX_test = X_test, fmodel = "selected", tau = ltau, rangevalY = rangeval_Y, rangevalX = rangeval_X,
                         fnbasisY = K_Y, fnbasisX = K_X, fncomp_Y = 5, fncomp_X =5)


# MSPE under true model
mean((Y_test - model_qr_full$pred)^2) # 0.4538092
# MSPE under true model
mean((Y_test - model_qr_true$pred)^2) # 0.3418706
# MSPE under true model
mean((Y_test - model_qr_selected$pred)^2) # 0.3121307

# Plot of the functional response in the test sample
ts.plot(t(Y_test))
# Plots of the predicted functional response in the test sample (selected model)
ts.plot(t(model_qr_selected$pred))


# To construct prediction intervals directly using the proposed method fit two different models at tau levels alpha/2 and (1-alpha/2)

model_qr_true_q1 = fpca(fY = Y, fX = X, fX_test = X_test, fmodel = "true", tau = alpha/2, rangevalY = rangeval_Y, rangevalX = rangeval_X,
                        mindex = true_index, fnbasisY = K_Y, fnbasisX = K_X, fncomp_Y = 5, fncomp_X =5)
model_qr_full_q1 = fpca(fY = Y, fX = X, fX_test = X_test, fmodel = "full", tau = alpha/2, rangevalY = rangeval_Y, rangevalX = rangeval_X,
                        fnbasisY = K_Y, fnbasisX = K_X, fncomp_Y = 5, fncomp_X =5)
model_qr_selected_q1 = fpca(fY = Y, fX = X, fX_test = X_test, fmodel = "selected", tau = alpha/2, rangevalY = rangeval_Y, rangevalX = rangeval_X,
                            fnbasisY = K_Y, fnbasisX = K_X, fncomp_Y = 5, fncomp_X =5)
model_qr_true_q2 = fpca(fY = Y, fX = X, fX_test = X_test, fmodel = "true", tau = (1-alpha/2), rangevalY = rangeval_Y, rangevalX = rangeval_X,
                        mindex = true_index, fnbasisY = K_Y, fnbasisX = K_X, fncomp_Y = 5, fncomp_X =5)
model_qr_full_q2 = fpca(fY = Y, fX = X, fX_test = X_test, fmodel = "full", tau = (1-alpha/2), rangevalY = rangeval_Y, rangevalX = rangeval_X,
                        fnbasisY = K_Y, fnbasisX = K_X, fncomp_Y = 5, fncomp_X =5)
model_qr_selected_q2 = fpca(fY = Y, fX = X, fX_test = X_test, fmodel = "selected", tau = (1-alpha/2), rangevalY = rangeval_Y, rangevalX = rangeval_X,
                            fnbasisY = K_Y, fnbasisX = K_X, fncomp_Y = 5, fncomp_X =5)


# Compute the score and CPD values
int_score_qrq_full = matrix(NA, nrow = nrow(Y_test), ncol = 2)
int_score_qrq_true = matrix(NA, nrow = nrow(Y_test), ncol = 2)
int_score_qrq_selected = matrix(NA, nrow = nrow(Y_test), ncol = 2)

for(is in 1:nrow(Y_test)){
  int_score_qrq_full[is,] = interval_score(Y_test[is,], model_qr_full_q1$pred[is,], model_qr_full_q2$pred[is,], alpha)
  int_score_qrq_true[is,] = interval_score(Y_test[is,], model_qr_true_q1$pred[is,], model_qr_true_q2$pred[is,], alpha)
  int_score_qrq_selected[is,] = interval_score(Y_test[is,], model_qr_selected_q1$pred[is,], model_qr_selected_q2$pred[is,], alpha)
}

# score and CPD values under full model
apply(int_score_qrq_full, 2, mean) # 6.7369230 0.1597667
# score and CPD values under full model
apply(int_score_qrq_true, 2, mean) # 3.90554459 0.07346667
# score and CPD values under full model
apply(int_score_qrq_selected, 2, mean) # 4.01211945 0.04986667
