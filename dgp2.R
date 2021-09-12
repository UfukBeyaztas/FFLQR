
# Daga generating process chi-square(1) errors 

# n.curves  : number of predictors
# ntest     : number of functions in testing phase
# ntrain    : number of functions in training phase
# lengthX   : number of discrete time points in X
# lengthY   : number of discrete time points in Y
# Lag       : Lag parameter that controls the correlation between the response variables
# sigma2.e  : noise level
DGP = function(n.curves, ntest, ntrain, lengthX, lengthY, Lag, sigma2.e){
  
  n.curves = n.curves
  ntest = ntest
  ntrain = ntrain
  lengthX = lengthX
  lengthY = lengthY
  Lag = Lag
  sigma2.e = sigma2.e

  ntot = ntest+ntrain
  t.x = seq(0,1,length=lengthX)
  t.y = seq(0,1,length=lengthY)
  total.main=n.curves
  
  
  gammaexpcov = function(n.point, t){ # the function used to calculate the covariance function of X_i(s)
    Sigma = matrix(0, n.point,n.point)
    
    for(i in 1:(n.point-1)){
      for(j in (i+1):n.point){
        Sigma[i,j] = exp(-(10*abs(t[j]-t[i]))^2)
      }
    }
  
    Sigma = (Sigma+t(Sigma))
    diag(Sigma) = 1
    return(Sigma)
  }
  
  true.main.index = c(2,4,5)
  
  f = function(x){
    exp(-5*x^2)
  }
  
  beta.fun = list() # the coefficient of the main effects in the true model
  beta.fun[[5]] = function(x1,t){
    sqrt(x1)*sqrt(t)
  }
  
  beta.fun[[4]] = function(x1,t){
    sin(1.5*pi*x1)*sin( pi*t)
  }
  
  beta.fun[[3]] = function(x1,t){
    f((x1-0.5))*f((t-0.5))+8*f((x1-1.5))*f((t-0.5))
  }
  
  beta.fun[[2]] = function(x1,t){
    exp(-3*(x1-1)^2) *exp(-5*(t-0.5)^2)
  }
  
  beta.fun[[1]] = function(x1,t){
      (1-x1)^2*(t-0.5)^2
  }
  
  beta.vals = lapply(1:length(beta.fun),function(k){outer(t.x,t.y,beta.fun[[k]])})
  
  
  Sigma.X = gammaexpcov(lengthX, t.x)
  
  eig = eigen(Sigma.X)
  eigval = eig$values
  b = cumsum(eigval)/sum(eigval)
  n.comp.sigma = which(b>0.999999999999)[1]
  eig0 = eigval[1:n.comp.sigma]
  sigma.5 =  eig$vectors[,1:n.comp.sigma] %*% sqrt(diag(eig0))
  range(abs(sigma.5 %*% t(sigma.5) - Sigma.X))
  
  
  V=list()
  
  for(i in 1:(n.curves+Lag)){
    z = matrix(rnorm(n.comp.sigma*ntot, 0,1), n.comp.sigma, ntot)
    V[[i]] = t(sigma.5 %*% z )
  }
  
  X = list() # the predictor curves
  for(j in 1:n.curves){
    X[[j]] = 0
    for(k in 0:Lag)
    {X[[j]]=X[[j]]+V[[j+k]]+10}
    X[[j]]=X[[j]]/sqrt(Lag+1)
  }
  
  library(goffda)
  
  E = r_ou(n=ntot, t = seq(0,1,len=lengthY), mu=0, alpha = 1, sigma = sigma2.e,
           x0=rchisq(n=ntot, df=1))$data
  
  tmp = lapply(1:length(true.main.index),  function(k){X[[true.main.index[k]]]%*%beta.vals[[true.main.index[k]]]/length(t.x)})
  Y.1 = Reduce("+",tmp)
  
  Y = Y.1+E
  t.x.all  =list()
  for(i in 1:n.curves){
    t.x.all[[i]]=t.x
  }
  
  X.train = list()
  X.test = list()
  for(j in 1:n.curves){
    X.train[[j]]=X[[j]][1:ntrain,]
    X.test[[j]]=X[[j]][-(1:ntrain),]
  }
  
  Y.train = Y[1:ntrain, ]
  Y.test = Y[-(1:ntrain), ]
  Y.f.train = Y.1[1:ntrain, ]
  Y.f.test = Y.1[-(1:ntrain), ]
  
  return(list("X_train" = X.train, "X_test" = X.test, "Y_train" = Y.train, "Y_test" = Y.test, "fun_Y_train" = Y.f.train, "fun_Y_test" = Y.f.test))
  
}
