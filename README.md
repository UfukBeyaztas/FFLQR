This repo contains an example R code used in the Monte Carlo experiments of the paper "Function-on-function linear quantile regression"
# Authors
Ufuk Beyaztas amd Han Lin Shang
# Description
1) auxiliary.R file contains all the necessary functions to perform the FFLQR method
2) dgp1.R file is used to generate data under Gaussian errors
3) dgp2.R file is used to generate data under chi-square(1) errors
4) run.R file contains a toy example
# Packages
library(fda) 
library(quantreg)
library(matrixStats)
library(goffda)
