This repo contains an example R code used in the Monte Carlo experiments of the paper "Function-on-function linear quantile regression"
# Authors
Ufuk Beyaztas amd Han Lin Shang
# Description
auxiliary.R file contains all the necessary functions to perform the FFLQR method \\
dgp1.R file is used to generate data under Gaussian errors\\
dgp2.R file is used to generate data under chi-square(1) errors\\
run.R file contains a toy example
# Packages
library(fda) 
library(quantreg)
library(matrixStats)
library(goffda)
