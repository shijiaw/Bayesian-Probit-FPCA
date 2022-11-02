rm(list=ls())
set.seed(029323)
library('mvtnorm')
library('fda')
library('orthogonalsplinebasis')
library('purrr')
library('truncnorm')
library("pracma")
library("LaplacesDemon")
library('xtable')
#library('ggplot2')
#library('Rmisc') 
#library('Sleuth2')
#library('ggpubr')
#library('stringr')

N <- 500
nncol <- 30
t <- seq(0, 1, 0.05)
T <- length(t)
K <- 3
sigma <- 1
d <- c(-Inf, 0.5, 1)
tau <- cumsum(exp(d))-3.5
L_tau <- 5
eps_tau <- 0.01
niter <- 50000
source('HMC_probit.R')
nburnin <- 0.5*niter

tau_meanRes <- matrix(NA, nr = nncol, nc = K)
tau_025Res <- matrix(NA, nr = nncol, nc = K)
tau_975Res <- matrix(NA, nr = nncol, nc = K)

Lambda_meanRes <- matrix(NA, nr = nncol, nc = K)
Lambda_025Res <- matrix(NA, nr = nncol, nc = K)
Lambda_975Res <- matrix(NA, nr = nncol, nc = K)

meancurve_RMSERes <- rep(NA, nr = nncol)
phi_RMSERes <- matrix(NA, nr = nncol, nc = K)

for(nsim in 1:nncol){
  print(nsim)
  source("mix_FPCA.R")
  
  Lambda_mean <- apply(Lambda_chain[,-(1:nburnin)], 1, mean)
  Theta_mean <- matrix(apply(Theta_chain[ ,-(1:nburnin)], 1, mean), nc = K)
  est_curve <- basismat%*%apply(c_chain[,-(1:nburnin)], 1, mean)
  phi_est <- t(obasis_eval%*%matrix(apply(Theta_chain[,-(1:nburnin)], 1, mean), nc = K))
  meancurve_RMSE <- sqrt(mean(((meancurve - est_curve)^2)))
  phi_RMSE <- rep(NA, K)
  for(k in 1:K){
    phi_RMSE[k] <- sqrt(mean(((phi[k,] - phi_est[k,])^2)))
  }
  
  gname_theta = paste("Figures/theta_chain", nsim,".eps",sep="")  
  gname_xi = paste("Figures/xi_chain", nsim,".eps",sep="")  
  gname_Lambda = paste("Figures/Lambda_chain", nsim,".eps",sep="") 
  gname_phi = paste("Figures/phi_est", nsim,".eps",sep="") 
  gname_mean = paste("Figures/mean_est", nsim,".eps",sep="")  
  gname_tau = paste("Figures/tau_chain", nsim,".eps",sep="")  
  
  tau_mean <- apply(tau_chain[,-(1:nburnin)], 1, mean)
  tau_025 <- apply(tau_chain[,-(1:nburnin)], 1, function(x)(quantile(x, 0.025)))
  tau_975 <- apply(tau_chain[,-(1:nburnin)], 1, function(x)(quantile(x, 0.975)))
  
  Lambda_mean <- apply(sqrt(Lambda_chain[,-(1:nburnin)]), 1, mean)
  Lambda_025 <- apply(sqrt(Lambda_chain[,-(1:nburnin)]), 1, function(x)(quantile(x, 0.025)))
  Lambda_975 <- apply(sqrt(Lambda_chain[,-(1:nburnin)]), 1, function(x)(quantile(x, 0.975)))
  
  tau_meanRes[nsim,] <- tau_mean
  tau_025Res[nsim,] <- tau_025
  tau_975Res[nsim,] <- tau_975
  
  Lambda_meanRes[nsim,] <- Lambda_mean
  Lambda_025Res[nsim,] <-  Lambda_025
  Lambda_975Res[nsim,] <-  Lambda_975
  
  meancurve_RMSERes[nsim] <- meancurve_RMSE
  phi_RMSERes[nsim,] <- phi_RMSE
  
  source("FPCA_plot.R")
}

write.table(matrix(round(tau_meanRes, digits = 5),nr=1), file = "tau_mean.txt", row.names = FALSE, col.names = FALSE)
write.table(matrix(round(tau_025Res, digits = 5),nr=1), file = "tau_025.txt", row.names = FALSE, col.names = FALSE)
write.table(matrix(round(tau_975Res, digits = 5),nr=1), file = "tau_975.txt", row.names = FALSE, col.names = FALSE)

write.table(matrix(round(Lambda_meanRes, digits = 5),nr=1), file = "Lambda_mean.txt", row.names = FALSE, col.names = FALSE)
write.table(matrix(round(Lambda_025Res, digits = 5),nr=1), file = "Lambda_025.txt", row.names = FALSE, col.names = FALSE)
write.table(matrix(round(Lambda_975Res, digits = 5),nr=1), file = "Lambda_975.txt", row.names = FALSE, col.names = FALSE)

write.table(round(meancurve_RMSERes, digits = 5), file = "meancurve_RMSERes.txt", row.names = FALSE, col.names = FALSE)
write.table(matrix(round(phi_RMSERes, digits = 5), nr=1), file = "phi_RMSERes.txt", row.names = FALSE, col.names = FALSE)

