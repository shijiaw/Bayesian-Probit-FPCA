rm(list=ls())
set.seed(1293)
library('mvtnorm')
library('fda')
library('orthogonalsplinebasis')
library('purrr')
library('truncnorm')
library("pracma")
library("LaplacesDemon")
library("xtable")

nrep <- 30
Klist <- c(2, 3, 4)
nncol <- length(Klist)
waic_Res <- matrix(NA, nr = nrep,nc = nncol)
epd_Res <- matrix(NA, nr = nrep,nc = nncol)
eaic_Res <- matrix(NA, nr = nrep,nc = nncol)
ebic_Res <- matrix(NA, nr = nrep,nc = nncol)
dic_Res <- matrix(NA, nr = nrep,nc = nncol)
N <- 300
t <- seq(0, 1, 0.05)
T <- length(t)
niter <- 50000
nburnin <- 0.6*niter

for(rep in 1:nrep){
  K <- 3
  #rep <- length(Klist)
  sigma <- 1
  d <- c(-Inf, 0.5, 1)
  tau <- cumsum(exp(d))-3.5
  
  #simulate curves
  curves <- matrix(NA, nr = N, nc = T)
  phi <- matrix(NA, nr = K, nc = T)
  xi <- matrix(NA, nr = N, nc = K)
  
  meancurve <- 3/(sqrt(2*pi*0.01))*exp(-(t-0.5)^2/(2*0.1^2))
  meancurve <- meancurve - mean(meancurve)
  
  # simulate age and gender
  X1 <- runif(N, 0, 1)
  X2 <- rdunif(N, b=0,a=1)
  X <- cbind(X1, X2)
  w <- c(0.6, -0.9)
  effects <- X%*%w
  
  lambda <- c(3, 0.6, 0.1)
  
  for(k in 1:K){
    phi[k,] <- sqrt(2)*sin(pi*k*t)
    xi[,k] <- rnorm(N ,0, lambda[k])
  }
  
  for(n in 1:N){
    curves[n,] <- effects[n] + meancurve + xi[n,]%*%phi + rnorm(T ,0, sigma)
  }
  
  y <- matrix(NA, nr = N, nc = T)
  y[which(curves < tau[1], arr.ind = T)] <- 1
  y[which((curves > tau[1])&(curves < tau[2]), arr.ind = T)] <- 2
  y[which((curves > tau[2])&(curves < tau[3]), arr.ind = T)] <- 3
  y[which(curves > tau[3], arr.ind = T)] <- 4
  #mu <- meancurve
  
  y_index <- list()
  for(j in 1:4){
    y_index[[j]] <- which(as.vector(y) == j)
  }
  
  #define orthonormal basis functions
  tobs = t
  norder   = 4
  knots = c(rep(0, norder-1), seq(0, 1, 0.1), rep(1, norder-1))
  nbasis   = length(knots) - norder
  
  obasis <- OrthogonalSplineBasis(knots, norder)
  obasis_eval <- evaluate(obasis, tobs)

  #define B-spline basis functions
  knots_B = seq(0, 1 ,0.1)
  nknots_B   = length(knots_B)
  nbasis_B   = length(knots_B) + norder - 2
  bsbasis = create.bspline.basis(range(knots_B),nbasis_B,norder,knots_B)
  basismat   = eval.basis(tobs, bsbasis)
  
  
  for(nsim in 1:nncol){
    K <- Klist[nsim]
    L_tau <- 5
    eps_tau <- 0.01
    print(nsim)

    xi_chain <- matrix(NA, nr = N*K, nc = niter)
    Theta_chain <- matrix(NA, nr = nbasis*K, nc = niter)
    Lambda_chain <- matrix(NA, nr = K, nc = niter)
    c_chain <- matrix(NA, nr = nbasis, nc = niter)
    Thetasample_chain <- matrix(NA, nr = nbasis*K, nc = niter)
    w_chain <- matrix(NA, nr = 2, nc = niter)
    tau_chain <- matrix(NA, nr = length(d), nc = niter)
    
    source('HMC_probit.R')
    source("mix_FPCA.R")
    
    source("WAIC.R")
    waic_index <- (niter - round(nburnin)):(niter-1)
    waic <- WAIC(y_index, X, xi_chain[ ,waic_index], Theta_chain[ ,waic_index], c_chain[ ,waic_index], w_chain[ ,waic_index], tau_chain[ ,waic_index+1], length(waic_index))
    epd <- EPD(y_index, X, xi_chain[ ,waic_index], Theta_chain[ ,waic_index], c_chain[ ,waic_index], w_chain[ ,waic_index], tau_chain[ ,waic_index+1], length(waic_index))
    eaic <- EAIC(y_index, X, xi_chain[ ,waic_index], Theta_chain[ ,waic_index], c_chain[ ,waic_index], w_chain[ ,waic_index], tau_chain[ ,waic_index+1], length(waic_index))
    ebic <- EBIC(y_index, X, xi_chain[ ,waic_index], Theta_chain[ ,waic_index], c_chain[ ,waic_index], w_chain[ ,waic_index], tau_chain[ ,waic_index+1], length(waic_index))
    dic <- DIC(y_index, X, xi_chain[ ,waic_index], Theta_chain[ ,waic_index], c_chain[ ,waic_index], w_chain[ ,waic_index], tau_chain[ ,waic_index+1], length(waic_index))
    waic_Res[rep, nsim] <- waic
    epd_Res[rep, nsim] <- epd
    eaic_Res[rep, nsim] <- eaic
    ebic_Res[rep, nsim] <- ebic
    dic_Res[rep, nsim] <- dic
  }
  
  
}

write.table(matrix(round(waic_Res, digits = 5), nr=1), file = "waic_Res.txt", row.names = FALSE, col.names = FALSE)
write.table(matrix(round(epd_Res, digits = 5), nr=1), file = "epd_Res.txt", row.names = FALSE, col.names = FALSE)
write.table(matrix(round(eaic_Res, digits = 5), nr=1), file = "eaic_Res.txt", row.names = FALSE, col.names = FALSE)
write.table(matrix(round(ebic_Res, digits = 5), nr=1), file = "ebic_Res.txt", row.names = FALSE, col.names = FALSE)
write.table(matrix(round(dic_Res, digits = 5), nr=1), file = "dic_Res.txt", row.names = FALSE, col.names = FALSE)


