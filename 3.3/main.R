rm(list=ls())
set.seed(029323)
library('mvtnorm')
library('fda')
library('orthogonalsplinebasis')
library('purrr')
library('truncnorm')
library("pracma")
library("LaplacesDemon")
library('ggplot2')
library('Rmisc') 
library('Sleuth2')
library('ggpubr')
library('stringr')
library('face')
library('refund')
#library('factoextra')

mode <- function(x){
  return(as.numeric(names(which.max(table(x)))))
}

Nrep <- 30
model1 <- rep(NA, Nrep)
model2 <- rep(NA, Nrep)
model3 <- rep(NA, Nrep)
model4 <- rep(NA, Nrep)

for(nrep in 1:Nrep){
  #simulate curves
  N <- 500
  t <- seq(0, 1, 0.05)
  T <- length(t)
  K <- 3
  sigma <- 1
  curves <- matrix(NA, nr = N, nc = T)
  phi <- matrix(NA, nr = K, nc = T)
  xi <- matrix(NA, nr = N, nc = K)
  
  meancurve <-  3/(sqrt(2*pi*0.01))*exp(-(t-0.5)^2/(2*0.1^2))
  meancurve <- meancurve - mean(meancurve)
  # simulate age and gender
  #X1 <- rdunif(N, b=1,a=100)
  X1 <- runif(N, 0, 1)
  X2 <- rdunif(N, b=0,a=1)
  X <- cbind(X1, X2)
  w <- c(0.6, -0.4)
  effects <- X%*%w
  
  #lambda <- 3/(2^((1:K)-1))
  p_label = c(0.27, 0.38, 0.35)
  # constraints:\sum pi*mu = 0
  dim <- 3
  mu <- matrix(0, nr = K, nc = dim)
  mu[1,1] <- 4
  mu[1,2] <- -2
  mu[,3] <- -(mu[,1]*p_label[1] + mu[,2]*p_label[2])/p_label[3]
  # 
  # sig_lambda <- matrix(c(rep(0.5, 3), rep(0.2, 3), rep(0.1, 3)), nr = 3)
  # sig_lambda[3, ] = c(0.8, 0.7, 0.5)

  sig_lambda <- matrix(0, nr = dim, nc = K)
  sig_lambda[, 1] = c(0.8, 0.2, 0.6)
  sig_lambda[, 2] = c(1.2, 1.5, 0)
  sig_lambda[, 3] = c(0.1, 0, 0.8)

  
  Z <- sample(1:dim, size = N, replace = TRUE, prob = p_label)
  #xi <- matrix(0, nr = K, nc = N)
  for(i in 1:N){
    xi[i, ] <- rnorm(K, mean = mu[,Z[i]], sd = sig_lambda[,Z[i]])
  }
  
  lambda <- 3/(2^((1:K)-1))
  for(k in 1:K){
    phi[k,] <- sqrt(2)*sin(pi*k*t)
    lambda[k] <- sqrt(sum(p_label*(sig_lambda[k,]^2 + mu[k,]^2)))
    #xi[,k] <- rnorm(N ,0, lambda[k])
  }
 
   
  for(n in 1:N){
    curves[n,] <- effects[n] + meancurve + xi[n,]%*%phi + rnorm(T ,0, sigma)
  }
  
  d <- c(-Inf, 0.5, 1)
  tau <- cumsum(exp(d))-3.5
  
  
  y <- matrix(NA, nr = N, nc = T)
  y[which(curves < tau[1], arr.ind = T)] <- 1
  y[which((curves > tau[1])&(curves < tau[2]), arr.ind = T)] <- 2
  y[which((curves > tau[2])&(curves < tau[3]), arr.ind = T)] <- 3
  y[which(curves > tau[3], arr.ind = T)] <- 4
  
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
  #deriv1 <- deriv(obasis)
  #obasis_deriv <- evaluate(deriv1, tobs)
  
  #define B-spline basis functions
  knots_B = seq(0, 1, 0.1)
  nknots_B   = length(knots_B)
  nbasis_B   = length(knots_B) + norder - 2
  bsbasis = create.bspline.basis(range(knots_B),nbasis_B,norder,knots_B)
  basismat   = eval.basis(tobs, bsbasis)
  
  
  beta_est <- matrix(NA, nr = N, nc = nbasis)
  xi_est <- matrix(NA, nr = N, nc = K)
  Theta_est <- matrix(NA, nr = nbasis, nc = K)
  Lambda_est <- 2/2.5^(0:(K-1))
  meancurve_est <- apply(curves, 2, mean)
  meancurve_est <- meancurve_est - mean(meancurve_est)
  
  xi_est <- xi  + matrix(rnorm(N*K, 0, 0.1), nr = N, nc = K)
  X_temp <- kronecker(xi_est, obasis_eval)
  Theta_est <- matrix(solve(t(X_temp)%*%X_temp + diag(10^(-8),nbasis*K))%*%t(X_temp)%*%as.vector(t(curves-matrix(rep(meancurve_est, N), nr = N, byrow = TRUE))), nr = nbasis)
  phi_est <- t(obasis_eval%*%Theta_est)
  beta_est <- xi_est%*%t(Theta_est)
  c_est <- solve(t(basismat)%*%basismat)%*%t(basismat)%*%(meancurve_est)
  w_est <- solve(t(X)%*%X)%*%t(X)%*%apply((curves - matrix(rep(basismat%*%c_est, N), nr = N, byrow = TRUE)), 1, mean)
  tau_est <- c(-3, 0, -1)
  
  # MCMC for FPCA 
  niter <- 20000
  nburnin <- 0.5*niter
  source('HMC_probit.R')
  source('Gibbs.R')
  source('mix_FPCA.R')
  
  
  ao <- 0.1
  bo <- 0.1
  muo <- 0
  sigmao <- 3
  aalpha <- 0.1 
  balpha <- 0.1
  ZResult1_est <- Z
  ZResult2_est <- Z
  muResult_est <- mu
  sigmaResult_est <- sigmaGibbs(ZResult1_est, muResult_est, ao, bo, t(xi_est))
  kappaResult_est <- 0.5 
  nn <- tabulate(ZResult1_est)
  pResult_est <-  pGibbs(kappaResult_est, nn)
  uResult_est <- uslice(pResult_est, Z)
  mu_temp <- rnorm(K, muo, sigmao)
  sigma_temp <- rep(3, K)
  
  L_tau <- 4
  eps_tau <- 0.008
  
  sig2e <- 1
  
  MCMCResults <- MCMCmove(ZResult1_est, ZResult2_est, muResult_est, kappaResult_est, sigmaResult_est, pResult_est, uResult_est, muo, sigmao, aalpha, balpha, ao, bo, xi_est,tau_est, Theta_est, Lambda_est, c_est, w_est)
  xi_chain <- MCMCResults$xi_chain 
  Theta_chain <- MCMCResults$Theta_chain
  Lambda_chain <- MCMCResults$Lambda_chain 
  c_chain <- MCMCResults$c_chain 
  w_chain <- MCMCResults$w_chain 
  tau_chain <- MCMCResults$tau_chain 
  sigmaResult_chain <- MCMCResults$sigmaResult_chain 
  muResult_chain <- MCMCResults$muResult_chain 
  pResult_chain <- MCMCResults$pResult_chain 
  uResult_chain <- MCMCResults$uResult_chain 
  ZResult_chain <- MCMCResults$ZResult_chain 
  KResult_chain <- MCMCResults$KResult_chain 
  kappaResult_chain <- MCMCResults$kappaResult_chain 
  ncluster_chain <- MCMCResults$ncluster_chain
  
  
  FinalK <- as.numeric(names(which.max(table(KResult_chain))))
  K_index <- which(KResult_chain == FinalK)
  final_mu <- matrix(NA, nr = K*FinalK, nc = length(K_index))
  final_sigma <- matrix(NA, nr = K*FinalK, nc = length(K_index))
  final_Z <- matrix(NA, nr = N, nc = length(K_index))
  final_p <- matrix(NA, nr = FinalK, nc = length(K_index))
  for(ind in 1:length(K_index)){
    muResult_est <- muResult_chain[[K_index[ind]]]
    order_mu <- order(apply(muResult_est, 2, function(x)(sum(abs(x)))), decreasing = TRUE)
    ZResult_temp1 <- ZResult_chain[[K_index[ind]]]
    ZResult_temp <- ZResult_chain[[K_index[ind]]]
    for(i in 1:length(order_mu)){
      ZResult_temp[which(ZResult_temp1 == i)] <- which(order_mu == i)
    }
    
    final_mu[ ,ind] <- matrix(muResult_est[,order_mu], nc = 1)
    final_sigma[ ,ind] <- matrix(sigmaResult_chain[[K_index[ind]]][,order_mu], nc = 1)
    final_Z[ ,ind] <- ZResult_temp
    final_p[ ,ind] <- pResult_chain[[K_index[ind]]][order_mu]
    
  }
  
  # Lambda_mean <- apply(Lambda_chain[,-(1:nburnin)], 1, mean)
  # Theta_mean <- matrix(apply(Theta_chain[ ,-(1:nburnin)], 1, mean), nc = K)
  # est_curve <- basismat%*%apply(c_chain[,-(1:nburnin)], 1, mean)
  # meancurve_RMSE <- sqrt(mean(((meancurve - est_curve)^2)))
  # phi_RMSE <- rep(NA, K)
  # for(k in 1:K){
  #   phi_RMSE[k] <- sqrt(mean(((phi[k,] - phi_est[k,])^2)))
  # }
  # 
  # gname_theta = paste("theta_chain", N,".eps",sep="")  
  # gname_xi = paste("xi_chain", N,".eps",sep="")  
  # gname_Lambda = paste("Lambda_chain", N,".eps",sep="") 
  # gname_phi = paste("phi_est", N,".eps",sep="") 
  # gname_mean = paste("mean_est", N,".eps",sep="")  
  # gname_tau = paste("tau_chain", N,".eps",sep="")  
  # 
  # tau_mean <- apply(tau_chain[,-(1:nburnin)], 1, mean)
  # tau_025 <- apply(tau_chain[,-(1:nburnin)], 1, function(x)(quantile(x, 0.025)))
  # tau_975 <- apply(tau_chain[,-(1:nburnin)], 1, function(x)(quantile(x, 0.975)))
  # 
  # Lambda_mean <- apply(sqrt(Lambda_chain[,-(1:nburnin)]), 1, mean)
  # Lambda_025 <- apply(sqrt(Lambda_chain[,-(1:nburnin)]), 1, function(x)(quantile(x, 0.025)))
  # Lambda_975 <- apply(sqrt(Lambda_chain[,-(1:nburnin)]), 1, function(x)(quantile(x, 0.975)))
  
  #gname_mu <- c("mu.eps",sep="")  
  #gname_sigma <- c("sigma.eps",sep="")  
  #gname_p <- c("p.eps",sep="")  
  estimate_label <- apply(final_Z, 1, function(x)(mode(x)))
  if(max(estimate_label)!= K){
    model1[nrep] <- NA
    model2[nrep] <- NA
    model3[nrep] <- NA
    model4[nrep] <- NA
  }else{
    model1[nrep] <- sum(estimate_label == Z)/N
    
    # method1, fpca + K-mean
    fd1.1 <- Data2fd(t(y), basisobj=bsbasis, lambda = 1)
    fpca <- pca.fd(fd1.1, nharm = 10, harmfdPar=fdPar(fd1.1), centerfns = TRUE)
    fpcavar <- fpca$varprop
    nfpca <- min(which(cumsum(fpcavar) > 0.99))
    score <- fpca$scores[ ,1:nfpca]
    
    km <- kmeans(score, centers = 3, nstart = 25)
    km$cluster
    km$centers
    km$withinss
    apply(km$centers, 1, function(x)(sum(x^2)))
    order_km <- order(apply(km$centers, 1, function(x)(sum(abs(x)))), decreasing = TRUE)
    #order_km <- order(km$withinss, decreasing = TRUE) 
    
    km_result <- rep(NA, N)
    for(i in 1:length(order_km)){
      km_result[which(km$cluster == i)] <- which(order_km == i)
    }
    
    model2[nrep] <- sum(km_result == Z)/N
    
    # FACE
    id <- rep(1:N,each = T)
    tt <- rep(1:T,times = N)
    res <- as.vector(t(y))
    #sel <- which(is.na(y))
    
    
    ## organize data and apply FACEs
    data <- data.frame(y=(res),
                       argvals = tt,
                       subj = id)
    #data <- data[data$y>4.5,]
    #fit_face <- face.sparse(data,argvals.new=1:T)
    
    ## set calculate.scores to TRUE if want to get scores
    fit_face <- face.sparse(data,argvals.new=1:T,calculate.scores=TRUE)
    face_scores <- fit_face$rand_eff$scores
    
    km_face <- kmeans(face_scores, centers = 3, nstart = 25)
    km_face$cluster
    km_face$centers
    order_face <- order(apply(km_face$centers, 1, function(x)(sum(abs(x)))), decreasing = TRUE)
    #order_km <- order(km$withinss, decreasing = TRUE) 
    
    face_result <- rep(NA, N)
    for(i in 1:length(order_km)){
      face_result[which(km_face$cluster == i)] <- which(order_face == i)
    }
    
    model3[nrep] <- sum(face_result == Z)/N
    
    # SC
    Fit.MM = fpca.sc(y, var = TRUE, simul = TRUE)
    sc_score <- Fit.MM$scores
    km_sc <- kmeans(sc_score, centers = 3, nstart = 25)
    order_sc <- order(apply(km_sc$centers, 1, function(x)(sum(abs(x)))), decreasing = TRUE)
    
    sc_result <- rep(NA, N)
    for(i in 1:length(order_km)){
      sc_result[which(km_sc$cluster == i)] <- which(order_sc == i)
    }
    
    model4[nrep] <- sum(sc_result == Z)/N
    
    
  }
}

write.table(model1, file = "model1.txt", row.names = FALSE, col.names = FALSE)
write.table(model2, file = "model2.txt", row.names = FALSE, col.names = FALSE)
write.table(model3, file = "model3.txt", row.names = FALSE, col.names = FALSE)
write.table(model4, file = "model4.txt", row.names = FALSE, col.names = FALSE)


