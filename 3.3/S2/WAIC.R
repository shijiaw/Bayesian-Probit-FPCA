log_like <- function(y_index, X, xi_est, Theta_est, c_est, w_est, tau_est){
  meancurve_est <- basismat%*%c_est
  beta_est <- xi_est%*%t(Theta_est)
  sig2e <- 1
  
  XXX <- as.vector(beta_est%*%t(obasis_eval) + matrix(rep(meancurve_est, N), nr = N, byrow = TRUE) + matrix(rep(X%*%w_est, T), nr = N))
  
  K <- length(tau_est)
  cutoff <- cumsum(exp(c(-Inf, tau_est[-1])))+tau_est[1]
  rt <-   log(pnorm(cutoff[1] - XXX[y_index[[1]]]))
  for(j in 2:K){
    rt <- c(rt, log(pnorm(cutoff[j]- XXX[y_index[[j]]]) - pnorm(cutoff[j-1]- XXX[y_index[[j]]])))
  }
  rt <- c(rt, log(1 - pnorm(cutoff[K] - XXX[y_index[[K+1]]])))
  sindex <- which(rt < - 180)
  rt[sindex] <- -180
  return(rt)
}

logSum <- function(logx){
  logxmax <- max(logx)
  rt <- log(sum(exp(logx - logxmax))) + logxmax
  return(rt)
}

WAIC <- function(y_index, X, xi_chain, Theta_chain, c_chain, w_chain, tau_chain, nsamp){
  logp <- matrix(NA, N*T, nsamp)
  #logp <- matrix(NA, N, nsamp)
  for(k in 1:nsamp){
    logp[ ,k] <- log_like(y_index, X, matrix(xi_chain[,k], nc = K), matrix(Theta_chain[,k], nc = K), c_chain[,k], w_chain[,k], tau_chain[,k])
  }
  Vlogp <- apply(logp, 1, var)
  print(sum(Vlogp))
  rt <- sum(Vlogp) + N*T*log(nsamp)
  #rt <- sum(Vlogp) + N*log(nsamp)
  for(i in 1:(N*T)){
    rt <- rt - logSum(logp[i,])
  }
  return(2*rt)
}

EPD <- function(y_index, X, xi_chain, Theta_chain, c_chain, w_chain, tau_chain, nsamp){
  logp <- matrix(NA, (N*T), nsamp)
  for(k in 1:nsamp){
    logp[ ,k] <- log_like(y_index, X, matrix(xi_chain[,k], nc = K), matrix(Theta_chain[,k], nc = K), c_chain[,k], w_chain[,k], tau_chain[,k])
  }
  #Vlogp <- apply(logp, 1, var)
  rt <- (N*T)*log(nsamp)
  for(i in 1:(N*T)){
    rt <- rt - logSum(logp[i,])
  }
  return(2*rt)
}

EAIC <- function(y_index, X, xi_chain, Theta_chain, c_chain, w_chain, tau_chain, nsamp){
  logp <- matrix(NA, (N*T), nsamp)
  for(k in 1:nsamp){
    logp[ ,k] <- log_like(y_index, X, matrix(xi_chain[,k], nc = K), matrix(Theta_chain[,k], nc = K), c_chain[,k], w_chain[,k], tau_chain[,k])
  }
  # double check xi_chain
  rt <- T*N*log(nsamp) + (nrow(tau_chain) + K + nrow(Theta_chain) + nrow(c_chain) + nrow(w_chain))
  for(i in 1:(N*T)){
    rt <- rt - logSum(logp[i,])
  }
  return(2*rt)
}

EBIC <- function(y_index, X, xi_chain, Theta_chain, c_chain, w_chain, tau_chain, nsamp){
  logp <- matrix(NA, (N*T), nsamp)
  for(k in 1:nsamp){
    logp[ ,k] <- log_like(y_index, X, matrix(xi_chain[,k], nc = K), matrix(Theta_chain[,k], nc = K), c_chain[,k], w_chain[,k], tau_chain[,k])
  }
  # double check xi_chain
  rt <- T*N*log(nsamp) + 0.5*log(N*T)*(nrow(tau_chain) + K + nrow(Theta_chain) + nrow(c_chain) + nrow(w_chain))
  for(i in 1:(N*T)){
    rt <- rt - logSum(logp[i,])
  }
  return(2*rt)
}


DIC <- function(y_index, X, xi_chain, Theta_chain, c_chain, w_chain, tau_chain, nsamp){
  logp <- matrix(NA, N*T, nsamp)
  for(k in 1:nsamp){
    logp[ ,k] <- log_like(y_index, X, matrix(xi_chain[,k], nc = K), matrix(Theta_chain[,k], nc = K), c_chain[,k], w_chain[,k], tau_chain[,k])
  }
  Vlogp <- var(apply(logp, 2, sum))
  xi_est <- matrix(apply(xi_chain, 1, mean), nc = K)
  Theta_est <- matrix(apply(Theta_chain, 1, mean), nc = K)
  c_est <- apply(c_chain, 1, mean)
  w_est <- apply(w_chain, 1, mean)
  tau_est <- apply(tau_chain, 1, mean)
  
  rt <- -2*sum(log_like(y_index, X, xi_est, Theta_est, c_est, w_est, tau_est)) + 4*Vlogp
  return(rt)
}




