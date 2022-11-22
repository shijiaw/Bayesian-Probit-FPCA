OneFPCAMCMCmove <- function(tau_est, xi_est, Theta_est, Lambda_est, c_est, w_est, N, T){
  meancurve_est <- basismat%*%c_est
  beta_est <- xi_est%*%t(Theta_est)
  sig2e <- 1
  
  XXX <- as.vector(beta_est%*%t(obasis_eval) + matrix(rep(meancurve_est, N), nr = N, byrow = TRUE) + matrix(rep(X%*%w_est, T), nr = N))
  tau_est <- HMC_tau(gradient_tau, energyH_tau, eps_tau, L_tau, tau_est, y_index, XXX)
  cutoff <- c(-Inf, cumsum(exp(c(-Inf, tau_est[-1])))+tau_est[1], Inf)
  inv_l <- cutoff[y]
  inv_u <- cutoff[y+1]  
  y_star <- rtruncnorm(length(XXX), a= inv_l, b=inv_u, mean = XXX, sd = 1)
  curves <- matrix(y_star, nr = N)
  
  sdobs <- curves - matrix(rep(meancurve_est, N), nr = N, byrow = TRUE) - matrix(rep(X%*%w_est, T), nc = T)
  
  # sample Theta
  Sig_hat <- rinvwishart(N + T + 1, t(sdobs)%*%sdobs + diag(rep(sig2e, T)) ) 
  Sig_temp <- solve(t(obasis_eval)%*%obasis_eval)%*%t(obasis_eval)
  mat <- Sig_temp%*%(Sig_hat - diag(sig2e, T))%*%t(Sig_temp)
  
  eigenmat <- eigen((mat + t(mat))/2)
  Theta_est <- eigenmat$vectors[,1:K]
  sign <- as.numeric((t(obasis_eval%*%Theta_est))[,4]>0)
  signindex <- which(sign == 0)
  if(length(signindex)>0){
    Theta_est[,signindex] <- -Theta_est[,signindex]
  }
  Lambda_est <- eigenmat$values[1:K]
  
  
  phi_est <- t(obasis_eval%*%Theta_est)
  
  # sample c
  Sigma_tildeinv <- solve(diag(sig2e, T) + t(phi_est)%*%diag(eigenmat$values[1:K])%*%phi_est)
  c_Sigma <- solve(N*t(basismat)%*%Sigma_tildeinv%*%basismat)
  cholcSig <- chol(c_Sigma)
  c_mean <-  c_Sigma%*%t(basismat)%*%Sigma_tildeinv%*%apply(curves, 2, sum)
  c_est <- cholcSig%*%rnorm(nbasis, 0, 1) + c_mean
  c_est <- c_est - mean(basismat%*%c_est)
  
  # sample w 
  meancurve_est <- basismat%*%c_est
  sdobs <- curves - matrix(rep(meancurve_est, N), nr = N, byrow = TRUE) - matrix(rep(X%*%w_est, T), nc = T)
  
  
  # sample xi
  Sigmaxi <- solve(phi_est%*%t(phi_est)/sig2e+diag(1/Lambda_est))
  cholSig <- chol(Sigmaxi)
  fix_mean <-  Sigmaxi%*%phi_est/sig2e
  for(n in 1:N){
    xi_est[n,] <- cholSig%*%rnorm(K, 0, 1) + fix_mean%*%sdobs[n,]
  }
  
  
  # sample w 
  w_Sigma <- solve(T*t(X)%*%X/sig2e + diag(c(1, 1)))
  cholwSig <- chol(w_Sigma)
  w_mean <- w_Sigma%*%(t(X)%*%apply((curves - matrix(rep(meancurve_est, N), nr = N, byrow = TRUE)- xi_est%*%t(Theta_est)%*%t(obasis_eval)), 1, sum)/sig2e)
  w_est <-  cholwSig%*%rnorm(2, 0, 1) + w_mean
 
  return(list(tau_est = tau_est, xi_est = xi_est, Theta_est = Theta_est, Lambda_est = Lambda_est, c_est = c_est, w_est = w_est))
}


gradient_tau <- function(y_index, tau, X){
  K <- length(tau)
  rt <- - diag(1/rep(100, K))%*%tau
  #n <- table(y)
  cut <- cumsum(exp(c(-Inf, tau[-1])))+tau[1]
  cutoff_u <- list()
  cutoff_l <- list()
  pnormcutoff_u <- list()
  pnormcutoff_l <- list()
  dnormcutoff_u <- list()
  dnormcutoff_l <- list()
  for(j in 1:K){
    cutoff_l[[j]] <- cut[j] - X[y_index[[j]]]
    cutoff_u[[j]] <- cut[j] - X[y_index[[j+1]]]
    pnormcutoff_u[[j]] <- pnorm(cutoff_u[[j]])
    pnormcutoff_l[[j]] <- pnorm(cutoff_l[[j]])
    dnormcutoff_u[[j]] <- dnorm(cutoff_u[[j]])
    dnormcutoff_l[[j]] <- dnorm(cutoff_l[[j]])
  }
  rt[K] <- exp(tau[K])*(sum(((dnormcutoff_l[[K]])/(pnormcutoff_l[[K]]-pnormcutoff_u[[K-1]])))-sum(dnormcutoff_u[[K]]/(1-pnormcutoff_u[[K]])))
 
  rt[1] <- sum(dnormcutoff_l[[1]]/pnormcutoff_l[[1]])-sum(dnormcutoff_u[[K]]/(1-pnormcutoff_u[[K]]))#+sum(((dnormcutoff_l[[2]]-dnormcutoff_u[[1]])/(pnormcutoff_l[[2]]-pnormcutoff_u[[1]])))+sum(((dnorm(cutoff_l[[3]])-dnorm(cutoff_u[[2]]))/(pnorm(cutoff_l[[3]])-pnorm(cutoff_u[[2]]))))
  
  for(i in 2:K){
    rt[1] <- rt[1] + sum(((dnormcutoff_l[[i]]-dnormcutoff_u[[i-1]])/(pnormcutoff_l[[i]]-pnormcutoff_u[[i-1]])))
   
  }
  for(j in 2:(K-1)){
    rt[j] <- sum(((dnormcutoff_l[[j]])/(pnormcutoff_l[[j]]-pnormcutoff_u[[j-1]]))) - sum(dnormcutoff_u[[K]]/(1-pnormcutoff_u[[K]]))
    
    for(k in (j+1):K){
      rt[j] <- rt[j] +sum(((dnormcutoff_l[[k]]-dnormcutoff_u[[k-1]])/(pnormcutoff_l[[k]] - pnormcutoff_u[[k-1]])))
    }
    rt[j] <- exp(tau[j])*rt[j]
  }
  return(rt)
}


energyH_tau <- function(y_index, tau, p, X){
  K <- length(tau)
  cutoff <- cumsum(exp(c(-Inf, tau[-1])))+tau[1]
  rt <- t(tau)%*%diag(1/rep(100, K))%*%tau/2 + 1/2*t(p)%*%p - sum(log(pnorm(cutoff[1] - X[y_index[[1]]])))
  for(j in 2:K){
    rt <- rt - sum(log(pnorm(cutoff[j]- X[y_index[[j]]]) - pnorm(cutoff[j-1]- X[y_index[[j]]])))
  }
  rt <- rt - sum(log(1 - pnorm(cutoff[K] - X[y_index[[K+1]]])))
  return(rt)
}

HMC_tau <- function(gradient, energyH, eps, L, current_tau, y, X){
  x <- current_tau
  p <- rnorm(length(x), 0, 1)
  current_p <- p
  p <- p + eps/2*gradient(y, x, X)
  
  for(i in 1:L){
    x <- x + eps*p
    if (i!=L) p <- p + eps*gradient(y, x, X)
  }
  
  p = p + eps*gradient(y, x, X)/2
  
  RH <- -energyH(y, x, p, X) + energyH(y, current_tau, current_p, X)
  if(!is.nan(RH)){
    if(log(runif(1)) <  RH){
      return(x)
    }else{
      return(current_tau)
    }
  }else{
    return(current_tau)
  }
  
}

