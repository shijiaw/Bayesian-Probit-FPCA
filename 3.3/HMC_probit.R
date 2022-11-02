
OneFPCAMCMCmove <- function(tau_est, xi_est, Theta_est, Lambda_est, c_est, w_est, Z, mu, sig_lambda){
  meancurve_est <- basismat%*%c_est
  beta_est <- xi_est%*%t(Theta_est)
  XXX <- as.vector(beta_est%*%t(obasis_eval) + matrix(rep(meancurve_est, N), nr = N, byrow = TRUE) + matrix(rep(X%*%w_est, T), nr = N))
  tau_est <- HMC_tau(gradient_tau, energyH_tau, eps_tau, L_tau, tau_est, y_index, XXX)
  cutoff <- c(-Inf, cumsum(exp(c(-Inf, tau_est[-1])))+tau_est[1], Inf)
  inv_l <- cutoff[y]
  inv_u <- cutoff[y+1]  
  y_star <- rtruncnorm(length(XXX), a= inv_l, b=inv_u, mean = XXX, sd = 1)
  curves <- matrix(y_star, nr = N)
  
  sdobs <- curves - matrix(rep(meancurve_est, N), nr = N, byrow = TRUE) - matrix(rep(X%*%w_est, T), nc = T)
  
  # sample Theta
  Sig_hat <- rinvwishart(N + T + 1, t(sdobs)%*%sdobs + diag(rep(0.1, T)) ) 
  Sig_temp <- solve(t(obasis_eval)%*%obasis_eval)%*%t(obasis_eval)
  mat <- Sig_temp%*%(Sig_hat - diag(sig2e, T))%*%t(Sig_temp)
  
  eigenmat <- eigen((mat + t(mat))/2)
  Theta_est <- eigenmat$vectors[,1:K]
  Lambda_est <- eigenmat$values[1:K]
  
  sign <- as.numeric((t(obasis_eval%*%Theta_est))[,4]>0)
  signindex <- which(sign == 0)
  if(length(signindex)>0){
    Theta_est[,signindex] <- -Theta_est[,signindex]
  }
  
  phi_est <- t(obasis_eval%*%Theta_est)
  
  # sample c
  Sigma_tildeinv <- solve(diag(sig2e, T) + t(phi_est)%*%diag(Lambda_est)%*%phi_est)
  c_Sigma <- solve(N*t(basismat)%*%Sigma_tildeinv%*%basismat)
  cholcSig <- chol(c_Sigma)
  c_mean <-  c_Sigma%*%t(basismat)%*%Sigma_tildeinv%*%apply(curves, 2, sum)
  c_est <- cholcSig%*%rnorm(nbasis, 0, 1) + c_mean
  c_est <- c_est - mean(basismat%*%c_est)
  
  meancurve_est <- basismat%*%c_est
  sdobs <- curves - matrix(rep(meancurve_est, N), nr = N, byrow = TRUE) - matrix(rep(X%*%w_est, T), nc = T)
  
  # sample xi
  ncluster <- length(table(Z))
  Sigmaxi <- list()
  cholSig <- list()
  for(nn in 1:ncluster){
    Sigmaxi[[nn]] <- solve(phi_est%*%t(phi_est)/sig2e+diag(1/sig_lambda[,nn]^2))
    cholSig[[nn]] <- chol(Sigmaxi[[nn]])
    #fix_mean <-  Sigmaxi%*%phi_est/sig2e
  }
  #Sigmaxi <- solve(phi_est%*%t(phi_est)/sig2e+diag(1/Lambda_est))
  #cholSig <- chol(Sigmaxi)
  #fix_mean <-  Sigmaxi%*%phi_est/sig2e
  for(n in 1:N){
    #xi_est[n,] <- cholSig[[Z[n]]]%*%rnorm(K, 0, 1) + fix_mean%*%sdobs[n,]
    xi_est[n,] <- cholSig[[Z[n]]]%*%rnorm(K, 0, 1) + Sigmaxi[[Z[n]]]%*%(phi_est%*%sdobs[n,]/sig2e + diag(1/sig_lambda[,Z[n]]^2)%*%mu[,Z[n]])
    #xi_est[n,] <- cholSig%*%rnorm(K, 0, 1) + Sigmaxi%*%(phi_est%*%sdobs[n,]/sig2e + diag(1/sig_lambda[,Z[n]]^2)%*%mu[,Z[n]])
  }
  # debug 
  
  
  
  # sample Lambda
  for(k in 1:K){
    Lambda_est[k] <- 1/rgamma(1, shape = 1+N/2, rate = 0.1+sum(xi_est[,k]^2)/2)
  }
  #print(Lambda_est)
  
  # sample w 
  w_Sigma <- solve(T*t(X)%*%X/sig2e + diag(c(1, 1)))
  cholwSig <- chol(w_Sigma)
  w_mean <- w_Sigma%*%(t(X)%*%apply((curves - matrix(rep(meancurve_est, N), nr = N, byrow = TRUE)- xi_est%*%t(Theta_est)%*%t(obasis_eval)), 1, sum)/sig2e)
  w_est <-  cholwSig%*%rnorm(2, 0, 1) + w_mean
  
  return(list(tau_est = tau_est, xi_est = xi_est, Theta_est = Theta_est, Lambda_est = Lambda_est, c_est = c_est, w_est = w_est, sig2e = sig2e))
}

OneMixMCMCmove <- function(ZResult, muResult, kappaResult, sigmaResult, pResult, uResult, muo, sigmao, aalpha, balpha, ao, bo, y){
  nn <- tabulate(ZResult)
  if(min(nn) == 0){
    nn_index <- which(nn == 0)
    ZResult[which(ZResult > nn_index)] <- ZResult[which(ZResult > nn_index)] - 1
    muResult <- muResult[,-nn_index]
    nn <- tabulate(ZResult)
  }
  
  sigmaResult_est <- sigmaGibbs(ZResult, muResult, ao, bo, y)
  
  pResult_est <- pGibbs(kappaResult, nn)
  
  muResult_est <- muGibbs(ZResult, sigmaResult_est, muo, sigmao, y, pResult_est)
  
  uResult_est <- uslice(pResult_est, ZResult)
  
  ZResult_est <- ZGibbs(uResult_est, pResult_est, muResult_est, sigmaResult_est, y)
  
  ncluster_est <- tabulate(ZResult_est)
  
  KResult_est <- length(pResult_est)
  
  kappaResult_est <- kappaGibbs(kappaResult, aalpha, balpha, KResult_est)
  
  return(list(sigmaResult_est = sigmaResult_est, nn = nn, pResult_est = pResult_est, muResult_est = muResult_est, uResult_est = uResult_est, ZResult_est = ZResult_est, ncluster_est = ncluster_est,  KResult_est =  KResult_est, kappaResult_est = kappaResult_est))
}

gradient_tau <- function(y_index, tau, X){
  K <- length(tau)
  rt <- - diag(1/rep(100, K))%*%tau
  #n <- table(y)
  cut <- cumsum(exp(c(-Inf, tau[-1])))+tau[1]
  cutoff_u <- list()
  cutoff_l <- list()
  dnormcutoff_u <- list()
  dnormcutoff_l <- list()
  pnormcutoff_u <- list()
  pnormcutoff_l <- list()
  for(j in 1:K){
    cutoff_l[[j]] <- cut[j] - X[y_index[[j]]]
    cutoff_u[[j]] <- cut[j] - X[y_index[[j+1]]]
    dnormcutoff_u[[j]] <- dnorm(cutoff_u[[j]])
    dnormcutoff_l[[j]] <- dnorm(cutoff_l[[j]])
    pnormcutoff_u[[j]] <- pnorm(cutoff_u[[j]])
    pnormcutoff_l[[j]] <- pnorm(cutoff_l[[j]])
  }
  rt[1] <- sum(dnormcutoff_l[[1]]/pnormcutoff_l[[1]])+sum(((dnormcutoff_l[[2]]-dnormcutoff_u[[1]])/(pnormcutoff_l[[2]]-pnormcutoff_u[[1]])))+sum(((dnormcutoff_l[[3]]-dnormcutoff_u[[2]])/(pnormcutoff_l[[3]]-pnormcutoff_u[[2]])))-sum(dnormcutoff_u[[3]]/(1-pnormcutoff_u[[3]]))
  rt[2] <- exp(tau[2])*(sum(((dnormcutoff_l[[2]])/(pnormcutoff_l[[2]]-pnormcutoff_u[[1]])))+sum(((dnormcutoff_l[[3]]-dnormcutoff_u[[2]])/(pnormcutoff_l[[3]]-pnormcutoff_u[[2]])))-sum(dnormcutoff_u[[3]]/(1-pnormcutoff_u[[3]])))
  rt[3] <- exp(tau[3])*(sum(((dnormcutoff_l[[3]])/(pnormcutoff_l[[3]]-pnormcutoff_u[[2]])))-sum(dnormcutoff_u[[3]]/(1-pnormcutoff_u[[3]])))
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
  
  if(log(runif(1)) <  -energyH(y, x, p, X) + energyH(y, current_tau, current_p, X)){
    return(x)
  }else{
    return(current_tau)
  }
}

