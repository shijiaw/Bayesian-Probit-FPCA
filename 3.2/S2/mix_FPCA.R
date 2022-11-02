
beta_est <- matrix(NA, nr = N, nc = nbasis)
xi_est <- matrix(NA, nr = N, nc = K)
Theta_est <- matrix(NA, nr = nbasis, nc = K)
Lambda_est <- 2/2.5^(0:(K-1))
meancurve_est <- apply(curves, 2, mean)
meancurve_est <- meancurve_est - mean(meancurve_est)

if(K <= 3){
  xi_est <- xi[,1:K]  + matrix(rnorm(N*K, 0, 0.1), nr = N, nc = K)
}else{
  xi_est <- cbind(xi, matrix(0, nr = nrow(xi), nc = K-3))  + matrix(rnorm(N*K, 0, 0.1), nr = N, nc = K)
}

X_temp <- kronecker(xi_est, obasis_eval)
Theta_est <- matrix(solve(t(X_temp)%*%X_temp)%*%t(X_temp)%*%as.vector(t(curves-matrix(rep(meancurve_est, N), nr = N, byrow = TRUE))), nr = nbasis)
phi_est <- t(obasis_eval%*%Theta_est)
beta_est <- xi_est%*%t(Theta_est)
c_est <- solve(t(basismat)%*%basismat)%*%t(basismat)%*%(meancurve_est)
w_est <- solve(t(X)%*%X)%*%t(X)%*%apply((curves - matrix(rep(basismat%*%c_est, N), nr = N, byrow = TRUE)), 1, mean)
tau_est <- c(-3.5, -1, -1)


#X <- as.vector(xi_est%*%t(obasis_eval%*%Theta_est))
XXX <- as.vector(beta_est%*%t(obasis_eval) + matrix(rep(meancurve_est, N), nr = N, byrow = TRUE) + matrix(rep(X%*%w_est, T), nr = N))
tau_est <- HMC_tau(gradient_tau, energyH_tau, eps_tau, L_tau, tau_est, y_index, XXX)
sig2e <- 1


for(iter in 1:niter){
  #print(iter)
  oneFPCAmove <- OneFPCAMCMCmove(tau_est, xi_est, Theta_est, Lambda_est, c_est, w_est, N, T)
  xi_est <- oneFPCAmove$xi_est
  Theta_est <- oneFPCAmove$Theta_est
  Lambda_est <- oneFPCAmove$Lambda_est
  c_est <- oneFPCAmove$c_est
  w_est <- oneFPCAmove$w_est
  tau_est <- oneFPCAmove$tau_est 
 
  xi_chain[,iter] <- matrix(xi_est, nc = 1)
  Theta_chain[,iter] <- matrix(Theta_est, nc = 1)
  Lambda_chain[,iter] <- Lambda_est
  c_chain[,iter] <- c_est
  w_chain[,iter] <- w_est
  tau_chain[,iter] <- tau_est
}


