
MCMCmove <- function(ZResult1_est, ZResult2_est, muResult_est, kappaResult_est, sigmaResult_est, pResult_est, uResult_est, muo, sigmao, aalpha, balpha, ao, bo, xi_est,tau_est, Theta_est, Lambda_est, c_est, w_est){
  xi_chain <- matrix(NA, nr = N*K, nc = niter)
  Theta_chain <- matrix(NA, nr = nbasis*K, nc = niter)
  Lambda_chain <- matrix(NA, nr = K, nc = niter)
  sig2e_chain <- rep(NA, niter)
  c_chain <- matrix(NA, nr = nbasis, nc = niter)
  w_chain <- matrix(NA, nr = 2, nc = niter)
  sigmaResult_chain <- list()
  muResult_chain <- list()
  pResult_chain <- list()
  uResult_chain <- list()
  ZResult_chain <- list()
  KResult_chain <- rep(NA, niter)
  kappaResult_chain <- rep(NA, niter)
  ncluster_chain <- list()
  tau_chain <- matrix(NA, nr = length(d), nc = niter)
  
  for(iter in 1:niter){
    #print(iter)
    OneMixmove <- OneMixMCMCmove(ZResult1_est, muResult_est, kappaResult_est, sigmaResult_est, pResult_est, uResult_est, muo, sigmao, aalpha, balpha, ao, bo, t(xi_est))
    sigmaResult_est <- OneMixmove$sigmaResult_est
    pResult_est <- OneMixmove$pResult_est
    muResult_est <- OneMixmove$muResult_est 
    uResult_est <- OneMixmove$uResult_est
    ZResult2_est <- OneMixmove$ZResult_est
    ncluster_est <- OneMixmove$ncluster_est  
    KResult_est <- OneMixmove$KResult_est
    kappaResult_est <- OneMixmove$kappaResult_est
    sigmaResult_chain[[iter]] <- sigmaResult_est
    muResult_chain[[iter]] <- muResult_est
    pResult_chain[[iter]] <- pResult_est
    uResult_chain[[iter]] <- uResult_est
    ZResult_chain[[iter]] <- ZResult2_est
    KResult_chain[iter] <- KResult_est
    kappaResult_chain[iter] <- kappaResult_est
    ncluster_chain[[iter]] <- ncluster_est
    
    
    oneFPCAmove <- OneFPCAMCMCmove(tau_est, xi_est, Theta_est, Lambda_est, c_est, w_est, ZResult1_est, muResult_est, sigmaResult_est)
    xi_est <- oneFPCAmove$xi_est
    Theta_est <- oneFPCAmove$Theta_est
    Lambda_est <- oneFPCAmove$Lambda_est
    c_est <- oneFPCAmove$c_est
    w_est <- oneFPCAmove$w_est
    tau_est <- oneFPCAmove$tau_est 
    ZResult1_est <- ZResult2_est
    
    xi_chain[,iter] <- matrix(xi_est, nc = 1)
    Theta_chain[,iter] <- matrix(Theta_est, nc = 1)
    Lambda_chain[,iter] <- Lambda_est
    c_chain[,iter] <- c_est
    w_chain[,iter] <- w_est
    #sig2e_chain[iter] <- sig2e
    tau_chain[,iter] <- tau_est
  }
  return(list(xi_chain = xi_chain, Theta_chain = Theta_chain, Lambda_chain = Lambda_chain, c_chain = c_chain, w_chain = w_chain, tau_chain = tau_chain, sigmaResult_chain = sigmaResult_chain, muResult_chain = muResult_chain, pResult_chain = pResult_chain, uResult_chain = uResult_chain, ZResult_chain = ZResult_chain, KResult_chain = KResult_chain, kappaResult_chain = kappaResult_chain, ncluster_chain = ncluster_chain))
}


