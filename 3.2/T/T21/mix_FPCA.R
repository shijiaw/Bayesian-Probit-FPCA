
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

lambda <- 3/(2^((1:K)-1))

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
#deriv1 <- deriv(obasis)
#obasis_deriv <- evaluate(deriv1, tobs)

#define B-spline basis functions
knots_B = seq(0, 1 ,0.1)
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
Theta_est <- matrix(solve(t(X_temp)%*%X_temp)%*%t(X_temp)%*%as.vector(t(curves-matrix(rep(meancurve_est, N), nr = N, byrow = TRUE))), nr = nbasis)
phi_est <- t(obasis_eval%*%Theta_est)
beta_est <- xi_est%*%t(Theta_est)
c_est <- solve(t(basismat)%*%basismat)%*%t(basismat)%*%(meancurve_est)
w_est <- solve(t(X)%*%X)%*%t(X)%*%apply((curves - matrix(rep(basismat%*%c_est, N), nr = N, byrow = TRUE)), 1, mean)
tau_est <- c(-2, -1, -1)

# MCMC for FPCA 
xi_chain <- matrix(NA, nr = N*K, nc = niter)
Theta_chain <- matrix(NA, nr = nbasis*K, nc = niter)
Lambda_chain <- matrix(NA, nr = K, nc = niter)
#sig2e_chain <- rep(NA, niter)
c_chain <- matrix(NA, nr = nbasis, nc = niter)
Thetasample_chain <- matrix(NA, nr = nbasis*K, nc = niter)
w_chain <- matrix(NA, nr = 2, nc = niter)
tau_chain <- matrix(NA, nr = length(d), nc = niter)

#source('HMC_probit.R')



#X <- as.vector(xi_est%*%t(obasis_eval%*%Theta_est))
XXX <- as.vector(beta_est%*%t(obasis_eval) + matrix(rep(meancurve_est, N), nr = N, byrow = TRUE) + matrix(rep(X%*%w_est, T), nr = N))
cutoff <- c(-Inf, cumsum(exp(c(-Inf, tau_est[-1])))+tau_est[1], Inf)
inv_l <- cutoff[y]
inv_u <- cutoff[y+1]  
y_star <- rtruncnorm(length(XXX), a= inv_l, b=inv_u, mean = X, sd = 1)
sdobs <- matrix(y_star, nr = N)
sig2e <- 1


for(iter in 1:niter){
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


