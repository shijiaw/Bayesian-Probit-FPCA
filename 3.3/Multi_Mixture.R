rm(list=ls())
library(invgamma)
library(gtools)
library(mvtnorm)
set.seed(120)

N = 500
K = 3
pi = c(0.2, 0.35, 0.45)
# constraints:\sum pi*mu = 0
dim <- 3
mu <- matrix(NA, nr = dim, nc = K)
mu[,1] <- c(3, -1, 0.2)
mu[,2] <- c(-1, 2, -0.2)
mu[,3] <- -(mu[,1]*pi[1] + mu[,2]*pi[2])/pi[3]
#mu = matrix(c(c(1, -1, 2), c(-3, -2, -1), c(-1, 1, -2)), nr = 3)
sigma = matrix(c(rep(0.1, 3), rep(0.3, 3), rep(0.3, 3)), nr = 3)
sig_lambda[3, ] = c(0.8, 0.7, 0.5)


### priors IG(ao, bo) and normal(muo, sigmao) ####
ao <- 0.1
bo <- 0.1

muo <- 0
sigmao <- 3

aalpha <- 0.1 
balpha <- 0.1

###simulate observations###
###simulate labels      ###
Z <- sample(1:K, size = N, replace = TRUE, prob = pi)
y <- matrix(0, nr = 3, nc = N)
for(i in 1:N){
  y[,i] <- rnorm(K, mean = mu[,Z[i]], sd = sigma[,Z[i]])
}

###choose initial value       ###
###true value plus small noise###
#temp <- abs(rnorm(K, 0, 1))
#pi_init <- rdirichlet(1, rep(1,4))
#remove(temp)
#mu_init <-  mu +rnorm(K, 0, 1)
#sigma_init <-  abs(sigma + rnorm(K, 0, 0.1))
#nn <- tabulate(Z)
### kappa: DP concentration parameter ###
kappa = 0.5

### create initial labels ###
### create initial value for parameters ###
sigma_init <- 0.1
Z = rep(NA, N)
Z[1] = 1

mu_init <- rnorm(dim, muo, sigmao)

nn <- c()
nn[1] <- 1

for(i in 2:N){
  prob <- c(nn, kappa)/sum(c(nn, kappa))
  Z[i] <- sample(1:(length(nn)+1), size = 1, prob = prob)
  if(Z[i] > length(nn)){
    nn[Z[i]] <- 1
    mu_init <- cbind(mu_init, rnorm(dim, muo, sigmao))
  }else{
    nn[Z[i]] <- nn[Z[i]] + 1
  }
}

K <- length(nn)

sigma_init <- rep(3, K)

#### Gibbs sampler ###
source("Gibbs.R")
niter <- 10000
sigmaResult <- list()
muResult <- list()
pResult <- list()
uResult <- list()
ZResult <- list()
KResult <- c()
kappaResult <- c()
ncluster <- list()

sigmaResult[[1]] <- sigmaGibbs(Z, mu_init, ao, bo, y)

pResult[[1]] <- pGibbs(kappa, nn)

muResult[[1]] <- muGibbs(Z, sigmaResult[[1]], muo, sigmao, y, pResult[[1]])

uResult[[1]] <- uslice(pResult[[1]], Z)

mu_temp <- rnorm(dim, muo, sigmao)

sigma_temp <- rep(3,dim)

ZResult[[1]] <- ZGibbs(uResult[[1]], pResult[[1]], muResult[[1]], sigmaResult[[1]], y)

ncluster[[1]] <- tabulate(ZResult[[1]])

KResult[1] <- length(pResult[[1]])

kappaResult[1] <- kappaGibbs(kappa, aalpha, balpha, KResult[1])

for(iter in 2:niter){
  sigmaResult[[iter]] <- sigmaGibbs(ZResult[[iter-1]], muResult[[iter-1]], ao, bo, y)
  
  nn <- tabulate(ZResult[[iter-1]])
  
  pResult[[iter]] <- pGibbs(kappaResult[iter-1], nn)
  
  muResult[[iter]] <- muGibbs(ZResult[[iter-1]], sigmaResult[[iter]], muo, sigmao, y, pResult[[iter]])
  
  uResult[[iter]] <- uslice(pResult[[iter]], ZResult[[iter-1]])
  
  ZResult[[iter]] <- ZGibbs(uResult[[iter]], pResult[[iter]], muResult[[iter]], sigmaResult[[iter]], y)
  
  ncluster[[iter]] <- tabulate(ZResult[[iter]])
  
  KResult[iter] <- length(pResult[[iter]])
  
  kappaResult[iter] <- kappaGibbs(kappaResult[iter-1], aalpha, balpha, KResult[iter])
}

