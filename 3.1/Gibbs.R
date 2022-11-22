###Gibbs function ####
###We want to sample parameters: ###
###mu, sigma, Z, p, u, beta###
###We first update mu, sigma ###
###Then update beta, p###
###updata u and label ###
###Finally update kappa, the concentration parameter ###
###Gibbs: sigma ###
sigmaGibbs <- function(Zvec, mu, ao, bo, y){
  nn <- tabulate(Zvec)
  K <- length(nn)
  a <- nn/2 + ao

  subY <- list()
  dim <- nrow(y)
  sigmaRt <- matrix(NA, nr = dim, nc = K)
  b <- matrix(NA, nr = dim, nc = K)
  for(i in 1:K){
    subY[[i]] <- as.matrix(y[ ,which(Zvec==i)])
#    if(length(subY[[i]])>1){
    if(i <= ncol(mu)){
      for(d in 1:dim){
        #b[d, i] <- bo + sum((subY[[i]]-mu[i])^2)/2
        b[d, i] <- bo + sum((subY[[i]][d,]-mu[d,i])^2)/2
        sigmaRt[d, i] <- sqrt(rinvgamma(1, a[i], b[d,i]))
      }
    }else{
      sigmaRt[,i] <- rep(3, dim)
    }
  }
  #sigmaRt[K+1] <- 3
  return(sigmaRt)
}

# constrained_normal
cons_sampling <- function(pi, mu, sigma){
  len <- length(mu)
  w <- pi*sigma
  c <- -sum(pi*mu)
  U <- c*w[-len]^2/(sum(w^2)) 
  if(len > 2){
    Sigma <- (diag(w[-len]^2)*sum(w^2) - w[-len]^2%*%t(w[-len]^2))/sum(w^2)
  }else{
    Sigma <- (w[-len]^2*sum(w^2) - w[-len]^2*(w[-len]^2))/sum(w^2)
  }
 
  z_tild <- rnorm(len-1, 0, 1)
  if(len > 2){
    eigen <- eigen(Sigma)
    x_tild <- U + eigen$vectors%*%diag(sqrt(eigen$values))%*%z_tild
  }else{
    x_tild <- U + sqrt(Sigma)*z_tild
  }
  
  x <- c(x_tild, c-sum(x_tild))
  output <- x/w*sigma + mu
  return(output)
}


###Gibbs: mu ###
muGibbs <- function(Zvec, sigma, muo, sigmao, y, pi){
  nn <- tabulate(Zvec)
  K <- length(nn)
  dim <- nrow(y)
  sigma_new <- matrix(NA, nr = dim, nc = K)
  for(i in 1:K){
    sigma_new[ ,i] <- sqrt(1/(1/sigmao^2+nn[i]/sigma[,i]^2)) # it's a vector
  }
  ## Fix the bug from here 
  subY <- list()
  munew <- matrix(NA, nr = dim, nc = K)
  muRt <- matrix(0, nr = dim, nc = K)
  for(i in 1:K){
    subY[[i]] <- y[,which(Zvec==i)]
      if(is.matrix(subY[[i]])){
        munew[,i] <- (1/sigma[,i]^2*muo+ 1/sigmao^2*apply(subY[[i]], 1, sum))/(1/sigma[,i]^2+1/sigmao^2*nn[i])
        #muRt[,i] <- rnorm(dim, munew[, i], sigma_new[, i])
      }else{
        munew[,i] <- (1/sigma[,i]^2*muo+1/sigmao^2*sum(subY[[i]]))/(1/sigma[,i]^2+1/sigmao^2*nn[i])
        #muRt[,i] <- rnorm(dim, munew[, i], sigma_new[, i])
      }
  }
  
  #for(j in 1:dim){
  for(j in 1:1){
    muRt[j, ] <- cons_sampling(pi, munew[j, ], sigma_new[j, ])
  }

  return(muRt)
}

###function for beta and p ###
pGibbs <- function(alpha, nn){
  K.max <- length(nn)
  beta <- rep(NA, K.max)
  p <- rep(NA, K.max)
  if(K.max > 1){
    for(i in 1:(K.max-1)){
      shape1 <- 1+nn[i]
      shape2 <- alpha+sum(nn[-(1:i)])
      beta[i] <- rbeta(1, shape1, shape2)
    }
    beta[K.max] <- rbeta(1, 1+nn[K.max], alpha)
    p[1] <- beta[1]
    for(i in 2:K.max){
      p[i] <- beta[i]*(1-sum(p[1:(i-1)]))
    }
  }else{
    shape1 <- 1+nn[1]
    shape2 <- alpha
    p[1] <- 1
  }
  
  return(p)
}

###slice sampler for u ###
uslice <- function(p, Zvec){
  u <- rep(NA, N)
  for(i in 1:N){
    u[i] <- runif(1, 0, p[Zvec[i]])
  }
  return(u)
}

###function for Z ###
uIndicator <- function(u, p){
  as.numeric(p > u)
}

ZGibbs <- function(u, p, mu, sigma, y){
  K <- length(p)
  yMatrix <- matrix(NA, nr = K+1, nc = N)
  posteriorZ <- matrix(NA, nr = K+1, nc = N)
  Zvec <- rep(NA, N)
  for(i in 1:K){
    yMatrix[i,] <- apply(dnorm(y, mu[,i], sigma[,i]),2,prod)
      #dnorm(y, mu[,i], sigma[,i])
    posteriorZ[i,] <- uIndicator(u, p[i])*yMatrix[i,] 
  }
  yMatrix[K+1,] <- apply(dnorm(y, mu_temp, sigma_temp), 2, prod) ##check here
  posteriorZ[K+1,] <- uIndicator(u, (1-sum(p)))*yMatrix[K+1,] 
  Znormalizer <- apply(posteriorZ, 2, sum)  
  for(i in 1:N){
    posteriorZ[,i] <- posteriorZ[,i]/Znormalizer[i]
    Zvec[i] <- sample(1:(K+1), size = 1, prob = posteriorZ[,i])
  }
  return(as.numeric(as.factor(Zvec)))
}

###update kappa ####
kappaGibbs <- function(kappa, aalpha, balpha, K){
  xi <- rbeta(1, kappa+1, N)
  eps <- (aalpha+K-1)/(aalpha+K-1+N*(balpha-log(xi)))
  ueps <- runif(1, 0, 1)
  if(ueps < eps){
    kappa <- rgamma(1, aalpha+K, balpha-log(xi))
  }else{
    kappa <- rgamma(1, aalpha+K-1, balpha-log(xi))
  }
  return(kappa)
}

