nyears <- 10
nLME <- 5
nunits <- 50
nsites <- 25

#-----------------------#
#-Multi-scale covariate-#
#-----------------------#

mean.cov <- runif(nLME, 3, 6) #mean covariate value in year 1 for each LME
sd.cov <- runif(1, 0, 0.5) #annual variation in covariate change
delta.cov <- runif(1, 0.9, 1.1) #mean annual percent change in covariate value
Y <- Y.e <- matrix(NA, nrow = nyears, ncol = nLME)
Y.e[1,] <- mean.cov
for(t in 2:nyears){
  Y.e[t,] <- rnorm(nLME, mean = Y.e[t-1,] * delta.cov, sd.cov)
}
X <- X.e <- array(NA, dim = c(nyears, nLME, nunits))
for(t in 1:nyears){
  for(l in 1:nLME){
    X.e[t,l,] <- rnorm(nunits, Y.e[t,l], 1)
  }
}
x <- x.e <- array(NA, dim = c(nyears, nLME, nunits, nsites))
for(t in 1:nyears){
  for(l in 1:nLME){
    for(i in 1:nunits){
      x.e[t,l,i,] <- rnorm(nsites, X.e[t,l,i], 1)
    }
  }
}

x.mean <- apply(x.e, c(1,2,3), mean)
x.sd <- apply(x.e, c(1,2,3), sd)
X.mean <- apply(x.mean, c(1,2), mean)
X.sd <- apply(x.mean, c(1,2), sd)
Y.mean <- mean(X.mean)
Y.sd <- sd(X.mean)

# Scale covariates    
for(t in 1:nyears){
  for(l in 1:nLME){
    Y[t,l] <- (X.mean[t,l] - Y.mean)/Y.sd
    for(i in 1:nunits){
      X[t,l,i] <- (x.mean[t,l,i] - X.mean[t,l])/X.sd[t,l]
      for(j in 1:nsites){
        x[t,l,i,j] <- (x.e[t,l,i,j] - x.mean[t,l,i])/x.sd[t,l,i]
      }
    }
  }
}

#------------------------#
#-Population growth rate-#
#------------------------#

gamma0 <- log(runif(1, 0.9, 1.1))
gamma1 <- runif(1, -0.25, 0.25)
gamma <- matrix(NA, nrow = nyears - 1, ncol = nLME)
for(k in 1:(nyears - 1)){
  for(l in 1:nLME){
    gamma[k,l] <- exp(gamma0 + gamma1 * Y[k,l])
  }
}

# LME mean abundance
Lambda.bar <- matrix(NA, nrow = nyears, ncol = nLME)
Lambda.bar[1,] <- runif(nLME, 100000, 200000)

for(t in 2:nyears){
  k <- t - 1
  for(l in 1:nLME){
    Lambda.bar[t,l] <- Lambda.bar[1,l] * prod(gamma[1:k,l])
  }
}

Area <- runif(nLME, 100, 200)


# Unit scale abundance
tau <- NULL
beta1 <- runif(1, -1, 1)
N <- Lambda <- Mu <- X <- array(NA, dim = c(nyears, nLME, nunits))
for(l in 1:nLME){
  #tau[l] <- abs(mean(Lambda.bar[,l])/nunits - runif(1, mean(Lambda.bar[,l])/nunits * 0.7, mean(Lambda.bar[,l])/nunits * 1.3))
  tau[l] <- runif(1, 0, 1)
  for(t in 1:nyears){
    for(i in 1:nunits){
      X[t,l,i] <- rnorm(1, 0, 1)
      Mu[t,l,i] <- rnorm(1, log(Lambda.bar[t,l]/nunits), tau[l])
      Lambda[t,l,i] <- exp(Mu[t,l,i] + beta1 * X[t,l,i])
      N[t,l,i] <- rpois(1, Lambda[t,l,i])
    }
  }
}


# Site scale abundance
COS <- sigma <- matrix(NA, nrow = nLME, ncol = nunits)
alpha1 <- runif(1, -1, 1)
n <- lambda <- mu <- x <- array(NA, dim = c(nyears, nLME, nunits, nsites))
for(l in 1:nLME){
  for(i in 1:nunits){
    sigma[l,i] <- runif(1, 0, 1)
    COS[l,i] <- runif(1, 50, 100)
    for(j in 1:nsites){
      for(t in 1:nyears){
        x[t,l,i,j] <- rnorm(1, 0, 1)
        mu[t,l,i,j] <- rnorm(1, log(Lambda[t,l,i]/COS[l,i]), sigma[l,i])
        lambda[t,l,i,j] <- exp(mu[t,l,i,j] + alpha1 * x[t,l,i,j])
        n[t,l,i,j] <- rpois(1, lambda[t,l,i,j])
      }
    }
  }
}




