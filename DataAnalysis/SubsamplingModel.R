#-------------#
#-Subsampling-#
#-------------#

#-----------#
#-Libraries-#
#-----------#

library(nimble)
library(coda)
library(parallel)

#------------#
#-Model code-#
#------------#

code <- nimbleCode({
  
  psi ~ dunif(0, 1)
  psi0 <- logit(psi)
  
  mu.gamma ~ dnorm(0, 0.001)
  mu.gamma.real <- mu.gamma
  tau.gamma ~ dgamma(0.001, 0.001)
  sd.gamma <- 1/sqrt(tau.gamma)
  
  mu.Lambda0 ~ dunif(0, 500000) #Expected mean abundance of a colony in year 1
  tau.Lambda0 ~ dgamma(0.1, 0.1) #Precision of colony variation in abundance in year 1
  sd.Lambda0 <- 1/sqrt(tau.Lambda0)
  
  for(i in 1:ncolonies){
    eps.gamma[i] ~ dnorm(0, tau.gamma)
    log(gamma0[i]) <- mu.gamma + eps.gamma[i]
    gamma0.real[i] <- gamma0[i]
    Lambda0[i] ~ dnorm(mu.Lambda0, tau.Lambda0)
    Lambda[1,i] <- Lambda0[i]
    
    for(t in 2:nyears){
      Lambda[t,i] <- Lambda[t-1,i] * gamma0[i]
    }
  }
  
  for(t in 1:nyears){
    for(i in 1:nunits){
      for(j in 1:nsites[i]){
        mu[t,i,j] <- (Lambda[t,units[t,i]]/(2 * COS[t,i])) * (1/psi)
        y.burrow[t,i,j] ~ dpois(mu[t,i,j])
        y.occ[t,i,j] ~ dbin(psi, y.burrow[t,i,j])
      }#end j
    }#end i
  }#end t
  
})

#--------------#
#-Compile data-#
#--------------#

# data <- list(y.burrow = y.burrow, y.occ = y.occ)
# 
# constants <- list(nyears = nyears, ncolonies = ncolonies, nunits = nunits, nsites = nsites, units = units, COS = COS)

data <- data

constants <- constants

#----------------#
#-Initial values-#
#----------------#

inits <- list(
  Lambda0 = Lambda0, #Stalls w/o
  mu.Lambda0 = mean(Lambda0), #Still runs w/o
  sd.Lambda0 = sd(Lambda0),
  tau.Lambda0 = 1/(sd(Lambda0))^2,
  mu.gamma = mu.gamma, #runif(1, min(gamma0 * 0.95, gamma0 * 1.05), max(gamma0 * 0.95, gamma0 * 1.05)),
  gamma0 = gamma0,
  sd.gamma = sd.gamma, #runif(1, 0.04, 0.06),
  tau.gamma = 1/sd.gamma^2,
  eps.gamma = gamma0 - mu.gamma,
  psi = expit(rnorm(1, 0.5539434, 0.019))
)

#--------------------#
#-Parameters to save-#
#--------------------#

params <- c(
  "Lambda0",
  "sd.gamma",
  # "theta.sd",
  #"p",
  "mu.gamma",
  "mu.gamma.real",
  "gamma0.real",
  "gamma0",
  "psi0"
)

#-------------------#
#-Parallel function-#
#-------------------#

par.fun <- function(seed, data, code, constants, params, inits){
  
  library(nimble)
  library(coda)
  
#---------------#
#-MCMC settings-#
#---------------#

model <- nimbleModel(code = code,
                     constants = constants,
                     data = data,
                     inits = inits)

MCMCconf <- configureMCMC(model, monitors = params)

#ncols.after1 <- (1:constants$ncolonies)[-constants$units[1,which(!(constants$units[1,] %in% unique(as.vector(constants$units[-1,]))))]]
ncols.after1 <- 1:constants$nunits

MCMCconf$removeSampler(target = c("mu.gamma", paste0("eps.gamma[", ncols.after1, "]"), paste0("Lambda0[", ncols.after1, "]")))

MCMCconf$addSampler(target = c("mu.gamma", paste0("eps.gamma[", ncols.after1, "]"), paste0("Lambda0[", ncols.after1, "]")),
                    type = "AF_slice")


MCMC <- buildMCMC(MCMCconf)

compiled.model <- compileNimble(model, MCMC)

ni <- 20000
nb <- 15000
nc <- 1
nt <- 5

out <- runMCMC(compiled.model$MCMC,
               niter = ni, nburnin = nb,
               nchains = nc, thin = nt,
               samplesAsCodaMCMC = TRUE)

return(out)
}

#-----------#
#-Run model-#
#-----------#

clustID <- makeCluster(3)

out <- parLapply(cl = clustID, X = 1:3,
                 fun = par.fun,
                 code = code,
                 data = data,
                 constants = constants,
                 inits = inits,
                 params = params)

stopCluster(clustID)

out <- mcmc.list(out)
