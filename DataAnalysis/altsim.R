#---------#
#-Library-#
#---------#

library(foreach)
library(parallel)
library(doParallel)
library(tidyverse)

#-------------#
#-Environment-#
#-------------#

envr <- environment()

#-----------#
#-Functions-#
#-----------#

source("monitoraltfun.R")

file2env <- function(filename, envr){
  load(filename)
  list2env(x = data, envir = envr)
}

#-------------------------#
#-Monitoring Alternatives-#
#-------------------------#

nsim <- 25 #Number of population process simulations
niter <- 25 #Number of observation process simulations per population
nscenario <- 7 #Number of population process scenarios
nalt <- 5
metadata <- expand_grid(scenario = 1:nscenario, simID = 1:nsim, iterID = 1:niter) #Metadata for simulation IDs 
ncombos <- nsim * niter * nscenario #Number of simulations in aggregate 

scenario <- 1
simID <- 1
iterID <- 1

clustID <- makeCluster(25)

registerDoParallel(clustID)

foreach(i = 1:ncombos) %dopar% {

#-Population process sim to load-#
path <- paste0("./scenario/scenario", metadata$scenario[i], "_", metadata$simID[i], ".Rds")

#-Alternative 1-#

#Colony selection: split panel; core: Aiktak, Buldir; rotating: random
#Site selection: core: preferential; rotating: random
#Sampling type: core: index plots; rotating: circular plots
#Sampling intensity: 10 years; 2 core colonies; 8 rotating colonies; Aiktak 10 sites; Buldir 1 site; rotating 35 sites; 1 replicate per site
#Data type: Burrow counts and occupancy

file2env(path, envr)
list2env(design.fun(colonies = colony.data$Island, 
                    core.units = c("Aiktak", "Buldir"),  
                    subsampling = colony.data$Subsampling, 
                    unit.selection = "splitpanel", 
                    scale.selection = rep("site", 10), 
                    site.selection = c(rep("preferential", 2), rep("random", 8)), 
                    sampling.type = c(rep("indexplot", 2), rep("circleplot", 8)), 
                    nunits = 10, 
                    ncore = 10, 
                    nsubsampled = 10, 
                    nsites = c(10, 1, rep(35, 8)), 
                    nreps = 1), envir = envr)

output <- list(
  data = list(y.burrow = y.burrow, y.occ = y.occ),
  constants = list(nyears = nyears, ncolonies = ncolonies, nunits = nunits, nsites = nsites, units = units, COS = COS),
  Lambda0 = Lambda0,
  mu.gamma = mu.gamma,
  sd.gamma = sd.gamma,
  gamma0 = gamma0,
  psi0 = psi0
)

save(output, file = paste0("./alternative/scenario", metadata$scenario[i], "_alt1_sim", metadata$simID[i], "_iter", metadata$iterID[i], ".Rds"))
     
#-Alternative 2-#

#Colony selection: split panel (only core units)
#Site selection: core:preferential
#Sampling type: index plots
#Sampling intensity: 10 years; 10 core colonies; 1 site; 1 replicate per site
#Data type: Burrow counts and occupancy

file2env(path, envr)
list2env(design.fun(colonies = colony.data$Island, 
                   core.units = NA,  
                   subsampling = colony.data$Subsampling, 
                   unit.selection = "splitpanel", 
                   scale.selection = rep("site", 10), 
                   site.selection = rep("preferential", 10), 
                   sampling.type = rep("indexplot", 10), 
                   nunits = 10, 
                   ncore = 10, 
                   nsubsampled = 10, 
                   nsites = rep(1, 10), 
                   nreps = 1), envir = envr)
output <- list(
  data = list(y.burrow = y.burrow, y.occ = y.occ),
  constants = list(nyears = nyears, ncolonies = ncolonies, nunits = nunits, nsites = nsites, units = units, COS = COS),
  Lambda0 = Lambda0,
  mu.gamma = mu.gamma,
  sd.gamma = sd.gamma,
  gamma0 = gamma0,
  psi0 = psi0
)

save(output, file = paste0("./alternative/scenario", metadata$scenario[i], "_alt2_sim", metadata$simID[i], "_iter", metadata$iterID[i], ".Rds"))

#-Alternative 3-#

#Colony selection: split panel (only core units)
#Site selection: core:random
#Sampling type: circle plots
#Sampling intensity: 10 years; 10 core colonies; 35 site; 1 replicate per site
#Data type: Burrow counts and occupancy

file2env(path, envr)
list2env(design.fun(colonies = colony.data$Island, 
                   core.units = NA,  
                   subsampling = colony.data$Subsampling, 
                   unit.selection = "splitpanel", 
                   scale.selection = rep("site", 10), 
                   site.selection = rep("random", 10), 
                   sampling.type = rep("circleplot", 10), 
                   nunits = 10, 
                   ncore = 10, 
                   nsubsampled = 10, 
                   nsites = rep(35, 10), 
                   nreps = 1), envir = envr)
output <- list(
  data = list(y.burrow = y.burrow, y.occ = y.occ),
  constants = list(nyears = nyears, ncolonies = ncolonies, nunits = nunits, nsites = nsites, units = units, COS = COS),
  Lambda0 = Lambda0,
  mu.gamma = mu.gamma,
  sd.gamma = sd.gamma,
  gamma0 = gamma0,
  psi0 = psi0
) 

save(output, file = paste0("./alternative/scenario", metadata$scenario[i], "_alt3_sim", metadata$simID[i], "_iter", metadata$iterID[i], ".Rds"))

#-Alternative 4-#

#Colony selection: split panel (only core units)
#Site selection: Aiktak & Buldir
#Sampling type: index plots
#Sampling intensity: 10 years; 2 core colonies; 1 site; 1 replicate per site
#Data type: Burrow counts and occupancy

file2env(path, envr)
list2env(design.fun(colonies = colony.data$Island, 
                   core.units = c("Aiktak", "Buldir"),  
                   subsampling = colony.data$Subsampling, 
                   unit.selection = "splitpanel", 
                   scale.selection = rep("site", 2), 
                   site.selection = rep("preferential", 2), 
                   sampling.type = rep("indexplot", 2), 
                   nunits = 2, 
                   ncore = 2, 
                   nsubsampled = 2, 
                   nsites = rep(1, 2), 
                   nreps = 1), envir = envr)
output <- list(
  data = list(y.burrow = y.burrow, y.occ = y.occ),
  constants = list(nyears = nyears, ncolonies = ncolonies, nunits = nunits, nsites = nsites, units = units, COS = COS),
  Lambda0 = Lambda0,
  mu.gamma = mu.gamma,
  sd.gamma = sd.gamma,
  gamma0 = gamma0,
  psi0 = psi0
)

save(output, file = paste0("./alternative/scenario", metadata$scenario[i], "_alt4_sim", metadata$simID[i], "_iter", metadata$iterID[i], ".Rds"))

#-Alternative 5-#

#Colony selection: split panel (only core units)
#Site selection: core:random
#Sampling type: circle plots
#Sampling intensity: 10 years; 15 core colonies; 35 site; 1 replicate per site
#Data type: Burrow counts and occupancy

file2env(path, envr)
list2env(design.fun(colonies = colony.data$Island, 
                    core.units = NA,  
                    subsampling = colony.data$Subsampling, 
                    unit.selection = "splitpanel", 
                    scale.selection = rep("site", 15), 
                    site.selection = rep("random", 15), 
                    sampling.type = rep("circleplot", 15), 
                    nunits = 15, 
                    ncore = 15, 
                    nsubsampled = 15, 
                    nsites = rep(35, 15), 
                    nreps = 1), envir = envr)

output <- list(
  data = list(y.burrow = y.burrow, y.occ = y.occ),
  constants = list(nyears = nyears, ncolonies = ncolonies, nunits = nunits, nsites = nsites, units = units, COS = COS),
  Lambda0 = Lambda0,
  mu.gamma = mu.gamma,
  sd.gamma = sd.gamma,
  gamma0 = gamma0,
  psi0 = psi0
)

save(output, file = paste0("./alternative/scenario", metadata$scenario[i], "_alt5_sim", metadata$simID[i], "_iter", metadata$iterID[i], ".Rds"))

# if(simID <= nsim){
#   scenario = scenario
#   simID = simID + 1
#   iterID = iterID
# }else{
#   if(iterID <= niter){
#     scenario = scenario
#     simID =  1
#     iterID = iterID + 1
#   }else{
#     if(scenario <= nscenario){
#      scenario = scenario + 1
#      simID = 1
#      iterID = 1
#     }
#   }
# }
  
}

stopCluster(clustID)

metadata <- expand_grid(scenario = 1:nscenario, altID = 1:5, simID = 1:nsim, iterID = 1:niter)

save(metadata, file = "metadata.Rds")
