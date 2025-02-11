#---------#
#-Library-#
#---------#

library(foreach)
library(parallel)
library(doParallel)

#-----------#
#-Functions-#
#-----------#

source("biosimfun.R")

#----------------------#
#-Biological Scenarios-#
#----------------------#

nsim <- 25

clustID <- makeCluster(25)

registerDoParallel(clustID)

foreach(i = 1:nsim) %dopar% {
  
  #-Scenario 1: Status quo-#
  
  data <- biosim.fun(scenario = "status quo", nyears = 10)
  save(data, file = paste0("./scenario/scenario1_", i, ".Rds")) 
  
  #-Scenario 2: IUCN Critical Endangered Decline-#
  
  data <- biosim.fun(scenario = "IUCN CR", nyears = 10)
  save(data, file = paste0("./scenario/scenario2_", i, ".Rds")) 
  
  #-Scenario 3: IUCN Vulnerable Decline-#

  data <- biosim.fun(scenario = "IUCN VU", nyears = 10)
  save(data, file = paste0("./scenario/scenario3_", i, ".Rds"))
  
  #-Scenario 4: Variable growth between colonies-#
  
  data <- biosim.fun(scenario = "Variable growth between colonies", nyears = 10)
  save(data, file = paste0("./scenario/scenario4_", i, ".Rds"))
  
  #-Scenario 5: Variable growth within colonies-#

  data <- biosim.fun(scenario = "Variable growth within colonies", nyears = 10)
  save(data, file = paste0("./scenario/scenario5_", i, ".Rds"))
  
  #-Scenario 6: Variable density within colonies-#
  
  data <- biosim.fun(scenario = "Variable density within colonies", nyears = 10)
  save(data, file = paste0("./scenario/scenario6_", i, ".Rds"))
  
  #-Scenario 7: Variable growth between and within colonies and variable density-#
  
  data <- biosim.fun(scenario = "Variation in everything", nyears = 10)
  save(data, file = paste0("./scenario/scenario7_", i, ".Rds"))

}

stopCluster(clustID)
