#-------------#
#-Environment-#
#-------------#

envr <- environment()

#----------#
#-Metadata-#
#----------#

load(file = "metadata2.Rds")

scenario <- head(metadata2$scenario, n = 1)
altID <- head(metadata2$altID, n = 1)
simID <- head(metadata2$simID, n = 1)
iterID <- head(metadata2$iterID, n = 1)

metadata2 <- metadata2[-1,]
save(metadata2, file = "metadata2.Rds")

#-----------#
#-Load data-#
#-----------#

path <- paste0("./alternative/scenario", scenario, "_alt", altID, "_sim",  simID, "_iter", iterID, ".Rds")
load(path)
list2env(output, envir = envr)

#-----------#
#-Run model-#
#-----------#

if(altID %in% c(1:7,11)){
  sys.source("SubsamplingModel.R", envir = envr, toplevel.env = envr)
}else{
  sys.source("IntegratedModel.R", envir = envr, toplevel.env = envr)
}


#-------------#
#-Save Output-#
#-------------#

output <- data.frame(matrix(NA, nrow = 1, ncol = 11))
colnames(output) <- c("scenario", "alternative", "simID", "iterID", "trend", "parameter", "truth", "mean", "sd", "coverage", "rhat")

for(k in 1:length(params)){
  if(length(get(params[k])) == 1){
    output <- rbind(output, 
                    data.frame("scenario" = scenario, "alternative" = altID, "simID" = simID, "iterID" = iterID, "trend" = NA,
                               "parameter" = params[k], "truth" = get(params[k]), 
                               "mean" = summary(out)[[1]][params[k],"Mean"],
                               "sd" = summary(out)[[1]][params[k],"SD"],
                               "coverage" = summary(out)[[2]][params[k],1] <= get(params[k]) & summary(out)[[2]][params[k],5] >= get(params[k]),
                               "rhat" = as.numeric(coda::gelman.diag(out[1:3][,params[k]])[[1]][,1])))
  }else{
    output <- rbind(output, 
                    data.frame("scenario" = scenario, "alternative" = altID, "simID" = simID, "iterID" = iterID, "trend" = "min",
                               "parameter" = params[k], "truth" = get(params[k])[which.min(get(params[k]))], 
                               "mean" = summary(out)[[1]][paste0(params[k], "[", which.min(get(params[k])), "]"),"Mean"],
                               "sd" = summary(out)[[1]][paste0(params[k], "[", which.min(get(params[k])), "]"),"SD"],
                               "coverage" = summary(out)[[2]][paste0(params[k], "[", which.min(get(params[k])), "]"),1] <= get(params[k])[which.min(get(params[k]))] & 
                                 summary(out)[[2]][paste0(params[k], "[", which.min(get(params[k])), "]"),5] >= get(params[k])[which.min(get(params[k]))],
                               "rhat" = as.numeric(coda::gelman.diag(out[1:3][,paste0(params[k], "[", which.min(get(params[k])), "]")])[[1]][,1])))
    output <- rbind(output, 
                    data.frame("scenario" = scenario, "alternative" = altID, "simID" = simID, "iterID" = iterID, "trend" = "max",
                               "parameter" = params[k], "truth" = get(params[k])[which.max(get(params[k]))], 
                               "mean" = summary(out)[[1]][paste0(params[k], "[", which.max(get(params[k])), "]"),"Mean"],
                               "sd" = summary(out)[[1]][paste0(params[k], "[", which.max(get(params[k])), "]"),"SD"],
                               "coverage" = summary(out)[[2]][paste0(params[k], "[", which.max(get(params[k])), "]"),1] <= get(params[k])[which.max(get(params[k]))] & 
                                 summary(out)[[2]][paste0(params[k], "[", which.max(get(params[k])), "]"),5] >= get(params[k])[which.max(get(params[k]))],
                               "rhat" = as.numeric(coda::gelman.diag(out[1:3][,paste0(params[k], "[", which.max(get(params[k])), "]")])[[1]][,1])))
  }
  
}# parameter loop

output <- rbind(output, 
                data.frame("scenario" = scenario, "alternative" = altID, "simID" = simID, "iterID" = iterID, "trend" = NA,
                           "parameter" = "mu.gamma.real", "truth" = get("mu.gamma.real"), 
                           "mean" = summary(out)[[1]]["mu.gamma","Mean"],
                           "sd" = summary(out)[[1]]["mu.gamma","SD"],
                           "coverage" = summary(out)[[2]]["mu.gamma",1] <= get("mu.gamma.real") & summary(out)[[2]]["mu.gamma",5] >= get("mu.gamma.real"),
                           "rhat" = as.numeric(coda::gelman.diag(out[1:3][,"mu.gamma"])[[1]][,1])))

output <- rbind(output, 
                data.frame("scenario" = scenario, "alternative" = altID, "simID" = simID, "iterID" = iterID, "trend" = "min",
                           "parameter" = "gamma0.real", "truth" = get("gamma0.real")[which.min(get("gamma0.real"))], 
                           "mean" = summary(out)[[1]][paste0("gamma0", "[", which.min(get("gamma0.real")), "]"),"Mean"],
                           "sd" = summary(out)[[1]][paste0("gamma0", "[", which.min(get("gamma0.real")), "]"),"SD"],
                           "coverage" = summary(out)[[2]][paste0("gamma0", "[", which.min(get("gamma0.real")), "]"),1] <= get("gamma0.real")[which.min(get("gamma0.real"))] & 
                             summary(out)[[2]][paste0("gamma0", "[", which.min(get("gamma0.real")), "]"),5] >= get("gamma0.real")[which.min(get("gamma0.real"))],
                           "rhat" = as.numeric(coda::gelman.diag(out[1:3][,paste0("gamma0", "[", which.min(get("gamma0.real")), "]")])[[1]][,1])))
output <- rbind(output, 
                data.frame("scenario" = scenario, "alternative" = altID, "simID" = simID, "iterID" = iterID, "trend" = "max",
                           "parameter" = "gamma0.real", "truth" = get("gamma0.real")[which.max(get("gamma0.real"))], 
                           "mean" = summary(out)[[1]][paste0("gamma0", "[", which.max(get("gamma0.real")), "]"),"Mean"],
                           "sd" = summary(out)[[1]][paste0("gamma0", "[", which.max(get("gamma0.real")), "]"),"SD"],
                           "coverage" = summary(out)[[2]][paste0("gamma0", "[", which.max(get("gamma0.real")), "]"),1] <= get("gamma0.real")[which.max(get("gamma0.real"))] & 
                             summary(out)[[2]][paste0("gamma0", "[", which.max(get("gamma0.real")), "]"),5] >= get("gamma0.real")[which.max(get("gamma0.real"))],
                           "rhat" = as.numeric(coda::gelman.diag(out[1:3][,paste0("gamma0", "[", which.max(get("gamma0.real")), "]")])[[1]][,1])))



output <- output[-1,]
save(output, file = paste0("./output/output_scenario", scenario, "_altID", altID, "_sim", simID, "_iter", iterID, ".Rds")) 
