#-------------#
#-Environment-#
#-------------#

envr <- environment()

#----------#
#-Metadata-#
#----------#

load(file = "metadata.Rds")

scenario <- head(metadata$scenario, n = 1)
altID <- head(metadata$altID, n = 1)
simID <- head(metadata$simID, n = 1)
iterID <- head(metadata$iterID, n = 1)

metadata <- metadata[-1,]
save(metadata, file = "metadata.Rds")

#-----------#
#-Load data-#
#-----------#

path <- paste0("./alternative/scenario", scenario, "_alt", altID, "_sim",  simID, "_iter", iterID, ".Rds")
load(path)
list2env(output, envir = envr)

#-----------#
#-Run model-#
#-----------#

sys.source("SubsamplingModel.R", envir = envr, toplevel.env = envr)

#-------------#
#-Save Output-#
#-------------#

output <- data.frame(matrix(NA, nrow = 1, ncol = 10))
colnames(output) <- c("scenario", "alternative", "simID", "iterID", "trend", "parameter", "truth", "mean", "sd", "rhat")

for(k in 1:length(params)){
  if(length(get(params[k])) == 1){
    output <- rbind(output, 
                    data.frame("scenario" = scenario, "alternative" = altID, "simID" = simID, "iterID" = iterID, "trend" = NA,
                               "parameter" = params[k], "truth" = get(params[k]), 
                               "mean" = summary(out)[[1]][params[k],"Mean"],
                               "sd" = summary(out)[[1]][params[k],"SD"],
                               "rhat" = as.numeric(coda::gelman.diag(out[1:3][,params[k]])[[1]][,1])))
  }else{
    output <- rbind(output, 
                    data.frame("scenario" = scenario, "alternative" = altID, "simID" = simID, "iterID" = iterID, "trend" = "min",
                               "parameter" = params[k], "truth" = get(params[k])[which.min(gamma0)], 
                               "mean" = summary(out)[[1]][paste0(params[k], "[", which.min(gamma0), "]"),"Mean"],
                               "sd" = summary(out)[[1]][paste0(params[k], "[", which.min(gamma0), "]"),"SD"],
                               "rhat" = as.numeric(coda::gelman.diag(out[1:3][,paste0(params[k], "[", which.min(gamma0), "]")])[[1]][,1])))
    output <- rbind(output, 
                    data.frame("scenario" = scenario, "alternative" = altID, "simID" = simID, "iterID" = iterID, "trend" = "max",
                               "parameter" = params[k], "truth" = get(params[k])[which.max(gamma0)], 
                               "mean" = summary(out)[[1]][paste0(params[k], "[", which.max(gamma0), "]"),"Mean"],
                               "sd" = summary(out)[[1]][paste0(params[k], "[", which.max(gamma0), "]"),"SD"],
                               "rhat" = as.numeric(coda::gelman.diag(out[1:3][,paste0(params[k], "[", which.max(gamma0), "]")])[[1]][,1])))
  }
  
}# parameter loop


output <- output[-1,]
save(output, file = paste0("./output/output_scenario", scenario, "_altID", altID, "_sim", simID, "_iter", iterID, ".Rds")) 
