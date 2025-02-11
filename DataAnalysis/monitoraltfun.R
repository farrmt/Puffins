#--------------------------#
#-Sampling design function-#
#--------------------------#

design.fun <- function(colonies,                       #Colony ID
                       core.units = "Aiktak",          #Core colony ID
                       subsampling,                    #Subsampling (yes or no)
                       #islandgroup,
                       unit.selection = "random",      #Colony selection
                       site.selection = "wholecolony", #Site selection
                       sampling.type,                  #Sampling type
                       scale.selection,                #Scale selection
                       nunits = 20,                    #Number of units (colonies) to sample
                       ncore = 0,                      #Number of core units
                       nsubsampled = 0,                #Number of subsampled colonies
                       nsites = 0,                     #Number of sites
                       nreps = 1                       #Sampling intensity
                       ){
  
  
  require(spatstat)
  require(sf)
  require(tidyverse)
  
  #-----------------------------------#
  #-Observation parameters and arrays-#
  #-----------------------------------#
  
  #Observation error around boat-based colony counts
  #theta.sd <- runif(1, 20, 30) #ADD realism
  colonycount.p <- 1
  
  #Detection probability of counting burrows
  #p <- runif(1, 0.7, 0.9)
  burrowcount.p <- 1
  
  #Detection probability of burrow occupancy
  burrowocc.p <- 1
  
  #Colony counts array
  Y <- array(NA, dim = c(nyears, nunits, nreps))
  
  #Burrow counts, occupied burrow counts arrays
  y.burrow <- y.occ <- array(NA, dim = c(nyears, nunits, max(nsites)))
  
  #Unit names, unit ID, change-of-support, and habitat area arrays
  unit.names <- units <- COS <- AREA.obs <- array(NA, dim = c(nyears, nunits))
  
  #Site location and arrays
  if(any(scale.selection == "site")){
    site.xy <- array(NA, dim = c(nunits, max(nsites)))
    site <- list()
  }
  
  
  #-------------------------#
  #-Unit (colony) selection-#
  #-------------------------#
  
  #Core unit selection
  if(unit.selection == "splitpanel"){
    # if(any(is.na(core.units))){
    if(length(core.units) < ncore){
      if(any(is.na(core.units))){
        core <- sample(which(subsampling == "Yes"), ncore)
      }else{
        needed <- ncore - length(core.units)
        core <- c(which(colonies %in% core.units),
                  sample(which(subsampling == "Yes")[-(which(colonies %in% core.units))], needed))
        #which(subsampling == "Yes")[order(N[1,which(subsampling == "Yes")], decreasing = T)][1:needed]
      }
    }else{
      core <- which(colonies %in% core.units)
    }
  }
  
  #Unit selection for each year
  for(t in 1:nyears){
    
    #Random selection
    if(unit.selection == "random"){
      if(all(scale.selection == "unit")){
        units[t,] <- sample(1:ncolonies, nunits)
      }
      else{
        if(all(scale.selection == "site")){
          units[t,] <- sample(which(subsampling == "Yes"), nunits) 
        }
        else{
          units.tmp <- sample(which(subsampling == "Yes"), nsubsampled)
          units[t,] <- c(units.tmp, sample((1:ncolonies)[-units.tmp], nunits - nsubsampled))
        }
      }
    }
    
    #Stratified selection
    if(unit.selection == "stratified"){
      if(all(scale.selection == "unit")){
        units[t,] <- c(sample(order(N[t,], decreasing = T)[1:round(ncolonies*0.25)], nunits - 3*round(nunits*0.25)),
                       sample(order(N[t,], decreasing = T)[(round(ncolonies*0.25) + 1):(round(ncolonies*0.25)*2)], round(nunits*0.25)),
                       sample(order(N[t,], decreasing = T)[(round(ncolonies*0.25)*2 + 1):(round(ncolonies*0.25)*3)], round(nunits*0.25)),
                       sample(order(N[t,], decreasing = T)[(round(ncolonies*0.25)*3 + 1):ncolonies], round(nunits*0.25)))
      }
      else{
        if(all(scale.selection == "site")){
          units[t,] <- c(sample(which(subsampling == "Yes")[order(N[t,which(subsampling == "Yes")], decreasing = T)][1:round(sum(subsampling == "Yes")*0.25)], nunits - 3*round(nunits*0.25)),
                         sample(which(subsampling == "Yes")[order(N[t,which(subsampling == "Yes")], decreasing = T)][(round(sum(subsampling == "Yes")*0.25) + 1):(round(sum(subsampling == "Yes")*0.25)*2)], round(nunits*0.25)),
                         sample(which(subsampling == "Yes")[order(N[t,which(subsampling == "Yes")], decreasing = T)][(round(sum(subsampling == "Yes")*0.25)*2 + 1):(round(sum(subsampling == "Yes")*0.25)*3)], round(nunits*0.25)),
                         sample(which(subsampling == "Yes")[order(N[t,which(subsampling == "Yes")], decreasing = T)][(round(sum(subsampling == "Yes")*0.25)*3 + 1):sum(subsampling == "Yes")], round(nunits*0.25)))
        }
        else{
          units[t,] <- c(sample(which(subsampling == "Yes")[order(N[t,which(subsampling == "Yes")], decreasing = T)][1:round(ncolonies*0.25)], nunits - 3*round(nunits*0.25)), #Update to not over or under subsample
                         sample(order(N[t,], decreasing = T)[(round(ncolonies*0.25) + 1):(round(ncolonies*0.25)*2)], round(nunits*0.25)),
                         sample(order(N[t,], decreasing = T)[(round(ncolonies*0.25)*2 + 1):(round(ncolonies*0.25)*3)], round(nunits*0.25)),
                         sample(order(N[t,], decreasing = T)[(round(ncolonies*0.25)*3 + 1):ncolonies], round(nunits*0.25)))
        }
      }
    }
    
    # #Island group
    # if(unit.selection == "islandgroup"){
    #   if(all(scale.election == "unit")){
    #     stop()
    #   }else{
    #     if(all(scale.selection == "site")){
    #       stop()
    #     }else{
    #       units[t,] <- sample(which(islandgroup == group[t] & subsampling == "Yes"))
    #     }
    #   }
    # }
    
    #Split-panel selection
    if(unit.selection == "splitpanel"){
      if(all(scale.selection == "unit")){
        if(ncore == nunits){
          units[t,] <- core
        }else{
          units[t,] <- c(core, sample((1:ncolonies)[-core], nunits - ncore))
        }
      }
      else{
        if(all(scale.selection == "site")){
          if(ncore == nunits){
            units[t,] <- core
          }else{
            units[t,] <- c(core, sample(which(subsampling == "Yes")[!which(subsampling == "Yes") %in% core], nunits - ncore))
          }
        }
        else{
          if(ncore == nunits){
            units[t,] <- core
          }else{
            units.tmp <- sample(which(subsampling == "Yes")[!which(subsampling == "Yes") %in% core], nsubsampled - ncore)
            units[t,] <- c(core, units.tmp, sample((1:ncolonies)[-c(core, units.tmp)], nunits - nsubsampled))
          }
        }
      }

    }
    
    #Names of selected colonies
    unit.names[t,] <- colonies[units[t,]]
    
  }
  
  #Reordering of IDs based on sampling
  if(sum(!(1:ncolonies %in% units)) > 0){
    
    #Number of sampled colonies out of all possible colonies
    nsampled <- sum(1:ncolonies %in% units)
    
    #Colony IDs of those not sampled
    notsampled <- which(!(1:ncolonies %in% units))
    
    #Remove unsampled colonies
    N <- N[,-notsampled]
    Lambda <- Lambda[,-notsampled]
    Lambda0 <- Lambda0[-notsampled]
    gamma0 <- gamma0[-notsampled]
    
    #Reorder IDs
    for(i in 1:nunits){
      for(t in 1:nyears){
        for(k in length(notsampled):1){
          if(units[t,i] > notsampled[k]){
            units[t,i] <- units[t,i] - 1
          }
        }
      }
    }
  }
  
  #---------------------#
  #-Observation process-#
  #---------------------#
  
  for(t in 1:nyears){
    for(i in 1:nunits){
      
      #Sampled colonies habitat area
      AREA.obs[t,i] <- AREA[units[t,i]]
      
      #Subsampling observation process
      #if(scale.selection == "site"|(scale.selection == "mixed" & i <= nsubsampled)){
        
      if(scale.selection[i] == "site"){
      
      #Density of observed colony
      lambda.tmp <- lambda[[unit.names[t,i]]][[t]]
      
      #True number of burrows of observed colony
      burrow.tmp <- burrow[[unit.names[t,i]]][[t]]
      
      #True number of occupied burrows of observed colony
      occ.burrow.tmp <- occ.burrow[[unit.names[t,i]]][[t]]
      
      #Burrow occupancy of observed colony
      psi.tmp <- psi[[unit.names[t,i]]]
      
      #Population groth of observed colony
      gamma.tmp <- gamma[[unit.names[t,i]]]
      
      #Preferential sampling: probability of site based on higher TUPU density
      #prob.site <- ((lambda.tmp/2) * (1/psi.tmp))/sum((lambda.tmp/2) * (1/psi.tmp))
      prob.site <- ((lambda.tmp * gamma.tmp^10/2) * (1/psi.tmp))/sum(((lambda.tmp * gamma.tmp^10/2) * (1/psi.tmp)))
      prob.site$v[is.na(prob.site$v)] <- 0
      
      #Random sampling: probability is random
      prob.ran <- prob.site
      prob.ran$v[prob.ran$v > 0] <- 1/sum(prob.ran$v > 0)
      
      #Index plot sampling
      #if(sampling.type == "indexplot"|sampling.type == "indexplot & colony.count"){
      if(sampling.type[i] == "indexplot"){
      
      #Change-of-support (315 m^2)
      COS[t,i] <- AREA[units[t,i]]/315
      
      #Preferential sampling site location (sites stay constant)
      # if(site.selection == "preferential" & t == 1){
      if(site.selection[i] == "preferential" & t == 1){
        site.xy[i,1:nsites[i]] <- sample(1:prod(dim(prob.site)), size = nsites[i], prob = prob.site$v)
      }
      
      #Random sampling site location (sites stay constant)
      # if(site.selection == "random" & t == 1){
      if(site.selection[i] == "random" & t == 1){
        site.xy[i,1:nsites[i]] <- sample(1:prod(dim(prob.ran)), size = nsites[i], prob = prob.ran$v) #This should be constant across years.
      }
      
      #Burrow counts and occupancy
      for(j in 1:nsites[i]){
        
        #Site location
        site[[paste(unit.names[t,i])]][[paste0("year", t)]][[paste0("indexplot", j)]] <- st_point(c(prob.site$xcol[which(prob.site$v == prob.site$v[site.xy[i,j]], arr.ind = T)[2]], 
                                                                                                    prob.site$yrow[which(prob.site$v == prob.site$v[site.xy[i,j]], arr.ind = T)[1]])) %>%
          st_sfc(.) %>% st_set_crs(3338) %>% st_buffer(., dist = sqrt(315)/2) %>% st_bbox(.) %>% st_as_sfc(.) %>% st_as_sf(.)
        
        #Burrow counts (perfect detection) 
        y.burrow[t,i,j] <- burrow.tmp[as.owin(site[[paste(unit.names[t,i])]][[paste0("year", t)]][[paste0("indexplot", j)]])]$n
        
        #Burrow occupancy (perfect detection)
        y.occ[t,i,j] <- occ.burrow.tmp[as.owin(site[[paste(unit.names[t,i])]][[paste0("year", t)]][[paste0("indexplot", j)]])]$n
        
        # for(k in 1:nreps){
        #   y[t,i,j,k] <- rbinom(1,burrow.tmp[as.owin(site[[paste(unit.names[t,i])]][[paste0("year", t)]][[paste0("indexplot", j)]])]$n, prob = p)
        #   occ[t,i,j,k] <- occ.burrow.tmp[as.owin(site[[paste(unit.names[t,i])]][[paste0("year", t)]][[paste0("indexplot", j)]])]$n
        # }
      }}
      
      #Circle plot sampling
      # if(sampling.type == "circleplot"|sampling.type == "colonycount & circleplot"){
      if(sampling.type[i] == "circleplot"){
        
        #Change-of-support (~19 m^2)
        COS[t,i] <- AREA[units[t,i]]/(pi * 2.5^2)
        
        # if(site.selection == "preferential"){
        #   site.xy[i,] <- sample(1:prod(dim(prob.site)), size = nsites, prob = prob.site$v)
        # }
        
        #Random sampling site location (sites rotate each year)
        # if(site.selection == "random"){
        if(site.selection[i] == "random"){
          site.xy[i,1:nsites[i]] <- sample(1:prod(dim(prob.site)), size = nsites[i], prob = NULL)
        }
        
        #Burrow counts and occupancy
        for(j in 1:nsites[i]){
          
          #Site location
          site[[paste(unit.names[t,i])]][[paste0("year", t)]][[paste0("circleplot", j)]] <- st_point(c(prob.site$xcol[which(prob.site$v == prob.site$v[site.xy[i,j]], arr.ind = T)[2]], 
                                                                                                      prob.site$yrow[which(prob.site$v == prob.site$v[site.xy[i,j]], arr.ind = T)[1]])) %>%
            st_sfc(.) %>% st_set_crs(3338) %>% st_buffer(., dist = 2.5) %>% st_as_sf(.)
          
          #Burrow counts (perfect detection)
          y.burrow[t,i,j] <- burrow.tmp[as.owin(site[[paste(unit.names[t,i])]][[paste0("year", t)]][[paste0("circleplot", j)]])]$n
          
          #Burrow occupancy (perfect detection)
          y.occ[t,i,j] <- occ.burrow.tmp[as.owin(site[[paste(unit.names[t,i])]][[paste0("year", t)]][[paste0("circleplot", j)]])]$n
          
          # for(k in 1:nreps){
          #   y[t,i,j,k] <- rbinom(1,burrow.tmp[as.owin(site[[paste(unit.names[t,i])]][[paste0("year", t)]][[paste0("circleplot", j)]])]$n, prob = p)
          #   occ[t,i,j,k] <- occ.burrow.tmp[as.owin(site[[paste(unit.names[t,i])]][[paste0("year", t)]][[paste0("circleplot", j)]])]$n
          # }
          
        }}#ADD QUADRAT PLOTS HERE
      }
      
      #Boat-based counts
      # if(scale.selection == "unit"|(scale.selection == "mixed" & i > nsubsampled)){
      if(scale.selection[i] == "unit"){
        # if(sampling.type == "colonycount"|sampling.type == "colonycount & circleplot"){
        if(sampling.type[i] == "colonycount"){
          for(k in 1:nreps){
            #Colony counts
            Y[t,i,k] <- rbinom(1,N[t,units[t,i]], prob = colonycount.p)
          }
        }
      }
    }
  }
  
  return(list(Y = Y,                 #Boat-based counts
              y.burrow = y.burrow,   #Observed burrow counts
              y.occ = y.occ,         #Observed occupied burrows
              Lambda0 = Lambda0,     #Initial expected colony abundance
              Lambda = Lambda,       #Expected colony abundance
              N = N,                 #True colony abundance
              gamma0 = gamma0, #Colony-specific population growth (log-scale)
              AREA.obs = AREA.obs,   #Colony habitat area
              COS = COS,             #Change-of-support
              units = units,         #Unit (colony) ID
              nyears = nyears,       #Number of years
              ncolonies = nsampled,  #Number of unique colonies sampled across years
              colonies = unit.names, #Colony names
              nunits = nunits,       #Number of colonies sampled per year
              nreps = nreps,         #Number of replicates
              nsites = nsites,       #Number of sites
              site = site,           #Site location
              #theta.sd = theta.sd,   #Boat-based count observation error
              colonycount.p = colonycount.p,
              #p = p                  #Detection probability of burrows
              burrowcount.p = burrowcount.p,
              burrowocc.p = burrowocc.p
              ))
}
  