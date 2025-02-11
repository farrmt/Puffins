#--------------------------------#
#-Biological simulation function-#
#--------------------------------#

biosim.fun <- function(scenario = "status quo", #Population growth scenario
                       nyears = 10,             #Number of years
                       LME = "Aleutians",        #Large marine ecosystem
                       #npix = 100
                       res.meter = 10           #10 meter resolution
                       ){

  #-Libraries-#
  require(tidyverse)
  require(sf)
  require(spatstat)
  require(nimble)
  
  #spatstat.options(npixel = npix)
  
  #--------------------------#
  #-Initial colony abundance-#
  #--------------------------#
  
  if(LME == "Aleutians"){
    
    #Read Aleutain colony data and generate colony counts from range of catalog values
    
    colony.data <- readxl::read_xlsx("../data/ColonyData.xlsx", sheet = "Aleutians")
    colony.data <- colony.data %>% 
      rowwise(.) %>%
      mutate(min.count = min(c(`AMNWR Colony Catalog`, `Beringian Colony Catalog`, `Byrd Colony Size`, `USFWS 2021`), na.rm = T),
             max.count = max(c(`AMNWR Colony Catalog`, `Beringian Colony Catalog`, `Byrd Colony Size`, `USFWS 2021`), na.rm = T))
    colony.data$min.count[is.infinite(colony.data$min.count)] <- 0
    colony.data$max.count[is.infinite(colony.data$max.count)] <- 0
    
    #Number of colonies in Aleutains
    ncolonies <- dim(colony.data)[1]
    
    #Simulate initial population abundance for each colony
    Lambda0 <- NULL
    for(i in 1:ncolonies){
      Lambda0[i] <- round(runif(1, min = colony.data$min.count[i], max = colony.data$max.count[i]))
    }
  }
  
  #LME true and expected abundance
  N <- Lambda <- array(NA, dim = c(nyears, ncolonies))
  
  #Expected abundance in year 1 
  #Lambda[1,] <- Lambda0
  
  #True abundance in year 1
  #N[1,] <- rpois(ncolonies, Lambda[1,])
  
  #------------------------#
  #-Population growth rate-#
  #------------------------#
  
  #Generation time for TUPU
  generation <- 11.3
  
  #IUCN time frame
  timeframe <- round(generation * 3)
  
  #Population growth and standard deviation parameters estimated in Pearson et al. 2023
  #Variation in population growth within colony estimated with Aiktak indexplot data
  if(scenario == "status quo"){
    #Mean Population growth across colonies within an LME
    mu.gamma <- 0.094
    #Variation in population growth across colonies within an LME
    sd.gamma <- 0.052
    #Variation in population growth within a colony
    var.plot.gamma <- 0.0004514045
    #Variation in density of TUPU
    var.lambda <- 0.2488428
  }
  
  #IUCN Critically Endangered: 80% decline in 3 generations or about 38% decline in 10 years
  if(scenario == "IUCN CR"){
    #Mean Population growth across colonies within an LME
    mu.gamma <- log((1 - 0.8) ^ (1/timeframe))
    #Variation in population growth across colonies within an LME
    sd.gamma <- 0.052
    #Variation in population growth within a colony
    var.plot.gamma <- 0.0004514045
    #Variation in density of TUPU
    var.lambda <- 0.2488428
  }
  
  #IUCN Vulnerable: 30% decline in 3 generations or about 10% decline in 10 years
  if(scenario == "IUCN VU"){
    #Mean Population growth across colonies within an LME
    mu.gamma <- log((1 - 0.3) ^ (1/timeframe))
    #Variation in population growth across colonies within an LME
    sd.gamma <- 0.052
    #Variation in population growth within a colony
    var.plot.gamma <- 0.0004514045
    #Variation in density of TUPU
    var.lambda <- 0.2488428
  }
  
  #Variable growth between colonies
  if(scenario == "Variable growth between colonies"){
    #Mean Population growth across colonies within an LME
    mu.gamma <- log((1 - 0.3) ^ (1/timeframe))
    #Variation in population growth across colonies within an LME
    sd.gamma <- 0.1
    #Variation in population growth within a colony
    var.plot.gamma <- 0.0004514045
    #Variation in density of TUPU
    var.lambda <- 0.2488428
  }
  
  #Variable growth within colonies
  if(scenario == "Variable growth within colonies"){
    #Mean Population growth across colonies within an LME
    mu.gamma <-  0.094
    #Variation in population growth across colonies within an LME
    sd.gamma <- 0.052
    #Variation in population growth within a colony
    var.plot.gamma <- 0.01
    #Variation in density of TUPU
    var.lambda <- 0.2488428
  }
  
  #Variable density within colonies
  if(scenario == "Variable density within colonies"){
    #Mean Population growth across colonies within an LME
    mu.gamma <-  0.094
    #Variation in population growth across colonies within an LME
    sd.gamma <- 0.052
    #Variation in population growth within a colony
    var.plot.gamma <-  0.0004514045
    #Variation in density of TUPU
    var.lambda <- 0.5
  }
    
    #Variable growth between and within colonies and variable density
    if(scenario == "Variation in everything"){
      #Mean Population growth across colonies within an LME
      mu.gamma <-  0.094
      #Variation in population growth across colonies within an LME
      sd.gamma <- 0.1
      #Variation in population growth within a colony
      var.plot.gamma <- 0.01
      #Variation in density of TUPU
      var.lambda <- 0.5
  }
  
  #Burrow occupancy
  psi0 <- logit(0.6383885)
  
  #Colony specific population growth
  gamma0 <- exp(rnorm(ncolonies, mu.gamma, sd.gamma))
  mu.gamma.real <- mean(log(gamma0))
  
  #-------------------#
  #-Population change-#
  #-------------------#
  
  #Colonization and persistence
  phi <- 1 #persistence
  zeta <- 0.5 #colonization
  Z <- NULL #presence/absence of a nesting birds
  
  #Population parameters
  win <- gamma <- lambda <- psi <- burrow <- occ.burrow <- list()
  AREA <- gamma0.real <- NULL
  
  #Island shapefiles
  islands <- st_read("../data/Spatial_domain/Aleutians/Colonies.shp")
  
  #Population growth, abundance, and density
  for(i in 1:ncolonies){
    if(colony.data$Subsampling[i]=="Yes"){
      
      #Island id to match with colony information
      island.id <- which(islands$Name %in% colony.data$Island[i])
      island <- islands[island.id,]
      
      #Create TUPU habitat as 50 m perimeter of island
      island.hab <- st_difference(island %>% st_geometry(.),
                                  island %>% st_geometry(.) %>%
                                    st_transform(., crs = 3338) %>% 
                                    st_buffer(., dist = -50) %>%
                                    st_transform(., crs = st_crs(island)))
      #Convert to spatstat owin
      win[[paste(colony.data$Island[i])]] <- as.owin(island.hab %>% st_transform(., crs = 3338))
      #Measure area of habitat
      AREA[i] <- st_area(island.hab)
      #X pixels
      xpix <- round(as.numeric(st_bbox(island %>% st_transform(., crs = 3338))[3] - st_bbox(island %>% st_transform(., crs = 3338))[1])/res.meter)
      #Y pixels
      ypix <- round(as.numeric(st_bbox(island %>% st_transform(., crs = 3338))[4] - st_bbox(island %>% st_transform(., crs = 3338))[2])/res.meter)
      #Pixel resolution
      spatstat.options(npixel = c(xpix, ypix))
      #TUPU population growth within a colony
      gamma[[paste(colony.data$Island[i])]] <-  exp(rGRFgauss(W = win[[paste(colony.data$Island[i])]], mu = log(gamma0[i]), var = var.plot.gamma, scale = 1))
      #TUPU realized population growth of a colony
      gamma0.real[i] <- exp(mean(log(gamma[[paste(colony.data$Island[i])]]$v), na.rm = T))
      #TUPU density within a colony #UPDATE TO 100 meter squared for mu and var
      lambda[[paste(colony.data$Island[i])]][[paste0("year", 1)]] <- exp(rGRFgauss(W = win[[paste(colony.data$Island[i])]], mu = log(Lambda0[i]/AREA[i]), var = var.lambda, scale = 1))
      #TUPU expected abundance in year 1
      Lambda[1,i] <- sum(lambda[[paste(colony.data$Island[i])]][[paste0("year", 1)]]$xstep *
                         lambda[[paste(colony.data$Island[i])]][[paste0("year", 1)]]$ystep *
                         lambda[[paste(colony.data$Island[i])]][[paste0("year", 1)]]$v, na.rm = T)
      #TUPU true abundance in year 1
      N[1,i] <- rpois(1, Lambda[1,i])
      #TUPU burrow occupancy
      #psi[[paste(colony.data$Island[i])]][[paste0("year", 1)]] <- expit(rGRFgauss(W = win[[paste(colony.data$Island[i])]], mu = logit(0.635050), var = 0.019, scale = 1))
      psi[[paste(colony.data$Island[i])]] <- expit(rGRFgauss(W = win[[paste(colony.data$Island[i])]], mu = psi0, var = 0.04185783, scale = 1))
      #Expected number of burrows
      mu.burrow <- (lambda[[paste(colony.data$Island[i])]][[paste0("year", 1)]]/2) * (1/psi[[paste(colony.data$Island[i])]])#[[paste0("year", 1)]])
      #True number of burrows
      burrow[[paste(colony.data$Island[i])]][[paste0("year", 1)]] <- rpoispp(mu.burrow)
      #True number of occupied burrows
      occ.burrow[[paste(colony.data$Island[i])]][[paste0("year", 1)]] <- rthin(burrow[[paste(colony.data$Island[i])]][[paste0("year", 1)]], psi[[paste(colony.data$Island[i])]])#[[paste0("year", 1)]])
    }else{
      #Radius of island
      radius.upper <- (colony.data$`Recorded Shoreline (km)`[i]/(2 * pi)) * 1000
      if(is.na(radius.upper)|radius.upper == 0){
        radius.upper <- sqrt((colony.data$`Recorded Island Area (Ha)`[i] * 10000)/pi)
      }
      if(is.na(radius.upper)|radius.upper == 0){
        radius.upper <- sqrt((colony.data$`Derived Island Area (Acres)`[i] * 40.4686)/pi)
      }
      if(radius.upper > 50){
        #Habitat band
        radius.lower <- radius.upper - 50
        #Convert to spatstat owin
        win[[paste(colony.data$Island[i])]] <- as.owin(st_difference(st_as_sf(disc(radius = radius.upper)), st_as_sf(disc(radius = radius.lower))))
        #Habitat area band of 50 meters
        AREA[i] <- st_area(st_as_sf(win[[paste(colony.data$Island[i])]]))
      }else{
        #Convert to spatstat owin
        win[[paste(colony.data$Island[i])]] <- as.owin(st_as_sf(disc(radius = radius.upper)))
        #Habitat area band of 50 meters
        AREA[i] <- st_area(st_as_sf(disc(radius = radius.upper)))
      }
      #X and Y pixels set to spatstat default
      ypix <- xpix <- 128
      #Pixel resolution
      spatstat.options(npixel = c(xpix, ypix))
      #TUPU density within a colony
      lambda.tmp <- exp(rGRFgauss(W = win[[paste(colony.data$Island[i])]], mu = log(Lambda0[i]/AREA[i]), var = var.lambda, scale = 1))
      #TUPU expected abundance in year 1
      Lambda[1,i] <- sum(lambda.tmp$xstep *
                         lambda.tmp$ystep *
                         lambda.tmp$v, na.rm = T)
      #TUPU true abundance in year 1
      N[1,i] <- rpois(1, Lambda[1,i])
      #Available to be colonized
      Z[i] <- ifelse(N[1,i] > 0, 1, 0) #presence/absence of a nesting birds
    }
    
    for(t in 2:nyears){
      
      if(colony.data$Subsampling[i] == "Yes"){
        #Probability of colony occupancy
        omega <- Z[i] * phi + (1 - Z[i]) * zeta
        #Latent occupancy of colony
        Z[i] <- rbinom(1, 1, omega)
        
        #Colonies w/ abundance >1 in previous year
        if(N[t-1,i] > 0){
          #TUPU density within a colony
          lambda[[paste(colony.data$Island[i])]][[paste0("year", t)]] <- lambda[[paste(colony.data$Island[i])]][[paste0("year", t-1)]] * gamma[[paste(colony.data$Island[i])]]
          #TUPU expected abundance in year 1
          Lambda[t,i] <- sum(lambda[[paste(colony.data$Island[i])]][[paste0("year", t)]]$xstep *
                               lambda[[paste(colony.data$Island[i])]][[paste0("year", t)]]$ystep *
                               lambda[[paste(colony.data$Island[i])]][[paste0("year", t)]]$v, na.rm = T)
          #TUPU true abundance in year 1
          N[t,i] <- rpois(1, Lambda[t,i])
          #TUPU burrow occupancy
          #psi[[paste(colony.data$Island[i])]][[paste0("year", t)]] <- expit(rGRFgauss(W = win[[paste(colony.data$Island[i])]], mu = logit(0.635050), var = 0.019, scale = 1))
          #Expected number of burrows
          mu.burrow <- (lambda[[paste(colony.data$Island[i])]][[paste0("year", t)]]/2) * (1/psi[[paste(colony.data$Island[i])]])#[[paste0("year", t)]])
          #True number of burrows
          burrow[[paste(colony.data$Island[i])]][[paste0("year", t)]] <- rpoispp(mu.burrow)
          #True number of occupied burrows
          occ.burrow[[paste(colony.data$Island[i])]][[paste0("year", t)]] <- rthin(burrow[[paste(colony.data$Island[i])]][[paste0("year", t)]], psi[[paste(colony.data$Island[i])]])#[[paste0("year", t)]])
        }
                #Colonies occupied in current year w/ abundance =0 in previous year MTF: avialable to be occupied
        if(N[t-1,i] == 0 & Z[i] == 1){
          #TUPU density within a colony
          lambda[[paste(colony.data$Island[i])]][[paste0("year", t)]] <- (lambda[[paste(colony.data$Island[i])]][[paste0("year", t-1)]] + 2/AREA[i]) * gamma[[paste(colony.data$Island[i])]]
          #TUPU expected abundance in year 1
          Lambda[t,i] <- sum(lambda[[paste(colony.data$Island[i])]][[paste0("year", t)]]$xstep *
                               lambda[[paste(colony.data$Island[i])]][[paste0("year", t)]]$ystep *
                               lambda[[paste(colony.data$Island[i])]][[paste0("year", t)]]$v, na.rm = T)
          #TUPU true abundance in year 1
          N[t,i] <- rpois(1, Lambda[t,i])
          #TUPU burrow occupancy
          #psi[[paste(colony.data$Island[i])]][[paste0("year", t)]] <- expit(rGRFgauss(W = win[[paste(colony.data$Island[i])]], mu = logit(0.635050), var = 0.019, scale = 1))
          #Expected number of burrows
          mu.burrow <- (lambda[[paste(colony.data$Island[i])]][[paste0("year", t)]]/2) * (1/psi[[paste(colony.data$Island[i])]])#[[paste0("year", t)]])
          #True number of burrows
          burrow[[paste(colony.data$Island[i])]][[paste0("year", t)]] <- rpoispp(mu.burrow)
          #True number of occupied burrows
          occ.burrow[[paste(colony.data$Island[i])]][[paste0("year", t)]] <- rthin(burrow[[paste(colony.data$Island[i])]][[paste0("year", t)]], psi[[paste(colony.data$Island[i])]])#[[paste0("year", t)]])
        }
        #Colonies unoccupied in current year w/ abundance =0 in previous year
        if(N[t-1,i] == 0 & Z[i] == 0){
          #TUPU density within a colony
          lambda[[paste(colony.data$Island[i])]][[paste0("year", t)]] <- lambda[[paste(colony.data$Island[i])]][[paste0("year", t-1)]] * gamma[[paste(colony.data$Island[i])]]
          #TUPU expected abundance in year 1
          Lambda[t,i] <- sum(lambda[[paste(colony.data$Island[i])]][[paste0("year", t)]]$xstep *
                               lambda[[paste(colony.data$Island[i])]][[paste0("year", t)]]$ystep *
                               lambda[[paste(colony.data$Island[i])]][[paste0("year", t)]]$v, na.rm = T)
          #TUPU true abundance in year 1
          N[t,i] <- 0
          #Expected number of burrows
          mu.burrow <- (lambda[[paste(colony.data$Island[i])]][[paste0("year", t)]]/2) * (1/psi[[paste(colony.data$Island[i])]])#[[paste0("year", t)]])
          #True number of burrows
          burrow[[paste(colony.data$Island[i])]][[paste0("year", t)]] <- rpoispp(mu.burrow)
          #True number of occupied burrows
          occ.burrow[[paste(colony.data$Island[i])]][[paste0("year", t)]] <- 0
        }}else{
        #Probability of colony occupancy
        omega <- Z[i] * phi + (1 - Z[i]) * zeta
        #Latent occupancy of colony
        Z[i] <- rbinom(1, 1, omega)
        
        #Colonies w/ abundance >1 in previous year
        if(N[t-1,i] > 0){
          Lambda[t,i] <- Lambda[t-1,i] * gamma0[i]
          N[t,i] <- rpois(1, Lambda[t,i])
          # N[t,i] <- rpois(1, lambda = Z[i] * N[t-1,i] * gamma[i])
        }
        #Colonies occupied in current year w/ abundance =0 in previous year MTF: avialable to be occupied
        if(N[t-1,i] == 0 & Z[i] == 1){
          Lambda[t,i] <- (Lambda[t-1,i] + 2) * gamma0[i]
          N[t,i] <- rpois(1, Lambda[t,i])
          #N[t,i] <- rpois(1, lambda = (Lambda0[i] + 1) * gamma[i]^(t-1))
        }
        #Colonies unoccupied in current year w/ abundance =0 in previous year
        if(N[t-1,i] == 0 & Z[i] == 0){
          Lambda[t,i] <- Lambda[t-1,i] * gamma0[i]
          N[t,i] <- 0
        }
      }
    }
  }
  
  #####
  # COS <- Lambda/rpois(ncolonies, lambda = 0.4437426 * 2.5^2 * pi) # Multiplicative change between colony and site
  # AREA <- COS * (2.5^2 * pi) 
  # maxsite <- floor(COS)
  # nsites <- ifelse(maxsite < 100, maxsite, 100)
  # n <- array(NA, dim = c(nyears, ncolonies, 100))
  # 
  # for(i in 1:ncolonies){
  #   for(t in 1:nyears){
  #     if(N[t,i] == 0){
  #       n[t,i,] <- 0
  #     }else{
  #       for(j in 1:nsites[i]){
  #         n[t,i,j] <- rpois(1, N[t,i]/COS[i])
  #         
  #         lambda0[t,i,j] <- log(Lambda[t,i]/AREA[i])
  #         lambda[t,i,j] <- exp(rGRFgauss(W = win, mu = lambda0, var = 0.260475, scale = 1))
  #         psi[t,i,j]
  #       }
  #       if(sum(n[t,i,1:nsites[i]]) > N[t,i]){
  #         N[t,i] <- sum(n[t,i,1:nsites[i]])
  #       }
  #     }
  #   }
  # }
  #####

return(list(Lambda0 = Lambda0,             #Initial expected abundance
            Lambda = Lambda,               #Expected abundance
            mu.gamma = mu.gamma,           #Mean expected population growth (log-scale)
            mu.gamma.real = mu.gamma.real, #Realized mean expected population growth (log-scale)
            gamma0 = gamma0,               #Colony-specific population growth
            gamma0.real = gamma0.real,     #Realized colony-specific population growth
            gamma = gamma,                 #Within colony population growth
            sd.gamma = sd.gamma,           #SD of expected population growth
            N = N,                         #True colony abundance
            AREA = AREA,                   #Habitat area of a colony
            lambda = lambda,               #Expected colony density
            psi0 = psi0,                   #Mean burrow occupancy logit-scale
            psi = psi,                     #Burrow occupancy
            burrow = burrow,               #True number of burrows
            occ.burrow = occ.burrow,       #True number of occupied burrows
            nyears = nyears,               #Number of years
            ncolonies = ncolonies,         #Number of colonies
            colony.data = colony.data      #Colony data
            ))

}#end function
