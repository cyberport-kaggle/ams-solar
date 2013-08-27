library(data.table)
library(reshape2)
library(foreach)
library(doMC)
library(profr)
registerDoMC(cores = 3)


buildMesonet <- function(train){
  data = list()
  stationInfo <- read.csv("data/station_info.csv", stringsAsFactors = FALSE)
  stationInfo$elon = stationInfo$elon + 360
  factors = c('dswrf_sfc','dlwrf_sfc','uswrf_sfc','ulwrf_sfc','ulwrf_tatm','pwat_eatm','tcdc_eatm','apcp_sfc',
              'pres_msl','spfh_2m','tcolc_eatm','tmax_2m','tmin_2m','tmp_2m','tmp_sfc')
  
  if(train){
    loadPath = file.path("data", "train")
    savePath = "data/cleaned/"
    print("Processing training set")
  }
  else{
    loadPath = file.path("data", "test")
    savePath = "data/cleanedTest/"
    print("Processing testing set")
  }
  
  for(stn in stationInfo$stid){
    cat("Loading data for", stn, "\n")
    data <- foreach (f = factors) %dopar% {
      #     for(f in factors){
      print(f)
      lon <- floor(stationInfo[stationInfo$stid == stn,]$elon)
      lat <- floor(stationInfo[stationInfo$stid == stn,]$nlat)
      
      lonidx <- (-1:2) + lon
      latidx <- (-1:2) + lat
      
      locs <- expand.grid(lon = lonidx, lat = latidx)
      latidx <- locs$lat
      lonidx <- locs$lon
      
      path <- Sys.glob(file.path(loadPath, paste0(f, "*", ".Rdata")))
      load(path)
      locs <- data.table(lon = lonidx, lat = latidx, key = c("lat","lon"))
      
      tbl <- tbl[locs, list(mean = mean(value), sd = sd(value)), by = list(date, lon, lat, hour)]
      tbl <- dcast(melt(tbl, id.vars=c("date", "lon", "lat", "hour")), 
                   date ~ variable + lon + lat + hour)
      
      names(tbl)[-1] <- paste(f, names(tbl)[-1], sep="_")
      data[[f]] <- data.table(tbl, key = "date")
    }
    
    tbl <- Reduce(merge, data)
    save(tbl, file = paste0(savePath, stn, ".RData"))
  }
}