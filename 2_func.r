stationFit <- function(stn, trainingPath=paste0(dataFolder, 'cleaned/')) {
    cat('Fitting model for', stn, '\n', sep=' ')
    f <- paste0(trainingPath, stn, '.RData')
    values <- get(load(f))
    values <- values[, -1]
    
    #perform PCA     
#     pcaFit = princomp(values)
#     values = pcaFit$scores[,1:10]
    
    y <- trainData[,stn]
    trainDf <- data.frame(y=y, values)
    fitCtrl <- trainControl(method = "cv",
                            number = 5)
    trainGrid = expand.grid(.mtry = seq(2,102, 20))
#     stnFit <- train(y ~ .,
#                     data=trainDf,
#                     method="rf",
#                     trControl=fitCtrl,
#                     tuneGrid = trainGrid,
#                     verbose=TRUE,
#                     ntree=500)

    stnFit <- randomForest(y~., data=trainDf, do.trace=TRUE, mtry=22, ntree=250)
    print(stnFit)
    save(stnFit, file=paste0('./data/RFmodels/', stn, '.Rdata'))
}

boostFit <- function(stn, trainingPath=paste0(dataFolder, 'cleaned/')) {
#   Caret tuning results min RMSE = 3050000
#   Tuning parameter 'shrinkage' was held constant at a value of 0.1
#   RMSE was used to select the optimal model using  the smallest value.
#   The final values used for the model were interaction.depth = 8, n.trees = 100 and shrinkage = 0.1. 
  
  cat('Fitting model for', stn, '\n', sep=' ')
  f <- paste0(trainingPath, stn, '.RData')
  values <- get(load(f))
  values <- values[, -1]
  
  #perform PCA     
  #     pcaFit = princomp(values)
  #     values = pcaFit$scores[,1:10]
  
  y <- trainData[,stn]
  trainDf <- data.frame(y=y, values)
  fitCtrl <- trainControl(method = "cv",
                          number = 2)
#   browser()
#   stnFit <- train(y ~ .,
#                   data=trainDf,
#                   method="gbm",
#                   trControl=fitCtrl,
#                   tuneLength = 10)
  
  stnFit <- gbm(y~., data=trainDf, shrinkage=0.1, interaction.depth=8, n.trees=100, verbose = TRUE)

  print(stnFit)
  save(stnFit, file=paste0('./data/GBMmodels/', stn, '.Rdata'))
}

svmFit <- function(stn, trainingPath=paste0(dataFolder, 'cleaned/')) {
  cat('Fitting model for', stn, '\n', sep=' ')
  f <- paste0(trainingPath, stn, '.RData')
  values <- get(load(f))
  values <- values[, -1]
  
  #perform PCA     
  pcaFit = princomp(values)
  #     values = pcaFit$scores[,1:10]
  
  y <- trainData[,stn]
  trainDf <- data.frame(y=y, values)
  fitCtrl <- trainControl(method = "cv",
                          number = 5)
  stnFit <- train(y ~ .,
                  data=trainDf,
                  method="svmRadial",
                  trControl=fitCtrl,
                  tuneLength = 5,
                  verbose=TRUE)
  
#   stnFit <- randomForest(y~., data=trainDf, do.trace=TRUE, mtry=10, ntree=500)
  
  print(stnFit)
  #     save(stnFit, pcaFit, file=paste0('../data/models/', stn, '.Rdata'))
}

predictStation <- function(stn, modelsFolder=modelsFolder) {
    # Given a station name stn, load the model, load the testing data, and use the model to predict.
    # returns an array of values, one item for each date in the testing data
    mdl <- get(load(paste0(dataFolder, modelsFolder, stn, '.Rdata')))
    rawTest <- get(load(paste0(dataFolder, testOutFolder, stn, '.RData')))
    dates <- rawTest$date
    testData <- rawTest[,c(-1)]
    
#     apply PCA
#     testData <- predict(mdl[[2]], testData)[,1:10]
    
    if(class(mdl) == "gbm") pred <- predict(mdl, testData, n.trees = 100)
    else pred <- predict(mdl, testData)
    
    res <- data.frame(date=dates)
    res[,stn] <- pred
    return(res)
}

selectiveRF <- function(stations=stationNames, outputFolder='selectiveRF/') {
  # Factor selection:
  # This is a more selective RF model, with only the factors that RF indicated were most important
  # Our cutoff was arbitrary, but we basically took the factors that ended up being in the top 50 most important factors
  # This was basically upSWFlux, dwnSWFlux, upLWFluxSurf, upLWFlux, tempSurf, cloud
  # 
  # Ensembles:
  # Some analysis showed that ensembles probably contained much more information than we originally thought
  # So instead of just taking the mean of the ensembles, this one also includes the min, max, and standard deviation
  # 
  # Hours:
  # It seems to make sense to combine hours for some factors, but not for others
  # For the flux measures, we'll scale (from W to a three hour period measure) and sum over the whole day
  # For temperature, we will not combine hours
  # For cloud cover, we will not combine hours
  #
  # Lat/Lon:
  # We'll use closes four points
  #
  # Total number of factors should be:
  # Flux measures: 4 factors * 4 points * 4 ensemble statistics = 64
  # Tempreature/cloud: 2 factors * 4 points * 5 hours * 4 ensemble statistics = 160 
  # Total = 224 factors per station
  
  #####
  # Subroutines
  #####
  processSingleStation <- function(thisStn, combine=TRUE) {
    cat('Processing for station', thisStn$stid, '\n')
    pts <- getPoints(thisStn$elon, thisStn$nlat)
    selector <- CJ(pts$lat,
                   pts$lon)
    # This subset contains all data points for the GEFS locations we're interested in
    stnTbl <- tbl[selector]
    
    # Further reduce the data as described above
    if (combine) {
      # If combine is true, then it means the value is a flux value, and we'll sum over hours and scale
      stnTbl <- stnTbl[, list(value=sum(value * 3600 * 3)), by='lat,lon,ens,date']
      stnTbl <- stnTbl[, list(min=min(value), max=max(value), sd=sd(value), mean=mean(value)), by='lat,lon,date']
      stnTbl <- data.table(cast(melt(stnTbl, id.var=c('lat', 'lon', 'date')),
                date ~ variable + lat + lon, value.var = 'value'))
    } else {
      stnTbl <- stnTbl[, list(min=min(value), max=max(value), sd=sd(value), mean=mean(value)), by='lat,lon,hour,date']
      stnTbl <- data.table(cast(melt(stnTbl, id.var=c('lat', 'lon', 'hour', 'date')),
                date ~ variable + lat + lon + hour, value.var = 'value'))
    }
    # Rename the variables to attach the variable name
    for (nm in names(stnTbl)) {
      if (nm != 'date') {
        setnames(stnTbl, nm, paste(thisVarName, nm, sep='_'))
      }
    }
    # Add a week of year variable
    stnTbl[, date:=date/100]
    return(stnTbl)  
  }
  
  mergeStationData <- function(stn, res, obj) {
    if (is.null(obj[[stn]])) {
      obj[[stn]] <- res
    } else {
      obj[[stn]] <- merge(obj[[stn]], res, by='date')
    }  
    return(obj)
  }
  
  ######
  # Execution
  ######
  
  # files we're interested in:
  fluxFiles <- trainRData[c(3, 13, 14, 15)]
  otherFiles <- trainRData[c(7, 12)]
  trainDts <- list()
  results <- list()
  # If a subset of stations are passed in, we have to get the indexes of the stations
  stationIndices <- which(stationNames %in% stations)
  
  # Load a flux file
  for (thisFile in fluxFiles) {
    cat('Loading data for factor', thisFile, '\n')
    load(file.path(dataFolder, trainFolder, thisFile))
    thisVarName <- trainFileToShortName[[thisFile]]
    
    # Process all stations
    stnDts <- foreach(i=stationIndices) %dopar% {
      thisStn <- stationInfo[i]
      return(processSingleStation(thisStn))
    }
    # Have to combine out here, since can't do it in the parallelized code
    for (i in 1:length(stnDts)) {
      thisStn <- stationInfo[stationIndices[i]]$stid
      trainDts <- mergeStationData(thisStn, stnDts[[i]], trainDts)
    }
    rm(tbl)
    gc()
  }
  
  # non-flux files
  for (thisFile in otherFiles) {
    cat('Loading data for factor', thisFile, '\n')
    load(file.path(dataFolder, trainFolder, thisFile))
    thisVarName <- trainFileToShortName[[thisFile]]
    
    # Process all stations
    stnDts <- foreach(i=stationIndices) %dopar% {
      thisStn <- stationInfo[i]
      return(processSingleStation(thisStn, combine=FALSE))
    }
    # Have to combine out here, since can't do it in the parallelized code
    for (i in 1:length(stnDts)) {
      thisStn <- stationInfo[stationIndices[i]]$stid
      trainDts <- mergeStationData(thisStn, stnDts[[i]], trainDts)
    }
    rm(tbl)
    gc()
  }
  # At this point, trainDts is a list with length == number of stations
  # Each list item is a stations' train dataset
  
  # Train each station
  for (i in 1:length(trainDts)) {
    thisStationName <- names(trainDts)[i]
    cat('Training model for station', thisStationName, '\n')
    
    # Annoyingly the date is formatted differently in the test data
    thisTrainDt <- merge(trainData[,c('date', thisStationName), with=FALSE], trainDts[[thisStationName]], by='date')
    setnames(thisTrainDt, thisStationName, 'y')
    
    # Add in a day of year kind of measure
    # However RF cannot take more than 32 levels in a categorical variable
    # Instead maybe calculate days from peak throughput (i.e. the summer solstice) so it can be continuous
    thisTrainDt[, daysFromMax := yday(as.IDate(as.character(date), format="%Y%m%d")) - yday(IDate('2013-06-20'))]
#     fitCtrl <- trainControl(method = "cv",
#                             number = 5)
#     stnFit <- train(y ~ .,
#                     data=thisTrainDt,
#                     method="rf",
#                     trControl=fitCtrl,
    t
#                     verbose=TRUE,
#                     ntree=500)
#     # Finds mtry=113 is optimal
    # Performance really stops improving at around 400 trees
    rf <- randomForest(y~., data=thisTrainDt, mtry=113, do.trace=TRUE, ntree=400)
    # Note that parallelized RF doesn't preserve the MSE
    rf <- foreach(i=rep(100, 4), .combine=combine, .multicombine=TRUE) %dopar% {
      randomForest(y~., data=thisTrainDt, mtry=113, do.trace=TRUE, ntree=i)
    }
    results[[thisStationName]] <- list(data=thisTrainDt, model=rf)
  }
  return(results)
  
}