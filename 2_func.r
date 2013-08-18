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

selectiveRF <- function(stations=stationNames,
                        outputFolder='selectiveRF/',
                        preprocess=TRUE,
                        train=TRUE,
                        predict=TRUE) {
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
  #
  # Additional Factors and cleaning:
  # - Added a "Days from summer solstice" factor.  0 indicates that the day is the summer solsitce (June 20)
  # and increases as you get farther from the summer solstice.  It maxes out at 365/2 = 183 days
  # - Bad training data where the value on day t for a station is equal to the value on day t-1 is removed
  
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
    # Annoyingly the date is formatted differently in the test data
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
  

  if (preprocess) {
    # files we're interested in:
    if (train) {
      fluxFiles <- trainRData[c(3, 13, 14, 15)]
      otherFiles <- trainRData[c(7, 12)]
    } else {
      fluxFiles <- testRData[c(3, 13, 14, 15)]
      otherFiles <- testRData[c(7, 12)]
    }
    cleanDts <- list()
    results <- list()
    # If a subset of stations are passed in, we have to get the indexes of the stations
    stationIndices <- which(stationNames %in% stations)
    
    # Load a flux file
    for (thisFile in fluxFiles) {
      cat('Loading data for factor', thisFile, '\n')
      if (train) {
        load(file.path(dataFolder, trainFolder, thisFile))
        thisVarName <- trainFileToShortName[[thisFile]]
      } else {
        load(file.path(dataFolder, testFolder, thisFile))
        thisVarName <- testFileToShortName[[thisFile]]
      }
      
      # Process all stations
      stnDts <- foreach(i=stationIndices) %dopar% {
        thisStn <- stationInfo[i]
        return(processSingleStation(thisStn))
      }
      # Have to combine out here, since can't do it in the parallelized code
      for (i in 1:length(stnDts)) {
        thisStn <- stationInfo[stationIndices[i]]$stid
        cleanDts <- mergeStationData(thisStn, stnDts[[i]], cleanDts)
      }
      rm(tbl)
      gc()
    }
    
    # non-flux files
    for (thisFile in otherFiles) {
      cat('Loading data for factor', thisFile, '\n')
      if (train) {
        load(file.path(dataFolder, trainFolder, thisFile))
        thisVarName <- trainFileToShortName[[thisFile]]
      } else {
        load(file.path(dataFolder, testFolder, thisFile))
        thisVarName <- testFileToShortName[[thisFile]]
      }
      
      # Process all stations
      stnDts <- foreach(i=stationIndices) %dopar% {
        thisStn <- stationInfo[i]
        return(processSingleStation(thisStn, combine=FALSE))
      }
      # Have to combine out here, since can't do it in the parallelized code
      for (i in 1:length(stnDts)) {
        thisStn <- stationInfo[stationIndices[i]]$stid
        cleanDts <- mergeStationData(thisStn, stnDts[[i]], cleanDts)
      }
      rm(tbl)
      gc()
    }
    # At this point, cleanDts is a list with length == number of stations
    if (train) {
      subfolder <- 'inputData/'
    } else {
      subfolder <- 'inputTestData/'
    }
    # Save each station individually
    dir.create(file.path(dataFolder, outputFolder))
    dir.create(file.path(dataFolder, outputFolder, subfolder))
    
    for (i in 1:length(cleanDts)) {
      thisStationName <- names(cleanDts)[i]
      cat('Saving training dataset for', thisStationName, '\n')
      thisDt <- cleanDts[[i]]
      save(thisDt, file=paste0(dataFolder, outputFolder, subfolder, thisStationName, '.RData'))
    }
    
    rm(cleanDts)
    gc()
  }
  
  if (train) {
    # Train each station
    cat('Training models\n')
    foreach(i=1:length(stations)) %dopar% {
      sink('log.txt', append=TRUE)
      thisStationName <- stations[i]
      cat('Training model for station', thisStationName, 'at', format(Sys.time(), "%a %b %d %X %Y"), '\n')
      
      load(paste0(dataFolder, outputFolder, 'inputData/', thisStationName, '.RData'))
      thisTrainDt <- merge(trainData[,c('date', thisStationName), with=FALSE], thisDt, by='date')
      setnames(thisTrainDt, thisStationName, 'y')
      
      # Some post-processing and factor generation
      # Remove bad training data
      badData <- which(diff(thisTrainDt$y, 1) == 0) + 1
      thisTrainDt <- thisTrainDt[!badData]
      cat('Dropped', length(badData), 'datapoints due to repeated values\n')
      # Add in days from summer solstice
      dayOfYear <- abs(yday(as.IDate(as.character(thisTrainDt$date), format="%Y%m%d")) - yday(as.IDate('2013-06-20')))
      daysFromSolstice <- 183-abs(183-dayOfYear) 
      thisTrainDt[, daysFromSolstice := daysFromSolstice]
      # The presence of the date column doens't seem to hurt the RF, but should probably
      # drop it, since it's not really part of the model
      thisTrainDt[, date:=NULL]
      #     fitCtrl <- trainControl(method = "cv",
      #                             number = 5)
      #     stnFit <- train(y ~ .,
#                           data=thisTrainDt,
      #                     method="rf",
      #                     trControl=fitCtrl,
      #                     verbose=TRUE,
      #                     ntree=500)
      #     # Finds mtry=113 is optimal
      # Performance really stops improving at around 400 trees
      # This is for checking the MSE of the models
      rf <- randomForest(y~., data=thisTrainDt, mtry=113, do.trace=TRUE, ntree=400)
      trainModel <- list(data=thisTrainDt, model=rf)
      save(trainModel, file=paste0(dataFolder, outputFolder, thisStationName, '.RData'))
      # Note that parallelized RF doesn't preserve the MSE
      #     rf <- foreach(i=rep(100, 4), .combine=combine, .multicombine=TRUE) %dopar% {
      #       randomForest(y~., data=thisTrainDt, mtry=113, do.trace=TRUE, ntree=i)
      #     }
      sink()
      return(NULL)
    }
  }
  
  if (predict) {
    cat('Predicting test data\n')
    # For each station, load in its model and its training data
    # You must have run the preprocess step while train=FALSE
    preds <- foreach(i=1:length(stations), .combine=data.table, .multicombine=TRUE) %dopar% {
      sink('log.txt', append=TRUE)
      thisStationName <- stations[i]
      cat('Predicting for station', thisStationName, '\n')
      # gives trainModel list where $data is the training data and $model is the model
      load(paste0(dataFolder, outputFolder, thisStationName, '.RData'))
      # gives thisDt object that is a data.table of the test data
      load(paste0(dataFolder, outputFolder, 'inputTestData/', thisStationName, '.RData'))
      
      # Add in days from summer solstice
      dayOfYear <- abs(yday(as.IDate(as.character(thisDt$date), format="%Y%m%d")) - yday(as.IDate('2013-06-20')))
      daysFromSolstice <- 183-abs(183-dayOfYear) 
      thisDt[, daysFromSolstice := daysFromSolstice]
      # drop date column
      thisDt[, date:=NULL]
      
      res <- predict(trainModel$model, thisDt)
      sink()
      return(res)
    }
    setnames(preds, names(preds), stations)
    preds <- data.table(Date=sampleSub$Date, preds)
    cat('Saving submission file\n')
    write.csv(preds, file=paste0(dataFolder, outputFolder, 'submission.csv'), row.names=FALSE)
  }
}