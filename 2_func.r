# Creates the data.tables for each Mesonet station. 
stationData <- function(factors, train){
  if(train) path = "data/train"
  else path = "data/test"
  res = list()
  
  #   Loop through factors
  for(f in factors){
    data = get(load(Sys.glob(file.path(path, paste0(f, '*', ".Rdata",collapse ='')))))

    #   For each factor, loop through stations
    for(s in 1:nrow(stationInfo)){
      locs = getPoints(stationInfo[s, elon], stationInfo[s, nlat], n = 2)
      return(locs)
    }
  }
}  




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

