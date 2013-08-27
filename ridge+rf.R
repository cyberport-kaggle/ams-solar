# Ridge regression for baseline 
# Get train and test data in format for ridge regression
# Process the train and test data by scaling and removing ZV values
# Fit the coefficients fro the ridge regresssion model
# Make predictions both from the train and test set. yridgetrain and yridgetest

library(caret)
library(data.table)
library(reshape2)
library(foreach)
library(doMC)
library(profr)
library(randomForest)
registerDoMC(cores = 4)

factors = c('dswrf_sfc','dlwrf_sfc','uswrf_sfc','ulwrf_sfc','ulwrf_tatm','pwat_eatm','tcdc_eatm','apcp_sfc',
            'pres_msl','spfh_2m','tcolc_eatm','tmax_2m','tmin_2m','tmp_2m','tmp_sfc')

ridgetrain = get(load("ridge/trainData.RData"))
ridgetest = get(load("ridge/testData.RData"))
ytrain = data.table(read.csv(file="data/train.csv"), key="Date")
stations = names(ytrain)[-1]

# Remove zero variance columns
zv = nearZeroVar(ridgetrain, saveMetrics=TRUE)$zeroVar
ridgetrain = ridgetrain[, which(!zv), with = FALSE]

# Scale and center results 
preProcValues <- preProcess(ridgetrain[, setdiff(colnames(ridgetrain), "date"), with = FALSE])
ridgetrain <- data.table(cbind(ridgetrain[, date], 
                        predict(preProcValues, ridgetrain[, setdiff(colnames(ridgetrain), "date"), with = FALSE])))
setnames(ridgetrain, colnames(ridgetrain)[1], "date")
setkey(ridgetrain, date)

# Convert response and predictors into matrices 
Xtrain = cbind(1, (as.matrix(ridgetrain[, setdiff(colnames(ridgetrain), c("date", "y")), with = FALSE])))
Ytrain = as.matrix(ytrain[,-1, with = FALSE])

# Convert into matrix format
Xtest = as.matrix(predict(preProcValues, ridgetest[, setdiff(colnames(ridgetest), "date"), with = FALSE]))
Xtest = cbind(1, Xtest)

# Calculate coefficients of ridge regression 
lambda = 128
beta = solve(t(Xtrain) %*% Xtrain + lambda * diag(ncol(Xtrain))) %*% t(Xtrain) %*% Ytrain

# Do ridge regression prediction 
yridgetrain = Xtrain %*% beta
yridgetest = Xtest %*% beta

# Get into right format for submission
yridgetrain = data.table(cbind(ridgetrain[, "date", with = FALSE], yridgetrain))
yridgetest = data.table(cbind(ridgetest[, "date", with = FALSE], yridgetest))

# RF for the residuals
# Calculate residuals for RF fitting residuals.train = ytrain - yridgetrain
# For each station 
#   Get data training data for random forest residual model (from data/cleaned/)
#   Fit a model for each Mesonet station to residuals.train (the residuals)
# Get the test data for the random forest residual model (from data/cleanedTest)
# Make a prediction of the residauls for each Mesonet station from test data
# and put into a data.table residuals.test
# The output is now the predicted resituals.  Add the ridge regression baseline to 
# form complete prediction final.ans = residuals.test + yridgetest

residuals.train = ytrain - yridgetrain

residualRFFit <- function(stn){
  trainingPath=paste0('data/cleaned/')
  cat('Fitting model for', stn, '\n', sep=' ')
  f <- paste0(trainingPath, stn, '.RData')
  values <- get(load(f))
  values <- values[, -1, with = FALSE]
  
  y <- residuals.train[,stn, with = FALSE]
  trainDf <- data.frame(y=y, values)
  colnames(trainDf)[1] = "y"
  
  fitCtrl <- trainControl(method = "cv",
                          number = 5)
  trainGrid = expand.grid(.mtry = seq(2,102, 20))
  # stnFit <- train(y ~ .,
  #                 data=trainDf,
  #                 method="rf",
  #                 trControl=fitCtrl,
  #                 tuneGrid = trainGrid,
  #                 verbose=TRUE,
  #                 do.trace=TRUE,
  #                 ntree=250)
  
  stnFit <- randomForest(y~., data=trainDf, do.trace=TRUE, mtry=44, ntree=250)
  
  save(stnFit, file=paste0('./data/RFresidualsmodels/', stn, '.Rdata'))
}

foreach(stn = stations) %dopar%{
# for(stn in stations[2:2]){
  residualRFFit(stn)
  return(stn)
}

pred <- list()
for(stn in stations){
  cat("Predicting for", stn, "\n")
  load(paste0("./data/RFresidualsmodels/", stn, ".Rdata"))
  load(paste0("./data/cleanedTest/", stn, ".RData"))
  pred[[stn]] <- predict(stnFit, tbl)
}

dates <- tbl[,date]
residuals.test <- do.call(cbind, pred)

# final.ans <- residuals.test + yridgetest[, -1, with = FALSE]

yridgetest <- read.csv("ridge/submission.csv")
final.ans <- yridgetest[,-1] + residuals.test

final.ans <- cbind(dates, final.ans)
colnames(final.ans)[1] <- "Date"
write.csv(final.ans, "rf+ridge submissions.csv", row.names=FALSE)

# mdl <- get(load("./data/RFresidualsmodels/ADAX.Rdata"))
# train <- get(load("./data/cleanedTest/ADAX.RData"))
# test <- get(load("./data/cleanedTest/ADAX.RData"))
# predict(mdl, test)