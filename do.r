library(ncdf4)
library(doMC)
library(foreach)
library(plyr)
library(reshape)
library(caret)
library(randomForest)
library(hash)

registerDoMC(cores=detectCores())

##########
# Config
##########


dataFolder <- 'data/'
testFolder <- 'test/' # Path to test files.  Set to '' if in root dataFolder
trainFolder <- 'train/' # Path to train files
trainOutFolder <- 'cleaned/' # output of reshaped training data
testOutFolder <- 'cleanedTest/' # output of reshaped test data
modelsFolder <- 'models/'

source('1_load.r')
source('2_func.r')

#########
# Run
#########


# Create the datasets
####################

## Only need to run this once!
buildDfs(train=TRUE)
buildDfs(train=FALSE)

# Use this for Caret
for (s in stationNames) {
    stationFit(s)
}

foreach(s=stationNames) %dopar% {
    stationFit(s)
    return(NULL)
}

#buildDfs('cleanedTest/', train=FALSE)


predDf <- list()
i <- 1
for (s  in stationNames) {
    cat('Prediction for station ', s, '\n', sep='')
    predDf[[i]] <- predictStation(s)
    i <- i + 1
}
res <- join_all(predDf, by="date")
write.csv(res, file = "submission.csv", row.names=FALSE)
