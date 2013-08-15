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

registerDoMC(cores=1)

# Create the datasets
####################

## Only need to run this once!
convertRawData <- function() {
    for (f in trainFiles) {
        ncdf2Rdata(paste0(dataFolder, trainFolder, f))
    }

    for (f in testFiles) {
        ncdf2Rdata(paste0(dataFolder, testFolder, f))
    }
}



##########
# Data exploration
##########

if (FALSE) {
    # Downard shortwave radiative flux
    fpath <- paste0(dataFolder, trainFolder, trainRData[3])
    load(fpath)
    # Get a station
    stn <- stationInfo[10]
    pts <- getPoints(stn$elon, stn$nlat)
    
    # Taking a look at the ensembles -- how much do they actually differ?
    # Key order is lat, lon, hour, ens, date
    tbl <- parseDate(tbl)
    tbl <- tbl[, ens:=as.factor(ens)]
    setkey(tbl, lat, lon, hour, ens, date)
    # If you look at the ensembles, there is actually pretty substantial separation
    # on some days
    selector <- CJ(pts$lat[1],
                     pts$lon[1],
                     18,
                     unique(tbl$ens),
                     seq(as.Date('1994-01-01'), as.Date('1994-12-31'), by='1 day'))
    ggplot(
           tbl[selector],
           aes(x=date, y=value, color=ens)
           ) + geom_point() + scale_colour_discrete()
    
    # Is there some correlation between the SD of the ensembles' predictions
    # and the correlation with the actual values?
    stnTrain <- parseDate(trainData)
    stnTrain <- stnTrain[selector[, 5, with=FALSE], c('date', stn$stid), with=FALSE]
 }

#########
# Run
#########


makeSubmission <- function() {
    predDf <- list()
    i <- 1
    for (s  in stationNames) {
        cat('Prediction for station ', s, '\n', sep='')
        predDf[[i]] <- predictStation(s)
        i <- i + 1
    }

    res <- join_all(predDf, by="date")

    write.csv(res, file = "submission.csv", row.names=FALSE)
}

#########
# Some attempts at ensembling the straightforward GBM and RF models
#########

if (FALSE) {
    # Ensembling GBM and RF by simple average
    rfSub <- data.table(read.csv('submissionRF.csv'))
    gbmSub <-data.table(read.csv('submissionGBM.csv'))
    rfMelt <- melt(rfSub, id.var='date')
    rfMelt$mdl <- 'rf'
    gbmMelt <- melt(gbmSub, id.var='date')
    gbmMelt$mdl <- 'gbm'
    melted <- data.table(rbind(gbmMelt, rfMelt))
    ensSub <- melted[,list(value=mean(value)), by='date,variable']
    newSub <- dcast(ensSub, date ~ variable)
    write.csv(newSub, file="ensSubmission.csv", row.names=FALSE)
}

if (FALSE) {
    # Weighted average of the RF and GBM models
    fpath <- paste0(dataFolder, 'GBMmodels/', stationNames[1], '.Rdata')
    load(fpath)
    errors <- data.table(stn=stationNames, gbm=NA, rf=NA)

    rfErrors <- foreach(stn=stationNames, .combine='c') %dopar% {
        fpath <- paste0(dataFolder, 'RFmodels/', stn, '.Rdata')
        cat('Opening ', fpath, '\n', sep='')
        load(fpath)
        err <- tail(stnFit$mse, 1)
        return(err)
    }

    gbmErrors <- foreach(stn=stationNames, .combine='c') %dopar% {
        fpath <- paste0(dataFolder, 'GBMmodels/', stn, '.Rdata')
        cat('Opening ', fpath, '\n', sep='')
        load(fpath)
        err <- tail(stnFit$train.error, 1)
        return(err)
    }

    errors <- data.table(station=stationNames, gbm=gbmErrors, rf=rfErrors)
    wts <- errors[, list(gbm=1-(gbm/sum(gbm, rf)), rf=1-(rf/sum(gbm, rf))), by=station]
    rfSub <- data.table(read.csv('submissionRF.csv'))
    gbmSub <-data.table(read.csv('submissionGBM.csv'))

    rfWeighted <- data.table(t(apply(rfSub, 1, function(x) {x * wts[, c(1, rf)]})))
    gbmWeighted <- data.table(t(apply(gbmSub, 1, function(x) {x * wts[, c(1, gbm)]})))
    newSub <- data.table(cbind(date=rfWeighted[,date], rfWeighted[,-1,with=FALSE] + gbmWeighted[,-1,with=FALSE]))
    write.csv(newSub, file="wtensSubmission.csv", row.names=FALSE)
}

##################
# Model exploration
################

if (FALSE) {
    # Some exploration of the models
    modelFolder <- paste0(dataFolder, 'RFmodels/')

    a <- foreach(stn=stationNames) %dopar% {
        cat(stn, '\n')
        fpath <- paste0(modelFolder, stn, '.Rdata')
        load(fpath) # gives stnFit object
        imp <- importance(stnFit)
        imp <- imp[order(imp, decreasing=TRUE),,drop=FALSE]
        return(imp)
    }

    # Try to figure out which factors in which dimensions to use
    b <- unlist(lapply(a, rownames))
    dim(b) <- c(1200, 98)

    impFactors <- data.table(data.frame(b, stringsAsFactors=FALSE))
    setnames(impFactors, stationNames)

    getFactorNames <- function(x) {
        sapply(strsplit(x, '_'), function(x){return(x[1])})
    }
    # Names only, excludes the location info
    impFactorNames <- impFactors[,lapply(.SD, getFactorNames)]
    # get unique factor names in top x factors for all stations
    uniqueFactors <- unique(unlist(apply(impFactorNames[1:100,], 2, unique)))

    # Hour 12 data consistently shows up at the bottom of the list
    # up and downward SW flux at later hours of the day (18+) generally are at the top
    # If you take top 400 factors for all models, then all of the 15 factors are being used
    # Number of important factors drops significantly at top 100, down to 9

    # Which models have the most error and why?
    a <- foreach(stn=stationNames) %dopar% {
        cat(stn, '\n')
        fpath <- paste0(modelFolder, stn, '.Rdata')
        load(fpath) # gives stnFit object
        return(tail(stnFit$mse, 1))
    }

    rfErrors <- data.table(stid=stationNames, error=unlist(a), key='stid')
    errorWithInfo <- stationInfo[rfErrors]
    ggplot(errorWithInfo, aes(x=nlat, y=error)) + geom_point()
    ggplot(errorWithInfo, aes(x=elon, y=error)) + geom_point()
    # This is interesting -- error seems to be correlated with longitude
    # stations with less negative longitude have higher errors
    ggplot(errorWithInfo, aes(x=elev, y=error)) + geom_point()
    # higher elevations mean lower errors

    gbmErrors <- foreach(stn=stationNames, .combine='c') %dopar% {
        fpath <- paste0(dataFolder, 'GBMmodels/', stn, '.Rdata')
        cat('Opening ', fpath, '\n', sep='')
        load(fpath)
        err <- tail(stnFit$train.error, 1)
        return(err)
    }    
    
    gbmErrors <- data.table(stid=stationNames, error=gbmErrors, key='stid')
    errorWithInfo <- errorWithInfo[gbmErrors]
    setnames(errorWithInfo, 'error', 'rfError')
    setnames(errorWithInfo, 'error.1', 'gbmError')
    save(errorWithInfo, file='modelErrors.RData')
    
    # GBM Models
    modelFolder <- paste0(dataFolder, 'GBMmodels/')
    fpath <- paste0(modelFolder, 'ACME', '.Rdata')
    load(fpath)
    a <- foreach(stn=stationNames) %dopar% {
        cat(stn, '\n')
        fpath <- paste0(modelFolder, stn, '.Rdata')
        load(fpath) # gives stnFit object
        return(summary(stnFit))
    }

  
}
