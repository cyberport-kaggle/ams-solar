##########
# Load basic info
##########
library(ncdf4)
library(doMC)
library(foreach)
library(plyr)
library(reshape)
library(caret)
library(randomForest)
library(hash)
library(data.table)
library(gbm)
library(ggplot2)


stationInfo <- data.table(read.csv(paste(dataFolder, 'station_info.csv', sep=''), stringsAsFactors = FALSE), key='stid')
stationNames <- stationInfo$stid
knownDims <- c('lat', 'lon', 'ens', 'fhour', 'intValidTime', 'intTime', 'time')
trainData <- data.table(read.csv(paste(dataFolder, 'train.csv', sep='')), key='Date')
setnames(trainData, 'Date', 'date')
sampleSub <- read.csv(paste0(dataFolder, 'sampleSubmission.csv'))

# Names of the netCDF files that contain the data
# The following are all of the variables provided
#Variable	Description	Units
#apcp_sfc	3-Hour accumulated precipitation at the surface	kg m-2
#dlwrf_sfc	Downward long-wave radiative flux average at the surface	W m-2
#dswrf_sfc	Downward short-wave radiative flux average at the surface	W m-2
#pres_msl	Air pressure at mean sea level	Pa
#pwat_eatm	Precipitable Water over the entire depth of the atmosphere	kg m-2
#spfh_2m	Specific Humidity at 2 m above ground	kg kg-1
#tcdc_eatm	Total cloud cover over the entire depth of the atmosphere	 %
#tcolc_eatm	Total column-integrated condensate over the entire atmos.	kg m-2
#tmax_2m	 Maximum Temperature over the past 3 hours at 2 m above the ground	 K
#tmin_2m	 Mininmum Temperature over the past 3 hours at 2 m above the ground	 K
#tmp_2m	 Current temperature at 2 m above the ground	 K
#tmp_sfc	 Temperature of the surface	 K
#ulwrf_sfc	 Upward long-wave radiation at the surface	 W m-2
#ulwrf_tatm	 Upward long-wave radiation at the top of the atmosphere	 W m-2
#uswrf_sfc	 Upward short-wave radiation at the surface	 W m-2
trainFiles <- read.csv('train_filenames.txt', header=FALSE, stringsAsFactors = FALSE)[[1]]
testFiles <- read.csv('test_filenames.txt', header=FALSE, stringsAsFactors = FALSE)[[1]]
# Resaved nc files, since reading Rdata is so much faster
trainRData <- sub('\\.nc', '.Rdata', trainFiles)
testRData <- sub('\\.nc', '.Rdata', testFiles)

# In R, each of the factors should have dimensions of (lon, lat, fhour, ens, time)
# This is reversed from netCDF, and what the python code will read


# Create a hash table of full variable names to shorter names
# For readability
varNames <- c("Total_precipitation",
              "Downward_Long-Wave_Rad_Flux",
              "Downward_Short-Wave_Rad_Flux",
              "Pressure",
              "Precipitable_water",
              "Specific_humidity_height_above_ground",
              "Total_cloud_cover",
              "Total_Column-Integrated_Condensate",
              "Maximum_temperature",
              "Minimum_temperature",
              "Temperature_height_above_ground",
              "Temperature_surface",
              "Upward_Long-Wave_Rad_Flux_surface",
              "Upward_Long-Wave_Rad_Flux",
              "Upward_Short-Wave_Rad_Flux")
shortVarNames <- c('totalPrecip',
                   'dwnLWFlux',
                   'dwnSWFlux',
                   'pressure',
                   'precipWater',
                   'specHumid',
                   'cloud',
                   'totalCondens',
                   'maxTemp',
                   'minTemp',
                   'temp2m',
                   'tempSurf',
                   'upLWFluxSurf',
                   'upLWFlux',
                   'upSWFlux')
names(varNames) <- shortVarNames
longNames <- hash(varNames)
shortNames <- invert(longNames)

# Since we are not using NC files, we need a map of filename to variable name
names(shortVarNames) <- trainRData
trainFileToShortName <- hash(shortVarNames)


if (FALSE) {
    # This works, but the size of a single variable's data in database is well over 1gig
    databaseName <- 'database.db'
    trainTableName <- 'solarTrain'
    testTableName <- 'solarTest'
    if (!exists('db')) {
        db <- dbDriver('SQLite')
    }
    if (!exists('dbc')) {
        dbc <- dbConnect(db, databaseName)
    }
    createDbTables <- function() {
        # Initialize the database tables
        # only need to run once.  Will wipe out all other data
        # Schema:
        # Single really tall table with columns:
        #  - variable
        #  - lat
        #  - lon
        #  - hour
        #  - ens
        #  - date
        #  - value
        for (tableName in c(trainTableName, testTableName)) {
            if (dbExistsTable(dbc, tableName)) {
                conf <- readline(paste0('Table ', tableName, ' already exists.  Drop?[y/n]'))
                if (conf == 'y') {
                    dbRemoveTable(dbc, tableName)
                } else {
                    cat('Did not create table', tableName, '.\n')
                    next
                }
            }

            emptydf <- data.frame(variable='foo',
                                  lat=1,
                                  lon=1,
                                  hour=1,
                                  ens=1,
                                  date=1,
                                  value=1)[-1,]
            sql <- dbBuildTableDefinition(dbc, tableName, value=emptydf, row.names=FALSE)
            dbGetQuery(dbc, sql)
            cat('Table', tableName, 'created.\n')
        }
    }
    loadDatabase <- function(fname, train=TRUE) {
        # Processes NCDF files and sticks them into an SQLite database
        # if train is true, then it uses the training database, else it puts it into the test database
        cat("Writing file ", fname, "to database.\n", sep="")
        nc <- nc_open(fname)
        dims <- getDimensions(nc)
        varName <- getVarName(nc)
        shortVarName <- shortNames[[varName]]
        values <- ncvar_get(nc, varName)
        dimnames(values) <- list(lon=dims$lon, lat=dims$lat, hour=dims$fhour, ens=dims$ens, date=dims$intTime)
        tbl <- data.table(melt(values))
        rm(values)
        gc()
        tbl$variable <- shortVarName
        # Need to chunk the writing to the table
        # Total 40million rows...
        if (train) {
            toTable <- trainTableName
        } else {
            toTable <- testTableName
        }
        chunksize <- 5000
        len <- nrow(tbl)
        i <- 1
        j <- chunksize
        while (j <= len) {
            cat('Writing to row ', j, '\n', sep='')
            dbWriteTable(dbc, toTable, tbl[i:j], append=TRUE, row.names=FALSE)
            i <- i+chunksize
            if (j == len) {
                j <- j + 1
            } else if((j + chunksize) > len) {
                j <- len
            } else {
                j <- j+chunksize
            }
        }
    }
}


getDimensions <- function(nc) {
    # Returns a list of the dimensions and their values from the netcdf file.
    # This is so that we don't have to continue to refer to dimensions just by their index
    # And for paranoid checking that the dimensions are in the right order in each file

    res <- list()
    res$lon <- ncvar_get(nc, 'lon')
    res$lat <- ncvar_get(nc, 'lat')
    res$fhour <- ncvar_get(nc, 'fhour')
    res$ens <- ncvar_get(nc, 'ens') # ensemble numbers; 0-10
    res$intTime <- ncvar_get(nc, 'intTime')
    return(res)
}

# Get main dimensions
nc <- nc_open(paste0(dataFolder, trainFolder, trainFiles[1]))
dataDims <- getDimensions(nc)
nc_close(nc)

getVarName <- function(nc) {
    # Variable names are different in each package, so we need to read it out.
    # Since we know all of the other variable names, we just eliminate them and the only
    # one left is the variable we want
    varNames <- c(names(nc$dim), names(nc$var))
    res <- varNames[!varNames %in% knownDims]
    if (length(res) == 1) {
        return(res)
    } else {
        warning(paste("Could not find the variable name in ", nc$filename, ".\nRemaining names are ", paste(res, collapse=" "), sep=""))
        return(NULL)
    }
}

getPoints <- function(lon, lat, dims=dataDims, n=1) {
    # Given a latitude and longitude, find the n nearest points in the GEFS grid 
    # in each direction
    # Parameters:
    # @lat: The latitude of the Mesonet location
    # @lon: the longitude of the Mesonet location
    # @n: number of points in each direction to find.  So n=1 should return four points
    # @dims: the list returned by getDimensions for the file
    # 
    # Returns the points as a start index in the netcdf file and a count
    lats <- dims$lat
    lons <- dims$lon
    # Since GEFS is spaced in 1 degree intervals, multiply by n to get the max distance in degrees from the lat and lon
    maxDist <- 1 * n
    # station values are negative values, so we add 360 to get it in the same scale as the GEFS coordinates
    if (lon < 0) {
        lon <- lon + 360
    }
    latIdx <- which(abs(lats - lat) <= maxDist)
    lonIdx <- which(abs(lons - lon) <= maxDist)
    if (length(latIdx) != (n * 2) && length(lonIdx) != (n * 2)) {
        warning("Something is wrong, didn't find the right number of coordinates")
        cat(paste("Lat: ", lat, " Lon: ", lon, sep=""))
        cat(paste("in Lats: ", lats, " Lons: ", lons, sep=""))
        return(NULL)
    } else {
        cat("Using GEFS points at lats: ", paste(lats[latIdx], collapse=" "), " and lons: ", paste(lons[lonIdx], collapse=" "), '\n', sep="")
        return(list(latStart=latIdx[1], latCnt=n*2, lonStart=lonIdx[1], lonCnt=n*2))
    }
}

ncdf2Rdata <- function(fname) {
    # Converts an NC file to an rdata file
    # Stores file as a single tall data table
    cat("Converting file ", fname, "\n", sep="")
    nc <- nc_open(fname)
    dims <- getDimensions(nc)
    varName <- getVarName(nc)
    shortVarName <- shortNames[[varName]]

    values <- ncvar_get(nc, varName)
    dimnames(values) <- list(lon=dims$lon, lat=dims$lat, hour=dims$fhour, ens=dims$ens, date=dims$intTime)
    tbl <- data.table(melt(values))
    # data.table(table(values, responseName='value')) also works, but not sure which is faster
    rm(values)
    gc()
    setkey(tbl, lat, lon, hour, ens, date)
    
    #assign(shortVarName, values)
    newName <- sub('\\.nc', '.Rdata', fname)
    save(tbl, file=newName)
    # manually removing the file
    rm(tbl)
    gc()
    nc_close(nc)
}

parseDate <- function(tbl, dates=NULL) {
  # Takes any data table who has a 'date' column and parses the date to the right type (instead of int)
  if (!('date' %in% names(tbl))) {
    stop('Date column not present')
  }
  oldKeys <- key(tbl)
  if (oldKeys != c('date')) {
    setkey(tbl, 'date')
  }
  if (is.null(dates)) {
    # Can be passed in, since dates are the same for all of the data files, so we don't have to recalculate
    # this every time, which actually takes a long time, since you have to coerce from integer to string to IDate
    # Gives a data table of unique dates
    dateCheck <- substr(tbl$date[1:100], 9, 10)
    if (all(dateCheck == "00")) {
      dates <- tbl[,as.IDate(as.character(date/100), format="%Y%m%d"), by='date']
    } else {
      dates <- tbl[,as.IDate(as.character(date), format="%Y%m%d"), by='date']
    }
  }
  # join it to data table
  tbl <- dates[tbl]
  tbl[, date:= NULL]
  tbl[, date:= V1]
  tbl[, V1:= NULL]
  setkeyv(tbl, oldKeys)
  return(tbl)
}