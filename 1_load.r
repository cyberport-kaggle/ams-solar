##########
# Load basic info
##########

stationInfo <- read.csv(paste(dataFolder, 'station_info.csv', sep=''), stringsAsFactors = FALSE)
stationNames <- stationInfo$stid
knownDims <- c('lat', 'lon', 'ens', 'fhour', 'intValidTime', 'intTime', 'time')
trainData <- read.csv(paste(dataFolder, 'train.csv', sep=''))
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
    cat("Converting file ", fname, "\n", sep="")
    nc <- nc_open(fname)
    dims <- getDimensions(nc)
    varName <- getVarName(nc)
    shortVarName <- shortNames[[varName]]

    values <- ncvar_get(nc, varName)
    dimnames(values) <- list(lon=dims$lon, lat=dims$lat, hour=dims$fhour, ens=dims$ens, date=dims$intTime)
    ncRData <- list(name=shortVarName, data=values)
    #assign(shortVarName, values)
    newName <- sub('\\.nc', '.Rdata', fname)
    save(ncRData, file=newName)
}

getVarData <- function(nc, lonIdx, latIdx, fhourIdx, ensIdx) {
    # Returns a data frame with one column.  Rows are dates, and the single column
    # is the variable data contained in the NetCDF file provided.  All other dimensions
    # must be specified
    # Parameters:
    # @nc: The open netCDF file
    # @latIdx: Index of the latitude dimension
    # @lonIdx: Index of the longitude dimension
    # @fhourIdx: Index of the forecast hour dimension
    # @ensIdx: Index of the ensemble dimension
    varName <- getVarName(nc)
    # append indexes to the name
    shortVarName <- paste(shortNames[[varName]], '_', paste(latIdx, lonIdx, fhourIdx, ensIdx, sep="."), sep="")
    # build inputs to ncvar_get
    startIdx <- c(lonIdx, latIdx, fhourIdx, ensIdx, 1)
    cnt <- c(1, 1, 1, 1, -1) # -1 on intTime to get all values
    values <- ncvar_get(nc, varName, start=startIdx, count=cnt)
    dates <- ncvar_get(nc, 'intTime')
    # build result data frame
    res <- data.frame(date=dates)
    res[shortVarName] <- values
    return(res)
}

getVarRData <- function(fpath, lonIdx, latIdx, fhourIdx, ensIdx) {
    # RData version of the above
    load(fpath) # creates an ncRData object
    shortVarName <- paste(ncRData[['name']], '_', paste(latIdx, lonIdx, fhourIdx, ensIdx, sep="."), sep="")
    values <- ncRData[['data']][lonIdx, latIdx, fhourIdx, ensIdx,]
    # TODO: take indexes of lat and lon, then use cast to reshape it instead of manually combining it
    res <- data.frame(date=names(values))
    res[shortVarName] <- values
}

getAllVarData <- function(nc, lonIdx, latIdx) {
    # For a given GEFS location and variable file, get and collapse all data
    # Returns a dataframe where columns are the variables at different combinations of
    # ensemble and forecast hour
    # so should have a total of 5 * 11 = 55 columns
    dims <- getDimensions(nc)
    tmp <- list()
    k <- 1
    for (i in 1:length(dims$fhour)) {
        for (j in 1:length(dims$ens)) {
            tmp[[k]] <- getVarData(nc, lonIdx, latIdx, i, j)
            k <- k + 1
        }
    }
    # Profiling shows that the join all is really slow, so instead we'll just copy the columns, since
    # in theory everything should be the same size.
    #dframe <- join_all(tmp, by='date')
    # Not using join_all reduces the processing time by 1/3
    dates <- data.frame(date=tmp[[1]][,1])
    colNames <- unlist(lapply(tmp, function(x) {return(colnames(x)[2])}))
    values <- data.frame(do.call(cbind, lapply(tmp, function(x) {return(x[,2])})))
    colnames(values) <- colNames
    dframe <- cbind(dates, values)
    if (sum(colnames(dframe) != 'date') != 55) {
        warning(paste("Flattened data does not contain 55 columns in file ", nc$filename, sep=""))
    }
    return(dframe)
}

getAllVarRData <- function(fpath, lonIdx, latIdx) {
    # rdata version of the above
    # uses global dataDims
    k <- 1
    tmp <- list()
    for (i in 1:length(dataDims$fhour)) {
        for (j in 1:length(dataDims$ens)) {
            tmp[[k]] <- getVarRData(fpath, lonIdx, latIdx, i, j)
            k <- k + 1
        }
    }

    dates <- data.frame(date=tmp[[1]][,1])
    colNames <- unlist(lapply(tmp, function(x) {return(colnames(x)[2])}))
    values <- data.frame(do.call(cbind, lapply(tmp, function(x) {return(x[,2])})))
    colnames(values) <- colNames
    dframe <- cbind(dates, values)
    if (sum(colnames(dframe) != 'date') != 55) {
        warning(paste("Flattened data does not contain 55 columns in file ", nc$filename, sep=""))
    }
    return(dframe)

}

getFullVarData <- function(nc) {
    # does not work, data is too large
    # Get *ALL* data out of the netcdf function.  getAllVarData only gets it for one latitude or longitude.
    dims <- getDimensions(nc)
    varName <- getVarName(nc)
    tmp <- list()
    k <- 1
    shortVarName <- paste(shortNames[[varName]], '_', paste(latIdx, lonIdx, fhourIdx, ensIdx, sep="."), sep="")
    values <- ncvar_get(nc, varName)
    dates <- dims$intTime
    res <- matrix(data=NA, nrow=(5113*11*5*9*16), ncol=6)
    colnames(res) <- c('lon', 'lat', 'hour', 'ens', 'date')
    i <- 1
    for (lo in dims$lon) {
        for (la in dims$la) {
            for (h in dims$h) {
                for (e in dims$ens) {
                    for (d in dims$intTime) {
                        rw <- c(lo, la, h, e, d, values[lo, la, l, h, e, d])
                        res[i,] <- rw
                        i <- i + 1
                        if (i %% 500) {
                            h(res)
                        }
                    }
                }
            }
        }
    }

    dimnames(values) <- list(lon=dims$lon, lat=dims$lat, hour=dims$fhour, ens=dims$ens, date=dims$intTime)

    # trying to melt/recast/adply, but run into memory barriers

}

combineHours <- function(df, method="mean") {
    # given a df that is a result of getAllVarData, reduces dimension of hours to 1
    #
    # For example, df will have column names with suffixes 1111, 1112, 1113 ... 1121, 1122... etc.
    # Forecast hour is the third number, so this function will combine columns:
    # 1111, 1121, 1131, 1141, 1151 to be 1101, and then the same with 1112, 1122, etc.
    #
    # Since getAllVarData has fixed lat and lon, we actually only need to flatten for the ensemble dimension

    if (!(method %in% c('mean', 'sum', 'max', 'min'))) {
        stop("Method is not valid")
    }

    # Get the suffixes of each column name
    # Split first by _, then split by .
    splitNames <- strsplit(colnames(df), '_')
    suffixes <- unlist(lapply(splitNames, function(x) {return(x[2])}))
    # This gives a list of numeric vectors, which indicate the dimensions of lat, lon, hour, and ensemble for each column of the df
    # Also element 1 of the list is NA since that was the date column
    idx <- strsplit(suffixes, '\\.')

    # Additionally, infer a few things about this data.frame from the names
    varName <- splitNames[[2]][1]
    lonIdx <- substr(splitNames[[2]][2], 1, 1)
    latIdx <- substr(splitNames[[2]][2], 3, 3)

    ens <- as.numeric(unlist(lapply(idx, function(x) {x[4]})))
    uniqueEns <- unique(ens[!is.na(ens)])
    res <- data.frame(date=df$date)
    for (e in uniqueEns) {
        # The columns to merge
        colIdx <- which(ens == e)
        newCol <- apply(df[,colIdx], 1, get(method))
        # add the new column to the res data frame
        newName <- paste(varName, '_', paste(lonIdx, latIdx, 0, e, sep="."), sep="")
        res[,newName] <- newCol
    }
    return(res)
}

combineEns <- function(df, method='mean') {
    # Basically exactly the same as combineHours, but for combining across ensembles
    if (!(method %in% c('mean', 'sum', 'max', 'min'))) {
        stop("Method is not valid")
    }

    # Get the suffixes of each column name
    # Split first by _, then split by .
    splitNames <- strsplit(colnames(df), '_')
    suffixes <- unlist(lapply(splitNames, function(x) {return(x[2])}))
    # This gives a list of numeric vectors, which indicate the dimensions of lat, lon, hour, and ensemble for each column of the df
    # Also element 1 of the list is NA since that was the date column
    idx <- strsplit(suffixes, '\\.')

    # Additionally, infer a few things about this data.frame from the names
    varName <- splitNames[[2]][1]
    lonIdx <- substr(splitNames[[2]][2], 1, 1)
    latIdx <- substr(splitNames[[2]][2], 3, 3)

    hrs <- as.numeric(unlist(lapply(idx, function(x) {x[3]})))
    uniqueHrs <- unique(hrs[!is.na(hrs)])
    res <- data.frame(date=df$date)
    for (hr in uniqueHrs) {
        # The columns to merge
        colIdx <- which(hrs == hr)
        newCol <- apply(df[,colIdx], 1, get(method))
        # add the new column to the res data frame
        newName <- paste(varName, '_', paste(lonIdx, latIdx, hr, 0, sep="."), sep="")
        res[,newName] <- newCol
    }
    return(res)

}


unlistData <- function(allData) {
    # Takes a list of data frames and combines them.  Expects all data frames to have a date column as the first column.
    # also all data frames should be the same number of rows.
    dates <- data.frame(date=allData[[1]][,1])
    colNames <- unlist(lapply(allData, function(x) {return(colnames(x)[-1])}))
    values <- data.frame(do.call(cbind, lapply(allData, function(x) {return(x[,-1])})))
    colnames(values) <- colNames
    allData <- cbind(dates, values)
    return(allData)
}

parSingleStationData <- function(stn, dims, subFolder=trainFolder, fileNames=trainFiles) {
    # Parallelized version of a single station data
    cat('Cleaning data for station', stn, '\n')

    stnInfo <- stationInfo[stationInfo$stid == stn,]
    # get four closest GEFS locations
    gefsLocs <- getPoints(stnInfo$elon, stnInfo$nlat, dims)
    latIdx <- seq(gefsLocs$latStart, gefsLocs$latStart + gefsLocs$latCnt - 1)
    lonIdx <- seq(gefsLocs$lonStart, gefsLocs$lonStart + gefsLocs$lonCnt - 1)
    allData <- foreach(f=fileNames) %dopar% {
        fileData <- list()
        i <- 1
        cat('Opening ', f, '\n', sep='')
        nc <- nc_open(paste0(dataFolder, subFolder, f))
        for (lat in latIdx) {
            for (lon in lonIdx) {
                dtmp <- getAllVarData(nc, lon, lat)
                # flatten ensembles
                fileData[[i]] <- combineEns(dtmp)
                i <- i+1
            }
        }
        nc_close(nc)
        return(unlistData(fileData))
    }
    return(unlistData(allData))
}

buildDfs <- function(train=TRUE) {
# For each Mesonet location:
## Find the four nearest GEFS Locations
## For each variable we want to include in the model
### For each GEFS location
#### Average over ensembles
#### average over prediction hours
#### Results in one column, add to dataset
## write df to folder
    if (train) {
        subFolder <- trainFolder
        fileNames <- trainFiles
        dfPath <- trainOutFolder
    } else {
        subFolder <- testFolder
        fileNames <- testFiles
        dfPath <- testOutFolder
    }
    
    # Read in base dimension information
    nc <- nc_open(paste0(dataFolder, subFolder, fileNames[1]))
    dims <- getDimensions(nc)
    nc_close(nc)

    for (stn in stationNames) {
        if (train) {
            allData <- parSingleStationData(stn, dims, subFolder=subFolder, fileNames=fileNames)
        } else {
            allData <- parSingleStationData(stn, dims, subFolder=subFolder, fileNames=fileNames)
        }
        #allData <- join_all(allData, by="date")
        cat('Writing data frame of ', ncol(allData) - 1, ' columns\n', sep="")
        fn <- paste0(dataFolder, dfPath, stn, '.RData')
        save(allData, file = fn)
    }
}


