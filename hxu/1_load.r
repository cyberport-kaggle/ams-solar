
##########
# Config
##########

library(plyr)
library(reshape)

dataFolder <- '../data/'
testFolder <- 'test/' # Path to test files.  Set to '' if in root dataFolder
trainFolder <- 'train/' # Path to train files

##########
# Load basic info
##########

stationInfo <- read.csv(paste(dataFolder, 'station_info.csv', sep=''))
stationNames <- stationInfo$stid
knownDims <- c('lat', 'lon', 'ens', 'fhour', 'intValidTime', 'intTime', 'time')
trainData <- read.csv(paste(dataFolder, 'train.csv', sep=''))

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
trainFiles <- read.csv('train_filenames.txt', header=FALSE)[[1]]
testFiles <- read.csv('test_filenames.txt', header=FALSE)[[1]]

# In R, each of the factors should have dimensions of (lon, lat, fhour, ens, time)
# This is reversed from netCDF, and what the python code will read


# Create a hash table of full variable names to shorter names
# For readability
library(hash)
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

getPoints <- function(lat, lon, dims, n=1) {
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
        cat(paste("Using GEFS points at lats: ", paste(lats[latIdx], collapse=" "), " and lons: ", paste(lons[lonIdx], collapse=" "), sep=""))
        return(list(latStart=latIdx[1], latCnt=n*2, lonStart=lonIdx[1], lonCnt=n*2))
    }
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
    shortVarName <- paste(shortNames[[varName]], '_', latIdx, lonIdx, fhourIdx, ensIdx, sep="")
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

combineHours <- function(nc, lonIdx, latIdx, ensIdx, method="mean") {
    # Combines the readings across each forecast hour in a given day using a particular method
    if (!(method %in% c('mean', 'sum'))) {
        stop("Method is not valid")
    }
    varName <- getVarName(nc)
    # 0 indicates that it's combined across all hours
    shortVarName <- paste(shortNames[[varName]], '_', latIdx, lonIdx, '0', ensIdx, sep="")
    # should be seq(12, 24, 3), but we're paranoid
    hours <- getDimensions(nc)$fhour
    # Build a list with all of the data at each hour
    tmp <- list()
    for (i in 1:length(dims$fhour)) {
        tmp[[i]] <- getVarData(nc, lonIdx, latIdx, i, ensIdx)
    }
    dframe <- join_all(tmp, by='date')
    dates <- dframe$date
    if (method == 'mean') {
        values <- apply(dframe[,-1], 1, mean)
    } else if(method == 'sum') {
        values <- apply(dframe[,-1], 1, sum)
    }
    
    res <- data.frame(date=dates)
    res[shortVarName] <- values
    return(res)
}


##########
# Load data
##########

# Libraries ncdf4 and RNetCDF have basically the same functionality
library(ncdf4)
nc <- nc_open(paste(dataFolder, trainFolder, trainFiles[1], sep=''))

# When extracting data from the ncdf file, it will be necessary to utilize the offsets
# for each dimension, otherwise it'll read the whole dataset out
# but this will require some manual determination of which values we want
