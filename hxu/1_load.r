
##########
# Config
##########

dataFolder <- '../data/'
testFolder <- 'test/' # Path to test files.  Set to '' if in root dataFolder
trainFolder <- 'train/' # Path to train files

##########
# Load basic info
##########

stationInfo <- read.csv(paste(dataFolder, 'station_info.csv', sep=''))
stationNames <- stationInfo$stid

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

##########
# Load data
##########

# Looks like these two libraries are basically the same.  I think ncdf4 might be a bit more user friendly?
# But neither of them really gives you the ability to introspect into the ncdf file, you have to just know it
# from the descriptive information that is printed out.
library(ncdf4)
nc <- nc_open(paste(dataFolder, trainFolder, trainFiles[1], sep=''))

library(RNetCDF)
nc <- open.nc(paste(dataFolder, trainFolder, trainFiles[1], sep=''))

# When extracting data from the ncdf file, it will be necessary to utilize the offsets
# for each dimension, otherwise it'll read the whole dataset out
# but this will require some manual determination of which values we want
