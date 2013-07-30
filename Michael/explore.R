source("Michael/02 func.R")

# take a look at the solar readings
readings = read.csv("data/train.csv")
str(readings)
summary(readings)

plot(readings$Date, readings$ACME)

# get familair with the NetCDF data format
library(ncdf)
nc = open.ncdf("data/train/apcp_sfc_latlon_subset_19940101_20071231.nc")
print(nc)
apcp = get.var.ncdf(nc, "Total_precipitation", start=c(1, 1, 1, 1, 1), count=c(1,1,-1,-1,-1))
time = get.var.ncdf(nc, "intTime")
b = get.var.ncdf(nc, "intValidTime")

# variables that are forecasted by GEFS 
varnames = c("apcp_sfc",
             "dlwrf_sfc",
             "dswrf_sfc",
             "pres_msl",
             "pwat_eatm",
             "spfh_2m",
             "tcdc_eatm",
             "tcolc_eatm",
             "tmax_2m",
             "tmin_2m",
             "tmp_2m",
             "tmp_sfc",
             "ulwrf_sfc",
             "ulwrf_tatm",
             "uswrf_sfc")

# Variable    NCDF Name	Description
# apcp_sfc	  Total_precipitation	3-Hour accumulated precipitation at the surface
# dlwrf_sfc	  Downward_Long-Wave_Rad_Flux	Downward long-wave radiative flux average at the surface
# dswrf_sfc	  Downward_Short-Wave_Rad_Flux	Downward short-wave radiative flux average at the surface
# pres_msl	  Pressure	Air pressure at mean sea level
# pwat_eatm	  Precipitable_water	Precipitable Water over the entire depth of the atmosphere
# spfh_2m	    Specific_humidity_height_above_ground	Specific Humidity at 2 m above ground
# tcdc_eatm	  Total_cloud_cover	Total cloud cover over the entire depth of the atmosphere
# tcolc_eatm	Total_Column-Integrated_Condensate	Total column-integrated condensate over the entire atmos.
# tmax_2m	    Maximum_temperature	 Maximum Temperature over the past 3 hours at 2 m above the ground
# tmin_2m	    Minimum_temperature	 Mininmum Temperature over the past 3 hours at 2 m above the ground
# tmp_2m	    Temperature_height_above_ground	 Current temperature at 2 m above the ground
# tmp_sfc	    Temperature_surface	 Temperature of the surface
# ulwrf_sfc	  Upward_Long-Wave_Rad_Flux_surface	 Upward long-wave radiation at the surface
# ulwrf_tatm	Upward_Long-Wave_Rad_Flux	 Upward long-wave radiation at the top of the atmosphere
# uswrf_sfc	  Upward_Short-Wave_Rad_Flux	 Upward short-wave radiation at the surface

# Longitude and Latitude of forecasts
# longitude -106 -105 -104 -103 -102 -101 -100 -99  -98  -97  -96  -95  -94  -93  -92  -91
# latitude 31 32 33 34 35 36 37 38 39

library(doMC)
registerDoMC(cores = 4)

# train a toy model for the MAYR station.  I picked MAYR it was the station located closest to a forecast grid point
# get data into a dataframe format.  each row is different time and each column represents difference ensembles, forecast
# time and variables 
MAYR = lapply(varnames, getVariable, lon = -99, lat = 37)
MAYR = scale.default(do.call(cbind, MAYR))
MAYR = as.data.frame(MAYR)

# load the response variable and add it to the MAYR data frame 
train = read.csv("data/train.csv")
MAYR$y = train$MAYR

# fit a linear model to MAYR 
fit = lm(y ~ ., MAYR)

# linear model seems to be a pretty good fit to the data. Rsquared of 0.86 
plot(fit$fitted.values, MAYR$y)
summary(fit)

# tried to train a SVM and RandomForest.  Took forever, never completed.
# fitControl = trainControl(method="cv", number=5)
# svmFit = train(y ~ ., data=MAYR, method="rf", 
#                trControl=fitControl, tuneLength=5)