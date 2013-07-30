getVariable <- function(lat, lon, var, path = "data/train"){
  
  filename = Sys.glob(file.path(path, paste0(var, "*", collapse="")))
  nc = open.ncdf(filename)
  
  varname = switch(var,
         apcp_sfc = "Total_precipitation",
         dlwrf_sfc = "Downward_Long-Wave_Rad_Flux",
         dswrf_sfc = "Downward_Short-Wave_Rad_Flux",
         pres_msl = "Pressure",
         pwat_eatm = "Precipitable_water",
         spfh_2m = "Specific_humidity_height_above_ground",
         tcdc_eatm = "Total_cloud_cover",
         tcolc_eatm = "Total_Column-Integrated_Condensate",
         tmax_2m = "Maximum_temperature",
         tmin_2m = "Minimum_temperature",
         tmp_2m = "Temperature_height_above_ground",
         tmp_sfc = "Temperature_surface",
         ulwrf_sfc = "Upward_Long-Wave_Rad_Flux_surface",
         ulwrf_tatm = "Upward_Long-Wave_Rad_Flux",
         uswrf_sfc = "Upward_Short-Wave_Rad_Flux",)
  
  lon = lon + 107
  lat = lat - 30
  
  data = get.var.ncdf(nc, varname, start = c(lon, lat, 1, 1, 1), count = c(1, 1, -1, -1, -1))
  data = aperm(data, c(3,2,1))
  dim(data) = c(5113, 55)
  rownames(data) = get.var.ncdf(nc, "intTime")
  
  return(data)
}

getVariableName <- function(var){
  filename = Sys.glob(file.path("data/train", paste0(var, "*", collapse="")))
  nc = open.ncdf(filename)
  
  print(nc)
}