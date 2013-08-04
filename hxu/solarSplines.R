require(ncdf4)


loadMesonetData <- function(filename, stationFilename = "station_info.csv") {
    data <- read.csv(file.path(data_path, filename))
    station_data <- read.csv(file.path(data_path, stationFilename))
    list(data = data[,-1], dates = data[,1], station_data = station_data)
}


loadData <- function(filename) {
  nc_open(filename)
}


getGrid <- function(data, date, fHour, eMember) {
    dateIdx <- which(ncvar_get(data, varid = "intTime") == date)
    fIdx <- which(ncvar_get(data, varid = "fhour") == fHour)
    eIdx <- which(ncvar_get(data, varid = "ens") == eMember)
    ncvar_get(data)[,,fIdx, eIdx, dateIdx]
}

getDailyMeanSumGrid <- function(data, date) {
    dateIdx <- which(ncvar_get(data, varid = "intTime") == date)
    res <- ncvar_get(data)[, , , , dateIdx]
    res  <- apply(res, c(1:2, 4), sum)
    apply(res, 1:2, mean) * 3600 * 3
}



buildSplines <- function(data, grid, station_data) {
    outdata <- vector(mode = "numeric",
                      length = nrow(station_data))
    cat(length(outdata), "\n")
    for ( i in seq_along(outdata) ) {
        slat <- station_data$nlat[i]
        slon <- station_data$elon[i]
        nearlat <- which(abs(ncvar_get(data, 'lat') - slat) < 2)
        nearlon <- which(abs(ncvar_get(data, 'lon') - slon - 360)  < 2)
        spline1 <- vector(mode = "numeric",
                          length = length(nearlon))
        for (l in seq_along(nearlat))
            spline1[l] <-  Spline(grid[nearlon, nearlat[l]],
                                  (slon - floor(slon)))
        outdata[i] <- Spline(spline1, (slat - floor(slat)))
    }
    outdata
}




Spline <- function(y, xi) {
    0.5*((2 * y[2]) + (y[3]-y[2]) * xi +
         (-y[4] + 4 * y[3] - 5 * y[2] + 2 * y[1]) * xi^2 +
         (y[4] - 3 * y[3] + 3 * y[2] - y[1]) * xi^3)
}

data_path <- "./"
data <- loadData(file.path(data_path, "dswrf_sfc_latlon_subset_20080101_20121130.nc")
listsubmit <- loadMesonetData(file.path(data_path, "sampleSubmission.csv"))
dates <- listsubmit$dates
header <- c("Date", names(listsubmit$data))

outdata <- lapply(dates, function(date) {
    cat(date, "\n")
    grid <- getDailyMeanSumGrid(data, date * 100)
    buildSplines(data, grid, listsubmit$station_data)
})

outdata <- do.call(rbind, outdata)
outdata <- cbind(dates, outdata)
colnames(outdata) <- header
write.csv(outdata, "submission.csv", row.names = FALSE)
