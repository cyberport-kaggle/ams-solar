
df <- read.csv('../data/train.csv')
colnames(df)

library(RNetCDF)

cdf <- open.nc('../data/test/spfh_2m_latlon_subset_20080101_20121130.nc')
var.inq.nc(cdf, 'time')
