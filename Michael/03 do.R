stations = read.csv("data/station_info.csv", stringsAsFactors = FALSE)
trainy = read.csv("data/train.csv")
pred = read.csv("data/test/sampleSubmission.csv", row.names = 1)

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

# Loop through each station in station_info.csv to train model and then make prediction
for(i in 1:nrow(stations)){
  #   Find the 4 closest GEFS pts to this station
  GEFS = closestGEFS(stations[i,1])
  
  #   Load the training data
  train = list()
  
  for(j in 1:nrow(GEFS)){
    train[[j]] = lapply(varnames, getVariable, lon = GEFS$lon[j], lat = GEFS$lat[j])
    train[[j]] = scale.default(do.call(cbind, train[[j]]))
    train[[j]] = as.data.frame(train[[j]])  
  }
  
  train = do.call(cbind, train)
  colnames(train) = 1:ncol(train)         #prevent duplicated variable names after merging from different GEFS pts
  train[is.na(train)] = 0                 # set missing values to zero
  
  #   filter out correlated variables
  trainCor = cor(train)
  highlyCor = findCorrelation(trainCor, cutoff=0.95)
  train = train[,-highlyCor]
  
  #   attach response variable
  train$y = trainy[, stations[i,1]]
  
  #   fit model to training set
  fit = lm(y ~ ., train)
  
  #   read in test set
  test = list()
  
  for(j in 1:nrow(GEFS)){
    test[[j]] = lapply(varnames, getVariable, lon = GEFS$lon[j], lat = GEFS$lat[j], path = "data/test")
    test[[j]] = scale.default(do.call(cbind, test[[j]]))
    test[[j]] = as.data.frame(test[[j]])  
  }
  
  test = do.call(cbind, test)
  colnames(test) = 1:ncol(test)         #prevent duplicated variable names after merging from different GEFS pts
  
  #   filter out correlated variables
  test = test[, -highlyCor]
  test[is.na(test)] = 0                # set NA fields to 0
  
  pred[, i] = predict(fit, test)
  
  print(i)
}

write.csv(pred, "submission.csv")

# Read in training data and fit model

# Read in test data and make prediction
# Store prediction in list with name of station  

# Merge predictions for eaach station and write to CSV 