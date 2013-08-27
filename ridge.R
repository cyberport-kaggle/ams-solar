library(caret)
library(data.table)
library(reshape2)
library(MASS)
library(parcor)
library(foreach)
library(doMC)
library(profr)
registerDoMC(cores = 4)

factors = c('dswrf_sfc','dlwrf_sfc','uswrf_sfc','ulwrf_sfc','ulwrf_tatm','pwat_eatm','tcdc_eatm','apcp_sfc',
            'pres_msl','spfh_2m','tcolc_eatm','tmax_2m','tmin_2m','tmp_2m','tmp_sfc')

loadRidgeGEFS = function(train = TRUE){
  ifelse(train, path <- "data/train", path <- "data/test")
  
  res = data.table()
  
  for (f in factors){
    cat("Loading", f, "\n")
    load(Sys.glob(file.path(path, paste0(f, "*", ".Rdata"))))
    tbl = tbl[, mean(value), by=list(date, lat, lon)]
    tbl = data.table(dcast(tbl, date~lat+lon, value.var = "V1"), key= "date")
    
    #   concatenate f with the column names first
    setnames(tbl, names(tbl)[-1], paste(f, names(tbl)[-1]))
    
    if (length(res) == 0) res = tbl
    else res = merge(res, tbl, by = "date")
  }
  return(res)
}

res = loadGEFS(train = TRUE)

# Remove zero variance columns
zv = nearZeroVar(res, saveMetrics=TRUE)$zeroVar
res = res[, which(!zv), with = FALSE]

# Scale and center results 
preProcValues <- preProcess(res[, setdiff(colnames(res), "date"), with = FALSE])
res <- data.table(cbind(res[, date], 
                       predict(preProcValues, res[, setdiff(colnames(res), "date"), with = FALSE])))
setnames(res, colnames(res)[1], "date")
setkey(res, date)

save(res, file = "ridge/trainData.RData")

# Load training labels 
ytrain = data.table(read.csv(file="data/train.csv"), key="Date")

# Convert response and predictors into matrices 
X = cbind(1, (as.matrix(res[, setdiff(colnames(res), c("date", "y")), with = FALSE])))
Y = as.matrix(ytrain[,-1, with = FALSE])
lambda = 128

# Calculate coefficients of ridge regression 
beta = solve(t(X) %*% X + lambda * diag(ncol(X))) %*% t(X) %*% Y

# load Test data
testData = loadGEFS(train = FALSE)
save(testData, file = "ridge/testData.RData")

# Convert into matrix format
Xtest = as.matrix(predict(preProcValues, testData[, setdiff(colnames(testData), "date"), with = FALSE]))
Xtest = cbind(1, Xtest)

# Do ridge regrassion prediction 
Ytest = Xtest %*% beta

# Get into right format for submission
ans = cbind(testData[,date], Ytest)
colnames(ans)[1] = "Date"
write.csv(ans, "ridge/submission.csv", row.names = FALSE)