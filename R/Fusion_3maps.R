# load libraries

if (!require(randomForest)) install.packages('randomForest')
if (!require(robust)) install.packages('robust')


######  INPUT DATA  ##########

#Gal <- raster(paste('./Input_Maps/Gal_1km_', cont, '.tif', sep=""))
Gal <- raster("./Maps/Gallaun/1km/bmAGB_Gallaun_crop.tif")
IIASA <- raster("./Maps/IIASA/1km/bmAg_IIASA2010_crop.tif")
Thur <- raster('./Maps/Thurner/1km/bmAg_Thurner_1km_crop.tif')  

#ref <- raster(paste('./Reference/Ref_', cont, '.tif', sep=""))
ref <- raster('./Maps/Ref/ref_ras.tif')

vcf <- raster("./Covariates/Landsat_VCF_2005/VCF_crop.tif")
#vcf <- raster(paste('./Strata/VCF_', cont, '.tif', sep=""))

#hei <- raster(paste('./Strata/HEI_', cont, '.tif', sep=""))
hei <- raster("./Covariates/Height/Height_crop.tif")

#cci <- raster(paste('./Strata/CCI_', cont, '.tif', sep=""))
cci <- raster('./Covariates/CCI_2005/CCI_crop.tif')

#water <- raster(paste('./Strata/WATER/Water_1km_', cont, '.tif', sep=""))
water <- raster("./Covariates/CCI_Water/Water_NL_Crop.tif") # AGGREGATE FIRST

###### WATER MASK APPLIED TO GalTCHI & ThurCINI ############

Gal[water == 2] <- NA
Thur[water == 2] <- NA
IIASA[water == 2] <- NA
rm(water)


##################################################################################
###########################   ERROR STRATA   #####################################
##################################################################################

Gal.er <- ref - Gal
Thur.er <- ref - Thur
IIASA.er <- ref - IIASA

#Gal.er <- Gal
#Thur.er <- Thur

error.map <- stack(Gal.er, Thur.er, IIASA.er, hei, vcf, Gal, Thur, IIASA, cci)
names(error.map) <- c("Gal.er", "Thur.er","IIASA.er", "hei", "vcf", "Gal", "Thur", "IIASA", "cci")
#rm(cci)

error <- as.data.frame(getValues(error.map))
error <- error[complete.cases(error),]
colnames(error) <- c("Gal.er", "Thur.er","IIASA.er", "hei", "vcf", "Gal", "Thur", "IIASA", "cci")
error$cci <- as.factor(error$cci)


### MODEL MAP ERROR with Random Forest  
# NOTE: mtry=2 to have comparable R2 without GLC2000, otherwise with 5 variables mtry=1 by default and the R2 drops)

Gal.rf.f <- formula(Gal.er ~ hei + cci + vcf + Thur + Gal + IIASA)
Thur.rf.f <- formula(Thur.er ~ hei + cci + vcf + Thur + Gal + IIASA)
IIASA.rf.f <- formula(IIASA.er ~ hei + cci + vcf + Thur + Gal + IIASA)
#Thur.rf.f <- formula(Thur.er ~ hei + vcf + Gal + Thur + cci)
set.seed(55)
Gal.rf <- randomForest(Gal.rf.f, error, importance=F, mtry=2)
set.seed(55)
Thur.rf <- randomForest(Thur.rf.f, error, importance=F, mtry=2)
set.seed(55)
IIASA.rf <- randomForest(IIASA.rf.f, error, importance=F, mtry=2)
# Gal.rf
# Thur.rf

rf.out <- matrix(0, nrow=3, ncol=2, dimnames=list(c("Predict_Error_Thurner", "Predict_Error_Gal", "Predict_Error_IIASA"), c("RMSE", "R2")))
rf.out[1,1] <- sqrt(Thur.rf$mse[500])
rf.out[1,2] <- mean(Thur.rf$rsq[500])
rf.out[2,1] <- sqrt(Gal.rf$mse[500])
rf.out[2,2] <- mean(Gal.rf$rsq[500])
rf.out[3,1] <- sqrt(IIASA.rf$mse[500])
rf.out[3,2] <- mean(IIASA.rf$rsq[500])

png(filename=paste("./Results/Strata/RF_Predict_Errors_VarImpPlot.png", sep=""), width=960, height=480, res=96)
opar <- par(mfrow=c(1,3))
varImpPlot(Gal.rf)
varImpPlot(Thur.rf)
varImpPlot(IIASA.rf)
dev.off()

varImp <- cbind(importance(Gal.rf), importance(Thur.rf), importance(IIASA.rf))
colnames(varImp) <- c('Gal_Imp', 'Thur_Imp', 'IIASA_Imp')
write.csv(varImp, paste("./Results/Strata/IncNodePur.csv", sep=""))
rm(varImp)

png(filename=paste("./Results/Strata/RF_Predict_Errors_CrossVal.png", sep=""), width=960, height=480, res=96)
opar <- par(mfrow=c(1,3))
plot(error$Gal.er,predict(Gal.rf, newdata=error), main="Cross Validation Error Gal", xlab="Error", ylab="Predicted error", xlim=c(min(error$Gal.er), max(error$Gal.er)), ylim=c(min(error$Gal.er),max(error$Gal.er)))
abline(0,1)
plot(error$Thur.er,predict(Thur.rf, newdata=error), main="Cross Validation Error Thur", xlab="Error", ylab="Predicted error", xlim=c(min(error$Thur.er), max(error$Thur.er)), ylim=c(min(error$Thur.er),max(error$Thur.er)))
abline(0,1)
plot(error$Thur.er,predict(IIASA.rf, newdata=error), main="Cross Validation Error IIASA", xlab="Error", ylab="Predicted error", xlim=c(min(error$IIASA.er), max(error$IIASA.er)), ylim=c(min(error$IIASA.er),max(error$IIASA.er)))
abline(0,1)
dev.off()

### PREDICT MAP ERROR

error.dat <- dropLayer(error.map, 1:3)
Gal.er <- predict(error.dat, Gal.rf, type= 'response', filename="./Results/Strata/Gal_Error_RF.tif", datatype='INT2S', overwrite=T, progress='text')
Thur.er <- predict(error.dat, Thur.rf, type= 'response', filename="./Results/Strata/Thur_Error_RF.tif", datatype='INT2S', overwrite=T, progress='text')
IIASA.er <- predict(error.dat, IIASA.rf, type= 'response', filename="./Results/Strata/IIASA_Error_RF.tif", datatype='INT2S', overwrite=T, progress='text')

rm(error.map)

err.all <- stack(Gal.er, Thur.er, IIASA.er)
mydata <- as.data.frame(getValues(err.all))
mydata.k <- mydata[complete.cases(mydata),]  # Remove NA, K-means cannot handle NAs
rm(list=ls(pattern='er'))
rm(list=ls(pattern='.rf'))

# K cluser
set.seed(55)
fit.k <- kmeans(mydata.k, 3)      # 8 cluster solution

# Convert to raster
mydata$cluster <- NA
mydata[rownames(mydata.k), "cluster"] <- fit.k$cluster
Strata <- Gal
values(Strata) <- mydata$cluster
dataType(Strata) <- 'INT1U'

rm(list=ls(pattern='mydata'))
rm(fit.k)
if (exists('wss')) rm(wss)
Sys.time()


########  PREDICT STRATA FOR NO DATA AREAS  #########

## Remove NA due to: 1. NA present in Thurcini -  2. Predicted Strata Error is NA for the categorical predictors (CCI) without reference data
## NA values in strata are now predicted using datasets without NA (Gal, vcf, hei) 
## Water (applied to Galtchi) was set to NA and remains NA (to have No strata in water)
## The RF model is trained using 10000 pixel randomly selected

str.map <- stack(Strata, Gal, vcf, hei)
names(str.map) <- c("strata", "Gal", "vcf", "hei")
#rm(hei)
#rm(vcf)

dat <- as.data.frame(getValues(str.map))
dat <- dat[complete.cases(dat),]
colnames(dat) <- c("strata", "Gal", "vcf", "hei")
dat$strata <- as.factor(dat$strata)
set.seed(30)
dat <- dat[sample(nrow(dat), 10000), ]

str.f <- formula(strata ~ Gal + vcf + hei)
set.seed(55)
str.rf <- randomForest(str.f, dat, importance=T)
varImpPlot(str.rf)

rf.out[2,1] <- 1-sum(diag(table(predict(str.rf), str.rf$y)))/sum(table(predict(str.rf), str.rf$y))
rf.out[2,2] <- NA
write.csv(rf.out, paste("./Results/Strata/RF_Predict_Errors_EU.csv", sep=""))
rm(rf.out)
rm(dat)

str.map <- dropLayer(str.map, 1)
str.map[[1]] <- mask(str.map[[1]], mask=Strata, inverse=TRUE, progress='text')    # to make predictions only where Strata is NA
strata.na <- predict(str.map, str.rf, type = 'response', progress='text')


## Apply predicted Strata in NA areas

strata <- stack(Strata, strata.na) 

strata.f <- function(x) {
  i <- is.na(x[,1])
  x[i,1] <- x[i,2]
  return(x[,1])
}
Strata.fn <- 'NL'
Strata <- calc(strata, fun=strata.f, filename=paste("./Results/Strata/Strata_", Strata.fn, ".tif", sep=""), datatype='INT1U', overwrite=T, progress='text')
names(Strata) <- paste("Strata_", Strata.fn, sep="")
rm(list=ls(pattern='str.'))

Sys.time()

-------------------------
  
Strata <- raster('./Results/Strata/Strata_NL.tif')


-------------------------
  
  ##################################################################################
  ###########################   REFERENCE DATA  ####################################
##################################################################################

########  CONSOLIDATE REFERENCE DATA #####  START HERE FOR STRATA AS GLC2000 OR VCF10
# Strata <- raster(paste("./Results/Strata/Strata_", Strata.fn, ".tif", sep=""))

Strata.s <- 3    
Strata.names <- 1:Strata.s

## Each NFI has a unique code, and there is a raster map where each plot of the NFI has the code value (e.g.: all plots in Spain have value 1, in France is 2, etc.)
codes <- raster(paste('./Reference/Ref_code_', cont, '.tif', sep=""))
code.names <- as.character(read.csv(paste("./Reference/Codes_CAM.csv", sep=""))[,2])
code.n <- maxValue(codes)

code.names <- 'NL'

# temp fix for testing
codes <- Gal
codes[codes] <- 1
codes[codes == 0] <- NA

## mask out Ref data in Water (NA in Strata)

ref[is.na(Strata)] <- NA


## Analyze distribution of Reference data by Strata

ref.str <- stack(ref, Strata, codes)
ref.dat <- as.data.frame(getValues(ref.str))
ref.dat <- ref.dat[complete.cases(ref.dat),]
colnames(ref.dat) <- c("AGB", "Strata", "Code")
str.code <- table(ref.dat$Code, ref.dat$Strata, dnn = c("Code", "Strata"))
rowSums(str.code)
colSums(str.code)

# adjust str.code in case ref data are missing for some strata (otherwise str.code will not have the Strata with no ref data)

if (ncol(str.code) < Strata.s) {
  str.code.new <- matrix(0, nrow=code.n, ncol=Strata.s)
  str.code <- as.matrix(str.code)
  for (i in 1:ncol(str.code)) {
    n.col <- as.numeric(colnames(str.code)[i])
    str.code.new[,n.col] <- str.code[,i]  }
  str.code <- str.code.new
  rm(list=c('str.code.new', 'n.col'))
}


# Compute % of reference data

str.code.x <- sweep(str.code, MARGIN = 2, colSums(str.code) * 0.01, '/')
str.code.x[is.nan(str.code.x)] <- 0
print(str.code.x, digits = 0)


# Plotting Ref data by Strata

Strata.fn <- 'NL' 
png(filename=paste("./Results/Reference/Reference_Data_", Strata.fn, "_Orig.png", sep=""))
plt <- barplot(str.code, beside=F, main=paste("Reference Data by ", Strata.fn, " strata - original", sep=""), xlab="Strata", legend = code.names, names.arg = Strata.names,
               cex.names = 0.7, col=c("red", "orange", "yellow", "green", "blue", "violet", "pink" , "cyan", "gray", "forestgreen", "orange3"), args.legend = list(x = "topright", cex=0.7))
text(plt, colSums(str.code), labels = colSums(str.code), pos = 3)
dev.off()

##################################################################################
###########################   BIOMASS FUSION   ###################################
##################################################################################

####### CALCULATE THE WEIGHTS & BIAS

## Input data

# saa <- raster(paste('./Input_Maps/saa_1km_', cont, '.tif', sep=""))
# bac <- raster(paste('./Input_Maps/bac_1km_', cont, '.tif', sep=""))  
# ref <- raster(paste('./Results/Reference/Ref_', Strata.fn, '_Cons.tif', sep=""))
# Strata <- raster(paste('./Results/Strata/Strata_', Strata.fn, ".tif", sep=""))
# Strata.s <- 8    
# Strata.n <- 9
# Strata.names <- 1:Strata.s
# water <- raster(paste('./Strata/WATER/Water_1km_', cont, '.tif', sep=""))
# saa[water == 1] <- NA
# bac[water == 1] <- NA
# rm(water)


# Error Maps

Gal.er <- ref - Gal
Thur.er <- ref - Thur
IIASA.er <- ref - IIASA

error.map <- stack(Strata, Gal.er, Thur.er, IIASA.er)

error <- as.data.frame(as.matrix(error.map, na.rm = TRUE))
error <- error[complete.cases(error),]
colnames(error) <- c("Strata", "Gal_er", "Thur_er", "IIASA_er")
error$Strata <- as.factor(error$Strata)
#rm(saa.er)
#rm(bac.er)
#rm(error.map)


# Plot Errors of IIASA and Thurner

png(filename="./Results/Reference/Errors of IIASA_Thur.png")
plot(error$IIASA_er, error$Thur_er, main=paste("Errors of input maps"), xlab="IIASA Error", ylab="Thur Error", xlim=c(min(error$IIASA_er, error$Thur_er), max(error$IIASA_er, error$Thur_er)), ylim=c(min(error$Thur_er, error$IIASA_er), max(error$Thur_er, error$IIASA_er)))
abline(0,1)
dev.off()


### Calculate the bias

Gal.bias <- tapply(error$Gal_er, error$Strata, mean, na.rm = TRUE)
Thur.bias <- tapply(error$Thur_er, error$Strata, mean, na.rm = TRUE)
IIASA.bias <- tapply(error$IIASA_er, error$Strata, mean, na.rm = TRUE)
bias <- cbind(Gal.bias, Thur.bias, IIASA.bias)
bias <- bias[1:Strata.s,]
rm(list=ls(pattern=".bias"))


### Calculate variance-covariance matrix and weight matrix using a ROBUST ESTIMATOR

X <- c(1,1,1)                      ##! adapt if more than 2 input maps. with 3 maps it should be X <- c(1,1,1)

# !! Wont Run becuase Gal dataset has a mean error of 0, due to the reference dataset being used in the map creation  !! 
cov <- vector('list',length = Strata.s)
for (i in 1:Strata.s){
  cov[[i]] <- covRob(error[error$Strata==i,2:4]) # set length (-:-) to the selected datasets
  cov[[i]] <- cov[[i]]$cov}

w <- vector('list',length = Strata.s)
for (i in 1:Strata.s){
  w[[i]] <- solve(t(X) %*% solve(cov[[i]]) %*% X) %*% t(X) %*% solve(cov[[i]]) }

weight <- matrix(1:(length(X)*Strata.s), ncol=(length(X)))
for (i in 1:Strata.s){ weight[i,] <- w[[i]] }
colnames(weight) <- c("Gal.weight", "Thur.weight", "IIASA.weight")
rownames(weight) <- c(1:Strata.s)


## Finalize Bias and Weights

weight[weight < 0] <- 0   # set weights to 0 - 1 limits
weight[weight > 1] <- 1   # set weights to 0 - 1 limits


#------------------------------------------------