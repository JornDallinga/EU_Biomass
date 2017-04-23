# load libraries

if (!require(randomForest)) install.packages('randomForest')
if (!require(robust)) install.packages('robust')
#if (!require(ranger)) install.packages('ranger')


######  INPUT DATA  ##########

#Gal <- raster(paste('./Input_Maps/Gal_1km_', cont, '.tif', sep=""))
Gal <- raster("./Maps/Gallaun/1km/bmAg_JR2000_ll_1km_eur.tif")
Gal <- raster("./Maps/Barredo/barredo_Alligned.tif")
Gal[Gal < 0] <- NA

# cant get Bar en Gal extents to match, fix required
#Bar <- projectRaster(Bar, Gal, method = 'ngb')
#writeRaster(Bar, filename = "./Maps/Barredo/barredo_reproj.tif")
#

#-------
IIASA <- raster("./Maps/IIASA/1km/bmAg_IIASA2010_Alligned.tif")
#IIASA <- crop(IIASA, Gal)

Thur <- raster('./Maps/Thurner/1km/bmAg_Thurner_1km_Alligned.tif')  
#Thur <- crop(Thur, Gal)
#ref <- raster(paste('./Reference/Ref_', cont, '.tif', sep=""))
ref <- raster('./Maps/Ref/ref_ras_EU2.tif')

#vcf <- raster("./Covariates/MODIS_VCF_2005/transformed/Mosaic/MODIS_VCF_Mosaic.tif")
#align_rasters("./Covariates/MODIS_VCF_2005/transformed/Mosaic/MODIS_VCF_Mosaic.tif", "./Maps/Gallaun/1km/bmAg_JR2000_ll_1km_eur.tif", dstfile = "./Covariates/MODIS_VCF_2005/transformed/Mosaic/MODIS_VCF_Mosaic_C.tif")
#vcf[vcf < 0] <- NA
#vcf[vcf > 100] <- NA
#writeRaster(vcf, filename = "./Covariates/MODIS_VCF_2005/transformed/Mosaic/MODIS_VCF_Mosaic_EU.tif")
vcf <- raster("./Covariates/Outputs/MODIS_VCF/Mosaic/Mosaic_aggre_align.tif")


#hei <- raster(paste('./Strata/HEI_', cont, '.tif', sep=""))
hei <- raster("./Covariates/Height/Height_align.tif")
#hei <- crop(hei, Gal)

#cci <- raster(paste('./Strata/CCI_', cont, '.tif', sep=""))
#cci <- raster('./Covariates/CCI_2005/CCI_2005_resample.tif')
#cci[cci < 0] <- NA
#writeRaster(cci, filename = './Covariates/CCI_2005/CCI_2005_resample_EU.tif')
cci <- raster("./Covariates/Outputs/CCI_2005/CCI_Mul_align.tif")



#water <- raster(paste('./Strata/WATER/Water_1km_', cont, '.tif', sep=""))
water <- raster("./Covariates/Outputs/CCI_Water/CCI_Water_aggr_Resam_R.tif")

###### WATER MASK APPLIED TO GalTCHI & ThurCINI ############

Gal[water == 2] <- NA
Thur[water == 2] <- NA
IIASA[water == 2] <- NA
#rm(water)

#c <- resample(x = Gal, Thur)

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
system.time(Gal.rf <- randomForest(Gal.rf.f, error, importance=F, mtry=2))
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

#rm(error.map)

err.all <- stack(Gal.er, Thur.er, IIASA.er)
mydata <- as.data.frame(getValues(err.all))
mydata.k <- mydata[complete.cases(mydata),]  # Remove NA, K-means cannot handle NAs
#rm(list=ls(pattern='er'))
#rm(list=ls(pattern='.rf'))

# K cluser
set.seed(55)
fit.k <- kmeans(mydata.k, 3)      # 8 cluster solution

# Convert to raster
mydata$cluster <- NA
mydata[rownames(mydata.k), "cluster"] <- fit.k$cluster
Strata <- Gal
values(Strata) <- mydata$cluster
dataType(Strata) <- 'INT1U'

#rm(list=ls(pattern='mydata'))
#rm(fit.k)
#if (exists('wss')) rm(wss)
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
#rm(rf.out)
#rm(dat)

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
Strata.fn <- 'EU'
Strata <- calc(strata, fun=strata.f, filename=paste("./Results/Strata/Strata_", Strata.fn, ".tif", sep=""), datatype='INT1U', overwrite=T, progress='text')
names(Strata) <- paste("Strata_", Strata.fn, sep="")
rm(list=ls(pattern='str.'))

Sys.time()

#-------------------------
  
Strata <- raster(paste("./Results/Strata/Strata_", Strata.fn, ".tif", sep=""))
#Strata <- crop(Strata, ex)

#-------------------------
  
  ##################################################################################
  ###########################   REFERENCE DATA  ####################################
##################################################################################

########  CONSOLIDATE REFERENCE DATA #####  START HERE FOR STRATA AS GLC2000 OR VCF10
# Strata <- raster(paste("./Results/Strata/Strata_", Strata.fn, ".tif", sep=""))

Strata.s <- 3    
Strata.names <- 1:Strata.s

## Each NFI has a unique code, and there is a raster map where each plot of the NFI has the code value (e.g.: all plots in Spain have value 1, in France is 2, etc.)
codes <- raster(paste('./Reference/Ref_code_EU.tif', sep=""))
codes <- crop(codes, Gal)
code.names <- as.character(read.csv(paste("./Reference/Codes_EU.csv", sep=""))[,2])
code.names <- code.names[unique(codes)]
code.n <- maxValue(codes)

#code.names <- 'NL' Single map testing

# temp fix for testing
#codes <- Gal
#codes[codes] <- 1
#codes[codes == 0] <- NA

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

Strata.fn <- 'EU' 
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

png(filename="./Results/Reference/Errors of IIASA_Thur.png") # Add Gal map and create 3 combi maps
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

# !! Wont Run because Gal dataset has a mean error of 0, due to the reference dataset being used in the map creation  !! 
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

##### for EUROPE: this section until line 561 is not essential and may be skipped

### calculate Error Variance of the Fused map

v <- vector('list', length = Strata.s)
for (i in 1:Strata.s){  v[[i]] <- solve(t(X) %*% solve(cov[[i]]) %*% X) }

v.err <- matrix(1:Strata.s)
for (i in 1:Strata.s){ v.err[i] <- v[[i]] }

colnames(v.err) <- c("Fus_Var")
rownames(v.err) <- c(1:Strata.s)
rm(list=c('cov', 'i', 'v', 'w', 'X'))


### Compute Error Variance of Input maps

Gal.var <- aggregate(error$Gal_er, by = list(error$Strata), FUN="var")
Thur.var <- aggregate(error$Thur_er, by = list(error$Strata), FUN="var")
IIASA.var <- aggregate(error$IIASA_er, by = list(error$Strata), FUN="var")

### Compute n. pixel per strata

n.pix <- matrix(1:Strata.s)
for (i in 1:Strata.s) { n.pix[i] <- freq(Strata, useNA='no', value = i) }
n.pix.w <- n.pix / (sum(n.pix))


### Compile Error Variances

err.var <- cbind(Gal.var, Thur.var[,2], IIASA.var[,2], v.err, n.pix, n.pix.w)
colnames(err.var) <- c("Strata", "Gal.var", "Thur.var","IIASA.var","Fus_Var", "N_Pix", "N_Pix_w")
err.var$Gal_Var_w <- err.var$Gal.var * err.var$N_Pix_w
err.var$Thur_Var_w <- err.var$Thur.var * err.var$N_Pix_w
err.var$IIASA_Var_w <- err.var$IIASA.var * err.var$N_Pix_w
err.var$Fus_Var_w <- err.var$Fus_Var * err.var$N_Pix_w


# Compile Bias, Weights, Error variance, and add Strata 9

fus.par <- aggregate(error$Strata, by = list(error$Strata), FUN="length")
colnames(fus.par) <- c("Strata", "N")
fus.par$Strata <- as.numeric(as.character(fus.par$Strata))
fus.par <- cbind(fus.par, bias, weight, err.var[,-1])   # adjust the number of classes
#fus.par[9,] <- c(9, 0, 0, 0, 0.5, 0.5, rep(0,8))
fus.par[Strata.s+1,] <- c(Strata.s+1, 0, 0, 0, 0, 0.5, 0.5, 0.5, rep(0,10)) 
# Should weights be set according to the number of input maps? e.g. 3 maps == weight 0.33? check this
bias <- fus.par[,3:5]                             # Add Strata 9 to bias and weight (for Fusion)
weight <- fus.par[,6:8]
dir.create("./Results/Fused_Map", showWarnings = F)
write.csv(fus.par, paste("./Results/Fused_Map/Bias_Weights_", Strata.fn, ".csv", sep=""), row.names = FALSE)


## Map Uncertainty (Standard deviation of Error of Fused map)

err.rcl <- cbind(fus.par$Strata[1:Strata.s], sqrt(fus.par$Fus_Var[1:Strata.s]))
uncer <- reclassify(Strata, err.rcl, filename=paste("./Results/Fused_Map/Uncertainty_", Strata.fn, ".tif", sep=""), datatype='FLT4S', overwrite=T)

rm(list=c('fus.par', 'i', 'n.pix', 'n.pix.w', 'uncer'))
rm(list=ls(pattern="err"))
rm(list=ls(pattern=".var"))



####### MAP FUSION  

## Set NA in Strata to a value 9 (Strata must have always values in Fusion, no NA)

Strata.n <- Strata.s + 1
Strata[is.na(Strata)] <- Strata.n  

#### FOR EUROPE: this function needs to be adapted to 3 maps instead of 2. After Line 593 there should be a "c <-..." and the map "c" should be added to line 594 with respective bias and weight 
###  FUSION           
IIASA[IIASA < 1] <- NA
Thur[Thur < 1] <- NA
Gal[Gal < 1] <- NA
maps <- stack(Gal, Thur, IIASA, Strata)


# # Biomass fusion function is broken...
# # the 'Strata.n' represents the number of stacked layers
# biomass.fusion <- function(x) {
#   result <- matrix(NA, dim(x)[1], 1) # create result matrix which stores the final biomass values
# 
#   for (n in 1:Strata.n) { # loop through number of strata's
#     ok <- !is.na(x[,Strata.n]) &  x[,Strata.n] == n # select strata
#     
#     #x[ok,1:(Strata.n-1)] <- sweep(x[ok,1:(Strata.n-1)],2,as.matrix(bias[n,]), FUN = "+") # add bias to cell values
#     x[ok,1:(Strata.n-1)] <- sweep(matrix(x[ok,1:(Strata.n-1)], ncol = Strata.n-1),2,as.matrix(bias[n,]), FUN = "+") # add bias to cell values
#     
#     #d <- matrix(NA, dim(x), (Strata.n-1)) # create matrix for # amount of maps to store the weights
#     #for (l in 1:(Strata.n-1)){
#       #d[ok,l] <- weight[n,l] # add weights to non NA values for the # of maps
#     #}
#     d <- matrix(NA, dim(x), (Strata.n-1)) # create matrix for # amount of maps to store the weights
#     d[ok,] <- as.matrix(weight[n,]) # add weights to non NA values for the # of maps
# 
#     #x[ok,][x[ok,1:(Strata.n-1) < 0]] <- NA # set weights to NA if the original map does not have a value
#     #x[ok,][is.na(x[ok,1:(Strata.n-1)]) | x[ok,1:(Strata.n-1)] < 0] <- NA # set weights to NA if the original map does not have a value
#     #[x[ok,1:(Strata.n-1)] == -10] <- NA
#     #x[ok,][x < 0] <- NA
#     #x[x[ok,] < 0] <- NA
#     d[ok,][is.na(x[ok,1:(Strata.n-1)])] <- NA # set weights to NA if the original map does not have a value
#     
#     #d <- d[!is.na(d)]
#     #d[ok,][!is.na(d)]
#     # -- recalculate weights --
#     p <- sum(d, na.rm = T) # calculate sum of weight values
#     pp <- d/p # divide weight values by sum to get the proportion to == 100
#     pp <- as.numeric(pp)
#     
#     g <- x[ok,1:(Strata.n-1)] # get biomass values
#     g[g < 0] <- 0
#     gg <- g[!is.na(g)] # select non Na biomass values
# 
#     #d[ok,][!is.na(d)]/sum(d[ok,][!is.na(d)], na.rm = T) test
#     # weight * bias
#     result[ok] <- as.matrix(as.integer(round(sum(gg * pp, na.rm = T)))) # multiply biomass values with weight and add as result
#     #result[ok] <- as.matrix(as.integer(round(sum(x[ok,][!is.na(x[ok,1:(Strata.n-1)])] * d[ok,][!is.na(d)]/sum(d[ok,][!is.na(d)], na.rm = T), na.rm = T)))) # multiply biomass values with weight and add as result
#     
#     #p <- rowSums(x[ok,1:(Strata.n-1)], na.rm = T) # sum biomass values accross rows for the # number of maps
#     #p <- rowSums(matrix(x[ok,1:(Strata.n-1)], ncol = Strata.n-1), na.rm = T) # sum biomass values accross rows for the # number of maps
#     #p <- sum(abs(x[ok,1:(Strata.n-1)]),na.rm = T)
#    # g <- mutate(Percent = matrix(x[ok,1:(Strata.n-1)] / sum(matrix(x[ok,1:(Strata.n-1)])), na.rm = T)
#     
#     #pp <- abs(x[ok,1:(Strata.n-1)])/p # divide sum of biomass through seperate biomass values to get proportional values
#     
#     #tot_mat <- x[ok,1:(Strata.n-1)] * pp # bias corrected biomass values * weight
#     #result[ok] <- as.matrix(as.integer(round(rowSums(tot_mat, na.rm = T)))) # sum weight applied biomass values and add to result
#   }
#   return(result)   
# }
# 
# 
# #old
# biomass.fusion <- function(x) {
#   result <- matrix(NA, dim(x)[1], 1)
#   for (n in 1:Strata.n) {
#     ok <- !is.na(x[,3]) &  x[,3] == n
#     a <- x[ok,1] + bias[n,1]              # for these pixels, take the values of map 1 and add the bias (output is a subset with only values for this Stratum)
#     b <- x[ok,2] + bias[n,2]
#     c <- x[ok,3] + bias[n,3]
#     
#     aa <- a * weight[n,1]
#     bb <- b * weight[n,2]
#     cc <- c * weight[n,3]
#     
#     
#     result[ok] <- a * weight[n,1] + b * weight[n,2] +  c * weight[n,2]  # compute fused biomass for the pixels belonging to this Stratum
#   }
#   return(result)
# }
# 
# 
# 
# 
# #old
# biomass.fusion <- function(x) {
#   result <- matrix(NA, dim(x)[1], 1)
#   for (n in 1:Strata.n) {
#     ok <- !is.na(x[,4]) &  x[,4] == n
#     g <- x[1:3] + as.matrix(bias[n,])
#     g <- x[ok,1:3] + as.matrix(bias[n,])
#     a <- x[ok,1] + bias[n,1]              # for these pixels, take the values of map 1 and add the bias (output is a subset with only values for this Stratum)
#     b <- x[ok,2] + bias[n,2]
#     c <- x[ok,3] + bias[n,3]
#     result[ok] <- g[1] * weight[n,1] + g[2] * weight[n,2] + g[3] * weight[n,3]  # compute fused biomass for the pixels belonging to this Stratum
#     #result[ok] <- a * weight[n,1] + b * weight[n,2] + c * weight[n,3]  # compute fused biomass for the pixels belonging to this Stratum
#   }
#   return(result)
# }

# working function. but very slow!
biomass.fusion <- function(x) {
  m <- matrix(x, nrow= 1, ncol=4)
  n <- m[,4]
  g <- m[1:(Strata.n-1)] + as.matrix(bias[n,])
  g[g < 0] <- 0
  w <- weight[n,1:(Strata.n-1)]
  w[is.na(g)]<- NA
  p <- sum(w, na.rm = T) # calculate sum of weight values
  pp <- w/p # divide weight values by sum to get the proportion to == 1
  pp <- as.numeric(pp)
  result <- as.integer(round(sum(pp*g, na.rm = T)))
  return(result)
}

m <- matrix(x, nrow= length(bands), ncol=length(dates))

ras <- Strata
ras[Strata == 1] <- Gal


x[1846]
x <- maps
x <- x[605]

r.samp <- sampleRandom(maps, size=(n+20), na.rm=TRUE, sp=FALSE, asRaster=FALSE) 
x <- r.samp[4,]
x <- matrix(x, ncol = 4)

plot(Gal)
e <- drawExtent()
maps1 <- crop(maps,e)


system.time(Fused.map <- calc(maps, fun = biomass.fusion, progress = 'text'))


Fused.map[Fused.map < 0] <- 0
Fused.map[Strata == 4] <- NA

#### for EUROPE: this section until line 690 aims to fill the NA in the fused map. the NA are due to a input biomass map with NA.
## another fused map is created using only the biomass map(s) with values, and then this additional map is mosaicked to the original fused map (of line 600)
## the lines 630-680 are aimed to blend the transition area to avoid artefacts: for Europe probably not appropriate and can be skipped - just merge the results
######## EXTEND TO SAATCHI AREA

### Create a SAATCHI BIAS-ADJ map (Full coverage)

# maps <- stack(Gal, Thur,IIASA, Strata)
# 
# # Should bias in the adj function be the mean of the input maps? since it will only take the values and bias from 1 map. Test below 
# #old
# adj <- function(x) {
#   result <- matrix(NA, dim(x)[1], 1)
#   for (n in 1:Strata.n) {
#     ok <- !is.na(x[,4]) &  x[,4] == n        # identify pixels belonging to a Stratum, without NA (logical: FALSE/TRUE vector for ALL pixels)
#     result[ok] <- x[ok,1] + bias[n,1]        # for these pixels, take the values of map 1 and add the bias (output is a subset with only values for this Stratum)
#   }
#   return(result)
# } 
# # new testing, may contain errors!!
# adj <- function(x) {
#   result <- matrix(NA, dim(x)[1], 1)
#   result1 <- matrix(NA, dim(x)[1], 1)
#   for (n in 1:Strata.n) {
#     ok <- !is.na(x[,4]) &  x[,4] == n        # identify pixels belonging to a Stratum, without NA (logical: FALSE/TRUE vector for ALL pixels)
#     result <- x[ok,1] + bias[n,1]        # for these pixels, take the values of map 1 and add the bias (output is a subset with only values for this Stratum)
#     result[2] <- x[ok,2] + bias[n,2] 
#     result[3] <- x[ok,3] + bias[n,3] 
#     result1[ok] <- mean(result, na.rm = T)
#       }
#   return(result1)
# } 
# 
# Gal.bias.adj <- calc(maps, fun = adj)
# Thur.bias.adj <- calc(maps1, fun = adj)
# 
# Gal.bias.adj[Gal.bias.adj < 0] <- 0
dir.create("Results/Fused_Map/Mosaic/", showWarnings = F)
writeRaster(Fused.map, filename = paste("Results/Fused_Map/Mosaic/Fused.map_", Strata.fn, ".tif", sep=""), datatype='INT1U', overwrite=T) # datatype = FLT4S
# 
# 
# ### FULL COVERAGE: MOSAIC the Baccini extent with the larger Saatchi extent 
# 
# overlap <- 1   # Define overlap area, in degree
# 
# ## Input data
# # Fused.map <- raster(paste("./Results/Fused_Map/Mosaic/Fused_map_", Strata.fn, ".tif", sep=""))
# # saa.bias.adj <- raster(paste("./Results/Fused_Map/Mosaic/Saa_Bias_Adj_", Strata.fn, ".tif", sep=""))
# 
# 
# ### Define Extents
# 
# ## Baccini Extent
# bac.ext <- read.csv('/media/DATA1/avita001/Baccini_extents.csv')
# # bac.ext <- read.csv("G:/GEOCARBON/workspace/Fusion/Baccini_extents.csv")
# 
# map1.ext <- extent(bac)
# if (cont=="AFR") { 
#   map1.ext@ymin <- bac.ext[1,3] }
# if (cont=='SAM') {
#   map1.ext@ymin <- bac.ext[2,3] }
# if (cont=='CAM') {
#   map1.ext@ymax <- bac.ext[3,2] }
# if (cont=='ASIA') {
#   map1.ext@ymax <- bac.ext[4,2] }
# 
# 
# ## Missing Extent (Saatchi - Baccini) + Overlap
# 
# map2.ext <- extent(saa)
# if (extent(saa)@ymin != map1.ext@ymin && map1.ext@ymin < 0) map2.ext@ymax <- map1.ext@ymin + overlap  # for AFR & SAM
# if (extent(saa)@ymin != map1.ext@ymin && map1.ext@ymin > 0) map2.ext@ymax <- map1.ext@ymin - overlap
# if (extent(saa)@ymax != map1.ext@ymax && map1.ext@ymax < 0) map2.ext@ymin <- map1.ext@ymax + overlap
# if (extent(saa)@ymax != map1.ext@ymax && map1.ext@ymax > 0) map2.ext@ymin <- map1.ext@ymax - overlap  # for CAM & ASIA
# 
# 
# #### Blend mosaic of Saatchi-bias.adj and Fused maps
# ## Compute distances and weights in the overlap area
# ## map1 is avg or Fused.map with smaller extent (baccini extent), map2 is saa or saa.bias.adj with larger extent (saatchi extent)
# 
# blend.f <- function(map1, map2, map.fn) {     
#   over.ext <- intersect(map1.ext, map2.ext)                       ## overlap extent (bac.ext is fixed, saa.ext depends on the continent)
#   overl <- crop(map1, over.ext)                                   ## overlap as raster
#   overl[] <- 0                                                    ## overlap values = 0 (NA using the R 'distance' function)
#   overl <- extend(overl, c(100, 0), value=1)                      ## extend to compute distance of overlap area (extend is both towards North and South)
#   overl <- crop(overl, map2.ext)                                  ## remove the extended part (North or South) not necessary
#   writeRaster(overl, filename='./temp/Overlap.tif', overwrite=T)  ## compute distance in GDAL
#   in.fn <- paste(getwd(),'/temp/Overlap.tif', sep="")
#   out.fn <- paste(getwd(),'/temp/Overlap_dist.tif', sep="")
#   dist.py <- paste("gdal_proximity.py ", in.fn, " ", out.fn, " -ot Float32 -distunits GEO", sep="")
#   system(dist.py)
#   wei.over <- raster(out.fn)
#   map1 <- crop(map1, map1.ext) 
#   map2 <- crop(map2, map2.ext)
#   wei.over <- calc(wei.over, fun = function(x) {x / overlap})             ## Compute weights from distances (overlap = maximum distance)
#   wei.over <- crop(wei.over, over.ext, filename='./temp/Weights_mosaic.tif', overwrite=T)  ## remove the extended part (North or South) not necessary
#   map1.over <- crop(map1, over.ext)
#   map2.over <- crop(map2, over.ext)
#   agb.over <- (map1.over * wei.over) + (map2.over * (1 - wei.over))
#   map.blend <- merge(agb.over, map1, map2, filename=paste('Results/', map.fn, Strata.fn, '.tif', sep=''), datatype='FLT4S', overwrite=T)
#   return(map.blend)
# }
# 
# Fused.final <- blend.f(map1 = Fused.map, map2 = saa.bias.adj, map.fn = 'FUSED_FINAL_')
# 
# rm(list=ls(pattern="fus."))
# rm(list=ls(pattern="saa."))
# rm(list=ls(pattern="bac."))
# rm(list=ls(pattern="dist"))
# rm(list=ls(pattern="map"))
# rm(list=c('blend.f', 'overlap'))
# Sys.time()
# 
# #------------------------------------------------
# 
# # Merging 2 output maps
# plot(Fused.map)
# 
# Fused.final <- merge(Fused.map, Gal.bias.adj)
writeRaster(Fused.map, "./Results/Fused_Map/Fused_Final.tif", overwrite = T)

##################################################################################
##############################   MODEL RESULTS   #################################
##################################################################################


########### TRAINING ERROR  (NB: using the Consolidated Reference dataset)

# Input data

# Strata <- raster(paste("./Results/Strata/Strata_", Strata.fn, ".tif", sep=""))
# ref <- raster(paste("Results/Reference/Ref_", Strata.fn, "_Cons.tif", sep=""))
# Fused.final <- raster(paste('Results/FUSED_FINAL_', Strata.fn, '.tif', sep=""))

Strata <- raster(paste("./Results/Strata/Strata_", Strata.fn, ".tif", sep=""))
ref <- raster('./Maps/Ref/ref_ras_EU2.tif')
ref <- crop(ref, Gal)
Fused.final <- raster(paste('Results/Fused_Map/Mosaic/Fused.map_final_EU.tif', sep=""))

fus.par <- read.csv(paste("./Results/Fused_Map/Bias_Weights_", Strata.fn, ".csv", sep=""))

dat.r <- stack(Gal, IIASA, Thur, Fused.final, ref, Strata)
dat <- as.data.frame(as.matrix(dat.r, na.rm = TRUE))
dat <- dat[complete.cases(dat),]
colnames(dat) <- c("Gal", "IIASA","Thur", "Fused", "Ref", "Strata")
dat$Mean <- (dat$Gal + dat$IIASA + dat$Thur) / 3
rm(dat.r)


# Bias, RMSE and SMSE (as Simple Average)

err <- dat - dat$Ref        # compute error for all maps
err$Strata <- dat$Strata    # rewrite the correct strata
bias.m <- apply(err, 2, mean)

rmse <- function(x){sqrt(mean(x^2))}
rmse.m <- apply(err, 2, FUN=rmse)

smse <- vector(length=Strata.s)
err.sd <- sqrt(fus.par$Fus_Var[1:Strata.s])

for (i in 1:Strata.s) {
  err.str <- err$Fused[which(err$Strata==i)]
  smse[i] <- mean((err.str/err.sd[i])^2)
}

err.all <- t(rbind(bias.m[-5:-6], rmse.m[-5:-6]))
n.plot <- nrow(err)
smse <- mean(smse)                     # simple Mean SMSE (not area-weighted)
err.all <- rbind(err.all, n.plot, smse)
colnames(err.all) <- c("Bias", "RMSE")
dir.create("./Results/Validation/", showWarnings = F)
dir.create("./Results/Validation/Training_Error", showWarnings = F)

write.csv(err.all, paste("Results/Validation/Training_Error/Training_Error_", Strata.fn, ".csv", sep=""))


# Plot Saatchi and Baccini errors

png(filename=paste("Results/Validation/Training_Error/Training_Error_Gal_IIASA_Thur_", Strata.fn, ".png", sep=""))
plot(dat$Ref, dat$Gal, xlim=c(0, max(dat)), ylim=c(0,max(dat)), type="p", pch=20, col='green', xlab="Reference data", ylab="Biomass map")
points(dat$Ref, dat$IIASA, pch=20, col='red')
points(dat$Ref, dat$Thur, pch=20, col='blue')
legend(x='topleft', legend=c("Gal", "IIASA", "Thur"), col=c('green', 'red', 'blue'), pch=20)
abline(0,1)
dev.off()


# Plot Fused errors

png(filename=paste("Results/Validation/Training_Error/Training_Error_Fused_", Strata.fn, ".png", sep=""))
plot(dat$Ref, dat$Fused, xlim=c(0, max(dat)), ylim=c(0,max(dat)), type="p", pch=20, col='blue', xlab="Reference data", ylab="Biomass map")
legend(x='topleft', legend="Fused map", col='blue', pch=20)
abline(0,1)
dev.off()

rm(list=c('dat', 'err', 'err.all', 'err.sd', 'err.str', 'Gal', 'bias.m', 'Fused.final', 'n.plot', 'ref', 'rmse', 'rmse.m', 'IIASA','Thur', 'smse'))










#### for EUROPE: compute Mean and total AGB in a way that is compatible with input maps (this step is not essential - to be considered)
########### TOTAL AND MEAN AGB, BY STRATA (ONLY FOR BACCINI EXTENT)


### PROJECT TO SIN (to calculate areas in equal-area prj. the results change only of 2-5% compared to results in WGS84)
library(plotKML)
library(gdalUtils)

## Testing, adjust SINU projection?
Gal <- raster("./Maps/Barredo/barredo_crop.tif")
IIASA <- raster("./Maps/IIASA/1km/bmAg_IIASA2010_crop.tif")
Thur <- raster('./Maps/Thurner/1km/bmAg_Thurner_1km_crop.tif')
Water <- raster('./Covariates/CCI_Water/Water_NL_Crop.tif')
Fused.final <- raster('./Results/Fused_Map/Fused_Final.tif')
Strata <- raster('./Results/Strata/Strata_NL.tif')

Gal_Sin <- raster("./Maps/Barredo/barredo_crop_SIN.tif")
IIASA_Sin <- raster("./Maps/IIASA/1km/bmAg_IIASA2010_crop_SIN.tif")
Thur_Sin <- raster('./Maps/Thurner/1km/bmAg_Thurner_1km_crop_SIN.tif')
Water_Sin <- raster('./Covariates/CCI_Water/Water_NL_Crop_SIN.tif')
Fused.final <- raster('./Results/Fused_Map/Fused_Final_SIN.tif')
Strata <- raster('./Results/Strata/Strata_NL_SIN.tif')

Gal_Sin <- reproject(Gal, CRS = '+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs', program = 'GDAL', method = 'near')
IIASA_Sin <- reproject(IIASA, CRS = '+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs', program = 'GDAL', method = 'near')
Thur_Sin <- reproject(Thur, CRS = '+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs', program = 'GDAL', method = 'near')
Water_Sin <- reproject(Water, CRS = '+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs', program = 'GDAL', method = 'near')
Fused.final  <- reproject(Fused.final , CRS = '+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs', program = 'GDAL', method = 'near')
Strata  <- reproject(Strata , CRS = '+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs', program = 'GDAL', method = 'near')

writeRaster(Gal_Sin, filename = "./Maps/Barredo/barredo_crop_SIN.tif", overwrite = T)
writeRaster(IIASA_Sin, filename = "./Maps/IIASA/1km/bmAg_IIASA2010_crop_SIN.tif", overwrite = T)
writeRaster(Thur_Sin, filename = "./Maps/Thurner/1km/bmAg_Thurner_1km_crop_SIN.tif", overwrite = T)
writeRaster(Water_Sin, filename = "./Covariates/CCI_Water/Water_NL_Crop_SIN.tif", overwrite = T)
writeRaster(Fused.final, filename = "./Results/Fused_Map/Fused_Final_SIN.tif", overwrite = T)
writeRaster(Strata, filename = "./Results/Strata/Strata_NL_SIN.tif", overwrite = T)


Gal.ext <- extent(Gal_Sin)
out.ext <- paste(Gal.ext@xmin, Gal.ext@ymin, Gal.ext@xmax, Gal.ext@ymax, sep=" ")

in.fn <- './Maps/Barredo/barredo_crop.tif'
out.fn <- './Maps/Barredo/barredo_crop_SIN.tif'
proj.py <- paste("gdalwarp -s_srs \'+proj=longlat +datum=WGS84\' -t_srs \'+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs\' -te ", out.ext, " -tr 926.6254331 926.6254331 -r near -srcnodata 255 -dstnodata 255 -ot Byte -overwrite ", in.fn, " ", out.fn, sep = "") 
system(proj.py)

## Project Water to SIN (if file is not available already)
if (!file.exists("./Covariates/CCI_Water/Water_NL_Crop_SIn.tif")) {
  in.fn <- "./Covariates/CCI_Water/Water_NL_Crop.tif"
  out.fn <- "./Covariates/CCI_Water/Water_NL_CropSIN.tif"
  proj.py <- paste("gdalwarp -s_srs \'+proj=longlat +datum=WGS84\' -t_srs \'+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs\' -te ", out.ext, " -tr 926.6254331 926.6254331 -r near -srcnodata 255 -dstnodata 255 -ot Byte -overwrite ", in.fn, " ", out.fn, sep = "") 
  system(proj.py)
} else Water_Sin <- raster(paste('./Covariates/CCI_Water/Water_NL_Crop_SIN.tif', sep=""))



Gal_Sin[Water_Sin == 2] <- NA
IIASA_Sin[Water_Sin == 2] <- NA
Thur_Sin[Water_Sin == 2] <- NA

## Project Strata (Byte) and Fused map (Float32) using GDAL

# in.fn <- paste(getwd(),'/Results/Strata/Strata_', Strata.fn, '.tif', sep="")
# out.fn <- paste(getwd(),'/Results/Strata/Strata_', Strata.fn, '_SIN.tif', sep="")
# proj.py <- paste("gdalwarp -s_srs \'+proj=longlat +datum=WGS84\' -t_srs \'+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs\' -te ", out.ext, " -tr 926.6254331 926.6254331 -r near -srcnodata 255 -dstnodata 255 -ot Byte -overwrite ", in.fn, " ", out.fn, sep = "") 
# system(proj.py)
# Strata <- raster(out.fn)
# 
# in.fn <- paste(getwd(),'/Results/FUSED_FINAL_', Strata.fn, '.tif', sep="")
# out.fn <- paste(getwd(),'/Results/FUSED_FINAL_', Strata.fn, '_SIN.tif', sep="")
# proj.py <- paste("gdalwarp -s_srs \'+proj=longlat +datum=WGS84\' -t_srs \'+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs\' -te ", out.ext, " -tr 926.6254331 926.6254331 -r near -srcnodata -3.39999995214e+038 -dstnodata -3.39999995214e+038 -ot Float32 -overwrite ", in.fn, " ", out.fn, sep = "") 
# system(proj.py)
# Fused.final <- raster(out.fn)



pix.area <- res(Gal_Sin)[1]
all.map <- stack(Gal_Sin, IIASA_Sin, Thur_Sin, Fused.final, Strata)

all.dat <- as.data.frame(as.matrix(all.map, na.rm = TRUE))
all.dat <- all.dat[complete.cases(all.dat),]                      ## here only data within the Baccini Extent is used
colnames(all.dat) <- c("Gal", "IIASA", "Thur", "Fused", "Strata")
all.dat$Strata <- as.factor(all.dat$Strata)

rm(list=c('all.map', 'bac', 'in.fn', 'out.fn', 'out.ext', 'proj.py', 'saa.ext', 'water'))


# MEAN AGB
Strata.s <- 3 
Strata.names <- 1:Strata.s

agb.mean.s <- aggregate(all.dat$Gal, by = list(all.dat$Strata), FUN="mean")
agb.mean.b <- aggregate(all.dat$IIASA, by = list(all.dat$Strata), FUN="mean")
agb.mean.t <- aggregate(all.dat$Thur, by = list(all.dat$Strata), FUN="mean")
agb.mean.f <- aggregate(all.dat$Fused, by = list(all.dat$Strata), FUN="mean")
agb.mean <- cbind(agb.mean.s[,2],agb.mean.b[,2], agb.mean.t[,2], agb.mean.f[,2])
agb.mean <- t(agb.mean)
mode(agb.mean) <- "numeric"
rownames(agb.mean) <- c("Gal", "IIASA", "Thur","Fused")
colnames(agb.mean) <- Strata.names
rm(list=ls(pattern="agb.mean."))


# TOTAL AGB  

agb.tot.s <- aggregate(all.dat$Gal, by = list(all.dat$Strata), FUN="sum")
agb.tot.b <- aggregate(all.dat$IIASA, by = list(all.dat$Strata), FUN="sum")
agb.tot.t <- aggregate(all.dat$Thur, by = list(all.dat$Strata), FUN="sum")
agb.tot.f <- aggregate(all.dat$Fused, by = list(all.dat$Strata), FUN="sum")
agb.tot <- cbind(agb.tot.s[,2],agb.tot.b[,2], agb.tot.t[,2], agb.tot.f[,2])
agb.tot <- t(agb.tot/10^9*(((pix.area^2)/10000))) # What is happening here? This is a hectare conversion?
mode(agb.tot) <- "numeric"
rownames(agb.tot) <- c("Gal", "IIASA", "Thur","Fused")
colnames(agb.tot) <- Strata.names
rowSums(agb.tot)  
rm(list=ls(pattern="agb.tot."))
rm(all.dat)


# Summary AGB stat

agb.all <- rbind(agb.mean, agb.tot)
rownames(agb.all) <- c("Gal_Mean", "IIASA_Mean", "Thur_Mean","Fus_Mean", "Gal_Tot", "IIASA_Tot", "Thur_Tot", "Fus_Tot")
colnames(agb.all) <- c(1:Strata.s)
colnames(agb.mean) <- c(1:Strata.s)
colnames(agb.tot) <- c(1:Strata.s)


# Overview Plot

png(filename=paste("Results/Fused_Map/Fusion_Summary_", Strata.fn, ".png", sep=""))
par(mfrow = c(2,2))
barplot(t(weight[1:Strata.s,]), beside = F, col=c("red", "green" ,"blue"), xlab="Strata", ylab="%", legend.text=c("Gal", "IIASA", "Thur"), main="Weigths")
barplot(t(bias[1:Strata.s,]), beside = T, col=c("red", "green" ,"blue"), xlab="Strata", ylab="Mg/ha", legend.text=c("Gal", "IIASA", "Thur"), main="Bias")
barplot(agb.mean[, 1:Strata.s], beside=T, col=c("red", "green" ,"blue", "darkgreen"), space=c(0.1, 2), xlab="Strata", ylab="Mg/ha", main="Mean AGB", legend.text=c("Gal", "IIASA", "Thur", "Fused"))
barplot(agb.tot[, 1:Strata.s], beside=T, col=c("red", "green" ,"blue","darkgreen"), space=c(0.1, 2), xlab="Strata", ylab="Pg", main="Total AGB", legend.text=c("Gal", "IIASA", "Thur", "Fused"))
dev.off()


### TOTAL AGB OF FUSED MAP FOR COMPLETE CONTINENT AREA (SAATCHI EXTENT)

maps <- stack(Gal_Sin, Fused.final)

agb.ha <- as.matrix(getValues(maps))
agb.ha <- agb.ha[complete.cases(agb.ha),]
agb.T <- agb.ha * ((pix.area^2)/10000)
Gal.cont <- (sum(agb.T[,1])) / 10^9 # Why the 10^9?
fus.cont <- (sum(agb.T[,2])) / 10^9
agb.all <- rbind(agb.all, Gal.cont, fus.cont)
#agb.all[7:8, 2:3] <- NA #


write.csv(agb.all, paste("Results/AGB_Stat_", Strata.fn, ".csv", sep=""))
rm(list=ls(pattern="agb"))
rm(Fused.final)
Sys.time()













#######################################################################################
########################## VALIDATION OF FUSED MAP  ###################################
#######################################################################################

## Method: Split consolidated Ref. data in Cal / Val sets (70% / 30%)
## make a fused map using Cal data and compute RMSE on Val data
## NB: directory are different from normal Fusion, outputs are in /Results/Validation/, and outputs are reduced

################## SPLIT REFERENCE DATA IN CAL / VAL DATASETS #########################
Sys.time()

## Reference dataset (Original)

ref <- raster('./Maps/Ref/ref_ras_EU2.tif') 
ref.fn <- names(ref)


## Select CAL / VAL datasets (as raster)

ref.p <- rasterToPoints(ref, spatial=T)
set.seed(10)
cal.p <- ref.p[sample(nrow(ref.p), size = nrow(ref.p)*0.70), ]
dir.create("./Results/Validation/Reference", showWarnings = F)
cal <- rasterize(cal.p, ref, field=ref.fn, filename=paste("./Results/Validation/Reference/Ref_", Strata.fn, "_Cal.tif", sep=""), overwrite=TRUE)
val <- mask(ref, cal, inverse=T, filename=paste("./Results/Validation/Reference/Ref_", Strata.fn, "_Val.tif", sep=""), overwrite=TRUE)
rm(list=ls(pattern='.p'))
rm(val)


##################################################################################
###############################   CALIBRATION   ##################################

ref <- cal
rm(cal)

######  INPUT DATA  ##########


#Gal <- raster(paste('./Input_Maps/Gal_1km_', cont, '.tif', sep=""))
Gal <- raster("./Maps/Gallaun/1km/bmAg_JR2000_ll_1km_eur.tif")
Gal <- raster("./Maps/Barredo/barredo_Alligned.tif")
Gal[Gal < 0] <- NA

# cant get Bar en Gal extents to match, fix required
#Bar <- projectRaster(Bar, Gal, method = 'ngb')
#writeRaster(Bar, filename = "./Maps/Barredo/barredo_reproj.tif")
#

#-------
IIASA <- raster("./Maps/IIASA/1km/bmAg_IIASA2010_Alligned.tif")
#IIASA <- crop(IIASA, Gal)

Thur <- raster('./Maps/Thurner/1km/bmAg_Thurner_1km_Alligned.tif')  
#Thur <- crop(Thur, Gal)
#ref <- raster(paste('./Reference/Ref_', cont, '.tif', sep=""))
ref <- raster('./Maps/Ref/ref_ras_EU2.tif')
ref <- crop(ref, Gal)

#vcf <- raster("./Covariates/MODIS_VCF_2005/transformed/Mosaic/MODIS_VCF_Mosaic.tif")
#align_rasters("./Covariates/MODIS_VCF_2005/transformed/Mosaic/MODIS_VCF_Mosaic.tif", "./Maps/Gallaun/1km/bmAg_JR2000_ll_1km_eur.tif", dstfile = "./Covariates/MODIS_VCF_2005/transformed/Mosaic/MODIS_VCF_Mosaic_C.tif")
#vcf[vcf < 0] <- NA
#vcf[vcf > 100] <- NA
#writeRaster(vcf, filename = "./Covariates/MODIS_VCF_2005/transformed/Mosaic/MODIS_VCF_Mosaic_EU.tif")
vcf <- raster("./Covariates/Outputs/MODIS_VCF/Mosaic/Mosaic_aggre_align.tif")


#hei <- raster(paste('./Strata/HEI_', cont, '.tif', sep=""))
hei <- raster("./Covariates/Height/Height_align.tif")
#hei <- crop(hei, Gal)

#cci <- raster(paste('./Strata/CCI_', cont, '.tif', sep=""))
#cci <- raster('./Covariates/CCI_2005/CCI_2005_resample.tif')
#cci[cci < 0] <- NA
#writeRaster(cci, filename = './Covariates/CCI_2005/CCI_2005_resample_EU.tif')
cci <- raster("./Covariates/Outputs/CCI_2005/CCI_Mul_align.tif")



#water <- raster(paste('./Strata/WATER/Water_1km_', cont, '.tif', sep=""))
water <- raster("./Covariates/Outputs/CCI_Water/CCI_Water_aggr_Resam_R.tif")

###### WATER MASK APPLIED TO GalTCHI & ThurCINI ############


###### WATER MASK APPLIED TO GalTCHI & ThurCINI ############




###### WATER MASK APPLIED TO SAATCHI ############

Gal[water == 2] <- NA
Thur[water == 2] <- NA
IIASA[water == 2] <- NA
rm(water)



##################################################################################
###########################   ERROR STRATA   #####################################
##################################################################################

Gal.er <- ref - Gal
IIASA.er <- ref - IIASA
Thur.er <- ref - Thur

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

dir.create("./Results/Validation/Strata/", showWarnings = F)
png(filename=paste("./Results/Validation/Strata/RF_Predict_Errors_VarImpPlot_", Strata.fn, ".png", sep=""), width=960, height=480, res=96)
opar <- par(mfrow=c(1,3))
varImpPlot(Gal.rf)
varImpPlot(Thur.rf)
varImpPlot(IIASA.rf)
dev.off()


png(filename=paste("./Results/Validation/Strata/RF_Predict_Errors_CrossVal_", Strata.fn, ".png", sep=""), width=960, height=480, res=96)
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
Gal.er <- predict(error.dat, Gal.rf, type= 'response', progress='text')
Thur.er <- predict(error.dat, Thur.rf, type= 'response', progress='text')
IIASA.er <- predict(error.dat, IIASA.rf, type= 'response', progress='text')

#rm(error.map)

err.all <- stack(Gal.er, Thur.er, IIASA.er)
mydata <- as.data.frame(getValues(err.all))
mydata.k <- mydata[complete.cases(mydata),]  # Remove NA, K-means cannot handle NAs
#rm(list=ls(pattern='er'))
#(list=ls(pattern='.rf'))

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
Sys.time()



########  PREDICT STRATA FOR NO DATA AREAS  #########

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
write.csv(rf.out, paste("./Results/Validation/Strata/RF_Predict_Errors_", Strata.fn, ".csv", sep=""))
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
Strata.fn <- 'EU'
Strata <- calc(strata, fun=strata.f, filename=paste("./Results/Validation/Strata/Strata_", Strata.fn, ".tif", sep=""), datatype='INT1U', overwrite=T, progress='text')
names(Strata) <- paste("Strata_", Strata.fn, sep="")
rm(list=ls(pattern='str.'))

Sys.time()



##################################################################################
###########################   REFERENCE DATA  ####################################
##################################################################################

########  CONSOLIDATE REFERENCE DATA #####  START HERE FOR STRATA AS GLC2000 OR VCF10
# Strata <- raster(paste("./Results/Strata/Strata_", Strata.fn, ".tif", sep=""))

Strata.s <- 3    
Strata.names <- 1:Strata.s

## Each NFI has a unique code, and there is a raster map where each plot of the NFI has the code value (e.g.: all plots in Spain have value 1, in France is 2, etc.)
#codes <- raster(paste('./Reference/Ref_code_', cont, '.tif', sep=""))
#code.names <- as.character(read.csv(paste("./Reference/Codes_CAM.csv", sep=""))[,2])
#code.n <- maxValue(codes)

code.names <- 'EU'

## Each NFI has a unique code, and there is a raster map where each plot of the NFI has the code value (e.g.: all plots in Spain have value 1, in France is 2, etc.)
codes <- raster(paste('./Reference/Ref_code_EU.tif', sep=""))
codes <- crop(codes, Gal)
code.names <- as.character(read.csv(paste("./Reference/Codes_EU.csv", sep=""))[,2])
code.names <- code.names[unique(codes)]
code.n <- maxValue(codes)

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

# adjust str.code in case ref data are missing for some code (datasets) due to random selection

if (nrow(str.code) < code.n) {
  str.code.new <- matrix(0, nrow=code.n, ncol=ncol(str.code)) # use ncol in case it is < Strata.s, which is adjusted in the next loop
  str.code <- as.matrix(str.code)
  for (i in 1:nrow(str.code)) {
    n.row <- as.numeric(rownames(str.code)[i])
    str.code.new[n.row,] <- str.code[i,]  }
  colnames(str.code.new) <- colnames(str.code)
  str.code <- str.code.new
  rm(list=c('str.code.new', 'n.row'))
}

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

Strata.fn <- 'EU' 
png(filename=paste("./Results/Validation/Reference/Reference_Data_", Strata.fn, "_Orig.png", sep=""))
plt <- barplot(str.code, beside=F, main=paste("Reference Data by ", Strata.fn, " strata - original", sep=""), xlab="Strata", legend = code.names, names.arg = Strata.names,
               cex.names = 0.7, col=c("red", "orange", "yellow", "green", "blue", "violet", "pink" , "cyan", "gray", "forestgreen", "orange3"), args.legend = list(x = "topright", cex=0.7))
text(plt, colSums(str.code), labels = colSums(str.code), pos = 3)
dev.off()


#### for EUROPE: consolidation is not necessary because there are not reference biomass maps, as indicated above

################ CONSOLIDATE: Extract Reference data from maps in under-represented strata

## Reference Maps: read the Maps (wall-to-wall) and the raster with only Reference pixel of the maps (Pixel with plots) in the Folders

ref.maps <- list.files('./Reference/Maps', pattern='.tif$')
N.map <- length(ref.maps)
if (cont == "SAM") N.map <- 0  # not needed since there is no consolidation needed with the current reference dataset

if (N.map > 0){      ### if there are no reference maps, the whole consolidation section will be skipped!
  
  map.list <- vector('list', length=N.map)
  for (i in 1:N.map) {  map.list[[i]] <- raster(paste('./Reference/Maps/', ref.maps[i], sep="")) }
  
  ref.plot <- list.files('./Reference/Maps/Map_Ref', pattern='.tif$')
  plot.list <- vector('list', length=length(ref.plot))
  for (i in 1:length(ref.plot)) { plot.list[[i]] <- raster(paste('./Reference/Maps/Map_Ref/', ref.plot[i], sep="")) }
  
  
  ## Mask out pixel in the Maps already sampled (by the plots)
  
  for (i in 1:N.map) { map.list[[i]] <- mask(map.list[[i]], plot.list[[i]], inverse=TRUE) }
  
  
  ## Harmonize Strata with maps
  
  strata.list <- vector('list', length=N.map)
  
  for (i in 1:N.map) {
    strata.list[[i]] <- crop(Strata, map.list[[i]])
    strata.list[[i]] <- mask(strata.list[[i]], map.list[[i]])
  }
  
  
  ## Add code of each map
  
  code.list <- vector('list', length=N.map)
  
  for (i in 1:N.map) {
    code.list[[i]] <- crop(codes, strata.list[[i]])
    code.list[[i]] <- mask(code.list[[i]], strata.list[[i]])
    strata.list[[i]] <- stack(strata.list[[i]], code.list[[i]])
  }
  
  
  ## Frequency table: Nr. of pixels of Reference Maps available for each strata
  
  str.maps <- matrix(0, code.n, Strata.s)
  rownames(str.maps) <- code.names     # Code
  colnames(str.maps) <- c(1:Strata.s)   # Strata
  
  for (i in 1:N.map){
    map <- raster(strata.list[[i]], layer=1)
    code <- raster(strata.list[[i]], layer=2)
    if (cont == "ASIA") code <- 7 else  code <- maxValue(code)  ## temp solution for ASIA
    for (n in 1:ncol(str.maps)) {
      str.maps[code,n] <- freq(map, useNA='no', value = n)
    }
  }
  
  
  ##### Extract random samples for under-represented strata
  
  N.min <- (sum(str.code) / ncol(str.code)) / 2    # the min number of ref data per strata is defined as half of average ref per strata
  N.perc <- (75/(N.map)) * 0.01
  
  sizes <- str.code
  sizes[,] <- 0
  for (n in 1:ncol(str.code)) {
    for (j in 1:nrow(str.code)) {
      if (str.maps[j,n] > 0)
        if (sum(str.code[,n]) < N.min && (str.code.x[j,n] < 75 || str.code.x[j,n] == 100))  # this condition may be simplified, removing the part after &&
          if (str.maps[j,n] > round((N.min*2 - sum(str.code[,n]))/N.map)) sizes[j,n] = round((N.min*2 - sum(str.code[,n]))/N.map) else sizes[j,n] = round(str.maps[j,n]/2)
          else if (any(str.code.x[,n] > 75) && str.code.x[j,n] < 75)
            if (str.maps[j,n] > sum(str.code[,n])*N.perc ) sizes[j,n] = round(sum(str.code[,n])*N.perc) else sizes[j,n] = round(str.maps[j,n]/2)
    }
  }
  # sizes
  
  
  ## Extract samples from each reference maps
  
  sample.list <- vector('list', length=N.map)
  
  for (i in 1:N.map){
    sample.list[[i]] <- map.list[[i]]
    sample.list[[i]] [] <- NA
    str.map <- raster(strata.list[[i]], layer=1)
    code <- raster(strata.list[[i]], layer=2)
    if (cont == "ASIA") code <- 7 else  code <- maxValue(code)  ## temp solution for ASIA
    for (n in 1:Strata.s) {
      if (sizes[code,n] > 0){
        sampl.str <- map.list[[i]]
        sampl.str[str.map != n | str.map == NA] <- NA     # mask out the pixels of other strata and the NA values in str.map
        set.seed(45)
        sampl.str <- sampleRandom(sampl.str, size = sizes[code,n], asRaster=TRUE)
        sample.list[[i]] <- merge(sample.list[[i]], sampl.str)
      }
    }
  }
  
  rm(str.map)
  rm(sampl.str)
  
  
  ## Merge in 1 Reference Dataset
  
  for (i in 1:N.map){ ref <- merge(sample.list[[i]], ref) }
  writeRaster(ref, filename=paste("./Results/Validation/Reference/Ref_", Strata.fn, "_Cons.tif", sep=""), overwrite=T)
  
  
  # Make consolidated Barplot to see new distribution
  
  ref.str.c <- stack(ref, Strata, codes)
  ref.dat.c <- as.data.frame(getValues(ref.str.c))
  ref.dat.c <- ref.dat.c[complete.cases(ref.dat.c),]
  colnames(ref.dat.c) <- c("AGB", "Strata", "Code")
  str.code.c <- table(ref.dat.c$Code, ref.dat.c$Strata, dnn = c("Code", "Strata"))
  
  png(filename=paste("./Results/Validation/Reference/Reference_Data_", Strata.fn, "_Cons.png", sep=""))
  plt.c <- barplot(str.code.c, beside=F, main=paste("Reference Data by ", Strata.fn, " - consolidated", sep=""), xlab="Strata", legend = code.names, names.arg = Strata.names, 
                   cex.names = 0.6, col=c("red", "orange", "yellow", "green", "blue", "violet", "pink" , "cyan", "gray", "forestgreen"), args.legend = list(x = "topright", cex=0.7))
  text(plt.c, colSums(str.code.c), labels = colSums(str.code.c), pos = 3)
  dev.off()
  
  ref.code <- t(rbind(table(ref.dat$Code), table(ref.dat.c$Code)))
  colnames(ref.code) <- c("Ref Original", "Ref Consolid")
  rownames(ref.code) <- code.names
  write.csv(ref.code, paste("./Results/Validation/Reference/Ref_codes_", Strata.fn, ".csv", sep=""))
  
  rm(list=c("sizes", "n", "map", "i", "j"))
} else {
  writeRaster(ref, filename=paste("./Results/Validation/Reference/Ref_", Strata.fn, "_Cons.tif", sep=""), overwrite=T)}
### End of if statement, to skip Consolidation -for Continents without Ref maps

rm(list=ls(pattern="ref."))
rm(list=ls(pattern="str."))
rm(list=ls(pattern="code"))
rm(list=ls(pattern=".list"))
rm(list=ls(pattern="N."))
rm(list=ls(pattern="plt"))


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

png(filename="./Results/Validation/Reference/Errors of IIASA_Thur.png")
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

# !! Wont Run because Gal dataset has a mean error of 0, due to the reference dataset being used in the map creation  !! 
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

### calculate weight matrix variance of the fused map

v <- vector('list', length = Strata.s)
for (i in 1:Strata.s){  v[[i]] <- solve(t(X) %*% solve(cov[[i]]) %*% X) }

v.err <- matrix(1:Strata.s)
for (i in 1:Strata.s){ v.err[i] <- v[[i]] }

colnames(v.err) <- c("Fus_Var")
rownames(v.err) <- c(1:Strata.s)
rm(list=c('cov', 'i', 'v', 'w', 'X'))

## Finalize Bias and Weights

weight[weight < 0] <- 0   # set weights to 0 - 1 limits
weight[weight > 1] <- 1   # set weights to 0 - 1 limits

### Compute Error Variance of Input maps

Gal.var <- aggregate(error$Gal_er, by = list(error$Strata), FUN="var")
Thur.var <- aggregate(error$Thur_er, by = list(error$Strata), FUN="var")
IIASA.var <- aggregate(error$IIASA_er, by = list(error$Strata), FUN="var")

### Compute n. pixel per strata

n.pix <- matrix(1:Strata.s)
for (i in 1:Strata.s) { n.pix[i] <- freq(Strata, useNA='no', value = i) }
n.pix.w <- n.pix / (sum(n.pix))


### Compile Error Variances

err.var <- cbind(Gal.var, Thur.var[,2], IIASA.var[,2], v.err, n.pix, n.pix.w)
colnames(err.var) <- c("Strata", "Gal.var", "Thur.var","IIASA.var","Fus_Var", "N_Pix", "N_Pix_w")
err.var$Gal_Var_w <- err.var$Gal.var * err.var$N_Pix_w
err.var$Thur_Var_w <- err.var$Thur.var * err.var$N_Pix_w
err.var$IIASA_Var_w <- err.var$IIASA.var * err.var$N_Pix_w
err.var$Fus_Var_w <- err.var$Fus_Var * err.var$N_Pix_w


# Compile Bias, Weights, Error variance, and add Strata 9

fus.par <- aggregate(error$Strata, by = list(error$Strata), FUN="length")
colnames(fus.par) <- c("Strata", "N")
fus.par$Strata <- as.numeric(as.character(fus.par$Strata))
fus.par <- cbind(fus.par, bias, weight, err.var[,-1])   # adjust the number of classes
#fus.par[9,] <- c(9, 0, 0, 0, 0.5, 0.5, rep(0,8))
fus.par[Strata.s+1,] <- c(Strata.s+1, 0, 0, 0, 0, 0.5, 0.5, 0.5, rep(0,10)) 
# Should weights be set according to the number of input maps? e.g. 3 maps == weight 0.33? check this
bias <- fus.par[,3:5]                             # Add Strata 9 to bias and weight (for Fusion)
weight <- fus.par[,6:8]
dir.create("./Results/Validation/Fused_Map/", showWarnings = F)
write.csv(fus.par, paste("./Results/Validation/Fused_Map/Bias_Weights_", Strata.fn, ".csv", sep=""), row.names = FALSE)


## Map Uncertainty (Standard deviation of Error of Fused map)

rm(list=c('fus.par', 'i', 'n.pix', 'n.pix.w', 'uncer'))
rm(list=ls(pattern="err"))


####### MAP FUSION  

## Set NA in Strata to a value (Strata must have always values in Fusion, no NA)

Strata.n <- Strata.s + 1
Strata[is.na(Strata)] <- Strata.n  

#### FOR EUROPE: this function needs to be adapted to 3 maps instead of 2. After Line 593 there should be a "c <-..." and the map "c" should be added to line 594 with respective bias and weight 
###  FUSION                          
# Double check, if 3 maps are used and 1 of them contains NA, the outcome will be NA. Adjusted to ignore NA?
maps <- stack(Gal, Thur, IIASA, Strata)

# the 'Strata.n' represents the number of stacked layers
biomass.fusion <- function(x) {
  m <- matrix(x, nrow= 1, ncol=4)
  n <- m[,4]
  g <- m[1:(Strata.n-1)] + as.matrix(bias[n,])
  g[g < 0] <- 0
  w <- weight[n,1:(Strata.n-1)]
  w[is.na(g)]<- NA
  p <- sum(w, na.rm = T) # calculate sum of weight values
  pp <- w/p # divide weight values by sum to get the proportion to == 1
  pp <- as.numeric(pp)
  result <- as.integer(round(sum(pp*g, na.rm = T)))
  return(result)
}

Fused.map <- calc(maps, fun = biomass.fusion, progress = 'text')
Fused.map[Fused.map < 0] <- 0


######## EXTEND TO SAATCHI AREA ################
#### NOT NEEDED FOR CALIBRATION MAP (ALL VALIDATION DATA ARE ON BACCINI EXTENT - CALIBRATION MAP)

Fused.final <- Fused.map
writeRaster(Fused.final, filename=paste('Results/Validation/Fused_Map/FUSED_FINAL_', Strata.fn, '.tif', sep=''), datatype='FLT4S', overwrite=T)
rm(list=c("Fused.map", "maps"))





##################################################################################
###############################   VALIDATION   ###################################
##################################################################################

Fused.final <- raster(paste('./Results/Validation/Fused_Map/FUSED_FINAL_', Strata.fn, '.tif', sep=''))

#Gal <- raster(paste('./Input_Maps/Gal_1km_', cont, '.tif', sep=""))
Gal <- raster("./Maps/Gallaun/1km/bmAg_JR2000_ll_1km_eur.tif")
Gal <- raster("./Maps/Barredo/barredo_Alligned.tif")
Gal[Gal < 0] <- NA
IIASA <- raster("./Maps/IIASA/1km/bmAg_IIASA2010_Alligned.tif")
Thur <- raster('./Maps/Thurner/1km/bmAg_Thurner_1km_Alligned.tif')   
Strata <- raster(paste("./Results/Validation/Strata/Strata_", Strata.fn, ".tif", sep=""))

## Validation dataset

ref <- raster(paste("./Results/Validation/Reference/Ref_", Strata.fn, "_Val.tif", sep=""))    

dat.r <- stack(Gal, IIASA, Thur, Fused.final, ref, Strata)
dat <- as.data.frame(as.matrix(dat.r, na.rm = TRUE))
dat <- dat[complete.cases(dat),]
colnames(dat) <- c("Gal", "IIASA","Thur", "Fused", "Ref", "Strata")
dat$Mean <- (dat$Gal + dat$IIASA + dat$Thur) / 3
rm(dat.r)


# Bias, RMSE and SMSE (as Simple Average)

err <- dat - dat$Ref        # compute error for all maps
err$Strata <- dat$Strata    # rewrite the correct strata
bias.m <- apply(err, 2, mean)

rmse <- function(x){sqrt(mean(x^2))}
rmse.m <- apply(err, 2, FUN=rmse)

fus.par <- read.csv(paste("./Results/Validation/Fused_Map/Bias_Weights_", Strata.fn, ".csv", sep=""))
smse <- vector(length=Strata.s)
err.sd <- sqrt(fus.par$Fus_Var[1:Strata.s])

for (i in 1:Strata.s) {
  err.str <- err$Fused[which(err$Strata==i)]
  smse[i] <- mean((err.str/err.sd[i])^2)
}

err.all <- t(rbind(bias.m[-5:-6], rmse.m[-5:-6]))
n.plot <- nrow(err)
smse <- mean(smse)                     # simple Mean SMSE (not area-weighted)
err.all <- rbind(err.all, n.plot, smse)
colnames(err.all) <- c("Bias", "RMSE")
write.csv(err.all, paste("./Results/Validation/Validation_Error_", Strata.fn, ".csv", sep=""))


# Plot Saatchi and Baccini errors

png(filename=paste("./Results/Validation/Validation_Error_Gal_IIASA_Thur_", Strata.fn, ".png", sep=""))
par(mfrow = c(1,1))
plot(dat$Ref, dat$Gal, xlim=c(0, max(dat)), ylim=c(0,max(dat)), type="p", tck= 0.02, pch=1, col='green', xlab="AGB Reference data (Mg/ha)", ylab="AGB Biomass map (Mg/ha)", main= paste0("Validation", Strata.fn))
points(dat$Ref, dat$IIASA, pch=1, col='red')
points(dat$Ref, dat$Thur, pch=1, col='blue')
legend(x='topleft', legend=c("Gal", "IIASA", "Thur"), col=c('green', 'red', 'blue'), pch=1, bty='n')
abline(0,1,lwd=0.02)
dev.off()

pdf(file=paste("./Results/Validation/Validation_Error_Gal_IIASA_Thur_", Strata.fn, ".pdf", sep=""))
par(mfrow = c(1,1))
plot(dat$Ref, dat$Gal, xlim=c(0, max(dat)), ylim=c(0,max(dat)), type="p", tck= 0.02, pch=1, col='green', xlab="AGB Reference data (Mg/ha)", ylab="AGB Biomass map (Mg/ha)", main= paste0("Validation", Strata.fn))
points(dat$Ref, dat$IIASA, pch=1, col='red')
points(dat$Ref, dat$Thur, pch=1, col='blue')
legend(x='topleft', legend=c("Gal", "IIASA", "Thur"), col=c('green', 'red', 'blue'), pch=1, bty='n')
abline(0,1,lwd=0.02)
dev.off()


# Plot Fused errors

png(filename=paste("./Results/Validation/Validation_Error_Fused_", Strata.fn, ".png", sep=""))
plot(dat$Ref, dat$Fused, xlim=c(0, max(dat)), ylim=c(0,max(dat)), type="p", pch=20, col='blue', xlab="Reference data", ylab="Biomass map", main="Validation")
legend(x='topleft', legend="Fused map", col='blue', pch=20)
abline(0,1)
dev.off()

rm(list=c('dat', 'err', 'err.all', 'err.sd', 'err.str', 'bac', 'bias', 'bias.m', 'fus.par', 'Fused.final', 'n.plot', 'pix.area', 'ref', 'rmse', 'rmse.m', 'smse', 'saa', 'weight'))

# }

Sys.time()
