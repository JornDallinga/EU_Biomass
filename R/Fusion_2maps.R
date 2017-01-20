# load libraries

if (!require(randomForest)) install.packages('randomForest')

######  INPUT DATA  ##########

#saa <- raster(paste('./Input_Maps/saa_1km_', cont, '.tif', sep=""))
saa <- raster("./Maps/Gallaun/1km/bmAGB_Gallaun_crop.tif")
IIASA <- raster("./Maps/IIASA/1km/bmAg_IIASA2010_crop.tif")
bac <- raster('./Maps/Thurner/1km/bmAg_Thurner_1km_crop.tif')  

#ref <- raster(paste('./Reference/Ref_', cont, '.tif', sep=""))
ref <- raster('./Maps/Ref/ref_ras.tif')

vcf <- raster("./Covariates/Landsat_VCF_2005/VCF_crop.tif")
#vcf <- raster(paste('./Strata/VCF_', cont, '.tif', sep=""))

#hei <- raster(paste('./Strata/HEI_', cont, '.tif', sep=""))
hei <- raster("./Covariates/Height/Height_reclass_crop.tif")

#cci <- raster(paste('./Strata/CCI_', cont, '.tif', sep=""))
cci <- raster('./Covariates/CCI_2005/CCI_crop.tif')

#water <- raster(paste('./Strata/WATER/Water_1km_', cont, '.tif', sep=""))
water <- raster("./Covariates/CCI_Water/Water_NL_Crop.tif") # AGGREGATE FIRST

###### WATER MASK APPLIED TO SAATCHI & BACCINI ############

saa[water == 2] <- NA
bac[water == 2] <- NA
IIASA[water == 2] <- NA
rm(water)


##################################################################################
###########################   ERROR STRATA   #####################################
##################################################################################

saa.er <- ref - saa
bac.er <- ref - bac
IIASA.er <- ref - IIASA

#saa.er <- saa
#bac.er <- bac

error.map <- stack(saa.er, bac.er, vcf, saa, bac, cci)
names(error.map) <- c("saa.er", "bac.er", "vcf", "saa", "bac", "cci")
#rm(cci)

error <- as.data.frame(getValues(error.map))
error <- error[complete.cases(error),]
colnames(error) <- c("saa.er", "bac.er", "vcf", "saa", "bac", "cci")
error$cci <- as.factor(error$cci)


### MODEL MAP ERROR with Random Forest  
# NOTE: mtry=2 to have comparable R2 without GLC2000, otherwise with 5 variables mtry=1 by default and the R2 drops)

saa.rf.f <- formula(saa.er ~ cci + vcf + bac + saa)
bac.rf.f <- formula(bac.er ~ cci + vcf + bac + saa)
#bac.rf.f <- formula(bac.er ~ hei + vcf + saa + bac + cci)
set.seed(55)
saa.rf <- randomForest(saa.rf.f, error, importance=F, mtry=2)
set.seed(55)
bac.rf <- randomForest(bac.rf.f, error, importance=F, mtry=2)
# saa.rf
# bac.rf

rf.out <- matrix(0, nrow=2, ncol=2, dimnames=list(c("Predict_Error_Saatchi", "Strata_Error_rate"), c("RMSE", "R2")))
rf.out[1,1] <- sqrt(bac.rf$mse[500])
rf.out[1,2] <- mean(bac.rf$rsq[500])
rf.out[2,1] <- sqrt(saa.rf$mse[500])
rf.out[2,2] <- mean(saa.rf$rsq[500])

png(filename=paste("./Results/Strata/RF_Predict_Errors_VarImpPlot.png", sep=""), width=960, height=480, res=96)
opar <- par(mfrow=c(1,2))
varImpPlot(saa.rf)
varImpPlot(bac.rf)
dev.off()

varImp <- cbind(importance(saa.rf), importance(bac.rf))
colnames(varImp) <- c('Saatchi_Imp', 'Baccini_Imp')
write.csv(varImp, paste("./Results/Strata/IncNodePur.csv", sep=""))
rm(varImp)

png(filename=paste("./Results/Strata/RF_Predict_Errors_CrossVal.png", sep=""), width=960, height=480, res=96)
opar <- par(mfrow=c(1,2))
plot(error$saa.er,predict(saa.rf, newdata=error), main="Cross Validation Error Saatchi", xlab="Error", ylab="Predicted error", xlim=c(min(error$saa.er), max(error$saa.er)), ylim=c(min(error$saa.er),max(error$saa.er)))
abline(0,1)
plot(error$bac.er,predict(bac.rf, newdata=error), main="Cross Validation Error Baccini", xlab="Error", ylab="Predicted error", xlim=c(min(error$bac.er), max(error$bac.er)), ylim=c(min(error$bac.er),max(error$bac.er)))
abline(0,1)
dev.off()

### PREDICT MAP ERROR

error.dat <- dropLayer(error.map, 1:2)
saa.er <- predict(error.dat, saa.rf, type= 'response', filename="./Results/Strata/Saatchi_Error_RF.tif", datatype='INT2S', overwrite=T, progress='text')
bac.er <- predict(error.dat, bac.rf, type= 'response', filename="./Results/Strata/Baccini_Error_RF.tif", datatype='INT2S', overwrite=T, progress='text')
rm(error.map)

err.all <- stack(saa.er, bac.er)
mydata <- as.data.frame(getValues(err.all))
mydata.k <- mydata[complete.cases(mydata),]  # Remove NA, K-means cannot handle NAs
rm(list=ls(pattern='er'))
rm(list=ls(pattern='.rf'))

# K cluser
set.seed(55)
fit.k <- kmeans(mydata.k, 8)      # 8 cluster solution

# Convert to raster
mydata$cluster <- NA
mydata[rownames(mydata.k), "cluster"] <- fit.k$cluster
Strata <- saa
values(Strata) <- mydata$cluster
dataType(Strata) <- 'INT1U'

rm(list=ls(pattern='mydata'))
rm(fit.k)
if (exists('wss')) rm(wss)
Sys.time()


########  PREDICT STRATA FOR NO DATA AREAS  #########

## Remove NA due to: 1. NA present in Baccini -  2. Predicted Strata Error is NA for the categorical predictors (CCI) without reference data
## NA values in strata are now predicted using datasets without NA (saa, vcf, hei) 
## Water (applied to Saatchi) was set to NA and remains NA (to have No strata in water)
## The RF model is trained using 10000 pixel randomly selected

str.map <- stack(Strata, saa, vcf)
names(str.map) <- c("strata", "saa", "vcf")
rm(hei)
rm(vcf)

dat <- as.data.frame(getValues(str.map))
dat <- dat[complete.cases(dat),]
colnames(dat) <- c("strata", "saa", "vcf")
dat$strata <- as.factor(dat$strata)
set.seed(30)
dat <- dat[sample(nrow(dat), 10000), ]

str.f <- formula(strata ~ saa + vcf)
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
Strata.fn <- 'EU'
Strata <- calc(strata, fun=strata.f, filename=paste("./Results/Strata/Strata_", Strata.fn, ".tif", sep=""), datatype='INT1U', overwrite=T, progress='text')
names(Strata) <- paste("Strata_", Strata.fn, sep="")
rm(list=ls(pattern='str.'))

Sys.time()

-------------------------
  
Strata <- raster('./Results/Strata/Strata_EU.tif')


-------------------------
  
  ##################################################################################
  ###########################   REFERENCE DATA  ####################################
##################################################################################

########  CONSOLIDATE REFERENCE DATA #####  START HERE FOR STRATA AS GLC2000 OR VCF10
# Strata <- raster(paste("./Results/Strata/Strata_", Strata.fn, ".tif", sep=""))

Strata.s <- 8    
Strata.names <- 1:Strata.s

## Each NFI has a unique code, and there is a raster map where each plot of the NFI has the code value (e.g.: all plots in Spain have value 1, in France is 2, etc.)
codes <- raster(paste('./Reference/Ref_code_', cont, '.tif', sep=""))
code.names <- as.character(read.csv(paste("./Reference/Codes_", cont, ".csv", sep=""))[,2])
code.n <- maxValue(codes)

# temp fix for testing
codes <- saa
codes[codes] <- 1

## mask out Ref data in Water (NA in Strata)

ref[is.na(Strata)] <- NA


## Analyze distribution of Reference data by Strata

ref.str <- stack(ref, Strata, codes)
ref.dat <- as.data.frame(getValues(ref.str))
ref.dat <- ref.dat[complete.cases(ref.dat),]
colnames(ref.dat) <- c("AGB", "Strata", "Code")
str.code <- table(ref.dat$Code, ref.dat$Strata, dnn = c("Code", "Strata"))
# rowSums(str.code)
# colSums(str.code)

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
# print(str.code.x, digits = 0)


# Plotting Ref data by Strata

Strata.fn <- 'NL' 
png(filename=paste("./Results/Reference/Reference_Data_", Strata.fn, "_Orig.png", sep=""))
plt <- barplot(str.code, beside=F, main=paste("Reference Data by ", Strata.fn, " strata - original", sep=""), xlab="Strata", legend = code.names, names.arg = Strata.names,
               cex.names = 0.7, col=c("red", "orange", "yellow", "green", "blue", "violet", "pink" , "cyan", "gray", "forestgreen", "orange3"), args.legend = list(x = "topright", cex=0.7))
text(plt, colSums(str.code), labels = colSums(str.code), pos = 3)
dev.off()
