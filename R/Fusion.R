# load libraries

require(randomForest)

######  INPUT DATA  ##########

saa <- raster(paste('./Input_Maps/saa_1km_', cont, '.tif', sep=""))
saa <- raster('C:/R_Projects/Biomass_Europe/Maps/Gallaun/1km/bmAg_JR2000_ll_1km_eur_Crop.tif')
bac <- raster(paste('./Input_Maps/bac_1km_', cont, '.tif', sep=""))  
ref <- raster(paste('./Reference/Ref_', cont, '.tif', sep=""))
vcf <- raster(paste('./Strata/VCF_', cont, '.tif', sep=""))
hei <- raster(paste('./Strata/HEI_', cont, '.tif', sep=""))
cci <- raster(paste('./Strata/CCI_', cont, '.tif', sep=""))
cci <- raster('Covariates/CCI_2005/CCI_crop.tif')
water <- raster(paste('./Strata/WATER/Water_1km_', cont, '.tif', sep=""))


###### WATER MASK APPLIED TO SAATCHI & BACCINI ############

saa[water == 1] <- NA
bac[water == 1] <- NA
rm(water)


##################################################################################
###########################   ERROR STRATA   #####################################
##################################################################################

saa.er <- ref - saa
bac.er <- ref - bac

error.map <- stack(saa.er, bac.er, hei, vcf, saa, bac, cci)
error.map <- stack(saa, cci)
names(error.map) <- c("saa","cci")
rm(cci)

error <- as.data.frame(getValues(error.map))
error <- error[complete.cases(error),]
colnames(error) <- c("saa","cci")
error$cci <- as.factor(error$cci)


### MODEL MAP ERROR with Random Forest  
# NOTE: mtry=2 to have comparable R2 without GLC2000, otherwise with 5 variables mtry=1 by default and the R2 drops)

saa.rf.f <- formula(saa ~ cci)
bac.rf.f <- formula(bac.er ~ hei + vcf + saa + bac + cci)
set.seed(55)
saa.rf <- randomForest(saa.rf.f, error, importance=F, mtry=2, ntree = 100)
set.seed(55)
bac.rf <- randomForest(bac.rf.f, error, importance=F, mtry=2)
# saa.rf
# bac.rf

rf.out <- matrix(0, nrow=2, ncol=2, dimnames=list(c("Predict_Error_Saatchi", "Strata_Error_rate"), c("RMSE", "R2")))
rf.out[1,1] <- sqrt(bac.rf$mse[500])
rf.out[1,2] <- mean(bac.rf$rsq[500])
rf.out[2,1] <- sqrt(saa.rf$mse[500])
rf.out[2,2] <- mean(saa.rf$rsq[500])

png(filename=paste("./Results/Strata/RF_Predict_Errors_VarImpPlot_", cont, ".png", sep=""), width=960, height=480, res=96)
opar <- par(mfrow=c(1,2))
varImpPlot(saa.rf)
varImpPlot(bac.rf)
dev.off()

varImp <- cbind(importance(saa.rf), importance(bac.rf))
colnames(varImp) <- c('Saatchi_Imp', 'Baccini_Imp')
write.csv(varImp, paste("./Results/Strata/IncNodePur_", cont, ".csv", sep=""))
rm(varImp)

png(filename=paste("./Results/Strata/RF_Predict_Errors_CrossVal_", cont, ".png", sep=""), width=960, height=480, res=96)
opar <- par(mfrow=c(1,2))
plot(error$saa.er,predict(saa.rf, newdata=error), main="Cross Validation Error Saatchi", xlab="Error", ylab="Predicted error", xlim=c(min(error$saa.er), max(error$saa.er)), ylim=c(min(error$saa.er),max(error$saa.er)))
abline(0,1)
plot(error$bac.er,predict(bac.rf, newdata=error), main="Cross Validation Error Baccini", xlab="Error", ylab="Predicted error", xlim=c(min(error$saa.er), max(error$saa.er)), ylim=c(min(error$saa.er),max(error$saa.er)))
abline(0,1)
dev.off()

