# load libraries

if (!require(randomForest)) install.packages('randomForest')
if (!require(robust)) install.packages('robust')


######  INPUT DATA  ##########


Gal <- raster("./Maps/Gallaun/1km/bmAg_JR2000_ll_1km_eur.tif")

#Bar <- raster("./Maps/Barredo/barredo_Alligned.tif")

IIASA <- raster("./Maps/IIASA/1km/bmAg_IIASA2010_Alligned.tif")


Thur <- raster('./Maps/Thurner/1km/bmAg_Thurner_1km_Alligned.tif')  

ref <- raster('./Maps/Ref/ref_ras_EU_final_round.tif')

vcf <- raster("./Covariates/Outputs/MODIS_VCF/Mosaic/Mosaic_aggre_align.tif")

hei <- raster("./Covariates/Height/Height_align.tif")

cci <- raster("./Covariates/Outputs/CCI_2005/CCI_Mul_align.tif")

water <- raster("./Covariates/Outputs/CCI_Water/CCI_Water_aggr_Resam_R.tif")

###### WATER MASK APPLIED TO GalTCHI & ThurCINI ############

# Set parameters

Strata.fn <- 'EU_3'# set file name
cluster_n <- 3 # set cluster solution or strata
n_maps <- 3 # set number of input maps
cov_names <- c("hei", "vcf", "cci")
Covariates <- stack(hei, vcf, cci)

code.names <- 'EU'


##################################################################################
###########################   ERROR STRATA   #####################################
##################################################################################
ls <- list.files("./Maps/all_maps/", pattern = ".tif$", full.names = T)
ls2 <- ls
ls2[1] <- ls[3]
ls2[3] <- ls[1]
ls <- ls2
# ls <- stack(Gal, Gal_2, Thur, IIASA)

e_ls <- list()
for (i in 1:n_maps){
   map <- ref - raster(ls[i])
   #map <- ref - raster(ls[i]) # for lists
   names(map) <- paste0("Map",i,"_er")
   e_ls[[i]] <- map
}

# remove water from datasets
Gal[water == 2] <- NA
Thur[water == 2] <- NA
IIASA[water == 2] <- NA
#rm(water)


##################################################################################
###########################   ERROR STRATA   #####################################
##################################################################################

# substract input maps from ref data
Gal.er <- ref - Gal
Thur.er <- ref - Thur
IIASA.er <- ref - IIASA

error.map <- stack(Gal.er, Thur.er, IIASA.er, hei, vcf, Gal, Thur, IIASA, cci)
names(error.map) <- c("Gal.er", "Thur.er","IIASA.er", "hei", "vcf", "Gal", "Thur", "IIASA", "cci")
#rm(cci)
error.map[water == 2] <- NA

error <- as.data.frame(getValues(error.map))
error <- error[complete.cases(error),]
colnames(error) <- c("Gal.er", "Thur.er","IIASA.er", "hei", "vcf", "Gal", "Thur", "IIASA", "cci")
error$cci <- as.factor(error$cci)


### MODEL MAP ERROR with Random Forest  
# NOTE: mtry=2 to have comparable R2 without GLC2000, otherwise with 5 variables mtry=1 by default and the R2 drops)

Gal.rf.f <- formula(Gal.er ~ hei + cci + vcf + Thur + Gal + IIASA)
Thur.rf.f <- formula(Thur.er ~ hei + cci + vcf + Thur + Gal + IIASA)
IIASA.rf.f <- formula(IIASA.er ~ hei + cci + vcf + Thur + Gal + IIASA)

set.seed(55)
Gal.rf <- randomForest(Gal.rf.f, error,  mtry=2)
set.seed(55)
Thur.rf <- randomForest(Thur.rf.f, error,  mtry=2)
set.seed(55)
IIASA.rf <- randomForest(IIASA.rf.f, error, mtry=2)

name_layers <- names(error.map)

# --------------Automation testing------------------------------------------------------------------------------------------------


# names(Covariates) <- c(cov_names)
# error.map  <- stack(error.map, Covariates, IIASA, Gal,  Thur)
# 
# #error.map <- stack(Gal.er, Thur.er, IIASA.er, hei, vcf, Gal, Thur, IIASA, cci)
# #rm(cci)
# error.map <- stack(IIASA.er, Gal.er, Thur.er, hei, vcf,cci, IIASA, Gal, Thur)
# 
# error <- as.data.frame(getValues(error.map))
# error <- error[complete.cases(error),]
# colnames(error) <- names(error.map)
# error$cci <- as.factor(error$cci)
# 
# 
# ### MODEL MAP ERROR with Random Forest  
# # NOTE: mtry=2 to have comparable R2 without GLC2000, otherwise with 5 variables mtry=1 by default and the R2 drops)
# 
# # fix this for automatization
# for (i in 1:length(n.maps)){
#   xnam <- paste0("Covariates[[", 1:6, "]]")
#   Gal.rf.f <- as.formula(paste("e_ls[[1]] ~ ", paste(xnam, collapse= "+")))
# }
# 
# 
# Map1_er <- error.map$Map1_er
# Map2_er <- error.map$Map2_er
# Map3_er <- error.map$Map3_er
# #Map4_er <- error.map$Map4_er
# 
# hei <- Covariates[[1]]
# vcf <- Covariates[[2]]
# cci <- Covariates[[3]]
# 
# # Gal, Thur, IIASA, Bar
# IIASA.rf.f <- formula(IIASA.er ~ hei + cci + vcf + IIASA + Gal + Thur)
# 
# Bar.rf.f <- formula(Map1_er  ~ hei + vcf + cci + Gal + IIASA + Thur)
# 
# IIASA.rf.f <- formula(IIASA.er  ~  hei + vcf + cci + Gal + IIASA + Thur)
# IIASA.rf.f <- formula(Map2_er  ~  hei + vcf + cci + Gal + IIASA + Thur)
# IIASA.rf.f <- formula(test2  ~  hei + vcf + cci + Gal + IIASA + Thur)
# IIASA.rf.f <- formula(test3  ~  hei + vcf + cci + Gal + IIASA + Thur)
# 
# Gal.rf.f <- formula(Map3_er  ~ Bar + Gal + IIASA + Thur + hei + vcf + cci)
# Thur.rf.f <- formula(Map4_er  ~  Bar + Gal + IIASA + Thur + hei + vcf + cci)
# 
# rf_models <- c(IIASA.rf.f,Gal.rf.f,Thur.rf.f)
# #rf_models <- c(Bar.rf.f, IIASA.rf.f,Gal.rf.f,Thur.rf.f)
# 
# Gal.rf <- randomForest(IIASA.rf.f, error,  mtry=2)

# --------------Automation testing------------------------------------------------------------------------------------------------


# bind models
rf_models <- c(IIASA.rf.f,Gal.rf.f,Thur.rf.f)

# load and run models in randomForest
dir.create("./Results/RF_Models", showWarnings = F)
for (i in 1:length(ls)){
  set.seed(55)
  Map.rf <- randomForest(rf_models[[i]], error, mtry=2)
  save(Map.rf,file = paste0("./Results/RF_Models/",name_layers[i], "_RF.RData"))
}
rm(Map.rf)

dimnam_t <- c()
for (i in 1:n_maps){
  dimnam <- paste0("Predict_Error_",name_layers[i])
  dimnam_t <- c(dimnam_t, dimnam) 
}

rf.out <- matrix(0, nrow=n_maps, ncol=2, dimnames=list(dimnam_t, c("RMSE", "R2")))

for (i in 1:n_maps){
  load(paste0("./Results/RF_Models/",name_layers[i], "_RF.RData"))
  rf.out[i,1] <- sqrt(Map.rf$mse[500])
  rf.out[i,2] <- mean(Map.rf$rsq[500])
}

png(filename=paste("./Results/Strata/RF_Predict_Errors_VarImpPlot_",Strata.fn,".png", sep=""), width=960, height=480, res=96)
opar <- par(mfrow=c(1,n_maps))
for (i in 1:n_maps){
  load(paste0("./Results/RF_Models/",name_layers[i], "_RF.RData"))
  varImpPlot(Map.rf, main = name_layers[i])
}
dev.off()

VarImp <- c()
for (i in 1:n_maps){
  load(paste0("./Results/RF_Models/",name_layers[i], "_RF.RData"))
  imp <- importance(Map.rf)
  VarImp <- cbind(VarImp,imp)
}

colnames(VarImp) <- name_layers
write.csv(VarImp, paste("./Results/Strata/IncNodePur_",Strata.fn,".csv", sep=""))
rm(VarImp)


png(filename=paste("./Results/Strata/RF_Predict_Errors_CrossVal_",Strata.fn,".png", sep=""), width=960, height=480, res=96)
opar <- par(mfrow=c(1,n_maps))
for (i in 1:length(ls)){
  load(paste0("./Results/RF_Models/",name_layers[i], "_RF.RData"))
  plot(error[[i]],predict(Map.rf, newdata=error), main= paste0("Cross Validation Error ",name_layers[i]), xlab="Error", ylab="Predicted error", xlim=c(min(error[[i]]), max(error[[i]])), ylim=c(min(error[[i]]),max(error[[i]])))
  abline(0,1)
}
dev.off()

### PREDICT MAP ERROR

error.dat <- dropLayer(error.map, 1:n_maps)

ls_ras <- list()
for (i in 1:length(ls)){
  load(paste0("./Results/RF_Models/",name_layers[i], "_RF.RData"))
  nam <- paste0(name_layers[i],"_Error_RF")
  assign(nam,predict(error.dat, Map.rf, type= 'response', filename=paste0("./Results/Strata/", name_layers[i],"_Error_RF_",Strata.fn,".tif"), datatype='INT2S', overwrite=T, progress='text'))
  ls_ras[[i]] <- paste0("./Results/Strata/", name_layers[i],"_Error_RF_",Strata.fn,".tif")
}

#rm(error.map)
err.all <- stack(unlist(ls_ras))
mydata <- as.data.frame(getValues(err.all))
mydata.k <- mydata[complete.cases(mydata),]  # Remove NA, K-means cannot handle NAs
#rm(list=ls(pattern='er'))
#rm(list=ls(pattern='.rf'))

# K cluser
set.seed(55)
fit.k <- kmeans(mydata.k, cluster_n)      
Strata.s <- cluster_n 


# Convert to raster
mydata$cluster <- NA
mydata[rownames(mydata.k), "cluster"] <- fit.k$cluster
Strata <- Covariates[[1]]
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

#str.map <- stack(e_ls)

# --------------Automation testing------------------------------------------------------------------------------------------------

str.map <- stack(Strata, Gal, vcf, hei)
names(str.map) <- c("strata", "Gal", "vcf", "hei")


#names(str.map) <- names(stack(Strata,e_ls))
#name_layers <- names(str.map)

# --------------Automation testing------------------------------------------------------------------------------------------------


#rm(hei)
#rm(vcf)

dat <- as.data.frame(getValues(str.map))
dat <- dat[complete.cases(dat),]
colnames(dat) <- names(str.map)
dat$strata <- as.factor(dat$strata)
set.seed(30)
dat <- dat[sample(nrow(dat), 10000), ]

str.f <- formula(strata ~ Gal + vcf + hei) 
set.seed(55)
str.rf <- randomForest(str.f, dat, importance=T)
varImpPlot(str.rf)

rf.out[2,1] <- 1-sum(diag(table(predict(str.rf), str.rf$y)))/sum(table(predict(str.rf), str.rf$y))
rf.out[2,2] <- NA
write.csv(rf.out, paste("./Results/Strata/RF_Predict_Errors_",Strata.fn,".csv", sep=""))
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

Strata.names <- 1:Strata.s

## Each NFI has a unique code, and there is a raster map where each plot of the NFI has the code value (e.g.: all plots in Spain have value 1, in France is 2, etc.)
codes <- raster(paste('./Reference/Ref_code_EU.tif', sep=""))
codes <- crop(codes, Covariates[[1]])
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

png(filename=paste("./Results/Reference/Reference_Data_", Strata.fn, "_Orig.png", sep=""))
plt <- barplot(str.code, beside=F, main=paste("Reference Data by ", Strata.fn, " strata - original", sep=""), xlab="Strata", legend = code.names, names.arg = Strata.names,
               cex.names = 0.7, col=c("red", "orange", "yellow", "green", "blue", "violet", "pink" , "cyan", "gray", "forestgreen", "orange3"), args.legend = list(x = "topright", cex=0.7))
text(plt, colSums(str.code), labels = colSums(str.code), pos = 3)
dev.off()

##################################################################################
###########################   BIOMASS FUSION   ###################################
##################################################################################

####### CALCULATE THE WEIGHTS & BIAS

# Error Maps
e_ls <- list()
for (i in 1:length(ls)){
  map <- ref - raster(ls[i])
  names(map) <- paste0("Map",i,"_er")
  e_ls[[i]] <- map
}

error.map <- stack(e_ls)
error.map  <- stack(Strata, error.map)


error.map[water == 2] <- NA
#error.map <- stack(Gal.er, Thur.er, IIASA.er, hei, vcf, Gal, Thur, IIASA, cci)
#rm(cci)

error <- as.data.frame(as.matrix(error.map, na.rm = TRUE))
error <- error[complete.cases(error),]
colnames(error) <- c(name_layers)
colnames(error)[1] <- c("Strata")
error$Strata <- as.factor(error$Strata)
#rm(saa.er)
#rm(bac.er)
#rm(error.map)


### Calculate the bias
rm(name_layers)
name_layers <- names(stack(e_ls))

bias <-c() 
for (i in 1:length(ls)){
  map.bias <- tapply(error[[i+1]], error$Strata, mean, na.rm = TRUE)
  bias <- cbind(bias, map.bias)
}
colnames(bias) <- c(paste0(name_layers, ".bias"))
#rm(list=ls(pattern=".bias"))


### Calculate variance-covariance matrix and weight matrix using a ROBUST ESTIMATOR

X <- rep.int(1,n_maps) 

cov <- vector('list',length = Strata.s)
for (i in 1:Strata.s){
  cov[[i]] <- covRob(error[error$Strata==i,2:ncol(error)]) # set length (-:-) to the selected datasets
  cov[[i]] <- cov[[i]]$cov}

w <- vector('list',length = Strata.s)
for (i in 1:Strata.s){
  w[[i]] <- solve(t(X) %*% solve(cov[[i]]) %*% X) %*% t(X) %*% solve(cov[[i]]) }

weight <- matrix(1:(length(X)*Strata.s), ncol=(length(X)))
for (i in 1:Strata.s){ weight[i,] <- w[[i]] }
colnames(weight) <- c(paste0(name_layers, ".weight"))
rownames(weight) <- c(1:Strata.s)


## Finalize Bias and Weights

weight[weight < 0] <- 0   # set weights to 0 - 1 limits
weight[weight > 1] <- 1   # set weights to 0 - 1 limits


#------------------------------------------------

### calculate Error Variance of the Fused map

v <- vector('list', length = Strata.s)
for (i in 1:Strata.s){  v[[i]] <- solve(t(X) %*% solve(cov[[i]]) %*% X) }

v.err <- matrix(1:Strata.s)
for (i in 1:Strata.s){ v.err[i] <- v[[i]] }

colnames(v.err) <- c("Fus_Var")
rownames(v.err) <- c(1:Strata.s)
rm(list=c('cov', 'i', 'v', 'w', 'X'))


### Compute Error Variance of Input maps

tot.var <-c()
tot.var <- aggregate(error[[2]], by = list(error$Strata), FUN="var")[1]
for (i in 1:length(ls)){
  map.var <- aggregate(error[[i+1]], by = list(error$Strata), FUN="var")[,2]
  tot.var <- cbind(tot.var,map.var)
}
colnames(tot.var) <- c("Strata",paste0(name_layers, ".var"))

### Compute n. pixel per strata

n.pix <- matrix(1:Strata.s)
for (i in 1:Strata.s) { n.pix[i] <- freq(Strata, useNA='no', value = i) }
n.pix.w <- n.pix / (sum(n.pix))


### Compile Error Variances

err.var <- cbind(tot.var,v.err, n.pix, n.pix.w)
err.var_c <- cbind(tot.var,v.err, n.pix, n.pix.w) # copy dataset for names

for (i in 1:(length(ls)+1)){
  err.var[length(err.var)+1] <- err.var[[i+1]] * err.var$n.pix.w
}
colnames(err.var) <- c(names(err.var_c), paste0(name_layers, "_var_w"), "Fus_Var_w")



fus.par <- aggregate(error$Strata, by = list(error$Strata), FUN="length")
colnames(fus.par) <- c("Strata", "N")
fus.par$Strata <- as.numeric(as.character(fus.par$Strata))
fus.par <- cbind(fus.par, bias, weight, err.var[,-1])   # adjust the number of classes
#fus.par[9,] <- c(9, 0, 0, 0, 0.5, 0.5, rep(0,8))
weight_fill <- 1/n_maps
fus.par[Strata.s+1,] <- c(Strata.s+1, 0, rep(0,n_maps), rep(weight_fill,n_maps), rep(0,((n_maps*2)+4))) # Should weights be set according to the number of input maps? e.g. 3 maps == weight 0.33? check this

bias <- fus.par[,3:(n_maps+2)]                             # Add Strata 9 to bias and weight (for Fusion)
weight <- fus.par[,(3+n_maps):(2+(n_maps*2))]
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

###  FUSION           
maps <- stack(Strata, stack(ls))
#maps[maps < 1] <- NA
n.layers <- nlayers(maps)

# set parameters for biomass fuction
bias_matrix = as.matrix(bias)
weight_matrix = as.matrix(weight)
Strata.minus1 = 1:(Strata.n-1)
Strata.plus1 = 2:(n_maps+1)

biomass.fusion3 <- function(x) {
  n <- x[1] # get the stratum. Stratum should be first in the raster stack!!
  g <- x[Strata.plus1] + bias_matrix[n,] # add bias to raster values
  #g[g < 0] <- 0 # set values below 0 to 0
  w <- weight_matrix[n,] # get correct strata weight values
  w[is.na(g)]<- NA # set weight to NA if (g) raster values are NA
  p <- sum(w, na.rm = T) # calculate sum of weight values
  pp <- w/p # divide weight values by sum to get the proportion to == 1
  result <- as.integer(sum(pp*g, na.rm = T)) # return raster value
  return(result)
}

system.time(Fused.map <- calc(maps, fun = biomass.fusion3, progress = 'text'))

dir.create("./Results/Fused_Map/New/", showWarnings = F)
Fused.map[Fused.map < 0] <- NA # set biomass values below 0 to 0
writeRaster(Fused.map, filename = paste("Results/Fused_Map/New/Fused.mapX_", Strata.fn, ".tif", sep=""), datatype='INT4U', overwrite=T) # datatype = FLT4S
Fused.map <- raster(paste("Results/Fused_Map/New/Fused.mapX_", Strata.fn, ".tif", sep=""))


Fused.final <- Fused.map


##################################################################################
##############################   MODEL RESULTS   #################################
##################################################################################


########### TRAINING ERROR  (NB: using the Consolidated Reference dataset)

# Input data

fus.par <- read.csv(paste("./Results/Fused_Map/Bias_Weights_", Strata.fn, ".csv", sep=""))
input_maps <- stack(ls)
input_maps[input_maps < 1] <- NA
dat.r <- stack(input_maps, Fused.final, ref, Strata)
rm(input_maps)


dat <- as.data.frame(as.matrix(dat.r, na.rm = TRUE))
dat <- dat[complete.cases(dat),]
colnames(dat) <- c(name_layers,"Fused", "Ref", "Strata")

a <- dat[[1]]
for (i in 1:(length(ls)-1)){
  sum_maps <- a + dat[[i+1]] 
}
dat$Mean <- sum_maps/n_maps


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

err.all <- t(rbind(bias.m[-(n.maps+2):-(n.maps+3)], rmse.m[-(n.maps+2):-(n.maps+3)]))
n.plot <- nrow(err)
smse <- mean(smse)                     # simple Mean SMSE (not area-weighted)
err.all <- rbind(err.all, n.plot, smse)
colnames(err.all) <- c("Bias", "RMSE")
dir.create("./Results/Validation/", showWarnings = F)
dir.create("./Results/Validation/Training_Error", showWarnings = F)

write.csv(err.all, paste("Results/Validation/Training_Error/Training_Error_", Strata.fn, ".csv", sep=""))


# Plot errors #[FIX]
# --------------Automation testing------------------------------------------------------------------------------------------------

png(filename=paste("Results/Validation/Training_Error/Training_Error_Gal_IIASA_Thur_", Strata.fn, ".png", sep=""))
plot(dat$Ref, dat$Map1_er, xlim=c(0, max(dat)), ylim=c(0,max(dat)), type="p", pch=20, col='green', xlab="Reference data", ylab="Biomass map")
points(dat$Ref, dat$Map2_er, pch=20, col='red')
points(dat$Ref, dat$Map3_er, pch=20, col='blue')
legend(x='topleft', legend=c("Gal", "IIASA", "Thur"), col=c('green', 'red', 'blue'), pch=20)
abline(0,1)
dev.off()


# Plot Fused errors #[FIX]
# --------------Automation testing------------------------------------------------------------------------------------------------

png(filename=paste("Results/Validation/Training_Error/Training_Error_Fused_", Strata.fn, ".png", sep=""))
plot(dat$Ref, dat$Fused, xlim=c(0, max(dat)), ylim=c(0,max(dat)), type="p", pch=20, col='blue', xlab="Reference data", ylab="Biomass map")
legend(x='topleft', legend="Fused map", col='blue', pch=20)
abline(0,1)
dev.off()

rm(list=c('dat', 'err', 'err.all', 'err.sd', 'err.str', 'Gal', 'bias.m', 'Fused.final', 'n.plot', 'ref', 'rmse', 'rmse.m', 'IIASA','Thur', 'smse'))



#######################################################################################
########################## VALIDATION OF FUSED MAP  ###################################
#######################################################################################

## Method: Split consolidated Ref. data in Cal / Val sets (70% / 30%)
## make a fused map using Cal data and compute RMSE on Val data
## NB: directory are different from normal Fusion, outputs are in /Results/Validation/, and outputs are reduced

################## SPLIT REFERENCE DATA IN CAL / VAL DATASETS #########################
Sys.time()

## Reference dataset (Original)


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



##################################################################################
###########################   ERROR STRATA   #####################################
##################################################################################

ls <- list.files("./Maps/all_maps/", pattern = ".tif$", full.names = T)

e_ls <- list()
for (i in 1:n_maps){
  map <- ref - raster(ls[[i]])
  names(map) <- paste0("Map",i,"_er")
  e_ls[[i]] <- map
}

###### WATER MASK APPLIED
error.map <- stack(e_ls)
error.map[water == 2] <- NA

names(error.map) <- names(stack(e_ls))
name_layers <- names(error.map)

names(Covariates) <- c(cov_names)
error.map  <- stack(error.map, Covariates)

#error.map <- stack(Gal.er, Thur.er, IIASA.er, hei, vcf, Gal, Thur, IIASA, cci)
#rm(cci)

error <- as.data.frame(getValues(error.map))
error <- error[complete.cases(error),]
colnames(error) <- names(error.map)
error$cci <- as.factor(error$cci) # FIX


### MODEL MAP ERROR with Random Forest  
# NOTE: mtry=2 to have comparable R2 without GLC2000, otherwise with 5 variables mtry=1 by default and the R2 drops)


# fix this for automatization
for (i in 1:length(n.maps)){
  xnam <- paste0("Covariates[[", 1:6, "]]")
  Gal.rf.f <- as.formula(paste("e_ls[[1]] ~ ", paste(xnam, collapse= "+")))
}

Map1_er <- error.map$Map1_er
Map2_er <- error.map$Map2_er
Map3_er <- error.map$Map3_er
Map4_er <- error.map$Map4_er

hei <- Covariates[[1]]
vcf <- Covariates[[2]]
cci <- Covariates[[3]]

Bar.rf.f <- formula(Map1_er  ~ Map1_er + Map2_er + Map3_er + Map4_er + hei + vcf + cci)
IIASA.rf.f <- formula(Map2_er  ~ Map1_er + Map2_er + Map3_er + Map4_er + hei + vcf + cci)
Gal.rf.f <- formula(Map3_er  ~ Map1_er + Map2_er + Map3_er + Map4_er + hei + vcf + cci)
Thur.rf.f <- formula(Map4_er  ~ Map1_er + Map2_er + Map3_er + Map4_er + hei + vcf + cci)

rf_models <- c(Bar.rf.f, IIASA.rf.f,Gal.rf.f,Thur.rf.f)

dir.create("./Results/RF_Models", showWarnings = F)
for (i in 1:length(ls)){
  set.seed(55)
  Map.rf <- randomForest(rf_models[[i]], error, importance=F, mtry=2)
  save(Map.rf,file = paste0("./Results/RF_Models/",name_layers[i], "_RF.RData"))
}
rm(Map.rf)

dimnam_t <- c()
for (i in 1:n_maps){
  dimnam <- paste0("Predict_Error_",name_layers[i])
  dimnam_t <- c(dimnam_t, dimnam) 
}

rf.out <- matrix(0, nrow=n_maps, ncol=2, dimnames=list(dimnam_t, c("RMSE", "R2")))

for (i in 1:n_maps){
  load(paste0("./Results/RF_Models/",name_layers[i], "_RF.RData"))
  rf.out[i,1] <- sqrt(Map.rf$mse[500])
  rf.out[i,2] <- mean(Map.rf$rsq[500])
}

png(filename=paste("./Results/Strata/RF_Predict_Errors_VarImpPlot_",Strata.fn,".png", sep=""), width=960, height=480, res=96)
opar <- par(mfrow=c(2,2))
for (i in 1:n_maps){
  load(paste0("./Results/RF_Models/",name_layers[i], "_RF.RData"))
  varImpPlot(Map.rf, main = name_layers[i])
}
dev.off()

VarImp <- c()
for (i in 1:n_maps){
  load(paste0("./Results/RF_Models/",name_layers[i], "_RF.RData"))
  imp <- importance(Map.rf)
  VarImp <- cbind(VarImp,imp)
}

colnames(VarImp) <- name_layers
write.csv(VarImp, paste("./Results/Strata/IncNodePur_",Strata.fn,".csv", sep=""))
rm(VarImp)


png(filename=paste("./Results/Strata/RF_Predict_Errors_CrossVal_",Strata.fn,".png", sep=""), width=960, height=480, res=96)
opar <- par(mfrow=c(1,n.maps))
for (i in 1:length(ls)){
  load(paste0("./Results/RF_Models/",name_layers[i], "_RF.RData"))
  plot(error[[i]],predict(Map.rf, newdata=error), main= paste0("Cross Validation Error ",name_layers[i]), xlab="Error", ylab="Predicted error", xlim=c(min(error[[i]]), max(error[[i]])), ylim=c(min(error[[i]]),max(error[[i]])))
  abline(0,1)
}
dev.off()



### PREDICT MAP ERROR # FIX


error.dat <- dropLayer(error.map, 1:n.maps)

dir.create("./Validation/RF_Models", showWarnings = F)
for (i in 1:length(ls)){
  set.seed(55)
  Map.rf <- randomForest(rf_models[[i]], error, importance=F, mtry=2)
  save(Map.rf,file = paste0("./Validation/RF_Models/",name_layers[i], "_RF.RData"))
}

ls_ras <- list()
for (i in 1:length(ls)){
  load(paste0("./Validation/RF_Models/",name_layers[i], "_RF.RData"))
  nam <- paste0(name_layers[i],"_Error_RF")
  assign(nam,predict(error.dat, Map.rf, type= 'response', filename=paste0("./Validation/Strata/", name_layers[i],"_Error_RF_",Strata.fn,".tif"), datatype='INT2S', overwrite=T, progress='text'))
  ls_ras[[i]] <- paste0("./Validation/Strata/", name_layers[i],"_Error_RF_",Strata.fn,".tif")
}

#rm(error.map)
err.all <- stack(unlist(ls_ras))
mydata <- as.data.frame(getValues(err.all))
mydata.k <- mydata[complete.cases(mydata),]  # Remove NA, K-means cannot handle NAs
#rm(list=ls(pattern='er'))
#rm(list=ls(pattern='.rf'))

# K cluser
set.seed(55)
fit.k <- kmeans(mydata.k, cluster_n)      
Strata.s <- cluster_n 


# Convert to raster
mydata$cluster <- NA
mydata[rownames(mydata.k), "cluster"] <- fit.k$cluster
Strata <- Covariates[[1]]
values(Strata) <- mydata$cluster
dataType(Strata) <- 'INT1U'

#rm(list=ls(pattern='mydata'))
#rm(fit.k)
#if (exists('wss')) rm(wss)
Sys.time()



########  PREDICT STRATA FOR NO DATA AREAS  #########

str.map <- stack(e_ls)
names(str.map) <- names(stack(e_ls))
name_layers <- names(str.map)

#rm(hei)
#rm(vcf)

dat <- as.data.frame(getValues(str.map))
dat <- dat[complete.cases(dat),]
colnames(dat) <- names(str.map)
dat$strata <- as.factor(dat$strata)
set.seed(30)
dat <- dat[sample(nrow(dat), 10000), ]

str.f <- formula(strata ~ Gal + vcf + hei) # Fix
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

Strata <- calc(strata, fun=strata.f, filename=paste("./Results/Validation/Strata/Strata_", Strata.fn, ".tif", sep=""), datatype='INT1U', overwrite=T, progress='text')
names(Strata) <- paste("Strata_", Strata.fn, sep="")
rm(list=ls(pattern='str.'))

Sys.time()



##################################################################################
###########################   REFERENCE DATA  ####################################
##################################################################################

########  CONSOLIDATE REFERENCE DATA #####  START HERE FOR STRATA AS GLC2000 OR VCF10

Strata.names <- 1:Strata.s


## Each NFI has a unique code, and there is a raster map where each plot of the NFI has the code value (e.g.: all plots in Spain have value 1, in France is 2, etc.)
codes <- raster(paste('./Reference/Ref_code_EU.tif', sep=""))
codes <- crop(codes, Gal) #Fix
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


png(filename=paste("./Results/Validation/Reference/Reference_Data_", Strata.fn, "_Orig.png", sep=""))
plt <- barplot(str.code, beside=F, main=paste("Reference Data by ", Strata.fn, " strata - original", sep=""), xlab="Strata", legend = code.names, names.arg = Strata.names,
               cex.names = 0.7, col=c("red", "orange", "yellow", "green", "blue", "violet", "pink" , "cyan", "gray", "forestgreen", "orange3"), args.legend = list(x = "topright", cex=0.7))
text(plt, colSums(str.code), labels = colSums(str.code), pos = 3)
dev.off()



##################################################################################
###########################   BIOMASS FUSION   ###################################
##################################################################################

####### CALCULATE THE WEIGHTS & BIAS


# Error Maps

# load input maps
error.map <- stack(e_ls)
error.map[water == 2] <- NA

names(error.map) <- names(stack(e_ls))
name_layers <- names(error.map)

names(Covariates) <- c(cov_names)
error.map  <- stack(error.map, Covariates)

error <- as.data.frame(getValues(error.map))
error <- error[complete.cases(error),]
colnames(error) <- names(error.map)
error$Strata <- as.factor(error$Strata)

#rm(error.map)


# Plot Errors of IIASA and Thurner # Fix

# png(filename= paste0("./Results/Validation/Reference/Errors of IIASA_Thur_", Strata.fn, ".png"))
# plot(error$IIASA_er, error$Thur_er, main=paste("Errors of input maps"), xlab="IIASA Error", ylab="Thur Error", xlim=c(min(error$IIASA_er, error$Thur_er), max(error$IIASA_er, error$Thur_er)), ylim=c(min(error$Thur_er, error$IIASA_er), max(error$Thur_er, error$IIASA_er)))
# abline(0,1)
# dev.off()


### Calculate the bias

bias <-c() 
for (i in 1:length(ls)){
  map.bias <- tapply(error[[i+1]], error$Strata, mean, na.rm = TRUE)
  bias <- cbind(bias, map.bias)
}
colnames(bias) <- c(paste0(name_layers, ".bias"))


### Calculate variance-covariance matrix and weight matrix using a ROBUST ESTIMATOR

X <- rep.int(1,n.maps)

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


tot.var <-c()
tot.var <- aggregate(error[[2]], by = list(error$Strata), FUN="var")[1]
for (i in 1:length(ls)){
  map.var <- aggregate(error[[i+1]], by = list(error$Strata), FUN="var")[,2]
  tot.var <- cbind(tot.var,map.var)
}
colnames(tot.var) <- c("Strata",paste0(name_layers, ".var"))

### Compute n. pixel per strata

n.pix <- matrix(1:Strata.s)
for (i in 1:Strata.s) { n.pix[i] <- freq(Strata, useNA='no', value = i) }
n.pix.w <- n.pix / (sum(n.pix))


### Compile Error Variances
err.var <- cbind(tot.var,v.err, n.pix, n.pix.w)
err.var_c <- cbind(tot.var,v.err, n.pix, n.pix.w) # copy dataset for names


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
weight_fill <- 1/n.maps
fus.par[Strata.s+1,] <- c(Strata.s+1, 0, rep(0,n.maps), rep(weight_fill,n.maps), rep(0,((n.maps*2)+3))) # Should weights be set according to the number of input maps? e.g. 3 maps == weight 0.33? check this

bias <- fus.par[,3:(n.maps+2)]                             # Add Strata 9 to bias and weight (for Fusion)
weight <- fus.par[,(3+n.maps):(2+(n.maps*2))]
dir.create("./Results/Validation/Fused_Map", showWarnings = F)
write.csv(fus.par, paste("./Results/Validation/Fused_Map/Bias_Weights_", Strata.fn, ".csv", sep=""), row.names = FALSE)


## Map Uncertainty (Standard deviation of Error of Fused map)

err.rcl <- cbind(fus.par$Strata[1:Strata.s], sqrt(fus.par$Fus_Var[1:Strata.s]))
uncer <- reclassify(Strata, err.rcl, filename=paste("./Results/Validation/Fused_Map/Uncertainty_", Strata.fn, ".tif", sep=""), datatype='FLT4S', overwrite=T)


rm(list=c('fus.par', 'i', 'n.pix', 'n.pix.w', 'uncer'))
rm(list=ls(pattern="err"))


####### MAP FUSION  

## Set NA in Strata to a value (Strata must have always values in Fusion, no NA)

Strata.n <- Strata.s + 1
Strata[is.na(Strata)] <- Strata.n  

###  FUSION                          

maps <- stack(ls, Strata)
maps[maps < 1] <- NA
n.layers <- nlayers(maps)

# set parameters for biomass fuction
bias_matrix = as.matrix(bias)
weight_matrix = as.matrix(weight)
Strata.minus1 = 1:(Strata.n-1)
Strata.plus1 = 2:(Strata.n)
Strata.plus1 = 2:(n.maps+1)

biomass.fusion3 <- function(x) {
  n <- x[1] # get the stratum. Stratum should be first in the raster stack!!
  g <- x[Strata.plus1] + bias_matrix[n,] # add bias to raster values
  #g[g < 0] <- 0 # set values below 0 to 0
  w <- weight_matrix[n,] # get correct strata weight values
  w[is.na(g)]<- NA # set weight to NA if (g) raster values are NA
  p <- sum(w, na.rm = T) # calculate sum of weight values
  pp <- w/p # divide weight values by sum to get the proportion to == 1
  result <- as.integer(sum(pp*g, na.rm = T)) # return raster value
  return(result)
}

system.time(Fused.map <- calc(maps, fun = biomass.fusion3, progress = 'text'))

dir.create("./Results/Validation/Fused_Map/", showWarnings = F)
Fused.map[Fused.map < 0] <- NA
writeRaster(Fused.map, filename = paste("Results/Validation/Fused_Map/Fused.map_", Strata.fn, ".tif", sep=""), datatype='INT4U', overwrite=T) # datatype = FLT4S
Fused.final <- raster(paste("Results/Validation/Fused_Map/Fused.map_", Strata.fn, ".tif", sep=""))

rm(list=c("Fused.final", "maps"))





##################################################################################
###############################   VALIDATION   ###################################
##################################################################################

Fused.final <- raster(paste('./Results/Validation/Fused_Map/Fused.map_', Strata.fn, '.tif', sep=''))
Strata <- raster(paste("./Results/Validation/Strata/Strata_", Strata.fn, ".tif", sep=""))

## Validation dataset

ref <- raster(paste("./Results/Validation/Reference/Ref_", Strata.fn, "_Val.tif", sep=""))    

dat.r <- stack(e_ls, Fused.final, ref, Strata)
dat <- as.data.frame(as.matrix(dat.r, na.rm = TRUE))
dat <- dat[complete.cases(dat),]
colnames(dat) <- c(name_layers,"Fused", "Ref", "Strata")

# dat$Mean <- (dat$Gal + dat$IIASA + dat$Thur) / 3
a <- dat[[1]]
for (i in 1:(length(ls)-1)){
  sum_maps <- a + dat[[i+1]] 
}
dat$Mean <- sum_maps/n_maps
rm(dat.r)




#test R2
#[FIX] for automation
fusion.lm = lm(dat$Fused ~ dat$Ref, data=dat)
Gal.lm = lm(dat$Gal ~ dat$Ref, data=dat)
IIASA.lm = lm(dat$IIASA ~ dat$Ref, data=dat)
Thur.lm = lm(dat$Thur~ dat$Ref, data=dat)

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


# Plot errors

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
