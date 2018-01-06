# Load libraries

if (!require(xlsx)) install.packages('xlsx')
if (!require(dplyr)) install.packages('dplyr')

dir.create(paste("./Results/Validation/BiomassPerCountry/", sep=""), showWarnings = F)

name.fn <- "EU" # set file name
Strata.fn <- name.fn # copy 
cluster_n <- 8 # set cluster solution
Strata.s <- cluster_n
Strata.fn <- paste0(name.fn,"_", cluster_n)
n_maps <- 3 # set number of input maps


# load ref raster datasets (Check which ref data to load!). Load either full raster or 30% validation raster.
#ref <- raster('./Maps/Ref/ref_ras_EU_final_round.tif') # all ref dataset
ref <- raster(paste("./Results/Validation/Reference/Ref_", Strata.fn, "_Val.tif", sep="")) # 30% validation ref dataset

# load Fused map
fused.map <- raster(paste('Results/Validation/Fused_Map/new/FUSED_FINAL_', Strata.fn, '.tif', sep=''))


# load biomass maps
Thur <- raster('./Maps/Thurner/1km/bmAg_Thurner_1km_Alligned.tif') 
IIASA <- raster("./Maps/IIASA/1km/bmAg_IIASA2010_Alligned.tif")
Gal <- raster("./Maps/all_maps/bmAg_JR2000_ll_1km_eur.tif")

# load strata raster
Strata <- raster(paste("./Results/Strata/Strata_", Strata.fn, ".tif", sep=""))

# load EU map
EU <- readOGR(dsn = "./data/Boundaries/Europe", layer = "Europe_in_Gal")

# load forest mask map
forest_m <- raster("./Forest_mask/Forest_mask_reproj.tif")

dir.create("./data/", showWarnings = F)
dir.create("./data/Boundaries", showWarnings = F)


#CountryShape <- raster::getData('GADM', country = Country, level=0, path = './data/Boundaries')
#coordsys <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")
#spTransform(CountryShape, coordsys)


Fused.final_Sin  <- reproject(fused.map , CRS = '+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs', program = 'GDAL', method = 'near')


CountryShape_Sin  <- reproject(EU , CRS = '+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs', program = 'GDAL', method = 'near')
Strata_Sin  <- reproject(Strata , CRS = '+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs', program = 'GDAL', method = 'near')
Gal.map_Sin  <- reproject(Gal, CRS = '+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs', program = 'GDAL', method = 'near')
IIASA.map_Sin <- reproject(IIASA, CRS = '+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs', program = 'GDAL', method = 'near')
Thur.map_Sin <- reproject(Thur, CRS = '+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs', program = 'GDAL', method = 'near')


# add forest mask here
all.map <- stack(Thur.map_Sin, IIASA.map_Sin, Gal.map_Sin ,Fused.final_Sin, Strata_Sin)
forest_m_Sin  <- reproject(forest_m, CRS = '+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs', program = 'GDAL', method = 'near')
forest_m_Sin[forest_m_Sin == 0] <- NA
all.map <- mask(all.map, forest_m_Sin, progress = 'text')

# Data prep
all.map <- crop(all.map, CountryShape_Sin)
all.map <- mask(all.map, CountryShape_Sin, progress = 'text')
all.dat <- as.data.frame(as.matrix(all.map, na.rm = TRUE))
#all.dat <- all.dat[complete.cases(all.dat),]                      ## here only data within the Baccini Extent is used
colnames(all.dat) <- c("Thur","IIASA","Gal","Fused" ,"Strata")
all.dat$Strata <- as.factor(all.dat$Strata)

n.rot <- length(names(all.map))-1 #rotation, looping amount
output <- matrix(data = NA, nrow =length(CountryShape_Sin), ncol = n.rot+1 )
names_col <- c("Thur","IIASA","Gal","Fused" ,"Strata")
colnames(output) <- c("Country",names_col[-length(names(all.map))])
pix.area <- res(Thur.map_Sin)[1]


# Compute AGB for the countries

for (b in 1:length(CountryShape_Sin)){
  #all.map <- crop(all.map, CountryShape_Sin[b,])
  
  flag <- tryCatch(sel.map <- crop(all.map, CountryShape_Sin[b,]), error = function(x) NULL)
  if (is.null(flag)){
    next
  } else {
    sel.map <- crop(sel.map , CountryShape_Sin[b,])
  }
  
  sel.map  <- mask(sel.map , CountryShape_Sin[b,])
  all.dat <- as.data.frame(as.matrix(sel.map, na.rm = TRUE))
  #all.dat <- all.dat[complete.cases(all.dat),]                      ## here only data within the Baccini Extent is used
  colnames(all.dat) <- names_col
  all.dat$Strata <- as.factor(all.dat$Strata)
  output[b,1] <- as.character(as.factor(CountryShape_Sin[b,]$NAME))
  for (i in 1:n.rot){
    agb.tot.s <- aggregate(all.dat[,i], by = list(all.dat$Strata), FUN="sum", na.rm=TRUE)
    if (nrow(agb.tot.s) == 0){
      next
    }
    agb.tot <- cbind(agb.tot.s[,2])
    agb.tot <- t(agb.tot*(((pix.area^2)/10000)))
    agb.sum <- rowSums(agb.tot) 
    output[b,i+1] <- agb.sum
  }
  print(paste0(b, " out of ", length(CountryShape_Sin))) 
}


# remove countries with all NA's
output_mill <- output[rowSums(is.na(output[,2:5]))!=4, ]
output_mill <- as.data.frame(output_mill)
# Set country column to character
output_mill$Country <- as.character(output_mill$Country)

# Identify factor columns
indx <- sapply(output_mill, is.factor)
# Replace factor columns by Numbers
output_mill[indx] <- lapply(output_mill[indx], function(x) as.numeric(as.character(x)))
# Recalculate biomass values to millions of tonnes
output_mill[2:5] <- output_mill[, 2:ncol(output_mill)]/1000000
# Round values
output_mill[2:5] <- round(output_mill[2:5])

# list of countries to remove
remove_country <- c("Iceland", "Austria", "Serbia", "San Marino","Cyprus") 

# remove names from dataframe
output2 <- output_mill[ ! output_mill$Country %in% remove_country, ]
# rename country
output2$Country[output2$Country == "The former Yugoslav Republic of Macedonia"] <- "Macedonia (FYROM)" 

# sort by name
sort.df <- with(output2,  output2[order(output2$Country) , ])
# write to excel file
write.xlsx(sort.df, file = paste("./Results/Validation/BiomassPerCountry/Biomass_",Strata.fn,"_new.xlsx", sep=""), sheetName = "Biomass", append = F)


# data prep for R2 computing
dat.r <- stack(Thur, IIASA, Gal, fused.map, ref, Strata)
forest_m[forest_m == 0] <- NA
dat.r <- mask(dat.r, forest_m, progress = 'text')
names(dat.r) <- c("Thur","IIASA","Gal","Fused" ,"ref", "Strata")
n.rot <- length(names(dat.r))-2
output_val <- matrix(data = NA, nrow =length(EU), ncol = n.rot+3 )
names_col <- c("Country","Thur","IIASA","Gal","Fused" ,"ref", "Strata")

colnames(output_val) <- names_col

# compute R2 for all EU countries
for (b in 1:length(EU)){
  #all.map <- crop(all.map, CountryShape_Sin[b,])
  
  flag <- tryCatch(sel.map <- crop(dat.r, EU[b,]), error = function(x) NULL)
  if (is.null(flag)){
    next
  } else {
    sel.map <- crop(sel.map , EU[b,])
  }
  
  sel.map  <- mask(sel.map , EU[b,])
  
  all.map <- as.data.frame(as.matrix(sel.map, na.rm = TRUE))
  all.map$Strata <- as.factor(all.map$Strata)
  output_val[b,"Country"] <- as.character(as.factor(EU[b,]$NAME))
 
  
  if (all(is.na(sel.map$ref@data@values)) == T){
    next
  }
  output_val[b,"ref"] <- sum(!is.na(sel.map$ref@data@values))
  
  resultsList <- list()
  
  for (i in 1:n.rot){
    
    sel <- subset(sel.map,i)
    sel.dat <-  getValues(sel)
    flag1 <- tryCatch(lmfit <- lm(sel.dat ~ all.map$ref,  data=all.map), error = function(x) NULL)
    if (is.null(flag1)){
      next
    } else {
      lmfit <- lm(sel.dat ~ all.map$ref)
    }
    resultsList[[i]] <- summary(lmfit)
    output_val[b,i+1] <- resultsList[[i]]$r.squared
  }
  print(paste0(b, " out of ", length(EU))) 
}

# remove countries with all NA's
output_val_df <- as.data.frame(output_val)
# Set country column to character
output_val_df$Country <- as.character(output_val_df$Country)

# Identify factor columns
indx <- sapply(output_val_df, is.factor)
# Replace factor columns by Numbers
output_val_df[indx] <- lapply(output_val_df[indx], function(x) as.numeric(as.character(x)))

# list of countries to remove
remove_country <- c("Iceland", "Austria", "Serbia", "San Marino","Cyprus", "Holy See (Vatican City)", "Guernsey", "Jersey", "Isle of Man") 

# remove names from dataframe
output3 <- output_val_df[ ! output_val_df$Country %in% remove_country, ]
# rename country
output3$Country[output3$Country == "The former Yugoslav Republic of Macedonia"] <- "Macedonia (FYROM)" 

# sort by name
sort.dff <- with(output3,  output3[order(output3$Country) , ])

sort.final <- sort.dff %>% mutate_if(is.numeric, funs(round(., 2)))


# write to excel file
write.xlsx(sort.final , file = paste("./Results/Validation/BiomassPerCountry/Biomass_",Strata.fn,"_R2_new.xlsx", sep=""), sheetName = "Biomass", append = F)


#model a linear model in R
#https://stackoverflow.com/questions/36716637/how-to-create-a-loop-for-a-linear-model-in-r



#summary(fusion.lm)
#summary(Gal.lm )
#summary(IIASA.lm)
#summary(Thur.lm)

#t <- extract(all.map$bmAg_IIASA2010_Alligned ,CountryShape_Sin1 , na.rm = T, fun = sum, df = T)


