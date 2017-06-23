Country <- 'NLD' # Netherlands

name.fn <- 'EU'# set file name
cluster_n <- 8 # set cluster solution
Strata.s <- cluster_n
Strata.fn <- paste0(name.fn,"_", cluster_n)
n_maps <- 3 # set number of input maps
n.maps <- 3 # set number of input maps


dir.create("./data/", showWarnings = F)
dir.create("./data/Boundaries", showWarnings = F)

EU <- readOGR(dsn = "./data/Boundaries/Europe", layer = "Europe_in_Gal")


CountryShape <- raster::getData('GADM', country = Country, level=0, path = './data/Boundaries')
coordsys <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")
spTransform(CountryShape, coordsys)

fused.map <- raster(paste("Results/Validation/Fused_map/FUSED_FINAL_", Strata.fn, ".tif", sep=""))
Fused.final_Sin  <- reproject(fused.map , CRS = '+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs', program = 'GDAL', method = 'near')

Thur <- raster('./Maps/Thurner/1km/bmAg_Thurner_1km_Alligned.tif') 
IIASA <- raster("./Maps/IIASA/1km/bmAg_IIASA2010_Alligned.tif")
Gal <- raster("./Maps/all_maps/bmAg_JR2000_ll_1km_eur.tif")
Gal <- raster("./Maps/Barredo/barredo_Alligned.tif")
Gal[Gal > 1000] <- NA
dataType(Gal) <- "INT4U"
Strata <- raster(paste("./Results/Strata/Strata_", Strata.fn, ".tif", sep=""))

CountryShape_Sin  <- reproject(EU , CRS = '+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs', program = 'GDAL', method = 'near')
CountryShape_Sin <- crop(CountryShape_Sin, Gal.map_Sin)
CountryShape_Sin1  <- reproject(CountryShape , CRS = '+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs', program = 'GDAL', method = 'near')
Strata_Sin  <- reproject(Strata , CRS = '+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs', program = 'GDAL', method = 'near')
Gal.map_Sin  <- reproject(Gal, CRS = '+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs', program = 'GDAL', method = 'near')
IIASA.map_Sin <- reproject(IIASA, CRS = '+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs', program = 'GDAL', method = 'near')
Thur.map_Sin <- reproject(Thur, CRS = '+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs', program = 'GDAL', method = 'near')


#fused_NL_c <- crop(Thur.map_Sin  , CountryShape_Sin1)
#fused_NL_c2 <- mask(fused_NL_c   , CountryShape_Sin1)
#pix.area <- res(Thur.map_Sin)[1]

# sum cell values
#cellStats(fused_NL_c2, 'sum', na.rm = T)

# copy
all.map <- stack(Thur.map_Sin, IIASA.map_Sin, Gal.map_Sin ,Fused.final_Sin, Strata_Sin)
all.map <- crop(all.map, CountryShape_Sin1)
all.map <- mask(all.map, CountryShape_Sin1)
all.dat <- as.data.frame(as.matrix(all.map, na.rm = TRUE))
#all.dat <- all.dat[complete.cases(all.dat),]                      ## here only data within the Baccini Extent is used
colnames(all.dat) <- c("Thur","IIASA","Gal","Fused" ,"Strata")
all.dat$Strata <- as.factor(all.dat$Strata)
agb.tot.s <- aggregate(all.dat$Fused, by = list(all.dat$Strata), FUN="sum", na.rm=TRUE)
length(names(all.map))-2


all.map <- stack(Thur.map_Sin, IIASA.map_Sin, Gal.map_Sin ,Fused.final_Sin, Strata_Sin)
n.rot <- length(names(all.map))-1
output <- matrix(data = NA, nrow =length(CountryShape_Sin), ncol = n.rot+1 )
names_col <- c("Thur","IIASA","Gal","Fused" ,"Strata")
colnames(output) <- c("Country",names_col[-length(names(all.map))])

#CountryShape_Sin <- CountryShape_Sin[-3,]
#

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


#t <- extract(all.map$bmAg_IIASA2010_Alligned ,CountryShape_Sin1 , na.rm = T, fun = sum, df = T)


