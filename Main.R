# Libaries

if (!require(rgdal)) install.packages('rgdal')
if (!require(rgeos)) install.packages('rgeos')
if (!require(raster)) install.packages('raster')
if (!require(gfcanalysis)) install.packages('gfcanalysis')
if (!require(SDMTools)) install.packages('SDMTools')
if (!require(devtools)) install.packages("devtools")
if (!require(MODIS)) install.packages("MODIS", repos="http://R-Forge.R-project.org")
if (!require(RCurl)) install.packages('RCurl')
#devtools::install_github('JornDallinga/VCF') # install this package atleast once, after that you wont have to run it everytime
if (!require(VCF)) install.packages('VCF')
if (!require(plotKML)) install.packages('plotKML')
if (!require(gdalUtils)) install.packages('gdalUtils')
if (!require(foreach)) install.packages('foreach')
if (!require(doParallel)) install.packages('doParallel')
if (!require(parallel)) install.packages('parallel')
if (!require(snow)) install.packages('snow')
if (!require(tcltk)) install.packages('tcltk')


# load functions
source("R/hansen.R")
source("R/SDM_function.R")
source("R/listing_files.R")
source("R/Unpack_VCF.R")
source("R/VCF.R")
source("R/Mosaic_Raster.R")

rasterOptions(tmpdir = "S:/R_Projects/temp_R")


# Read biomass points
shape <- readOGR(dsn = "./Ref_datasets/", layer = "Netherlands_NFI")
# Transform coordinate system
BufferWGS <- spTransform(shape, CRS("+proj=longlat +datum=WGS84"))
# country of analysis
ccodes()
Country <- 'NLD' # Netherlands
Country <- 'DEU'

# Read Raster
getwd()
list.files('./Maps/Gallaun/1km')
ras <- raster('./Maps/Gallaun/1km/bmAg_JR2000_ll_1km_eur_Crop.tif')

# -------------------------------------------------------------------------------------------------------------- #

# set year threshold
year <- 2000

# Set year values based on threshold to NA
BufferWGS[BufferWGS$YEAR < year ,] <- NA

# delete years with NA values
BufferWGS <- BufferWGS[!is.na(BufferWGS$YEAR),]


# -------------------------------------------------------------------------------------------------------------- #


# Exctract cell values, but keep it as a dataframe
o <- extract(ras, BufferWGS, method = 'simple', cellnumbers = T, sp = T) 
head(o)

# Unique
duplicated(o$cells)

# subset unique
n_occur <- data.frame(table(o$cells))
n_occur[n_occur$Freq > 1,]
df <- o[o$cells %in% n_occur$Var1[n_occur$Freq > 1],]


# replace spatial points AGB values with mean AGB values
lis <- unique(df$cells)
for (i in 1:length(lis)){
  df_sel <- df[df$cells == lis[i], ]
  df_sel$AGB <- mean(df_sel$AGB) # calc mean value of AGB values
  df_mean <- df_sel[1,] # select 1 spatial point out of multiple. Does not matter which one, since they all overlap the same pixel. Here we select the first
  df_2 <- df[!df$cells %in% df_sel$cells,]
  df <- rbind(df_2, df_mean)
}

# delete duplicates from dataframe
df_del <- o[!o$cells %in% df$cells,]

# merge new mean values with orginal dataframe
df_final <- rbind(df_del, df)

# check if any duplicates are remaining
n_occur <- data.frame(table(df_final$cells))
n_occur[n_occur$Freq > 1,]

# -------------------------------------------------------------------------------------------------------------- #

# crop raster to spatialpointdataframe, improve computation time
c <- crop(ras, df_final)

# mask to only select pixels/cells that overlap with spatial points
m <- mask(c, df_final)

# raster cells to spatial polygons
rastopol <- rasterToPolygons(m)

# Union features
union_points <- gUnaryUnion(rastopol)

# write to output file for verifying in other software (e.g. ArcGIS)
writeOGR(union_points, dsn = paste0(getwd(), '/data'), layer = "rastopolygon", driver = "ESRI Shapefile")

# unlink existing file (if it exists). Known bugs arise if you simply overwrite the .rds file.
unlink("data/BufferWGS.rds", recursive = FALSE)
# write to RDS file
saveRDS(union_points, file = "data/BufferWGS.rds", ascii = FALSE, version = NULL,
        compress = TRUE, refhook = NULL)

# -------------------------------------------------------------------------------------------------------------- #

# set parameters
input <- paste0(getwd(),"/data/BufferWGS.rds")
download_loc <- paste0(getwd(),"/data/Hansen_download/")
output <- paste0(getwd(), '/data/Hansen_output/')
Threshold <- 30
Year <- 2000


# downloading tree cover data from Hansen (Global Forest Watch)
# can take some time!!
# skip for large datasets, see Ref_data.R
Hansen(input, Threshold, Year, download_loc, output)


# lets test a method that select plots per tile
list_f <- list.files("./data/Hansen_download/", pattern = "treecover2000", full.names = T, recursive = T)

#------------------------------------------------------------------------------------------------------
#-------- Applies threshold and SDM function for plots that intersect with a tile ----------
#------------------------------------------------------------------------------------------------------

df <- rastopol[0,]
for (i in 1:length(list_f)){
  lr <- raster(list_f[i])
  c <- crop(rastopol, lr)
  
  for (b in 1:length(c)){
    pb <- winProgressBar(title = "progress bar", min = 0,
                         max = length(c), width = 300)
    Sys.sleep(0.1)
    setWinProgressBar(pb, b, title=paste(round(b/length(c)*100, 0),
                                         "% done"))
    
    pol <- crop(lr,c[b,])
    pol[pol < Threshold] <- 0 
    pol[pol >= Threshold] <- 1
    SDM <- SDM_function(pol)
    
    if (SDM$lanscape.division.index > .15 | is.na(SDM$lanscape.division.index)){
      c$bmAg_JR2000_ll_1km_eur[b] <- NA
    }
    if (b == length(c))
      close(pb)
  }
  cat(" | Tile ", i , " Out of ", length(list_f) , " Done ")
  df <- rbind(df, c)
}
# delete features with NA values
pol_sel <- df[!is.na(df$bmAg_JR2000_ll_1km_eur),]

# safe output 
writeOGR(pol_sel, dsn = paste0(getwd(), '/data/Output'), layer = "pol_sel_EU", driver = "ESRI Shapefile")

#------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------

### old
# plot
plot(H[[1]])
plot(input, add = T)

# set features to NA if below/above a certain threshold (here landscape.division.index > .15)
for (i in 1:length(rastopol)){
  # create progress bar
  pb <- winProgressBar(title = "progress bar", min = 0,
                       max = length(rastopol), width = 300)
  Sys.sleep(0.1)
  setWinProgressBar(pb, i, title=paste(round(i/length(rastopol)*100, 0),
                                       "% done"))
                    
  c <- crop(H[[1]],rastopol[i,])
  SDM <- SDM_function(c)
  if (SDM$lanscape.division.index > .15)
    rastopol$bmAg_JR2000_ll_1km_eur[i] <- NA
  if (i == length(rastopol))
    close(pb)
}

# delete features with NA values
pol_sel <- rastopol[!is.na(rastopol$bmAg_JR2000_ll_1km_eur),]

# safe output 
writeOGR(pol_sel, dsn = paste0(getwd(), '/data/Output'), layer = "pol_sel_NL", driver = "ESRI Shapefile")


# -------------------------------------------------------------------------------------------------------------- #

# Read biomass polygons
ref_pol <- readOGR(dsn = "./data/Output", layer = "pol_sel_NL")

# select intersecting points from rasterized polygons
ref_data <- crop(df_final, ref_pol)

# -------------------------------------------------------------------------------------------------------------- #

# Rasterize points 
# The whole upper section of this script can be simplified by running the rasterize function,
# instead of keep working with spatialpointsdataframes
ref_ras <- rasterize(ref_data, ras, ref_data$bmAg_JR2000_ll_1km_eur_Crop, fun=mean) 
writeRaster(ref_ras, filename = './Maps/Ref/ref_ras.tif') 

# -------------------------------------------------------------------------------------------------------------- #
# Download country boundaries
Country <- 'DEU' # run ccodes() to check for country codes
dir.create("./data/")
dir.create("./data/Boundaries")
CountryShape <- getData('GADM', country = Country, level=0, path = './data/Boundaries')
coordsys <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")
spTransform(CountryShape, coordsys)

# VCF
# Set parameters VCF
dataFolder <- './data/VCF/'
Year <- 2005
CountryShape <- CountryShape

# Downloading VCF data
VCF(dataFolder, Year, CountryShape = CountryShape, mosaic = F) # Mosaic will take a long time using R
# Best to load the downloaded data in Arcmap and transform/mosaic from there

# Reprojection
#library(plotKML)
library(gdalUtils)

# Reproject all the VCF files selected in your list
ls <- list.files("./data/VCF/extract_sexton/", pattern = ".tif", full.names = T)
ls2 <- list.files("./data/VCF/extract_sexton/", pattern = ".tif", full.names = F)
dir.create('./Covariates/Landsat_VCF_2005/transformed', showWarnings = F)
# Reprojection loop
# use gdalwarp from the gdalUtils package (significantly faster than the projectRaster::Raster package)
for (i in 1:length(ls)){
  lr <- raster(ls[i])
  VCF_WGS <- reproject(lr, CRS = '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0', program = 'GDAL', method = 'near')
  writeRaster(VCF_WGS, filename = paste0('./Covariates/Landsat_VCF_2005/transformed/', 'WGS_',ls2[i]))
}

#####################################################################################################
#####################################################################################################
#####################################################################################################

# Broken, does not work properly

# Testing MODIS VCF download
# download manually for now, use url below
dirs <- "ftp://ftp.glcf.umd.edu/glcf/Global_VCF/Collection_5/2005"


fls <- character()
library(RCurl)
getContent <- function(dirs) {
  urls <- paste(dirs, "/", sep="")
  fls <- strsplit(getURL(urls, dirlistonly=TRUE), "\r?\n")
  ok <- sapply(fls, length) > 0
  unlist(mapply(paste, urls[ok], fls[ok], sep="", SIMPLIFY=FALSE),
         use.names=FALSE)
}

while (length(dirs)) {
  message(length(dirs))
  urls <- getContent(dirs)
  isgz <- grepl("gz$", urls)
  fls <- append(fls, urls[isgz])
  dirs <- urls[!isgz]
}

urls <- getContent(dirs)
isgz <- grepl("gz$", urls)
fls <- append(fls, urls[isgz])

# Download and unpack MODIS VCF.
# Downloads ALL the tiles. no documentation found on tile numbers for the VCF product.
dir.create("./Covariates/MODIS_VCF_2005", showWarnings = F)
dir.create("./Covariates/", showWarnings = F)
url_ls <- fls

for (i in 1:length(fls)){
  url_s <- url_ls[i]
  downloadname <- strsplit(url_s, "/")
  downloadname <- unlist(downloadname)
  downloadname <- downloadname[length(downloadname)]
  file_name <- paste0("D:/R_Projects/EU_Biomass/Covariates/MODIS_VCF_2005/", downloadname)
  file_name2 <- substr(file_name, 1, nchar(file_name)-3)
  if (!file.exists(file_name2)){
    download.file(url = url_s, destfile = file_name)
    R.utils::gunzip(filename=file_name, remove=T, overwrite = F)
  }
}

#####################################################################################################
#####################################################################################################
#####################################################################################################


### Reprojecting MODIS (UTM) to MODIS (WGS-84)

dir.create('./Covariates/MODIS_VCF_2005/transformed/', showWarnings = F)
dir.create('./Covariates/MODIS_VCF_2005/transformed/WGS_MODIS/', showWarnings = F)

ls_files <- list.files("./Covariates/MODIS_VCF_2005/MODIS/", pattern = '\\.tif$', recursive = T,full.names = T)
ls2 <- list.files("./Covariates/MODIS_VCF_2005/MODIS/", pattern = "\\.tif$", full.names = F, recursive = T)
ls2 <- sub('.*/', '', ls2)
# use gdalwarp from the gdalUtils package (significantly faster than the projectRaster::Raster package)
for (i in 1:length(ls_files)){
  filename_M <- paste0('./Covariates/MODIS_VCF_2005/transformed/WGS_MODIS/', 'WGS_',ls2[i])
  if (!file.exists(filename_M)){
    print('starting reprojection')
    lr <- raster(ls_files[i])
    VCF_WGS <- reproject(lr, CRS = '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0', program = 'GDAL', method = 'near', overwrite= F)
    writeRaster(VCF_WGS, filename = paste0('./Covariates/MODIS_VCF_2005/transformed/WGS_MODIS/', 'WGS_',ls2[i]), overwrite = F)
  } else {
    print('file exists')
  }
}

### Resampling
ls_files <- list.files("./Covariates/MODIS_VCF_2005/transformed/WGS_MODIS/", pattern = '\\.tif$', recursive = T,full.names = T)
ls2 <- list.files("./Covariates/MODIS_VCF_2005/transformed/WGS_MODIS/", pattern = "\\.tif$", full.names = F, recursive = T)
ls2 <- sub('.*/', '', ls2)

# Set resolution
Gal <- raster("./Maps/Gallaun/1km/bmAg_JR2000_ll_1km_eur.tif")
xsize <- res(Gal)[1]
ysize <- res(Gal)[2]

# Use Gdal library (python)
dir.create('./Covariates/MODIS_VCF_2005/transformed/WGS_MODIS_resampled/', showWarnings = F)
# use gdal_translate to resample using the method 'mode' in R
for (i in 1:length(ls_files)){
  filename_S <- paste0('./Covariates/MODIS_VCF_2005/transformed/WGS_MODIS_resampled/', 'resample_',ls2[i])
  if (!file.exists(filename_S)){
    gdal_translate(src_dataset = ls_files[i], ot = "Int16" ,dst_dataset = paste0('./Covariates/MODIS_VCF_2005/transformed/WGS_MODIS_resampled/', 'resample_',ls2[i]),tr = c(xsize,ysize), r = "mode")
  }
}


# Reclassify MODIS 
dir.create('./Covariates/MODIS_VCF_2005/transformed/WGS_MODIS_resampled_Reclassified/', showWarnings = F)

ls_files <- list.files("./Covariates/MODIS_VCF_2005/transformed/WGS_MODIS_resampled/", pattern = '\\.tif$', recursive = T,full.names = T)
ls2 <- list.files("./Covariates/MODIS_VCF_2005/transformed/WGS_MODIS_resampled/", pattern = "\\.tif$", full.names = F, recursive = T)
ls2 <- sub('.*/', '', ls2)

# Reclassifying cells
for (i in 1:length(ls_files)){
  filename_S <- paste0('./Covariates/MODIS_VCF_2005/transformed/WGS_MODIS_resampled_Reclassified/', 'reclass_',ls2[i])
  if (!file.exists(filename_S)){
    lr <- raster(ls_files[i])
    lr[lr > 100] <- NA # remove water and other thematics from the MODIS data
    writeRaster(lr, filename = paste0('./Covariates/MODIS_VCF_2005/transformed/WGS_MODIS_resampled_Reclassified/', 'reclass_',ls2[i]), overwrite = F)
  }
}


### Mosaic using Gdal
# Mosaic all reprojected and resampled raster tiles
dir.create('./Covariates/MODIS_VCF_2005/transformed/Mosaic/', showWarnings = F)
ls_files_m <- list.files("./Covariates/MODIS_VCF_2005/transformed/WGS_MODIS_resampled_Reclassified/", pattern = '\\.tif$', recursive = T,full.names = T)
mosaic_rasters(gdalfile = ls_files_m, dst_dataset = "./Covariates/MODIS_VCF_2005/transformed/Mosaic/Mosaic.tif", overwrite = T, ot = "Int16")

### plotting
t <- raster("./Covariates/MODIS_VCF_2005/transformed/Mosaic/Mosaic.tif")
t
plot(t)

### cropping
t <- raster("./Covariates/MODIS_VCF_2005/transformed/Mosaic/Mosaic.tif")
crop(t,Gal, filename = "./Covariates/MODIS_VCF_2005/transformed/Mosaic/MODIS_VCF_Mosaic.tif", overwrite = T)


### plotting

t <- raster("./Covariates/MODIS_VCF_2005/transformed/Mosaic/MODIS_VCF_Mosaic.tif")
plot(t)



#-----------------CCI---------------------------------------------------------------------

cci <- raster("./ESACCI-LC-L4-LCCS-Map-300m-P5Y-2005-v1.6.1.tif/ESACCI-LC-L4-LCCS-Map-300m-P5Y-2005-v1.6.1.tif")
ex <- extent(Gal)
cci <- crop(cci2, Gal)
writeRaster(cci, filename = './Covariates/CCI_2005/clip_cci_EU.tif')

Gal <- raster("./Maps/Gallaun/1km/bmAg_JR2000_ll_1km_eur.tif")

ls_files <- list.files('./Covariates/CCI_2005/', pattern = 'clip_cci_EU.tif', full.names = T)
xsize <- res(Gal)[1]
ysize <- res(Gal)[2]
gdal_translate(src_dataset = ls_files[1], ot = "Int16" ,dst_dataset = './Covariates/CCI_2005/CCI_2005_resample.tif',tr = c(xsize,ysize), r = "mode")
cci <- raster('./Covariates/CCI_2005/CCI_2005_resample.tif')

writeRaster(cci, filename = './Covariates/CCI_2005/CCI_2005_resample_EU.tif', overwrite = F)

#-----------------Water---------------------------------------------------------------------

water <- raster("./Covariates/CCI_Water/ESACCI-LC-L4-WB-Map-150m-P13Y-2000-v4.0.tif")

ex <- extent(Gal)
water <- crop(water, Gal)
writeRaster(water, filename = './Covariates/CCI_Water/CCI_Water_EU.tif', overwrite = T)

Gal <- raster("./Maps/Gallaun/1km/bmAg_JR2000_ll_1km_eur.tif")

ls_files <- list.files('./Covariates/CCI_Water/', pattern = 'CCI_Water_EU.tif', full.names = T)
xsize <- res(Gal)[1]
ysize <- res(Gal)[2]
gdal_translate(src_dataset = ls_files[1], ot = "Int16" ,dst_dataset = './Covariates/CCI_Water/CCI_Water_resample.tif',tr = c(xsize,ysize), r = "mode")
water <- raster('./Covariates/CCI_Water/CCI_Water_resample.tif')

writeRaster(water, filename = "./Covariates/CCI_Water/CCI_Water_resample_EU.tif", overwrite = T)
water <- raster("./Covariates/CCI_Water/CCI_Water_resample_EU.tif")

#-----------------Height---------------------------------------------------------------------

height <- raster("./Covariates/Height/Height_large_extent.tif")
Gal <- raster("./Maps/Gallaun/1km/bmAg_JR2000_ll_1km_eur.tif")

ex <- extent(Gal)
height <- crop(height, Gal)
writeRaster(height, filename = './Covariates/Height/Height_EU.tif', overwrite = T)



ls_files <- list.files('./Covariates/Height/', pattern = 'Height_EU.tif', full.names = T)
xsize <- res(Gal)[1]
ysize <- res(Gal)[2]
gdal_translate(src_dataset = ls_files[1], ot = "Int16" ,dst_dataset = './Covariates/Height/Height_resample.tif',tr = c(xsize,ysize), r = "mode")
height <- raster('./Covariates/Height/Height_resample.tif')

writeRaster(height, filename = "./Covariates/Height/Height_resample_EU.tif", overwrite = T)
height <- raster("./Covariates/Height/Height_resample_EU.tif")


