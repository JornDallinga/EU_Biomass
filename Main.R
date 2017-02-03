# Libaries

if (!require(rgdal)) install.packages('rgdal')
if (!require(rgeos)) install.packages('rgeos')
if (!require(raster)) install.packages('raster')
if (!require(gfcanalysis)) install.packages('gfcanalysis')
if (!require(SDMTools)) install.packages('SDMTools')
if (!require(devtools)) install.packages("devtools")
if (!require(MODIS)) install.packages("MODIS", repos="http://R-Forge.R-project.org")
if (!require(RCurl)) install.packages('RCurl')
#devtools::install_github('JornDallinga/VCF')
if (!require(VCF)) install.packages('VCF')
if (!require(plotKML)) install.packages('plotKML')
if (!require(gdalUtils)) install.packages('gdalUtils')

# load functions
source("R/hansen.R")
source("R/SDM_function.R")
source("R/listing_files.R")
source("R/Unpack_VCF.R")
source("R/VCF.R")
source("R/Mosaic_Raster.R")


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
H <- Hansen(input, Threshold, Year, download_loc, output)

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



