#### Reference data preperation
# Make all reference data compatible

# -------------------------------------------------------------------------------------------------------------- #


# Read all biomass points
shape <- readOGR(dsn = "./Ref_datasets/Netherlands", layer = "Netherlands_NFI")
shape1 <- readOGR(dsn = "./Ref_datasets/Croatia", layer = "Croatia_NFI")
shape2 <- readOGR(dsn = "./Ref_datasets/Germany", layer = "Germany_NFI")
shape3 <- readOGR(dsn = "./Ref_datasets/France", layer = "France_NFI_10km")
shape4 <- readOGR(dsn = "./Ref_datasets/Italy", layer = "Italy_NFI_1km")
shape5 <- readOGR(dsn = "./Ref_datasets/Spain", layer = "Spain_NFI_2000_2008")

# Transform coordinate system
#shape <- spTransform(shape, CRS("+proj=longlat +datum=WGS84")) # if needed

# select AGB columns
shape <- shape[,c("AGB")]
shape1 <- shape1[,c("AGB")]
shape2 <- shape2[,c("AGB")]
shape3 <- shape3[,c("SumOfAGB")]
shape4 <- shape4[,c("W4apv_ha")]
shape5 <- shape5[,c("SumOfAGB")]

# rename AGB columns if needed
names(shape3) <- "AGB"
names(shape4) <- "AGB"
names(shape5) <- "AGB"

# add unique ID to country/continent etc
shape$country_ID <- 1
shape1$country_ID <- 2
shape2$country_ID <- 3
shape3$country_ID <- 4
shape4$country_ID <- 5
shape5$country_ID <- 6

# merge the dataframes together
dat <- rbind(shape, shape1, shape2, shape3, shape4, shape5)

# Create country code raster
# -------------------------------------------------------------------------------------------------------------- #


# rasterize country ID's to create a raster grid with country codes
ras <- raster('./Maps/Gallaun/1km/bmAg_JR2000_ll_1km_eur.tif')

## create NA raster to fill with values later
ras_NA <- ras
ras_NA[ras] <- NA
ref_ras <- rasterize(dat, ras_NA, dat$country_ID, filename = './Reference/Ref_code_EU.tif', overwrite = T) 

## Each NFI has a unique code, and there is a raster map where each plot of the NFI has the code value (e.g.: all plots in Spain have value 1, in France is 2, etc.)
codes <- raster(paste('./Reference/Ref_code_EU.tif', sep=""))
code.names <- as.character(read.csv(paste("./Reference/Codes_EU.csv", sep=""))[,2])
code.n <- maxValue(codes)

# -------------------------------------------------------------------------------------------------------------- #


# Read biomass points
BufferWGS <- dat
# set plot year threshold
year <- 2000

# Set year values based on threshold to NA
#BufferWGS[BufferWGS$YEAR < year ,] <- NA

# delete years with NA values
#BufferWGS <- BufferWGS[!is.na(BufferWGS$YEAR),]

# country/extent of analysis
Country  <- readOGR(dsn = "./data/Boundaries/Europe", layer = "Europe_in_Gal")

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
# can take some time!! Below a faster solution
H <- Hansen(input, Threshold, Year, download_loc, output)


#------------------------------------------------------------------------------------------------------

# list Hansen files that you wish to process
list_f <- list.files("./data/Hansen_download/", pattern = "treecover2000", full.names = T, recursive = T)
# Reclassify Hansen
dir.create('./data/Hansen_Reclass', showWarnings = F)

list_f <- list.files("./data/Hansen_download/", pattern = "treecover2000", full.names = T, recursive = T)
ls2 <- list.files("./data/Hansen_download/", pattern = "treecover2000", full.names = F, recursive = T)


# Works, but will take quite some GB and time
mosaic_rasters(list_f, dst_dataset = "./data/Hansen_download/mosaic.tif", overwrite = T, ot = "Int16")

# lets test something else
list_f
lr <- raster(list_f[7])
# Work in progress
#------------------------------------------------------------------------------------------------------
# tryCatch(!is.null(crop(r,extent(rect))), error=function(e) return(FALSE)) 

for (i in 1:length(rastopol)){
  pb <- winProgressBar(title = "progress bar", min = 0,
                       max = length(rastopol), width = 300)
  Sys.sleep(0.1)
  setWinProgressBar(pb, i, title=paste(round(i/length(rastopol)*100, 0),
                                       "% done"))
  for (b in 1:length(list_f)){
    lr <- raster(list_f[b])
    if (tryCatch(!is.null(crop(lr,extent(rastopol[i,]))), error=function(e) return(FALSE)) == T){
      c <- crop(lr,rastopol[i,])
      c[c < Threshold] <- 0 
      c[c >= Threshold] <- 1
      SDM <- SDM_function(c)
      if (SDM$lanscape.division.index > .15 | is.na(SDM$lanscape.division.index))
        rastopol$bmAg_JR2000_ll_1km_eur[i] <- NA
    }
  if (i == length(rastopol))
    close(pb)
  }
}



#------------------------------------------------------------------------------------------------------


# Reclassify
for (i in 1:length(list_f)){
  filename_S <- paste0('./data/Hansen_Reclass/', 'reclass_',ls2[i])
  if (!file.exists(filename_S)){
    lr <- raster(list_f[i])
    lr[lr >= Threshold] <- 1 
    lr[lr < Threshold] <- 0 
    writeRaster(lr, filename = paste0('./data/Hansen_Reclass/', 'reclass_',ls2[i]), overwrite = F)
  }
}




#------------------------------------------------------------------------------------------------------

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
s