#### Reference data preperation
# Make all reference data compatible

if (!require(velox)) install.packages('velox')
if (!require(SpaDES)) install.packages('SpaDES')


# -------------------------------------------------------------------------------------------------------------- #


# Read all biomass points
shape <- readOGR(dsn = "./Corrected_data/NL_AGB", layer = "NL_AGB")
shape1 <- readOGR(dsn = "./Ref_datasets/Croatia", layer = "Croatia_NFI")
shape2 <- readOGR(dsn = "./Ref_datasets/Germany", layer = "Germany_NFI")
shape3 <- readOGR(dsn = "./Ref_datasets/France", layer = "France_NFI_10km")
shape4 <- readOGR(dsn = "./Ref_datasets/Italy", layer = "Italy_NFI_1km")
shape5 <- readOGR(dsn = "./Corrected_data/Spain_AGB", layer = "Spain_AGB")
shape6 <- readOGR(dsn = "./Ref_datasets/Sweden/Sweden_NFI_Plots/Sweden_NFI_Plots", layer = "Sweden_NFI_2007_2016")

# Transform coordinate system
#shape <- spTransform(shape, CRS("+proj=longlat +datum=WGS84")) # if needed

# select AGB columns
shape <- shape[,c("AGB")]
shape1 <- shape1[,c("AGB")]
shape2 <- shape2[,c("AGB")]
shape3 <- shape3[,c("SumOfAGB")]
shape4 <- shape4[,c("W4apv_ha")]
shape5 <- shape5[,c("AGB_T_HA")]
shape6 <- shape6[,c("AGB")]

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
shape6$country_ID <- 7

# merge the dataframes together
dat <- rbind(shape, shape1, shape2, shape3, shape4, shape5, shape6)
rm(list=ls(pattern='shape'))

# Create country code raster
# -------------------------------------------------------------------------------------------------------------- #


# rasterize country ID's to create a raster grid with country codes
ras <- raster('./Maps/Gallaun/1km/bmAg_JR2000_ll_1km_eur.tif')
#ras <- raster('./Maps/Barredo/barredo_Alligned.tif')

## create NA raster to fill with values later
ras_NA <- ras
ras_NA[ras] <- NA
ref_ras <- rasterize(dat, ras_NA, dat$country_ID, filename = './Reference/Ref_code_EU.tif', overwrite = T) 

## Each NFI has a unique code, and there is a raster map where each plot of the NFI has the code value (e.g.: all plots in Spain have value 1, in France is 2, etc.)
codes <- raster('./Reference/Ref_code_EU.tif')
code.names <- as.character(read.csv(paste("./Reference/Codes_EU.csv", sep=""))[,2])
code.n <- maxValue(codes)

# -------------------------------------------------------------------------------------------------------------- #



# Read biomass points
#BufferWGS <- dat
# set plot year threshold
year <- 2000

# Set year values based on threshold to NA
#BufferWGS[BufferWGS$YEAR < year ,] <- NA

# delete years with NA values
#BufferWGS <- BufferWGS[!is.na(BufferWGS$YEAR),]

# country/extent of analysis
#Country  <- readOGR(dsn = "./data/Boundaries/Europe", layer = "Europe_in_Gal")

# -------------------------------------------------------------------------------------------------------------- #


# Exctract cell values, but keep it as a dataframe
#df <- extract(ras, BufferWGS, method = 'simple', cellnumbers = T, sp = T) 
#head(df)


# -------------------------------------------------------------------------------------------------------------- #
#df_final <- df

# crop raster to spatialpointdataframe, improve computation time
#df_final <- crop(df_final, ras)

ref_bio <- rasterize(dat, ras_NA, dat$AGB, fun = mean, filename = './Reference/Ref_bio_EU.tif', overwrite = T) 

# raster cells to spatial polygons
rastopol <- rasterToPolygons(ref_bio)

# mask to only select pixels/cells that overlap with spatial points
# m <- mask(c, df_final)
#
# # Union features
#union_points <- gUnaryUnion(rastopol)
#
# # write to output file for verifying in other software (e.g. ArcGIS)
#writeOGR(union_points, dsn = paste0(getwd(), '/data'), layer = "rastopolygon", driver = "ESRI Shapefile")
#
# # unlink existing file (if it exists). Known bugs arise if you simply overwrite the .rds file.
#unlink("data/BufferWGS.rds", recursive = FALSE)
# # write to RDS file
#saveRDS(union_points, file = "data/BufferWGS.rds", ascii = FALSE, version = NULL,
#        compress = TRUE, refhook = NULL)

# -------------------------------------------------------------------------------------------------------------- #

# set parameters
input <- paste0(getwd(),"/data/BufferWGS.rds")
download_loc <- paste0(getwd(),"/data/Hansen_download/")
output <- paste0(getwd(), '/data/Hansen_output/')
Threshold <- 30
Year <- 2000


# downloading tree cover data from Hansen (Global Forest Watch)
# interupt this function after downloading!
Hansen(input, Threshold, Year, download_loc, output)

#------------------------------------------------------------------------------------------------------

# list Hansen files that you wish to process
list_f <- list.files("./data/Hansen_download/", pattern = "treecover2000", full.names = T, recursive = T)


#------------- testing large scale extraction. Fastest method so far.------------------------------#

# lets test a method that select plots per tile
list_f <- list.files("./data/Hansen_download/", pattern = "treecover2000", full.names = T, recursive = T)

# create dataframe to write final results
df3 <- rastopol[0,]
# set cores
detectCores()
beginCluster( detectCores() -2) # use all but one core

# following function uses vortex and SpaDES, improved raster computing performance. 
system.time(for (i in 1:length(list_f)){
  
  pb <- winProgressBar(title = "progress bar", min = 0,
                       max = length(list_f), width = 300)
  Sys.sleep(0.1)
  setWinProgressBar(pb, i, title=paste(round(i/length(list_f)*100, 0),
                                       "% done"))
  
  lr <- raster(list_f[i])
  
  print(paste0(" | Tile ", i , " Out of ", length(list_f) , ": cropping"))
  
  # Split into sections - saves automatically to path
  sections = splitRaster(lr, nx=2, ny=2, path="./data/temp_data/")
  
  
  count <- 1
  for (b in sections){
    c <- crop(rastopol, b, progress = T)
    vx <- velox(b)
    if (!is.null(c)){
      lr_ex_velox <- vx$extract(sp = c, fun = NULL)
      lr_m <- lapply(lr_ex_velox ,mean)
      lr_sd  <- lapply(lr_ex_velox ,sd)
      c$Mean <- unlist(lr_m)
      c$SD <- unlist(lr_sd)
      df3 <- rbind(df3, c)
      print(paste0(" | section ", count , " Out of ", length(sections) , ": extracted"))
      count <- count + 1
    } else {
      print(paste0(" | section ", count , " Out of ", length(sections) , ": empty"))
      count <- count + 1
    }

  }
  # ----------------------- OLD SECTION, WORKS BUT IS MUCH SLOWER
  # c <- crop(rastopol, lr, progress = T)
  # 
  # print(paste0(" | Tile ", i , " Out of ", length(list_f) , ": extracting raster values of ", length(c), " features "))
  # 
  # lr_ex <- extract(x = lr,y = c, na.rm = T, progress = 'text')
  # 
  # print(paste0(" | Tile ", i , " Out of ", length(list_f) , ": calculating mean and sd "))
  # 
  # lr_m <- lapply(lr_ex ,mean)
  # lr_sd  <- lapply(lr_ex ,sd)
  # 
  # c$Mean <- unlist(lr_m)
  # c$SD <- unlist(lr_sd)
  
  if (i == length(list_f)) {
    close(pb)
    endCluster()
    print("Writing output to disk...")
    plot.sel3 <- rbind(df3[df3$SD <= 15 & df3$Mean >= 50 & df3$Ref_bio_EU > (df3$Mean / 2), ], df3[df3$SD <= 15 & df3$Mean < 50 & df3$Ref_bio_EU < (2 * df3$Mean), ])
    writeOGR(plot.sel3, dsn = paste0(getwd(), '/data/Output'), layer = "pol_sel_EU_final", driver = "ESRI Shapefile")
    
  }
  print(paste0(" | Tile ", i , " Out of ", length(list_f) , " Done "))
  
  # df3 <- rbind(df3, c)
  do.call(file.remove, list(list.files("./data/temp_data/", full.names = TRUE)))
})


# -------------------------------------------------------------------------------------------------------------- #

# Read biomass polygons
ref_pol <- readOGR(dsn = "./data/Output", layer = "pol_sel_EU_final")

# Rasterize points 
rasterize(ref_pol, ras, ref_pol$Ref_bio_EU, fun=mean, filename = './Maps/Ref/ref_ras_EU_final.tif', progress = 'text', overwrite = T) 

#rm(ref_pol, df, df_2, df_del, df_final, df_mean, df_sel, i, lis, n_occur)

# -------------------------------------------------------------------------------------------------------------- #
raster_ref <- raster('./Maps/Ref/ref_ras_EU_final.tif')
writeRaster(round(raster_ref), filename = './Maps/Ref/ref_ras_EU_final_round.tif', overwrite = T)
plot(raster_ref)

