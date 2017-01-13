# Load CCI
lc <- list.files('./Covariates/CCI_2005/',pattern='.tif$', full.names = T)
CCI <- raster(lc[1])

# Load Landsat VCF
lc2 <- list.files('C:/R_Projects/Biomass_Europe/Covariates/Landsat_VCF_2005/',pattern='Mosaic_TC_2005.tif$', full.names = T, include.dirs = T)
VCF_2005 <- raster(lc2)

# Load Biomass map
ras <- raster('./Maps/Gallaun/1km/bmAg_JR2000_ll_1km_eur.tif')

# -------------------------------------------------------------------------------------------------------------- #

# temp crop to NL
CCI <- crop(CCI, VCF_2005)

# check frequency of raster values
freq(CCI)

# copy dataset
CCI_new <- CCI

# reclassify CCI
# Evergreen forest, Deciduous forest, Woodland, Mosaic Vegetation, Shrubland, Grassland, Cropland, Other land

CCI_new[CCI == 60] <- 1 # Evergreen
CCI_new[CCI == 70] <- 2# Deciduous forest
CCI_new[CCI == 90] <- 3 # Woodland
CCI_new[CCI == 40 | CCI == 150 | CCI == 110 | CCI == 100] <- 4 # Mosaic vegetation
CCI_new[CCI >= 120 & CCI <= 122 | CCI == 152 | CCI == 180] <- 5 # Shrubland
CCI_new[CCI == 130] <- 6 # Grassland
CCI_new[CCI >= 10 & CCI <= 30] <- 7 # Cropland
CCI_new[CCI == 190 | CCI >= 200 & CCI <= 202] <- 8 # Other land
CCI_new[CCI == 210] <- NA # Water

# check frequency of raster values
freq(CCI_new)

# -------------------------------------------------------------------------------------------------------------- #

# Resample CCI
aggr <- aggregate(CCI_new, 2, fun = mode, na.rm = T)
CCI_resample <- resample(aggr, ras, method = 'ngb')

# -------------------------------------------------------------------------------------------------------------- #

# crop to area
CCI_crop <- crop(CCI_resample, CCI)
ras_crop <- crop(ras, CCI)

# write rasters
writeRaster(CCI_crop, filename = 'Covariates/CCI_2005/CCI_crop.tif')
writeRaster(ras_crop, filename = 'C:/R_Projects/Biomass_Europe/Maps/Gallaun/1km/bmAg_JR2000_ll_1km_eur_Crop.tif')
CCI <- raster("Covariates/CCI_2005/CCI_crop.tif")
# -------------------------------------------------------------------------------------------------------------- #
# -------------------------------------------------------------------------------------------------------------- #
# -------------------------------------------------------------------------------------------------------------- #

# LOAD VCF

VCF <- raster("./Covariates/CCI_Water/VCF_aggr_6.tif")
plot(VCF)

# check frequency of raster values
freq(VCF)

# reclassify VCF
VCF_new <- VCF
VCF_new[VCF1 >= 0 & VCF1 <= 12.5] <- 1
VCF_new[VCF >= 12.6 & VCF <= 25] <- 2
VCF_new[VCF >= 25.1 & VCF <= 37.5] <- 3
VCF_new[VCF >= 37.6 & VCF <= 50] <- 4
VCF_new[VCF >= 50.1 & VCF <= 62.5] <- 5
VCF_new[VCF >= 62.6 & VCF <= 75] <- 6
VCF_new[VCF >= 75.1 & VCF <= 87.5] <- 7
VCF_new[VCF >= 87.6 & VCF <= 100] <- 8

# Resample VCF
aggr <- aggregate(VCF_new, 2, fun = mode, na.rm = TRUE)
VCF_resample <- resample(aggr, ras, method = 'ngb')

# -------------------------------------------------------------------------------------------------------------- #

# crop to area
ras_crop <- raster("./Maps/Gallaun/1km/bmAg_JR2000_ll_1km_eur_Crop.tif")
VCF_crop <- crop(VCF_resample, ras_crop)

# write rasters
writeRaster(VCF_crop, filename = 'Covariates/Landsat_VCF_2005/VCF_crop.tif', overwrite = T)

# -------------------------------------------------------------------------------------------------------------- #
# -------------------------------------------------------------------------------------------------------------- #
# -------------------------------------------------------------------------------------------------------------- #

bac_resample <- resample(bac, saa, method = 'ngb')
