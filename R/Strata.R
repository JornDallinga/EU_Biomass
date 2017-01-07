# Load CCI
lc <- list.files('C:/R_Projects/Biomass_Europe/Covariates/CCI_2005/',pattern='.tif$', full.names = T)
CCI <- raster(lc[1])

# Load Landsat VCF
lc2 <- list.files('C:/R_Projects/Biomass_Europe/Covariates/Landsat_VCF_2005/',pattern='Mosaic_TC_2005.tif$', full.names = T, include.dirs = T)
VCF_2005 <- raster(lc2)

# Load Biomass map
ras <- raster('C:/R_Projects/Biomass_Europe/Maps/Gallaun/1km/bmAg_JR2000_ll_1km_eur.tif')

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

