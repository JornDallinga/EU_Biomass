ras <- raster('./Maps/Gallaun/1km/bmAg_JR2000_ll_1km_eur.tif')
ras1 <- raster('./Maps/Barredo/Barredo_wgs1.tif')
ras2 <- raster('./Maps/IIASA/1km/bmAg_IIASA2010.tif')
ras3 <- raster('./Maps/Thurner/1km/bmAg_Thurner_1km.tif')

# Extent of raster
e <- drawExtent()
#e <- extent(ras2)
p <- as(e, 'SpatialPolygons') 
#pp <- rasterToPolygons(ras2, dissolve=TRUE)
BufferWGS <- spTransform(p, CRS("+proj=longlat +datum=WGS84"))
p@proj4string <- CRS("+proj=longlat +datum=WGS84")
