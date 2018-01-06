

############################### Load libraries ###################################################

require(raster)
require(dplyr)
if (!require(rgdal)) install.packages('rgdal')
if (!require(rgeos)) install.packages('rgeos')

############################### Read csv 2 shapefile ###################################################

mydf <- read.csv("./Corrected_data/AGB_DATA_EU2_OCT17_Spain.csv", header = T, sep = ",")
xy <- mydf[,c(2,3)]

spdf <- SpatialPointsDataFrame(coords = xy, data = mydf,
                               proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))

writeOGR(obj=spdf, dsn="./Corrected_data/Spain_AGB", layer="Spain_AGB", driver="ESRI Shapefile") # this is in geographical projection

###############################  replace values by ID  ###################################################

Mydata <- read.csv("./Corrected_data/AGB_DATA_EU3_OCT17_NL.csv.csv", header = T, sep = ",")
#raster::projection(Mydata) = "+init=espg:4326" # WGS84 coords
#shapefile(Mydata, "Mydata.shp")

# load shapefile
shape_NL <- readOGR(dsn = "./Ref_datasets/Netherlands", layer = "Netherlands_NFI")


# copy file
shape_new <- shape_nl

# replace values
shape_new$plot_id <- as.integer(shape_new$plot_id)
shape_new$AGB[match(Mydata$PLOT_ID, shape_new$plot_id)] <- Mydata$AGB_T_HA


###############################  check if values are replaced correctly  ###################################################


# Check if ID's match
all(Mydata$PLOT_ID %in% shape_new$plot_id)

# order/sort shapefile ID to compare
shape_new <- shape_new[order(shape_new$plot_id),]
Mydata <- arrange(Mydata,Mydata$PLOT_ID)

# round AGB values
shape_new$AGB <- round(shape_new$AGB)
Mydata$AGB_T_HA <- round(Mydata$AGB_T_HA)

# test if the shapefile is corrected with AGB values
identical(as.numeric(round(shape_new$AGB)),as.numeric(round(Mydata$AGB_T_HA)))

writeOGR(obj=shape_new, dsn="./Corrected_data/NL_AGB", layer="NL_AGB", driver="ESRI Shapefile") # this is in geographical projection

# write csv
#library(xlsx) #load the package
#write.xlsx(x = test, file = "./Corrected_data/test.xlsx",
#           sheetName = "TestSheet", row.names = FALSE)
#write.xlsx(x = Mydata, file = "./Corrected_data/Mydata.xlsx",
#           sheetName = "TestSheet", row.names = FALSE)

