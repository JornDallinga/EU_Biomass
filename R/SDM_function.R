SDM_function <- function(S){
  # create asc file
  asc.from.raster(S)

  #### Creating statistics from forest cover map
  c1.data <- ClassStat(S, cellsize = 30, latlon = T)
  
  #### Set all values to 0 if the area has no forest cover
  if (nrow(c1.data) == 1 & (c1.data)[1,1] == 0) {
    c1.data <- c1.data[NA,]
  } else if (nrow(c1.data) > 1 & (c1.data)[1,1] == 0) {
    c1.data <- c1.data[-1,]
  } else {   
  }

  return (c1.data)
}