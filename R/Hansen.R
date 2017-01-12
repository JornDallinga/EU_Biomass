Hansen <- function(input, Threshold, year = Year, download_loc ,output){
  ## Create variable Area Of Interest (aio)
  aoi <- readRDS(file = input, refhook = NULL)
  
  ## Calculate tiles needed to cover the AOI
  tiles <- calc_gfc_tiles(aoi)
  print(length(tiles))
  
  ## Download GFC data
  download_tiles(tiles, download_loc)
  
  cat("---- extraction ----")
  ## Extract data from tiles
  gfc_extract <- extract_gfc(aoi, download_loc, filename= paste0(output,"GFC_extract.tif"), overwrite=TRUE, data_year = 2015, progress = 'text')
  cat("---- extraction done ----")
  
  cat("---- applying threshold ----")
  ## Apply threshold to extracted data 
  gfc_thresholded <- threshold_gfc(gfc_extract, Threshold=Threshold, 
                                   filename= paste0(output,"GFC_extract_thresholded.tif"), overwrite=TRUE, progress = 'text')
  cat("---- applying threshold done ----")
  ## Masking for water perc calculations
  #mask_water <- mask(gfc_thresholded, aoi)
  names(gfc_thresholded) <- c('forest2000','lossyear','gain','lossgain','datamask')
  #brick(paste0(output,"GFC_extract_thresholded.tif"))
  
  ## Masking gfc data to aoi
  # Masking appears to be slower than just loading the entire aoi. Might require further testing.
  #mask_gfc <- mask(gfc_thresholded, aoi)
  mask_gfc <- gfc_thresholded
  
  ## anual stack of years (old, takes too long)
  #annual_Hansen <- annual_stack(mask_gfc, data_year = 2015)
  
  # Select a year from the annual Hansen dataset
  #select_y <- sprintf("y%s",year)
  #subset_Year <- subset(annual_Hansen, subset = select_y)
  
  #-----
  cat("---- extraction ----")
  
  #-----
  
  # Create binary forest cover map and figure output
  
  if (year == 2000){
    ## create forest cover mask year 2000
    Figure_output <- mask_gfc
    
    #Figure_output$forest2000[Figure_output$forest2000 == 1] <- 2 # Forest
    #Figure_output$forest2000[Figure_output$forest2000 == 0] <- 1 # Non-Forest
    Figure_output$forest2000[Figure_output$datamask == 2] <- 2 # water
    #Figure_output$datamask[Figure_output$datamask == 0] <- 4 # No data
    #suppressWarnings(Figure_output$datamask[Figure_output$datamask < 3] <- NA) # Nodata for merging
    #Figure_output <- merge(Figure_output$datamask, Figure_output$forest2000, overlap = T)
    Figure_output <- Figure_output$forest2000
    names(Figure_output) <- "Hansen"
    
    mask_gfc <- mask_gfc$forest2000
    
  } else if (year >= 2001 & year <= 2014){
    
    # faster than annual_stack::gfc_analysis
    n <- substr(year, 3, 4)
    n <- gsub("(?<![0-9])0+", "", n, perl = TRUE)
    n <- as.numeric(n)
    # forest loss
    mask_gfc$lossyear[mask_gfc$lossyear <= n & mask_gfc$lossyear != 0] <- 1 # loss
    mask_gfc$lossyear[mask_gfc$lossyear > n] <- 0 # no loss
    mask_gfc$forest2000[mask_gfc$lossyear == 1] <- 0 # remove loss from forest2000
    # forest gain
    mask_gfc$gain[mask_gfc$lossgain == 1] <- 1 # 1 = gain, 0 = empty
    mask_gfc$forest2000[mask_gfc$gain == 1] <- 1 # add gain to forest2000
    
    ## create forest cover mask of years 2001 till 2014
    Figure_output <- mask_gfc$forest2000
    
    #Figure_output[Figure_output == 2] <- 8 # Non-forest
    #Figure_output[Figure_output == 1] <- 2 # Forest
    #Figure_output[Figure_output == 8] <- 1 # Non-forest
    
    #Figure_output[Figure_output == 3] <- 1 # Non-forest
    #Figure_output[Figure_output == 4 ] <- 2 # Forest gain to Forest
    #Figure_output[Figure_output == 5 ] <- 2 # Forest loss and gain to Forest
    Figure_output[mask_gfc$datamask == 2] <- 2 # Water
    names(Figure_output) <- "Hansen"
    
    #subset_Year[subset_Year == 2 ] <- 0
    #subset_Year[subset_Year == 3 ] <- 0
    #subset_Year[subset_Year == 4 ] <- 1
    #subset_Year[subset_Year == 5 ] <- 1
    #subset_Year[subset_Year == 6 ] <- 0
    mask_gfc <- mask_gfc$forest2000
  } else {
    warning("invalid year")
  }
  
  names(mask_gfc) <- "Hansen"
  
  #  <- list(mask_gfc, Water_perc, Figure_output)
  return_list <- list(mask_gfc, Figure_output)
  
  return (return_list)
}