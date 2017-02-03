VCF <- function(dataFolder, Year, CountryShape, mosaic = F){
  # Get buffer locations
  
  buffer <- CountryShape
  pr <- getPR(buffer)
  pr_copy <- pr
  dir <- dataFolder #"data/" # If doesnt work add "./"
  extract <- 'extract_sexton/'
  dir.create(paste0(dir,'extract_sexton'), showWarnings = F) 

  # check for Path/Row existance
  if (is.list(pr_copy)) { # Assuming the list provided is the variable returned by getPR() function
    pr_copy <- pr_copy$PR
  }
  
  pr_copy <- sprintf('%06d', pr_copy) # Make the pr individual objects always 6 digits
  
  p <- substr(pr_copy,1,3)
  r <- substr(pr_copy,4,6)
  
  
  ## checking for patch/row existance
  start_list <- list()
  t <- 1
  for (i in 1:length(pr_copy)){
    urlP <- sprintf('LandsatTreecover/WRS2/p%s/r%s/p%sr%s_TC_%d/', p[t], r[t], p[t], r[t], Year) #Path part of the url
    urlF <- sprintf('p%sr%s_TC_%d.tif.gz', p[t], r[t], Year) # Filename part of the url
    url <- sprintf('%s%s%s', baseURL = 'ftp://ftp.glcf.umd.edu/glcf/', urlP, urlF)
    Check <- url.exists(url = url)
    start_list[t] <- Check
    t <- t + 1
  }

  if (!('TRUE' %in% start_list)){
      return_list <- 'NA'
      print("Missing Path/Rows within selected buffer")
      return(return_list)

  } else {
    pr$map <- subset(pr$map, start_list != F)
    pr$PR <- subset(pr$PR, start_list != F)
    pr$PATH <- subset(pr$PATH, start_list != F)
    pr$ROW <- subset(pr$ROW, start_list != F)
    # Download data
    
    downloadPR(pr = pr, Year, paste0(dir, 'download'), log = NULL,
               baseURL = "ftp://ftp.glcf.umd.edu/glcf/")
    
    
    
    # Create names for unpacking
    p_filename <- sprintf("%03d", pr$PATH)
    r_filename <- sprintf("%03d", pr$ROW)
    pr_filename <- sprintf("p%sr%s_TC_%s.tif", p_filename, r_filename, Year)
    
    # list all files in folder
    list_file <- list.files(sprintf('%s/%s',dir,extract), full.names=FALSE)
    
    # Listing the pr_filenames in the list
    x <- listing_files(list_file, pr_filename)
    
    # Unpack VCF data
    print("Unpacking downloaded files")
    Unpack_VCF(pr_filename, x, extract, Year, pr, dir)
    
    # list raster files
    list_file <- list.files(sprintf('%s/%s',dir,extract), full.names=FALSE)
    x_list <- listing_files(list_file, pr_filename)
    
    # Mosaicing if multiple data sets are listed, else it takes a single raster
    if (mosaic == T){
      Masked_Raster <- Mosaic_Raster(x_list, dir, extract, buffer, pr_filename)
    } else {
      return(x_list)
    }

    names(Masked_Raster) <- "Sexton"
    
    return_r <- Masked_Raster
    
    return (return_r)
  }
  
}