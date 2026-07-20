get_sp_data <- function(sitecode,
                        startdate= as.Date(startdate),
                        enddate = as.Date(enddate)) {
  
  # create the data directory if it doesn't exist
  data_dir <- 'data/'
  if(!dir.exists(data_dir))
    dir.create(data_dir)
  
  # and the raw sub-directory to save these data
  save_dir <- glue::glue(data_dir,'raw/')
  if(!dir.exists(save_dir))
    dir.create(save_dir)
  
  # for loop around each site
  for(i in 1:length(sitecode)){
    
    # define the sitecode from the list of sites
    site <- sitecode[i]
    
    # call streampulse portal for the data
    d <- try(StreamPULSE::request_data(sitecode = site, startdate = startdate, enddate = enddate))
    
    # save the list that streampulse returns in 'data/raw'
    saveRDS(d, 
            glue::glue(save_dir, '{site}_raw.rds'))
    
    # show the 3 items in the list
    return(d)
    
  } # end for loop
} # end function


