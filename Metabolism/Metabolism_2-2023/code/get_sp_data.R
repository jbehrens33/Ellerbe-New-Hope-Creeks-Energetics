
get_sp_data <- function(sitecode) {
  
  data_dir <- 'data/'
  if(!dir.exists(data_dir))
    dir.create(data_dir)
  
  save_dir <- glue::glue(data_dir,'raw/')
  if(!dir.exists(save_dir))
    dir.create(save_dir)
  
  if(length(unique(sitecode)) > 1){
    for(i in 1:length(sitecode)){
      
      site <- sitecode[i]
      
      d <- try(StreamPULSE::request_data(sitecode = site))
      
      saveRDS(d, 
              glue::glue(save_dir, '{site}_raw.rds'))
      
    } # end for loop
  } else { # end if statement 
    site <- sitecode[1]
    
    d <- try(StreamPULSE::request_data(sitecode = sitecode))
    
    saveRDS(d, 
            glue::glue(save_dir, '{site}_raw.rds'))
  }
  
} # end function
