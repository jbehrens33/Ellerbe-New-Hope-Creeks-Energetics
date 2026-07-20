prep_sp_data_v2 <- function(raw_dir = 'data/raw/',                # define where the sensor data is saved
                            rc_dir = 'data/rating_curves',        # if you have rating curves, save them here
                            rm_flagged = list('Bad Data', 
                                              'Interesting',     # remove flagged data as Bad Data, Questionable, Interesting
                                              'Questionable')){
  # list of the sites to work with
  sites <- list.files(raw_dir,
                      full.names = TRUE)
  
  # for loop around each site
  for(i in 1:length(unique(sites))){
    
    # which site are we working on
    raw <- readRDS(sites[i])
    
    # extract the site name from the file name
    site <- gsub('_raw.rds','',
                 stringr::str_split(sites[i], '/')[[1]][3])
    
    # use prep_metabolism
    prep <- try(StreamPULSE::prep_metabolism(raw,                          # which data; raw is a list and contains more than just the time-series
                                             type = 'bayes',               # we want to run a Bayesian model
                                             estimate_areal_depth = TRUE,  # yes estimate areal depth
                                             rm_flagged = rm_flagged,      # what data are we removing
                                             fillgaps = 'interpolation'    # interpolate any gaps
    )
    )
    
    
    # if prep_metabolism() fails, we get sent here
    # this is an example of how we process data from NC_NHC where there is no USGS gage station
    # but there is externally available discharge and depth data
    if(inherits(prep, 'try-error') & site == 'NC_NHC'){
      
      # read in NHC rating curve data, from Brooke and Steve
      rc <- readr::read_csv(glue::glue(rc_dir, '/hydro_nhc-q_5-4-2023.csv')) %>% 
        # for high and low flows, we'll approximate Q via Manning's equation
        mutate(Q_high = ifelse(Q_ratingcurve > 6.3,
                              Q_mannings,
                              Q_ratingcurve),
               USGSDischarge_m3s = ifelse(Q_high < 0.08, # call this so that prep_metabolism recognizes it
                              Q_mannings,
                              Q_high)) %>% 
        rename(Depth_m = depth) %>% 
        dplyr::distinct(DateTime_UTC, .keep_all = TRUE)
      # re-run the prep function, being specific for this site
      # define the depth and discharge data
      Z_data <- rc$Depth_m
      Q_data <- rc$USGSDischarge_m3s
      
      # Set Q as "USGS Discharge m3s" so that prep_metabolism can run
      raw2<-raw

      raw2$data <- raw2$data %>% 
        pivot_wider(names_from = "variable", values_from = "value") %>% 
        left_join(select(rc, c(DateTime_UTC, USGSDischarge_m3s, Depth_m))) %>% # left-join the discharge data
        pivot_longer(names_to = "variable", values_to = "value", -c(DateTime_UTC, region, site, flagtype, flagcomment)) %>% 
        filter(!is.na(value)) %>% 
        filter(!flagtype %in% rm_flagged) %>% 
        select(DateTime_UTC, region, site, value, variable, flagtype, flagcomment) %>%  # put in order expected by prep_metabolism
        as.data.frame()
      
      prep <- prep_metabolism(raw2,                            # which data; raw is a list and contains more than just the time-series
                              type = 'bayes',                 # we want to run a Bayesian model
                              estimate_areal_depth = TRUE,    # yes estimate areal depth
                              rm_flagged = rm_flagged,
                              fillgaps = 'interpolation', 

                             # populate this argument
                              zq_curve = list(Z = Z_data,
                                              Q = Q_data,
                                              fit = 'power',
                                              ignore_oob_Z = FALSE))
    }
    
    # prepare to save outputs as a list
    prep_dir <- 'data/prep/'
    if(!dir.exists(prep_dir))
      dir.create(prep_dir)
    
    # and as a csv of data to go into the model
    sm_input_dir <- 'data/sm_input/'
    if(!dir.exists(sm_input_dir))
      dir.create(sm_input_dir)
    
    # and save the data
    saveRDS(prep,
            glue::glue(prep_dir, '{site}_prep.rds'))
    
    readr::write_csv(prep$data, 
                     glue::glue(sm_input_dir, '{site}_sm_input.csv'))
    
    return(head(prep$data))
    
  } # end for loop
} # end function
