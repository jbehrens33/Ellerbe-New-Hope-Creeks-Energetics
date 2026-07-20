prep_sp_data <- function(raw_dir = 'data/raw/',                # define where the sensor data is saved
                         rc_dir = 'data/rating_curves',        # if you have rating curves, save them here
                         rm_flagged = list('Bad Data',         # remove flagged data as Bad Data, Questionable, Interesting
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
      rc <- readr::read_csv(glue::glue(rc_dir, '/hydro_nhc-q_5-4-2023')) %>% 
        
        # for high and low flows, we'll approximate Q via Manning's equation
        mutate(Q_use = ifelse(Q_ratingcurve > 6.3,
                              Q_mannings,
                              Q_ratingcurve),
               Q_use = ifelse(Q_use < 0.08,
                              Q_mannings,
                              Q_use)) 
      # re-run the prep function, being specific for this site
      # define the depth and discharge data
      Z_data <- rc$depth
      Q_data <- rc$Q_use
      
      prep <- prep_metabolism(raw,                            # which data; raw is a list and contains more than just the time-series
                              type = 'bayes',                 # we want to run a Bayesian model
                              estimate_areal_depth = TRUE,    # yes estimate areal depth
                              rm_flagged = list('Bad Data',   # remove flagged data
                                                'Questionable', 
                                                'Interesting'),
                              fillgaps = 'interpolation',      # interpolate any gaps
                              
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
