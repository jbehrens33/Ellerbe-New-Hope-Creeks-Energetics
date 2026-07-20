prep_sp_data <- function(raw_dir = 'data/raw/'){
  
  sites <- list.files(raw_dir,
                      full.names = TRUE)
  
  for(i in 1:length(unique(sites))){
    
    raw <- readRDS(sites[i])
    
    site <- gsub('_raw.rds','',
                 stringr::str_split(sites[i], '/')[[1]][3])
    
    prep <- try(prep_metabolism(raw,                            # which data; raw is a list and contains more than just the time-series
                                type = 'bayes',                 # we want to run a Bayesian model
                                estimate_areal_depth = TRUE,    # yes estimate areal depth
                                rm_flagged = list('Bad Data',   # remove flagged data
                                                  'Questionable', 
                                                  'Interesting'),
                                fillgaps = 'interpolation'      # interpolate any gaps
    )
    )
    
    # this error will come out for NHC:
    # Error: Missing discharge and depth data.
    # Not enough information to proceed.
    # Might parameter "zq_curve" be of service?
    
    # if prep_metabolism() fails, we get sent here
    # we will load in external discharge, depth, or rating curve coefficients and re-run prep_metabolism
    if(inherits(prep, 'try-error')){
      
      # read in NHC rating curve data, from Brooke and Steve
      rc <- readr::read_csv('data/rating_curves/hydro_nhc-q_12212021.csv') %>% 
        
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
    
    prep_dir <- 'data/prep/'
    if(!dir.exists(prep_dir))
      dir.create(prep_dir)
    
    sm_input_dir <- 'data/sm_input/'
    if(!dir.exists(sm_input_dir))
      dir.create(sm_input_dir)
    
    saveRDS(prep,
            glue::glue(prep_dir, '{site}_prep.rds'))
    
    readr::write_csv(prep$data, 
                     glue::glue(sm_input_dir, '{site}_sm_input.csv'))
    
  } # end for loop
} # end function
