model_mle <- function(sm_input_dir = 'data/sm_input/',
                      years = NULL,
                      calc_water_year = TRUE) {
  
  sites_long <- list.files(sm_input_dir)

  sites <- gsub('_sm_input.csv','',
                sites_long)
  
  site_years <- data.frame(
    site = rep(sites, each = length(years)),
    year = rep(years, times = length(sites)))
  
  tracker <- matrix(nrow = length(site_years$site),
                    ncol = 5)
  colnames(tracker) <- c('Site', 'Year', 'Has Data', 'Fit Error', 'Fit Time')
  
  # define the model
  init_name = streamMetabolizer::mm_name(
    type = 'mle',
    pool_K600 = 'none',
    err_obs_iid = TRUE,
    err_proc_iid = FALSE,
    ode_method = 'trapezoid',
    deficit_src = 'DO_mod',
    engine = 'nlm',
    GPP_fun = 'linlight',
    ER_fun = 'constant')
  
  # set all specs, including priors
  init_specs = streamMetabolizer::specs(model_name = init_name)
  
  for(i in 1:length(site_years$site)){
    
    site <- site_years$site[i]
    year <- site_years$year[i]
    
    tracker[i, 'Site'] <- site
    tracker[i, 'Year'] <- year
    
    d <- readr::read_csv(glue::glue(sm_input_dir, site,'_sm_input.csv'))

    if(calc_water_year) {
      d_use <- d %>% 
        mutate(water_year = ifelse(lubridate::month(solar.time) > 9,
                                   lubridate::year(solar.time) + 1,
                                   lubridate::year(solar.time))) %>% 
        filter(water_year %in% year) %>% 
        select(-water_year)
    } else {
      d_use <- d %>% 
        filter(lubridate::year(d$solar.time)%in% year)
    }
    
    if(nrow(d_use) == 0){
      tracker[i, 'Has Data'] <- 'No Data'
    } else {
      tracker[i, 'Has Data'] <- 'Has Data'
    }
    
    start_date <- dplyr::first(lubridate::date(d_use$solar.time))
    end_date <- dplyr::last(lubridate::date(d_use$solar.time))
    
    fit_err <- FALSE
    tryCatch(
      
      # tries to run the model with specs and data defined above
      {dat_metab = streamMetabolizer::metab(init_specs,
                                            data = d_use)
      
      # extract the outputs from the model
      dat_fit = streamMetabolizer::get_fit(dat_metab)
      },
      
      # evalute fault tolerance: did the model fit?
      error = function(e) {
        fit_err <- TRUE
      }#,
      
      # warning = function(w) {
      #   fit_err <<- TRUE
      # }
    )
    
    # did the model fit or not?
    if(fit_err){
      tracker[i, 'Fit Error'] <- 'Fit Error'
      tracker[i, 'Fit Time'] <- 'NA'
    } else {
      tracker[i, 'Fit Error'] <- 'Fit Successful'
      tracker[i, 'Fit Time'] <- streamMetabolizer::get_fitting_time(dat_metab)[3]
    }
    
    # save model output
    # where are we saving data to
    write_dir = glue::glue('data/model_runs/')
    
    if(!dir.exists(write_dir)){
      dir.create(write_dir)
    }
    
    if(!dir.exists(glue::glue(write_dir,'/MLE')))
      dir.create(glue::glue(write_dir, '/MLE'))
    if(!dir.exists(glue::glue(write_dir,'/MLE/dat_fit')))
      dir.create(glue::glue(write_dir, '/MLE/dat_fit'))
    if(!dir.exists(glue::glue(write_dir,'/MLE/specs')))
      dir.create(glue::glue(write_dir, '/MLE/specs'))
    if(!dir.exists(glue::glue(write_dir,'/MLE/mod_obs_DO')))
      dir.create(glue::glue(write_dir, '/MLE/mod_obs_DO'))
    
    # write the fit data
    fn_prefix <- paste0(write_dir, 'MLE/dat_fit/', site, '_', start_date, '_', end_date, '_')
    readr::write_csv(dat_fit, paste0(fn_prefix, 'dat_fit.csv'))
    
    # write the specs
    specs_out <- data.frame(unlist(streamMetabolizer::get_specs(dat_metab)))
    fn_prefix <- paste0(write_dir, 'MLE/specs/', site, '_', start_date, '_', end_date, '_')
    readr::write_csv(specs_out, paste0(fn_prefix, 'specs.csv'))
    
    # write the output data
    data_out = streamMetabolizer::get_data(dat_metab)
    fn_prefix = paste0(write_dir, 'MLE/mod_obs_DO/', site, '_', start_date, '_', end_date, '_')
    readr::write_csv(data_out, paste0(fn_prefix, 'mod_and_obs_DO.csv'))
    
  } # end for loop
  
  # return the matrix summarizing the model fits
  readr::write_csv(data.frame(tracker),
                   'data/model_runs/MLE/MLE_results.csv')
  
  return(data.frame(tracker))
  
} # end function
