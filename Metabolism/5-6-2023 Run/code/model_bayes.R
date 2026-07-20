model_bayes <- function(sm_input_dir = 'data/sm_input/',
                        years = NULL,
                        calc_water_year = TRUE) {
  
  sites_long <- list.files(sm_input_dir)
  
  sites <- gsub('_sm_input.csv','',
                sites_long)
  
  site_years <- data.frame(
    site = rep(sites, each = length(years)),
    year = rep(years, times = length(sites)))
  
  tracker <- matrix(nrow = length(site_years$site),
                    ncol = 7)
  # we will track each site-year and record if there: ,
  colnames(tracker) = c('Site',
                        'Year',
                        'Has Data',       # is data (T/F)
                        'Setup Error',    # a K-Q relationship can be defined (T/F)
                        'Fit Error',      # there was an error fitting the model (T/F)
                        'Fit Error Text', # record the text from a model error (text)
                        'Fit Time')       # and if the model run, how long did it take (numeric)
  
  # Define the model name, which establishes core hyperparameters. These would only be
  # changed as a last resort if the model just refuses to converge. i.e., if
  # we have to settle for a simpler model
  bayes_name <- streamMetabolizer::mm_name(
    type = 'bayes',
    pool_K600 = 'binned',
    err_obs_iid = TRUE,
    err_proc_iid = TRUE,
    ode_method = 'trapezoid',
    deficit_src = 'DO_mod',
    engine = 'stan')
  
  # Define fault tolerances
  has_data <- TRUE
  setup_err <- FALSE
  fit_err <- FALSE
  
  # define write dir
  write_dir <- glue::glue('data/model_runs/Bayes')
  
  # and create the directory if it doesn't exist
  if(!dir.exists(write_dir))
    dir.create(write_dir)
  
  # set for loop
  for(i in 1:length(site_years$site)){
    
    site_id <- site_years[i, 'site']
    run_year <- site_years[i, 'year']
    
    # keep track of our results by site and by year
    tracker[i, 'Site'] <- site_id
    tracker[i, 'Year'] <- run_year
    
    # read in the data
    input_dat <- tryCatch({
      
      # read in site model input data
      readr::read_csv(glue::glue(sm_input_dir,'{site_id}_sm_input.csv'))
    },
    # fault tolerance if the file couldn't be read in; not, jump to the next site
    error = function(e) {
      cat('No data')
      next
    }
    )
    
    # Step 2: filter data to the run year or water year
    if(calc_water_year) {
      input_dat_sub <- input_dat %>% 
        mutate(water_year = ifelse(lubridate::month(solar.time) > 9,
                                   lubridate::year(solar.time) + 1,
                                   lubridate::year(solar.time))) %>% 
        filter(water_year %in% run_year) %>% 
        select(-water_year)
    } else {
      input_dat_sub <- input_dat %>% 
        filter(lubridate::year(d$solar.time)%in% run_year)
    }
    
    # if there are no rows in the dataframe, there is no data from that site-year and jump to the next
    if(nrow(input_dat_sub) == 0){
      tracker[i, 'Has Data'] <- FALSE
    } else {
      tracker[i, 'Has Data'] <- TRUE
    }
    
    # Step 4: define KQ priors and model specifications
    # pull the earliest and latest date in the modeled run for file naming purposes
    tryCatch({
      # define start and end dates from the filtered dataset
      start_date <- lubridate::date(dplyr::first(input_dat_sub$solar.time))
      end_date <- lubridate::date(dplyr::last(input_dat_sub$solar.time))
      
      # read in results from MLE model run
      mle_priors <- readr::read_csv(glue::glue('data/model_runs/MLE/MLE_diagnostics.csv'))
      
      # subset to pull gas exchange results
      median_K600_mle <- mle_priors %>%
        dplyr::filter(site == site_id) %>%
        dplyr::select(site, year, K_median)
      
      # create an object that sorts discharge in log space
      Q_by_day <- tapply(log(input_dat_sub$discharge),
                         substr(input_dat_sub$solar.time, 1, 10),
                         mean)
      
      # define each node of log discharge by 0.2 units
      K600_lnQ_node_centers <- seq(from = min(Q_by_day, na.rm = TRUE),
                                   to = max(Q_by_day, na.rm = TRUE),
                                   by = 0.2)
      
      # and how many nodes are there
      n_nodes <- length(K600_lnQ_node_centers)
      
      # define the model specs and other hyperparameters
      bayes_specs <- streamMetabolizer::specs(
        model_name = bayes_name,
        K600_lnQ_nodes_centers = K600_lnQ_node_centers,
        K600_lnQ_nodediffs_sdlog = 0.1,
        K600_lnQ_nodes_meanlog = rep(log(12), length(K600_lnQ_node_centers)),  # default; see ?specs
        K600_lnQ_nodes_sdlog = rep(0.1, n_nodes),
        K600_daily_sigma_sigma = (dplyr::filter(median_K600_mle,
                                                year == run_year) %>%
                                    dplyr::pull(K_median))* 0.02,
        burnin_steps = 1000,
        saved_steps = 500)
    },
    # if we can't define the K-Q relationship or model specs, there model setup fails
    error = function(e){
      setup_err <- TRUE
    }
    )
    
    # record in the tracker if the model setup is successful, and if not, jump to the next run
    if(setup_err == TRUE) {
      tracker[i, 'Setup Error'] <- TRUE
    } else {
      tracker[i, 'Setup Error'] <- FALSE
    }
    
    # Step 5: run the model
    tryCatch({
      # run the model
      dat_metab <- streamMetabolizer::metab(specs = bayes_specs,
                                            data = input_dat_sub)
      # if the model returns a warning, fit error is TRUE and we move
      # to the next site-year
      model_error <- dat_metab@fit$errors
      
      # if an error exists, model_error object will contain some character object(s), of some length
      if(length(model_error) > 0) {
        fit_err <- TRUE
      } else {
        fit_err <- FALSE
        model_error <- NA
        
        fit_time <- streamMetabolizer::get_fitting_time(dat_metab)
        tracker[i, 'Fit Time'] <- fit_time[3]
        dat_fit <- streamMetabolizer::get_fit(dat_metab)
      }
    },
    # if an error (of literal class 'error') this code will run instead:
    error = function(e){
      fit_err <- TRUE
      model_error <- "streamMetabolizer::metab() function ERROR"
    },
    # this will run at end of tryCatch no matter what
    finally = {
      tracker[i, 'Fit Error'] <- fit_err
      tracker[i, 'Fit Error Text'] <- model_error
    }
    )
    
    # Step 6: write model output
    # define fault tolerance
    write_err <- FALSE
    
    # create, if necessary, the directory structure of the output files
    tryCatch({
      if(!dir.exists(glue::glue(write_dir, '/daily')))
        dir.create(glue::glue(write_dir, '/daily'))
      if(!dir.exists(glue::glue(write_dir, '/overall')))
        dir.create(glue::glue(write_dir, '/overall'))
      if(!dir.exists(glue::glue(write_dir, '/KQ_overall')))
        dir.create(glue::glue(write_dir, '/KQ_overall'))
      if(!dir.exists(glue::glue(write_dir, '/specs')))
        dir.create(glue::glue(write_dir, '/specs'))
      if(!dir.exists(glue::glue(write_dir, '/datadaily')))
        dir.create(glue::glue(write_dir, '/datadaily'))
      if(!dir.exists(glue::glue(write_dir, '/mod_and_obs_DO')))
        dir.create(glue::glue(write_dir, '/mod_and_obs_DO'))
      
      # define file name prefix
      fn_prefix = paste0('/', site_id, '_', start_date, '_', end_date, '_')
      
      # create objects to save from the model outputs
      specs_out = data.frame(unlist(streamMetabolizer::get_specs(dat_metab)))
      daily_out = data.frame(streamMetabolizer::get_data_daily(dat_metab))
      data_out = data.frame(streamMetabolizer::get_data(dat_metab))
      
      # write the model output objects to the directory needed
      readr::write_csv(data.frame(dat_fit$daily),
                       paste0(write_dir, '/daily', fn_prefix, 'daily.csv'))
      readr::write_csv(data.frame(dat_fit$overall),
                       paste0(write_dir, '/overall',fn_prefix, 'overall.csv'))
      readr::write_csv(data.frame(dat_fit$KQ_overall),
                       paste0(write_dir,'/KQ_overall',fn_prefix, 'KQ_overall.csv'))
      readr::write_csv(specs_out,
                       paste0(write_dir,'/specs', fn_prefix, 'specs.csv'))
      readr::write_csv(daily_out,
                       paste0(write_dir,'/datadaily',fn_prefix, 'datadaily.csv'))
      readr::write_csv(data_out,
                       paste0(write_dir,'/mod_and_obs_DO',fn_prefix, 'mod_and_obs_DO.csv'))
    },
    # fault tolerance: if there was an error writing these files, flag the error
    error = function(e) {
      write_err <- TRUE
    }
    )
    
  }# end for loop
  
  readr::write_csv(data.frame(tracker),
                   glue::glue(write_dir, '/neon_bayes_row.csv')) #save run results if you want to
  
  # write the site-year tracker information
  return(tracker)
  
} # end function