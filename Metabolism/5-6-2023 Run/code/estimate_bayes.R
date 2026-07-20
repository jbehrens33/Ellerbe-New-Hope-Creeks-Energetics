estimate_bayes <- function(file_dir = 'data/model_runs/Bayes/',
                           save_output = TRUE) {
  
  # find the models
  mod_dir <- glue::glue(file_dir, 'daily')
  
  # how many models are there? This will be the object we loop around
  n_mods <- length(list.files(mod_dir,
                              full.names = TRUE))
  
  for(i in 1:n_mods) {
    
    # find an individual file
    file <- list.files(mod_dir,
                       full.names = TRUE)[i]
    
    # get the name
    nms <- stringr::str_split(file, '/')[[1]][5] %>% 
      stringr::str_split(.,'_') %>% 
      unlist() 
    site <- paste(nms[1], nms[2], sep = '_')
    
    # and date info
    start_date <- stringr::str_split(file, c('_','/'))[[1]][4] 
    
    end_date <- stringr::str_split(file, c('_','/'))[[1]][5] 
    
    year <- substr(start_date, 1,4)
    
    # read the file in 
    d <- readr::read_csv(file)
    
    if('GPP_50pct' %in% names(d)) {
      
      # remove the impossible results and rename column headers for easier use
      d_clean <- d %>%
        dplyr::filter(GPP_50pct > -0.5,
                      ER_50pct < 0.5) %>%
        dplyr::select(date,
                      GPP = GPP_50pct,
                      GPP.lower = GPP_2.5pct,
                      GPP.upper = GPP_97.5pct,
                      GPP.n_eff = GPP_n_eff,
                      GPP.Rhat = GPP_Rhat,
                      ER = ER_50pct,
                      ER.lower = ER_2.5pct,
                      ER.upper = ER_97.5pct,
                      ER.n_eff = ER_n_eff,
                      ER.Rhat = ER_Rhat,
                      K600 = K600_daily_50pct,
                      K600.lower = K600_daily_2.5pct,
                      K600.upper = K600_daily_97.5pct,
                      K600.n_eff = K600_daily_n_eff,
                      K600.Rhat = K600_daily_Rhat
        )
    }
    
    # find the DO model and observation file
    dir_DO <- glue::glue(file_dir, 'mod_and_obs_DO')
    
    # read it in
    DO <- try(readr::read_csv(glue::glue(dir_DO, '/{site}_{start_date}_{end_date}_mod_and_obs_DO.csv')))
    
    if(inherits(DO, 'try-error')){
      next
    }
    
    # find the resolution of the sensor data
    if(!'date' %in% names(DO)){
      DO <- DO %>%
        dplyr::mutate(date = lubridate::date(solar.time))
    }
    
    res <- glue::glue(diff(DO$solar.time) %>%
                        dplyr::first(),
                      'min')
    
    # summarise the input DO data (mean, min/max/range) and other input data (Q, light, depth)
    DO_clean <- DO %>%
      dplyr::group_by(date) %>%
      dplyr::summarise(temp.water = mean(temp.water, na.rm = TRUE),
                       DO.obs = mean(DO.obs, na.rm = TRUE),
                       DO.sat = mean(DO.sat, na.rm = TRUE),
                       DO.min = min(DO.obs, na.rm = TRUE),
                       DO.max = max(DO.obs, na.rm = TRUE),
                       depth = mean(depth, na.rm = TRUE),
                       discharge = mean(discharge, na.rm = TRUE),
                       shortwave = mean(light, na.rm = TRUE)) %>%
      dplyr::mutate(DO.amp = DO.max - DO.min,
                    DO.psat = (DO.obs/DO.sat)*100)
    
    # merge by date and get what data we want
    merged <- dplyr::left_join(d_clean,
                               DO_clean,
                               'date') %>%
      dplyr::mutate(resolution = res,
                    site = site) %>%
      dplyr::select(site, resolution, date, GPP, GPP.lower, GPP.upper, GPP.n_eff, GPP.Rhat,
                    ER, ER.lower, ER.upper, ER.n_eff, ER.Rhat, K600, K600.lower, K600.upper,
                    K600.n_eff, K600.Rhat, DO.obs, DO.sat, DO.amp, DO.psat, depth, temp.water,
                    discharge, shortwave)
    
    # save the site-year specific estimates to a new directory
    write_dir <- glue::glue(file_dir, 'estimates/')
    
    if(!dir.exists(write_dir))
      dir.create(write_dir)
    
    readr::write_csv(merged,
                     glue::glue(write_dir, 'neon_metab_{site}_{year}.csv'))
    
  } # end for loop
  
  # once all estimates are made, compile them into a single file
  estimates <- purrr::map_dfr(list.files(write_dir,
                                         full.names = TRUE),
                              readr::read_csv)
  # and save it if desired
  if(save_output){
    readr::write_csv(estimates,
                     glue::glue(file_dir, 'metab_estimates.csv'))
  }
  
  # and print the output compiled estimates
  return(estimates)
}
