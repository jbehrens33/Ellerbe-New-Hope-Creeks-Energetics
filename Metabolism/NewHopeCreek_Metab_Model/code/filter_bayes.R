filter_bayes <- function(mod_dir = 'data/model_runs',
                         save_output = TRUE) {
  
  # where are the estimate files located
  est_dir <- glue::glue(mod_dir, '/Bayes/')
  
  # read in the estimates and diagnostic file
  est <- readr::read_csv(glue::glue(est_dir, 'metab_estimates.csv'))
  diags <- readr::read_csv(glue::glue(est_dir, 'diag_metab.csv'))
  
  # identify the good site years based on the following criteria from the diagnostic file
  good_site_years <- diags %>%
    dplyr::filter(n_days > 365*0.6,
                  err_obs_iid_sigma_Rhat < 1.2,
                  err_proc_iid_sigma_Rhat < 1.2,
                  K600_daily_sigma_Rhat < 1.2,
                  K_median < 100,
                  ER_K_r2 < 0.6,
                  neg_GPP < 60,
                  pos_ER < 60) %>%
    dplyr::mutate(site_year = paste(site, year, sep = '-')) %>%
    dplyr::pull(site_year)
  
  # find the good site years in the estimate file and keep them
  good_est <- est %>%
    dplyr::mutate(year = lubridate::year(date),
                  site_year = paste(site, year, sep = '-')) %>%
    dplyr::filter(site_year %in% good_site_years) %>%
    dplyr::select(-site_year)
  
  # save the good data as a file if desired
  if(save_output){
    write_dir <- 'data/filtered_estimates/'
    
    if(!dir.exists(write_dir))
      dir.create(write_dir)
    
    readr::write_csv(good_est,
                     glue::glue(write_dir,'/filtered_metab_estimates.csv'))
  }
  
  # return the object
  return(good_est)
}