diag_bayes <- function(mod_dir = 'data/model_runs/Bayes',
                       save_output = TRUE) {
  
  # find the daily estimate output files
  mod_dir_fits <- glue::glue(mod_dir,'/daily')
  
  # and print them
  n_mods <- list.files(mod_dir_fits,
                       full.names = TRUE)
  
  # prepare the output data frame
  out <- data.frame(site = character(),
                    year = character(),
                    n_days = numeric(),
                    f_days = numeric(),
                    resolution = character(),
                    K600_daily_sigma_Rhat = numeric(),
                    err_obs_iid_sigma_Rhat = numeric(),
                    err_proc_iid_sigma_Rhat = numeric(),
                    K_median = numeric(),
                    K_range = numeric(),
                    neg_GPP = numeric(),
                    pos_ER = numeric(),
                    ER_K_r2 = numeric())
  
  # loop around each file
  for(i in 1:length(n_mods)) {
    
    # identify which file
    file <- list.files(mod_dir_fits,
                       full.names = TRUE)[i]
    
    # and read it in
    d <- readr::read_csv(file)
    
    # get the site name
    nms <- stringr::str_split(file, '/')[[1]][5] %>% 
      stringr::str_split(.,'_') %>% 
      unlist() 
    site <- paste(nms[1], nms[2], sep = '_')
    
    # and date
    start_date <- stringr::str_split(file, c('_','/'))[[1]][4] 
    
    end_date <- stringr::str_split(file, c('_','/'))[[1]][5] 
    
    year <- substr(start_date, 1,4)
    
    # start pulling data from this file
    if('GPP_50pct' %in% names(d)) {
      
      # how many days did we get estimates from?
      days <- d %>%
        dplyr::summarise(n = length(GPP_50pct)) %>%
        dplyr::pull()
      
      # how many impossible GPP estimates are there?
      GPP_neg <- d %>%
        dplyr::filter(GPP_50pct < -0.5) %>%
        dplyr::summarise(n = length(GPP_50pct),
                         GPP_neg_per = (n/days)*100) %>%
        dplyr::pull(GPP_neg_per)
      
      # how many positive ER estimates are there?
      ER_pos <- d %>%
        dplyr::filter(ER_50pct > 0.5) %>%
        dplyr::summarise(n = length(ER_50pct),
                         ER_pos_per = (n/days)*100) %>%
        dplyr::pull(ER_pos_per)
      
      # summarise the modeled gas exchange
      K <- d %>%
        dplyr::summarise(K_median = median(K600_daily_50pct, na.rm = TRUE),
                         K_min = min(K600_daily_50pct, na.rm = TRUE),
                         K_max = max(K600_daily_50pct, na.rm = TRUE)) %>%
        dplyr::mutate(K_range = K_max - K_min)
      
      # check equifinality as the R2 between ER and K600
      ER_K_r2 <- summary(lm(data = d,
                            ER_50pct ~ K600_daily_50pct))$adj.r.squared
    } else  # end if statement
      next
    
    # find the modeled and estimate DO file
    dir_DO <- glue::glue(mod_dir,'/mod_and_obs_DO')
    
    # and read it in
    DO <- try(readr::read_csv(glue::glue(dir_DO, '/{site}_{start_date}_{end_date}_mod_and_obs_DO.csv')))
    
    if(inherits(DO, 'try-error')){
      next
    } # end if statement
    
    # what is the resolution of the sensor data?
    res <- glue::glue(diff(DO$solar.time) %>%
                        dplyr::first(),
                      'min')
    
    # find the overall gas exchange model fit file
    dir_KQ <- glue::glue(mod_dir,'/KQ_overall')
    
    KQ <- try(readr::read_csv(glue::glue(dir_KQ, '/{site}_{start_date}_{end_date}_KQ_overall.csv')))
    
    if(inherits(KQ, 'try-error')){
      next
    } # end if statement
    
    # pull the convergence statistic for K600
    K600_daily_sigma_Rhat <- KQ %>%
      dplyr::pull(K600_daily_sigma_Rhat)
    
    # find the overall model fit file
    dir_overall <- glue::glue(mod_dir,'/overall')
    
    # and read it in
    overall <- try(readr::read_csv(glue::glue(dir_overall, '/{site}_{start_date}_{end_date}_overall.csv')))
    
    if(inherits(overall, 'try-error')){
      next
    } # end if statement
    
    # extract the model convergence statistics for process and observation errors
    err_proc_iid_sigma_Rhat <- overall %>%
      dplyr::pull(err_proc_iid_sigma_Rhat)
    
    err_obs_iid_sigma_Rhat <- overall %>%
      dplyr::pull(err_obs_iid_sigma_Rhat)
    
    # compile the output diagnostic file
    out <- out %>%
      dplyr::add_row(site = site,
                     year = year,
                     n_days = days,
                     f_days = (days/365)*100,
                     resolution = res,
                     K600_daily_sigma_Rhat = K600_daily_sigma_Rhat,
                     err_obs_iid_sigma_Rhat = err_obs_iid_sigma_Rhat,
                     err_proc_iid_sigma_Rhat = err_proc_iid_sigma_Rhat,
                     K_median = K$K_median,
                     K_range = K$K_range,
                     neg_GPP = GPP_neg,
                     pos_ER = ER_pos,
                     ER_K_r2 = ER_K_r2)
  } # end for loop
  
  # save the diagnostic file
  if(save_output){
    readr::write_csv(out,
                     glue::glue(mod_dir,'/diag_metab.csv'))
  }
  
  # print to the console; can be saved as an object
  return(out)
  
} # end function
