eval_mle <- function(mle_dir = 'data/model_runs/MLE'){
  
  dir_fits <- glue::glue(mle_dir, '/dat_fit')
  
  fits <- list.files(dir_fits,
                     full.names = TRUE)
  
  mle_diagnostic <- data.frame(
    site = character(),
    year = character(),
    GPP_neg = numeric(),
    ER_pos = numeric(),
    K_median = numeric(),
    K_sd = numeric(),
    no_dat_num = numeric(),
    no_dat_dates = character()
  )
  
  for(i in 1:length(fits)) {
    
    fit <- readr::read_csv(fits[i])
    
    nms <- stringr::str_split(fits[i], '/')[[1]][5] %>% 
      stringr::str_split(.,'_') %>% 
      unlist() 
    site <- paste(nms[1], nms[2], sep = '_')
      
    year <- stringr::str_split(fits[i], c('/','_'))[[2]][5] %>% 
      substr(.,1,4)
    
    GPP_neg <- fit %>%
      dplyr::filter(GPP.daily < -0.5) %>%
      dplyr::summarise(n = n()) %>%
      dplyr::pull(n)
    
    ER_pos <- fit %>%
      dplyr::filter(ER.daily > 0.5) %>%
      dplyr::summarise(n = n()) %>%
      dplyr::pull(n)
    
    K_sum <- fit %>%
      dplyr::filter(K600.daily > 0) %>%
      dplyr::summarise(medianK = median(K600.daily, na.rm = TRUE),
                       sdK = sd(K600.daily, na.rm = TRUE))
    
    no_dat_num <- sum(is.na(fit$GPP.daily))
    
    no_dat_dates <- fit %>%
      dplyr::filter(is.na(GPP.daily)) %>%
      dplyr::pull(date)
    
    mle_diagnostic <- mle_diagnostic %>%
      dplyr::add_row(
        site = site,
        year = year,
        GPP_neg = GPP_neg,
        ER_pos = ER_pos,
        K_median = dplyr::pull(K_sum, medianK),
        K_sd = dplyr::pull(K_sum,sdK),
        no_dat_num = no_dat_num,
        no_dat_dates = stringr::str_c(as.character(no_dat_dates), 
                                      collapse = ' ')
      )
  } # end for loop
  
  readr::write_csv(mle_diagnostic,
                   glue::glue(mle_dir, '/MLE_diagnostics.csv'))
  
  return(mle_diagnostic)
}
