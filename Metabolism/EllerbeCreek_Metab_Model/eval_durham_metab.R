## Header ----
## Script name: eval_durham_metab.R
##
## Purpose of script: Evaluate streampulse models from Durham, NC streams
##
## Author: Nick Marzolf
## Date Created: 
## Date Modified: 2022-04-28
## Email: nicholas.marzolf@duke.edu
##
## load packages:  
library(tidyverse)
library(dplyr)
library(ggplot2)
library(streamMetabolizer)
library(StreamPULSE)
library(lubridate)
library(glue)
library(cowplot)
##
## clear the environment if needed
rm(list = ls())
##

# assign sites
run_sites <- c('NC_EllerbeGlenn', 'NC_EllerbeClub')

# loop to access the large model outputs into usable chunks
for(i in run_sites){
  
  # which site
  sitecode <- run_sites[i]
  
  # read in the RDS model output, which is large
  d <- readRDS(glue('data/model_runs/raw_outputs/{sitecode}_mod.rds'))
  
  # break that object into more usable csv files and save them
  # be sure to recreate the folder directory 
  # DO predictions
  write_csv(d$predictions, 
            glue('data/model_runs/mod_and_obs_DO//{sitecode}_mod_obs_DO.csv'))
  # model specs
  write_csv(d$details, 
            glue('data/model_runs/specs/{sitecode}_specs.csv'))
  
  # extract the fitted data
  fit <-get_fit(d$fit)
  write_csv(fit$daily, 
            glue('data/model_runs/daily/{sitecode}_daily.csv'))
  write_csv(fit$KQ_overall, 
            glue('data/model_runs/KQ_overall/{sitecode}_KQ_overall.csv'))
  write_csv(fit$overall, 
            glue('data/model_runs/KQ_overall/{sitecode}_overall.csv'))
}

# loop to make and save some plots
for(j in 1:3) {
  sitecode <- run_sites[j]
  
  d <- read_csv(glue('data/model_runs/daily/{sitecode}_daily.csv'))
  
  preds <- ggplot(d %>% 
                    filter(GPP_mean > -0.5,
                           ER_mean < 0.5), 
                  aes(x = date))+
    geom_point(aes(y = GPP_mean, color = 'GPP'))+
    geom_ribbon(aes(ymin = GPP_2.5pct, ymax = GPP_97.5pct, fill = 'GPP'), alpha = 0.5)+
    geom_point(aes(y = ER_mean, color = 'ER'))+
    geom_ribbon(aes(ymin = ER_2.5pct, ymax = ER_97.5pct, fill = 'ER'), alpha = 0.5)+
    scale_color_manual(name = element_blank(),
                       values = c('darkgoldenrod4', 'darkgreen'))+
    scale_fill_manual(name = element_blank(),
                      values = c('darkgoldenrod4', 'darkgreen'))+
    ylab(expression(paste('g ', O[2], m^-2, d^-1)))+
    ggtitle(sitecode)+
    facet_wrap(lubridate::year(date))
  
  profile <- ggplot(d %>% 
                      filter(GPP_mean > -0.5,
                             ER_mean < 0.5), 
                    aes(x = GPP_mean, y = ER_mean))+
    geom_point(alpha = 0.5)+
    geom_hline(yintercept = 0, linetype = 'dashed')+
    geom_vline(xintercept = 0, linetype = 'dashed')+
    ylab(expression(paste('ER (g ', O[2], m^-2, d^-1,')')))+
    xlab(expression(paste('GPP (g ', O[2], m^-2, d^-1,')')))
  
  plot_grid(preds, profile,
            ncol = 2, align = 'hv')
  
  ggsave(glue('figures/preds/{sitecode}_preds.png'),
         dpi = 300, width = 9, height = 4)
}


# evaluate the model fits
rating_tbl <- data.frame(
  site = character(),
  year = character(),
  K600_daily_sigma_Rhat = numeric(),
  err_proc_iid_sigma_Rhat = numeric()
)

for(k in 1:length(run_sites)){
  sitecode <- run_sites[k]
  
  daily <- read_csv(glue('data/model_runs/daily/{sitecode}_daily.csv'))
  KQ_overall <- read_csv(glue('data/model_runs/KQ_overall/{sitecode}_KQ_overall.csv'))
  overall <- read_csv(glue('data/model_runs/KQ_overall/{sitecode}_overall.csv'))
  
  length = length(daily$GPP_mean)
  GPP_neg <- daily %>% 
    filter(GPP_mean < -0.5) %>% 
    summarise(GPP_imp_per = n()/length) %>% 
    pull()
  
  ER_pos <- daily %>% 
    filter(ER_mean < 0.5) %>% 
    summarise(ER_imp_per = n()/length) %>% 
    pull()
  
  K600_sum <- daily %>% 
    summarise(meanK = mean(K600_daily_mean, na.rm = TRUE),
              minK600 = min(K600_daily_mean, na.rm = TRUE),
              maxK600 = max(K600_daily_mean, na.rm = TRUE)) %>% 
    mutate(rangeK = maxK600 - minK600)
  
  equifin_r2 <- summary(lm(data = daily,
                           ER_mean ~ K600_daily_mean))[8]
  
  rating_tbl <- rating_tbl %>% 
    add_row(
      site = sitecode,
      K600_daily_sigma_Rhat = KQ_overall$K600_daily_sigma_Rhat,
      err_proc_iid_sigma_Rhat = overall$err_proc_iid_sigma_Rhat) %>% 
    mutate(
      RhatL = K600_daily_sigma_Rhat > 1.2 | err_proc_iid_sigma_Rhat > 1.2,
      RhatH = !RhatL,
      KrangeL = K600_sum$rangeK > 50,
      KrangeM = K600_sum$rangeK <= 50 & K600_sum$rangeK > 15,
      KrangeH = K600_sum$rangeK < 15,
      GPP_L = GPP_neg > 50,
      GPP_M = GPP_neg <= 50 & GPP_neg > 25,
      GPP_H = GPP_neg <= 25,
      ER_L = ER_pos > 50,
      ER_M = ER_pos <= 50 & ER_pos > 25,
      ER_H = ER_pos <= 25
      )
}





