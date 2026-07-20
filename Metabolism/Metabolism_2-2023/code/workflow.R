library(dplyr)
library(StreamPULSE)

prep_bayes_model()

sp_sites <- StreamPULSE::query_available_data(region = 'NC') 

my_sites <- c('NC_NHC', 'NC_EllerbeGlenn', 'NC_EllerbeClub')


get_sp_data(sitecode = my_sites)

prep_sp_data()

model_mle(years = 2021:2023)

eval_mle()

#This next one takes some time
model_bayes()

diag_bayes()

estimate_bayes()

filter_bayes()