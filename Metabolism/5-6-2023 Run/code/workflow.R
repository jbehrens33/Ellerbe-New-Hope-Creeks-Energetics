# Step 1) this will source in all functions (except this workflow) into the Global Env ----
scripts <- list.files('code/', 
                      pattern = '*.R', 
                      full.names = TRUE)
# exclude this script!
scripts <- scripts[!scripts %in% 'code/workflow.R']

# source them in 
sapply(scripts, source, .GlobalEnv)

# Step 2) Before you model  ----
## this fun installs streamMetabolizer (and dependencies) and rstan (and dependencies) from github
prep_bayes_model()

# if this returns an error, running this below will suffice
# install.packages("StanHeaders", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
# install.packages("rstan", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
# remotes::install_github("https://github.com/appling/unitted")
# remotes::install_github("https://github.com/USGS-R/streamMetabolizer")
library('rstan')
library('StanHeaders')
library('unitted')
library('streamMetabolizer')
library(StreamPULSE)

#install.packages(c('cli', 'fansi', 'vctrs'))
lapply(c('cli', 'fansi', 'vctrs'), library, character.only = TRUE)

# Step 3) Find your sites ----
# this will show you all the stream pulse sites for a given region, in this case North Carolina
sp_sites <- StreamPULSE::query_available_data(region = 'NC') 

# you can define your sites from sp_sites or define your sites if you know them
my_sites <- c('NC_NHC')

# Step 4) Get the data ----
# get and save the data from streampulse
get_sp_data(sitecode = 'NC_NHC')


# format for modeling
## this is a wrapper around StreamPULSE::prep_metabolism
## there are more arguments to modify within this function, check them out ?prep_metabolism
prep_sp_data()

# Step 5) model and evaluate using MLE ----
# the user will have to define which years (calendar) they are interested in
model_mle(years = 2021:2023)

# evaluate the MLE outputs
eval_mle()

# Step 6) model in Bayes ----
#model_bayes(years = 2021:2023)

model_bayes_nhc(years = 2021:2023)

# Step 7) Prepare the usable outputs ----
# Find bad data in annual outputs
diag_bayes()

# compile outputs
estimate_bayes()

# filter out bad data from outputs
filter_bayes()
