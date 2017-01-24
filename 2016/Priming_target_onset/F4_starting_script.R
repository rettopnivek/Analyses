#--------------------#
# Starting script    #
# Kevin Potter       #
# Updated 11/30/2016 #
#--------------------#

# Check if the variable 'type' already exists
if ( exists('type') ) {
  sel = NULL;
  sel = which( ls() == 'type' )
  # Clear workspace
  rm( list = ls()[-sel] )
} else {
  # Clear workspace
  rm( list = ls() )
}

# Save current directory
orig_dir = getwd()

###
### Load in useful packages
###

# For geting github packages
# install.packages(devtools)
# library(devtools)

# Integration of C++ and R
# install.packages(Rcpp)
library(Rcpp)

# Miscellanous functions for modeling/plotting
# install_github("rettopnivek/utilityf")
library(utilityf)

# Distribution functions for assorted sequential sampling models
# install_github("rettopnivek/seqmodels")
library(seqmodels)

# Functions for plotting response time and choice data
# install_github("rettopnivek/rtplots")
library(rtplots)

# Functions for preprocessing response time and choice data
# install_github("rettopnivek/rtclean")
library(rtclean)

# Bayesian estimation software
# install.packages(rstan)
library(rstan)
# To run chains in parallel
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# Define useful functions
source('F2_useful_functions.R')
