#--------------------#
# Starting script    #
# Kevin Potter       #
# Updated 11/21/2016 #
#--------------------#

# Clear workspace
rm( list = ls() )

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

###
### Create scripts for Stan
###

createScriptOS = F
if ( createScriptOS ) {
  
  setwd('Stan_scripts')
  fileName = "Wald_race_stan_functions.txt"
  wr_f = readChar(fileName, file.info(fileName)$size)
  fileName = "WR_one_subject_ex.txt"
  os = readChar(fileName, file.info(fileName)$size)
  model_script = paste( wr_f, os, sep = "" )
  writeChar( model_script, "WR_OS.stan" )
  
  # Clean up workspace
  rm( wr_f, os, model_script )
  
}

createScriptMS = F
if ( createScriptMS ) {
  
  setwd('Stan_scripts')
  fileName = "Wald_race_stan_functions.txt"
  wr_f = readChar(fileName, file.info(fileName)$size)
  fileName = "WR_multi_subject.txt"
  ms = readChar(fileName, file.info(fileName)$size)
  model_script = paste( wr_f, ms, sep = "" )
  writeChar( model_script, "WR_MS.stan" )
  
  # Clean up workspace
  rm( wr_f, ms, model_script )
  
}

# Clean up workspace
rm( createScriptMS, createScriptOS )

setwd( orig_dir )
