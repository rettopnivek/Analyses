#--------------------#
# Starting script    #
# Kevin Potter       #
# Updated 07/13/2017 #
#--------------------#

# Check if the variable 'model_type' already exists
if ( exists('model_type') ) {
  sel = NULL;
  sel = which( ls() == 'model_type' )
  # Clear workspace
  rm( list = ls()[-sel] )
} else {
  # Clear workspace
  rm( list = ls() )
}

# Save current directory
orig_dir = getwd()

# Index
# Lookup - 01:  Load in useful packages and functions
# Lookup - 02:  Load in data

###
### Load in useful packages and functions
###
# Lookup - 01

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

# Convenience functions for MLE
# install_github("rettopnivek/mle")
library(mle)

# Functions for numerical derivation
# install.packages( 'numDeriv' )
library( numDeriv )

# Optimization routines
# install.packages( 'optimx' )
library(optimx)

# Define useful functions
source('F2_model_construction.R')

# Define useful color palette
# http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", 
                "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

###
### Load in data
###
# Lookup - 02

# Navigate to folder with data
setwd( 'Data' )
load( 'SD_v_FC.RData' ) # Load in data
setwd( orig_dir )

# Copy data-frame and relabel columns for ease of
# navigation
d = allDat

if ( !exists( 'nROUSE_res' ) ) {
  colnames( d ) = c( "S", "Cnd", "Ac", "RT", "excl", "PD", 
                   "PDL", "Ta", "TaL", "Co", "CoL", "PT", 
                   "PTL", "Ch", "ChL" )
} else {
  colnames( d ) = c( "S", "Cnd", "Ac", "RT", "excl", "PD", 
                     "PDL", "Ta", "TaL", "Co", "CoL", "PT", 
                     "PTL", "Ch", "ChL", "TL", "FL", "zTL", "zFL" )
}
