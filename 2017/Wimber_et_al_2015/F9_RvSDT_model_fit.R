#--------------------#
# RvSDT model fit    #
# Kevin Potter       #
# Updated 01/20/2017 #
#--------------------#

# Clear workspace
rm( list = ls() )

# Save current directory
orig_dir = getwd()

# Define code segments to run
runCode = c( T )

###
### Load in useful packages and data
###

# For geting github packages
# install.packages(devtools)
# library(devtools)

# Miscellanous functions for modeling/plotting
# install_github("rettopnivek/utilityf")
library(utilityf)

# Load in package for Bayesian estimation
library(rstan)
# For parallel processing
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# Load in useful functions
source( 'F1_useful_functions.R' )

# Load in data
setwd( 'Data' )
load( 'Recog_mem.RData' )
setwd( orig_dir )

###
### Recall versus comparison process model
###

if ( runCode[1] ) {
  
  setwd( 'Stan_scripts' )
  setwd( 'RvSDT' )
  
  # Create design matrices
  Xd = matrix( 1, nrow( d ), 2 )
  Xd[ d$Cnd != 3, 2 ] = 0
  Xc = cbind( rep( 1, nrow(d) ) )
  
  # Define priors
  Priors = rbind(
    c( 1.537, 1 ),  # beta_dp
    c( 0, 1 ),
    c( 0, 1 ),  # beta_c
    c( 1, 1 ),  # phi
    c( 8, 2 ),  # nu
    c( 2, 4 ),  # tau
    c( 2, 0 )   # L_Omega
  )
  
  # Create input for Stan
  stanDat = list(
    No = nrow( d ),
    Ns = length( unique( d$S ) ),
    Ni = length( unique( d$I ) ),
    Nd = 2,
    Nc = 1,
    Y = d$Ch, 
    Co = d$Co, 
    subjIndex = d$S,
    itemIndex = d$I,
    Xd = Xd,
    Xc = Xc,
    Priors = Priors
  )
  
  warm = 500 # Warm-up
  niter = 1250 # Number of samples to approximate posterior
  chains = 8 # Number of chains to run
  
  startTime = Sys.time() # To assess run-time
  
  # Compile model
  sm = stan_model(stanc_ret = stanc_builder("RvSDT_model.stan"))
  
  # Draw samples
  fit = sampling( sm, data = stanDat, 
                  warmup = warm, iter = warm+niter, 
                  chains = chains,
                  control = list( adapt_delta = .995,
                                  max_treedepth = 14 ) )
  
  post = extract(fit)
  runTime = Sys.time() - startTime # To assess run-time
  
  setwd( orig_dir )
  
  # Save results to a .RData file
  fName = "C:/Users/Kevin/Documents/Posteriors from Stan"
  setwd(fName)
  setwd( "Wimber_et_al" )
  save( stanDat, fit, post, file = 'RvSDT_model_post.RData' )
  setwd( orig_dir )
  
}