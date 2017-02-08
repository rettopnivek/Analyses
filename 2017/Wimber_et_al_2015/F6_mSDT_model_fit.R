#--------------------#
# mSDT model fit     #
# Kevin Potter       #
# Updated 01/24/2017 #
#--------------------#

# Clear workspace
rm( list = ls() )

# Save current directory
orig_dir = getwd()

# Define code segments to run
runCode = c( T, F, F )

# Index
# Lookup - 01:  Load in useful packages and data
# Lookup - 02:  Model 1 (4 d', one mixture probability)

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
### Fit model using Stan
###

if ( runCode[1] ) {
  
  setwd( 'Stan_scripts' )
  setwd( 'mSDT' )
  
  # Create design matrices
  Xd = matrix( 0, nrow( d ), 4 )
  for ( i in 1:4 ) Xd[ d$Cnd == i, i ] = 1
  Xc = cbind( rep( 1, nrow(d) ) )
  Xm = Xc
  
  # Define priors
  Priors = rbind(
    c( 1.537, 1 ),  # beta_d
    c( 1.537, 1 ),
    c( 1.537, 1 ),
    c( 1.537, 1 ),
    c( 0, 1 ),  # beta_c
    c( logit( .9 ), .3 ),  # beta_m
    c( 2, 4 ),  # sigma_lambda
    c( 2, 4 ),  # tau
    c( 2, 0 )   # L_Omega
  )
  
  # Create input for Stan
  stanDat = list(
    No = nrow( d ),
    Ns = length( unique( d$S ) ),
    Nd = 4,
    Nc = 1,
    Nm = 1, 
    Y = d$Ch, 
    Co = d$Co, 
    subjIndex = d$S,
    Xd = Xd,
    Xc = Xc,
    Xm = Xm, 
    Priors = Priors
  )
  
  warm = 750 # Warm-up
  niter = 1250*3 # Number of samples to approximate posterior
  chains = 8 # Number of chains to run
  
  startTime = Sys.time() # To assess run-time
  
  # Compile model
  sm = stan_model(stanc_ret = stanc_builder("mSDT_model_no_item.stan"))
  
  # Draw samples
  fit = sampling( sm, data = stanDat, 
                  warmup = warm, iter = warm+niter, 
                  chains = chains,
                  thin = 3, 
                  control = list( adapt_delta = .995,
                                  max_treedepth = 14 ) )
  
  post = extract(fit)
  runTime = Sys.time() - startTime # To assess run-time
  
  setwd( orig_dir )
  
  # Save results to a .RData file
  fName = "C:/Users/Kevin/Documents/Posteriors from Stan"
  setwd(fName)
  setwd( "Wimber_et_al" )
  save( stanDat, fit, post, file = 'mSDT_model_M1_post.RData' )
  setwd( orig_dir )
  
}

