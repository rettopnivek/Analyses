#--------------------#
# WR model check     #
# Kevin Potter       #
# Updated 02/04/2017 #
#--------------------#

# Clear workspace
rm( list = ls() )

# Save current directory
orig_dir = getwd()

# Define code segments to run
runCode = c( T, F )

# Index
# Lookup - 01:  Load in useful packages and data
# Lookup - 02:  Pre-process data
# Lookup - 03:  Model 1

###
### Load in useful packages and data
###
# Lookup - 01

# For geting github packages
# install.packages(devtools)
# library(devtools)

# Miscellanous functions for modeling/plotting
# install_github("rettopnivek/utilityf")
library(utilityf)

# Miscellanous functions for modeling/plotting
# install_github("rettopnivek/rtplots")
library(rtplots)

# Package for sequential sampling models
# install_github("rettopnivek/seqmodels")
library(seqmodels)

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
### Pre-process data 
###
# Lookup - 02

# Remove missing responses
sel = which( is.na( d$RT ) )
d = d[ -sel, ]

# Determine sample size
N = length( unique( d$S ) )

# Remove overly fast responses
sel = which( d$RT < .2 )
d = d[ -sel, ]

###
### Model 1
###
# Lookup - 01

if ( runCode[1] ) {
  
  setwd( 'Stan_scripts' )
  setwd( 'WR' )
  
  # Define matrix of condition indices
  condIndex = cbind( 1, 2, 
                     d$Cnd*d$Co + (d$Cnd+4)*(1-d$Co),
                     (d$Cnd+4)*d$Co + d$Cnd*(1-d$Co) )
  
  # Define priors
  Priors = rbind(
    c( log(2.7), .4 ), # Drift rate
    c( log(2.7), .4 ),
    c( log(2.7), .4 ),
    c( log(2.7), .4 ),
    c( log(1.8), .4 ),
    c( log(1.8), .4 ),
    c( log(1.8), .4 ),
    c( log(1.8), .4 ),
    c( log(3.4), .4 ), # Thresholds
    c( log(3.4), .4 ), 
    c( 0.7, .1 ), # phi
    c( 8, 2 ), # lambda
    c( 1, 1 ),  # tau
    c( 2, 0 )   # L_Omega
  )
  
  # Create input for Stan
  stanDat = list(
    No = nrow( d ),
    Ns = length( unique( d$S ) ),
    Ni = length( unique( d$I ) ),
    Nx = 8,
    Nk = 2,
    Y = cbind( d$RT, d$Ch ), 
    subjIndex = d$S,
    condIndex = condIndex,
    min_RT = aggregate( d$RT, list( d$S ), min )$x, 
    Priors = Priors
  )
  
  warm = 500 # Warm-up
  niter = 1250 # Number of samples to approximate posterior
  chains = 8 # Number of chains to run
  
  # Generate starting values
  st_val = c()
  tmp = matrix( runif( 100, -.3, .3 ), 10, 10 )
  diag( tmp ) = 1
  while( det( tmp ) < 0 ) {
    tmp = matrix( runif( 100, -.3, .3 ), 10, 10 )
    diag( tmp ) = 1
  }
  for ( i in 1:chains) {
    st_val = c( st_val, list( list(
      xi = matrix( runif( stanDat$Ns * stanDat$Nx,
                          1, 3 ), 
        stanDat$Ns, stanDat$Nx ),
      kappa = matrix( runif( stanDat$Ns * stanDat$Nk,
                          1, 3 ), 
                   stanDat$Ns, stanDat$Nk ),
      theta = runif( stanDat$Ns, .4, .9 ),
      Mu = runif( (stanDat$Nx + stanDat$Nk ), 1, 3 ),
      Tau = runif( (stanDat$Nx + stanDat$Nk ), .1, .5 ),
      L_Omega = t(chol(tmp)),
      phi = runif( 1, .4, .9 ),
      lambda = runif( 1, 6, 12 )
    ) ) )
  }
  
  fitModel = T
  
  if (fitModel) {
    
    startTime = Sys.time() # To assess run-time
    
    # Compile model
    sm = stan_model(stanc_ret = stanc_builder("WR_model.stan"))
    
    # Draw samples
    fit = sampling( sm, data = stanDat, 
                    warmup = warm, iter = warm+niter,
                    chains = chains,
                    seed = 19493, 
                    init = st_val, 
                    control = list( adapt_delta = .95,
                                    max_treedepth = 14 ) )
    
    post = extract(fit)
    runTime = Sys.time() - startTime # To assess run-time
    print( runTime )
    
    setwd( orig_dir )
    
    # Save results to a .RData file
    fName = "C:/Users/Kevin/Documents/Posteriors from Stan"
    setwd(fName)
    setwd( "Wimber_et_al" )
    save( stanDat, fit, post, file = 'WR_model_M1_post.RData' )
    setwd( orig_dir )
    
  }
  
}