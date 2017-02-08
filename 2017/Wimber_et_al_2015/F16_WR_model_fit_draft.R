#--------------------#
# WR model fit       #
# Kevin Potter       #
# Updated 02/04/2017 #
#--------------------#

# Clear workspace
rm( list = ls() )

# Save current directory
orig_dir = getwd()

# Define code segments to run
runCode = c( T, F, F )

# Index
# Lookup - 01:  Load in useful packages and data
# Lookup - 02:  Pre-process data
# Lookup - 02:  Model 1

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
# Lookup - 03

if ( runCode[1] ) {
  
  setwd( 'Stan_scripts' )
  setwd( 'WR' )
  
  # Design matrices for threshold
  X_k_R = cbind( 1, d$Co )
  X_k_R[ d$Co == 0 , 2 ] = -1
  X_k_L = X_k_R
  X_k_L[,2] = -X_k_L[,2]
  
  # Design matrices for drift rate
  X_x_R = matrix( 0, nrow(d), 8 )
  X_x_L = X_x_R
  for ( i in 1:4 ) {
    
    tmp = d$Co;
    tmp[ tmp == 0 ] = -1
    
    X_x_R[ d$Cnd == i, 1:2 + 2*(i - 1) ] = 
      cbind( 1, tmp[ d$Cnd == i ] )
    X_x_L[ d$Cnd == i, 1:2 + 2*(i - 1) ] = 
      cbind( 1, -tmp[ d$Cnd == i ] )
    
  }
  
  # Use effects coding for correct position
  Co_ef = d$Co
  Co_ef[ Co_ef == 0 ] = -1
  
  # Define priors
  Priors = rbind(
    c( 2.0, .2 ), # beta_x
    c( 0.2, .2 ),
    c( 2.0, .2 ), 
    c( 0.2, .2 ),
    c( 2.0, .2 ), 
    c( 0.2, .2 ),
    c( 2.0, .2 ), 
    c( 0.2, .2 ),
    c( 3.0, .5 ), # beta_c
    c( 0.0, .2 ), 
    c( 0.7, .1 ), # phi
    c( 8, 2 ), # lambda
    c( 1, 1 ),  # sigma_zeta
    c( 1, 1 ),  # tau
    c( 2, 0 )   # L_Omega
  )
  
  # Create input for Stan
  stanDat = list(
    No = nrow( d ),
    Ns = length( unique( d$S ) ),
    Ni = length( unique( d$I ) ),
    Nx = ncol( X_x_L),
    Nk = ncol( X_k_L),
    Y = cbind( d$RT, d$Ch ), 
    Co = Co_ef, 
    subjIndex = d$S,
    itemIndex = d$I,
    X_x_R = X_x_R,
    X_x_L = X_x_L,
    X_k_R = X_k_R,
    X_k_L = X_k_L,
    min_RT = aggregate( d$RT, list( d$S ), min )$x, 
    Priors = Priors
  )
  
  warm = 10 # Warm-up
  niter = 10 # Number of samples to approximate posterior
  chains = 1 # Number of chains to run
  
  # Generating starting values
  st_val = c()
  tmp = matrix( runif( 9, -.3, .3 ), 3, 3 )
  diag( tmp ) = 1
  for ( i in 1:chains) {
    st_val = c( st_val, list( list(
      beta_x = runif( ncol(X_x_L),
                      c( 1, -.2, 1, -.2, 1, -.2, 1, -.2 ),
                      c( 3, .2, 3, .2, 3, .2, 3, .2 ) 
                      ),
      beta_k = runif( ncol(X_k_L), 
                      c(.8, -.2 ), 
                      c(2, .2 ) ),
      eta = matrix( runif( stanDat$Ns*3, -.2, .2 ), 
                    stanDat$Ns, 3 ),
      L_Omega = t(chol(tmp)),
      tau = runif( 3, .1, .5 ),
      zeta = runif( stanDat$Ni, -.2, .2 ),
      sigma_zeta = runif( 1, .1, .5 ),
      theta = runif( stanDat$Ns, .4, .9 ),
      phi = runif( 1, .4, .9 ),
      lambda = runif( 1, 6, 12 )
    ) ) )
  }
  
  # Check starting values
  stanDatcheck = c( stanDat, st_val[[1]] )
  sm = stan_model(stanc_ret = stanc_builder(
    "WR_items_model_likelihood_check.stan"))
  fit = sampling( sm, data = stanDatcheck, 
                  warmup = 0, iter = 1,
                  chains = 1, 
                  algorithm = "Fixed_param" )
  
  fitModel = F
  if (fitModel) {
    
    startTime = Sys.time() # To assess run-time
    
    # Compile model
    sm = stan_model(stanc_ret = stanc_builder("WR_items_model.stan"))
    
    # Draw samples
    fit = sampling( sm, data = stanDat, 
                    warmup = warm, iter = warm+niter,
                    init = st_val, 
                    chains = chains,
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