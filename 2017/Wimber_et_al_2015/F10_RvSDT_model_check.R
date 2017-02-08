#--------------------#
# RvSDT model check  #
# Kevin Potter       #
# Updated 01/26/2017 #
#--------------------#

# Clear workspace
rm( list = ls() )

# Save current directory
orig_dir = getwd()

# Define code segments to run
runCode = c( T )

# Index
# Lookup - 01:  Load in useful packages and data
# Lookup - 02:  Simulation example (Hierarchical)

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
### Simulation example (Hierarchical)
###
# Lookup - 02

if ( runCode[1] ) {
  
  # Define design
  Ni = 144 # Number of items
  Ns = 10 # Number of subjects
  No = Ni*Ns # Total number of observations
  Nd = 2; # Total number of conditions for d'
  Nc = 1; # Total number of conditions for criterion
  Co = numeric(No) # Position of correct choice
  for (ns in 1:Ns) Co[ 1:Ni + Ni*(ns-1) ] = sample( rep(0:1,each=72) )
  # Condition indices
  IT = rep( 2:1, each = 72 ); IT = rep(IT,Ns)
  B = rep( c( rep(1,54),rep(2,18) ), 2 ); B = rep(B,Ns)
  Cnd = c( rep(3,54), rep(4,18), rep(1,54), rep(2,18) )
  Cnd = rep(Cnd,Ns)
  
  # Define generating parameters
  beta_dp = c( 1.537, -.4 )
  beta_c = 0.1
  Omega = rbind( c(1, .3), c(.3, 1) )
  tau = c(.71,.71)
  mu_lambda = logit(.2)
  sigma_lambda = .25
  
  # Define design matrices
  Xd = matrix( 0, No, 2 );
  Xd[ Cnd == 3, 2 ] = 1
  Xc = matrix( 1, No, 1 );
  
  subjIndex = rep( 1:Ns, each = Ni )
  itemIndex = numeric( No )
  for (ns in 1:Ns) 
    itemIndex[ 1:Ni + Ni*(ns-1) ] = sample( 1:Ni )
  
  # Simulate data using Stan
  
  setwd('Stan_scripts')
  setwd('RvSDT')
  
  # Create input for Stan
  stanDat = list(
    No = No,
    Ns = Ns,
    Ni = Ni,
    Nd = Nd,
    Nc = Nc,
    Co = Co,
    subjIndex = subjIndex,
    itemIndex = itemIndex,
    Xd = Xd,
    Xc = Xc,
    beta_dp = beta_dp,
    beta_c = array( beta_c, dim = 1 ), 
    Omega = Omega,
    tau = tau,
    mu_lambda = mu_lambda,
    sigma_lambda = sigma_lambda
  )
  
  startTime = Sys.time() # To assess run-time
  
  # Compile model
  sm = stan_model(stanc_ret = stanc_builder("RvSDT_model_rng.stan"))
  
  # Draw samples
  fit = sampling( sm, data = stanDat, 
                  warmup = 0, iter = 1, 
                  chains = 1, algorithm = 'Fixed_param' )
  
  sim = extract(fit)
  resp = as.vector( sim$Y ) # Extract simulated responses
  
  runTime = Sys.time() - startTime
  print( runTime ); rm( startTime )
  
  # Parameter recovery in Stan
  
  # Define priors
  Priors = rbind(
    c( 1.537, 1 ),  # beta_dp
    c( 0, 1 ),
    c( 0, 1 ),  # beta_c
    c( -1.386, 1 ), # mu_lambda
    c( 2, 4 ),  # sigma_lambda
    c( 2, 4 ),  # tau
    c( 2, 0 )   # L_Omega
  )
  
  # Create input for Stan
  stanDat = list(
    No = No,
    Ns = Ns,
    Ni = Ni,
    Nd = Nd,
    Nc = Nc,
    Y = resp, 
    Co = Co, 
    subjIndex = subjIndex,
    itemIndex = itemIndex,
    Xd = Xd,
    Xc = Xc,
    Priors = Priors
  )
  
  warm = 250 # Warm-up
  niter = 1250 # Number of samples to approximate posterior
  chains = 8 # Number of chains to run
  
  startTime = Sys.time() # To assess run-time
  
  # Compile model
  sm = stan_model(stanc_ret = stanc_builder("RvSDT_model.stan"))
  
  # Draw samples
  # Draw samples
  fit = sampling( sm, data = stanDat, 
                  warmup = warm, iter = warm+niter, 
                  chains = chains,
                  control = list( adapt_delta = .995,
                                  max_treedepth = 14 ) )
  
  post = extract(fit)
  
  runTime = Sys.time() - startTime
  print( runTime ); rm( startTime )
  
}