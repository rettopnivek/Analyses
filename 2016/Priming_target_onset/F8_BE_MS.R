#---------------------------#
# Sequential sampling model #
# using hierarchical BE     #
# Kevin Potter              #
# 11/21/2016                #
#---------------------------#

# Initialize script
source('F5_starting_script.R')

# Load in data
setwd( 'Data' )
load( 'Priming_offset.RData' )
setwd( orig_dir )

# For easy manipulation
d = allDat
colnames( d ) = c('S','PD','PT','O','Co','RT',
                  'Ch','Ac','OT','DT','CDT','Cnd')

# Indicate whether model estimation should be carried out 
# or if previous estimates should be loaded
modelFit = F

###
### Create inputs for model estimation
###

# Define version of model to fit (unless previously defined)
if ( !exists( 'type' ) ) type = 1

# Define priors

# Previous drift rate values
#     T    F
# STP 3.6  1.1
# LTP 3.3  1.4
# SFP 2.0  2.6
# LFP 2.9  1.5

if ( type == 1 ) {
  Priors = cbind(
    c( c(1.2,1.2), # Kappa
       c( 3.6, 3.3, 2.0, 2.9, 1.1, 1.4, 2.6, 1.5 ), # xi
       rep( 8, 10 ) ), # Theta (Prop. of tau)
    c( rep( 2, 2 ), # Kappa
       rep( 2, 8 ), # Xi 
       rep( 2, 10 ) ), # Theta (Prop. of tau)
    c( rep( 5, 2 ), # Kappa
       rep( 5, 8 ), # Xi 
       rep( 2, 10 ) ), # Theta (Prop. of tau)
    c( rep( 8, 2 ), # Kappa
       rep( 8, 8 ), # Xi 
       rep( 2, 10 ) ) # Theta (Prop. of tau)
  )
}

# Define folder location and file name to save output
folderName = "C:/Users/Kevin/Documents/Posteriors from Stan/Priming_target_onset"
outName = paste("Posterior_estimates_",type,".RData",sep="")

if (modelFit) {
  
  # Extract data and covariates
  input = model_structure_create(type, d, Priors)
  
  warm = 625 # Warm-up period
  niter = 625 # Number of samples to approximate posterior per chain
  chains = 8 # Number of chains to run
  
  setwd('Stan_scripts')
  
  stan_seed = 8195
  
  startTime = Sys.time() # To assess run-time
  fit = stan( file = 'WR_MS.stan', data = input, 
              warmup = warm, iter = warm+niter, 
              chains = chains,
              control = list( adapt_delta = .9 ),
              seed = stan_seed, init_r = .5 )
  
  post = extract(fit)
  # Report run time
  runTime = Sys.time() - startTime
  print( runTime )
  rm( startTime )
  
  # Extract convergence diagnostics
  conv = convergence_extract( fit )
  
  # Save posterior estimates
  setwd( folderName )
  save( post, conv, input, file = outName )
  
} else {
  # Load in posterior estimates
  setwd( folderName )
  load( outName )
}

setwd( orig_dir )

# Clean up workspace
rm( folderName, outName )