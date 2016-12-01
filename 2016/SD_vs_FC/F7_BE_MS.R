#---------------------------#
# Sequential sampling model #
# using hierarchical BE     #
# Kevin Potter              #
# 12/01/2016                #
#---------------------------#

# Initialize script
source('F4_starting_script.R')

# Load in data
setwd( 'Data' )
load( 'SD_v_FC.RData' )
setwd( orig_dir )

# For easy manipulation
d = allDat
colnames( d ) = c('S','Cnd','Ac','RT','E','PD',
                  'Ta','Co','PT','Ch','TL','FL','zTL','zFL')

# Indicate whether model estimation should be carried out 
# or if previous estimates should be loaded
modelFit = F
# Indicate if posterior samples should be saved
saveFit = F

# Index
# Lookup - 01:  Create inputs for model estimation
#          01a: Model 1
#          01b: Model 2
#          01c: Model 3
# Lookup - 02:  Estimate model parameters for single subject
# Lookup - 03:  Plot model results


###
### Create inputs for model estimation
###
# Lookup - 01

# Define version of model to fit (unless previously defined)
if ( !exists( 'type' ) ) type = 1

# Model 1
# Lookup - 01a
if ( type == 1 ) {
  
  Priors = cbind(
    c( c(1,1,0,0,0), # kappa
       rep( c( 2.5, 1.4 ), each = 2 ), # xi
       7 ), # theta (Prop. of tau)
    c( c(.25,.25,.1,.1,.1), # kappa
       rep( .5, 4 ), # Xi 
       2 ), # theta (Prop. of tau)
    c( rep( 4, 5 ), # kappa
       rep( 4, 4 ), # xi 
       3 ), # theta (Prop. of tau)
    c( rep( 4, 5 ), # kappa
       rep( 4, 4 ), # xi
       2 ) # theta (Prop. of tau)
  )
  
  # Define function to initialize values
  init_f = function( input, nCh ) {
    
    # Initialize output
    out = c()
    for (nc in 1:nCh) out = c( out, list() )
    
    # Generate initial values
    for (nc in 1:nCh) {
      
      out[[nc]] = list(
        kappa_S = runif( input$Ns, .5, 2 ),
        kappa_L = runif( input$Ns, .5, 2 ),
        kappa_SP = runif( input$Ns, -.2, .2 ),
        kappa_LP = runif( input$Ns, -.2, .2 ),
        kappa_B = runif( input$Ns, -.2, .2 ),
        xi = matrix( runif( 4*input$Ns, .5, 4 ), 
                     input$Ns, 4 ),
        theta = runif( input$Ns, .5, .9 ),
        kappa_mu = as.array( runif( 4, c(.5,.5,-.2,-.2),c(2,2,.2,.2) ) ),
        kappa_sig = as.array( runif( 5, .2, 2 ) ),
        xi_mu = as.array( runif( 4, .5, 4 ) ),
        xi_sig = as.array( runif( 4, .2, 2 ) ),
        theta_alpha = runif( 1, 2, 15 ),
        theta_beta = runif( 1, 2, 15 )
      )
      
    }
    
    return( out )
  }
  
}

# Model 2
# Lookup - 01b
if ( type == 2 ) {
  
  Priors = cbind(
    c( 1, # kappa
       rep( c( 2.5, 1.4 ), each = 4 ), # xi
       7 ), # theta (Prop. of tau)
    c( .25, # kappa
       rep( .5, 8 ), # Xi 
       2 ), # theta (Prop. of tau)
    c( rep( 4, 1 ), # kappa
       rep( 4, 8 ), # xi 
       3 ), # theta (Prop. of tau)
    c( rep( 4, 1 ), # kappa
       rep( 4, 8 ), # xi
       2 ) # theta (Prop. of tau)
  )
  
  # Define function to initialize values
  init_f = function( input, nCh ) {
    
    # Initialize output
    out = c()
    for (nc in 1:nCh) out = c( out, list() )
    
    # Generate initial values
    for (nc in 1:nCh) {
      
      out[[nc]] = list(
        kappa = runif( input$Ns, .5, 2 ),
        xi = matrix( runif( 8*input$Ns, .5, 4 ), 
                     input$Ns, 8 ),
        theta = runif( input$Ns, .5, .9 ),
        kappa_mu = runif( 1, .5, 2 ),
        kappa_sig = runif( 1, .2, 2 ),
        xi_mu = as.array( runif( 8, .5, 4 ) ),
        xi_sig = as.array( runif( 8, .2, 2 ) ),
        theta_alpha = runif( 1, 2, 15 ),
        theta_beta = runif( 1, 2, 15 )
      )
      
    }
    
    return( out )
  }
  
}

# Model 3
# Lookup - 01c
if ( type == 3 ) {
  
  Priors = cbind(
    c( c(1,1,0,0,0), # kappa
       rep( c( 2.5, 1.4 ), each = 4 ), # xi
       7 ), # theta (Prop. of tau)
    c( c(.25,.25,.1,.1,.1), # kappa
       rep( .5, 8 ), # Xi 
       2 ), # theta (Prop. of tau)
    c( rep( 4, 5 ), # kappa
       rep( 4, 8 ), # xi 
       3 ), # theta (Prop. of tau)
    c( rep( 4, 5 ), # kappa
       rep( 4, 8 ), # xi
       2 ) # theta (Prop. of tau)
  )
  
  # Define function to initialize values
  init_f = function( input, nCh ) {
    
    # Initialize output
    out = c()
    for (nc in 1:nCh) out = c( out, list() )
    
    # Generate initial values
    for (nc in 1:nCh) {
      
      out[[nc]] = list(
        kappa_S = runif( input$Ns, .5, 2 ),
        kappa_L = runif( input$Ns, .5, 2 ),
        kappa_SP = runif( input$Ns, -.2, .2 ),
        kappa_LP = runif( input$Ns, -.2, .2 ),
        kappa_B = runif( input$Ns, -.2, .2 ),
        xi = matrix( runif( 8*input$Ns, .5, 4 ), 
                     input$Ns, 8 ),
        theta = runif( input$Ns, .5, .9 ),
        kappa_mu = as.array( runif( 3, c(.5,.5,-.2),c(2,2,.2) ) ),
        kappa_sig = as.array( runif( 5, .2, 2 ) ),
        xi_mu = as.array( runif( 8, .5, 4 ) ),
        xi_sig = as.array( runif( 8, .2, 2 ) ),
        theta_alpha = runif( 1, 2, 15 ),
        theta_beta = runif( 1, 2, 15 )
      )
      
    }
    
    return( out )
  }
  
}

# Define folder location and file name to save output
folderName = "C:/Users/Kevin/Documents/Posteriors from Stan/SD_vs_FC"
outName = paste("Posterior_estimates_",type,".RData",sep="")

###
### Carry out model estimation
###

if (modelFit) {
  
  # Extract data and covariates
  input = model_structure_create(type, d, Priors)
  
  warm = 300 # Warm-up period
  niter = 1250 # Number of samples to approximate posterior per chain
  chains = 8 # Number of chains to run
  
  setwd('Stan_scripts')
  
  if ( type == 1 ) stan_seed = 2032
  if ( type == 2 ) stan_seed = 7654
  if ( type == 3 ) stan_seed = 4281
  
  startVal = init_f( input, chains )
  
  startTime = Sys.time() # To assess run-time
  fit = stan( file = input$mName, data = input[1:16], 
              warmup = warm, iter = warm+niter, 
              chains = chains,
              init = startVal,
              control = list( adapt_delta = .9 ),
              seed = stan_seed )
  
  post = extract(fit)
  # Report run time
  runTime = Sys.time() - startTime
  print( runTime )
  rm( startTime )
  
  # Extract convergence diagnostics
  conv = convergence_extract( fit )
  
  # Save posterior estimates
  setwd( folderName )
  if (saveFit) save( post, conv, input, file = outName )
  
} else {
  # Load in posterior estimates
  setwd( folderName )
  load( outName )
}

setwd( orig_dir )

# Clean up workspace
rm( folderName, outName )