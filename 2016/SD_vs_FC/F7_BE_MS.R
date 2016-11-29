#---------------------------#
# Sequential sampling model #
# using hierarchical BE     #
# Kevin Potter              #
# 11/28/2016                #
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
modelFit = T

###
### Create inputs for model estimation
###

# Define version of model to fit (unless previously defined)
if ( !exists( 'type' ) ) type = 1

# Define priors

if ( type == 1 ) {
  
  Priors = cbind(
    c( rep(1,4), # kappa
       rep(c(2.5,1.4),each=4), # xi
       7 ), # theta (Prop. of tau)
    c( rep( 1, 4 ), # kappa
       rep( 1, 8 ), # Xi 
       2 ), # Theta (Prop. of tau)
    c( rep( 4, 4 ), # kappa
       rep( 4, 8 ), # xi 
       3 ), # theta (Prop. of tau)
    c( rep( 4, 4 ), # kappa
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
        kappa = matrix( runif( input$C[1]*input$Ns, .5, 2 ), 
                        input$Ns, input$C[1] ),
        xi = matrix( runif( input$C[2]*input$Ns, .5, 4 ), 
                     input$Ns, input$C[2] ),
        theta = matrix( runif( input$C[3]*input$Ns, .5, .9 ), 
                        input$Ns, input$C[3] ),
        kappa_mu = as.array( runif( input$C[1], .5, 2 ) ),
        kappa_sig = as.array( runif( input$C[1], .2, 2 ) ),
        xi_mu = as.array( runif( input$C[2], .5, 4 ) ),
        xi_sig = as.array( runif( input$C[2], .2, 2 ) ),
        theta_alpha = as.array( runif( input$C[3], 2, 15 ) ),
        theta_beta = as.array( runif( input$C[3], 2, 15 ) )
      )
      
    }
    
    return( out )
  }
  
}

# Define folder location and file name to save output
folderName = "C:/Users/Kevin/Documents/Posteriors from Stan/SD_vs_FC"
outName = paste("Posterior_estimates_",type,".RData",sep="")

if (modelFit) {
  
  # Extract data and covariates
  input = model_structure_create(type, d, Priors)
  
  warm = 100 # Warm-up period
  niter = 100 # Number of samples to approximate posterior per chain
  chains = 1 # Number of chains to run
  
  setwd('Stan_scripts')
  
  stan_seed = 3265
  startVal = init_f( input, chains )
  
  startTime = Sys.time() # To assess run-time
  fit = stan( file = 'WR_MS_fixed.stan', data = input, 
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
  # setwd( folderName )
  # save( post, conv, input, file = outName )
  
} else {
  # Load in posterior estimates
  setwd( folderName )
  load( outName )
}

setwd( orig_dir )

# Clean up workspace
rm( folderName, outName )