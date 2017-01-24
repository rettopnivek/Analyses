#--------------------#
# mSDT model check   #
# Kevin Potter       #
# Updated 01/19/2017 #
#--------------------#

# Clear workspace
rm( list = ls() )

# Save current directory
orig_dir = getwd()

# Define code segments to run
runCode = c( F, F, F )

# Index
# Lookup - 01:  Load in useful packages and data

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
### Parameter recovery (MLE)
###
# Lookup - 02

if ( runCode[1] ) {
  
  # Define generating parameters
  gp = c( .45, .4, 0, .98, .93 )
  
  # Define design
  No = 144
  IT = rep( 2:1, each = 72 )
  B = rep( c( rep(1,54),rep(2,18) ), 2 )
  Co = sample( rep(0:1,each=72) )
  
  # Define probability of picking right
  dp = gp[1:2];
  crt = gp[3]
  lmb = gp[4:5]
  
  # Define function to determine sum of the log-likelihoods
  mle_f = function( prm, dat, priors = NULL ) {
    
    # Extract data
    Y = dat[,1]
    IT = dat[,2]
    B = dat[,3]
    Co = dat[,4]
    
    # Extract and transform parameters
    dp = prm[1:2]
    crt = prm[3]
    lmb = logistic( prm[4:5] )
    
    # Calculate probability of picking right
    theta = mSDT_prob( dp[IT], crt, lmb[B], Co )
    
    # Calculate sum of the log-likelihood
    sll = sum( dbinom( Y, 1, theta, log = T ) )
    if ( is.na( sll ) ) sll = -Inf
    
    return( sll )
  }
  
  # Define function to generate random starting values
  st_f = function() {
    
    prm = runif( 5,
                 c(-1,-1,-.5,1,1),
                 c(1,1,.5,4,4) )
    
    return( prm )
  }
  
  # Define number of repetitions
  nRep = 20
  # Define matrix to store difference between estimates 
  # and generating parameters
  allPar = matrix( NA, nRep, 5 )
  diffsc = matrix( NA, nRep, 5 )
  
  # Create a progress bar using a base R function
  pb = txtProgressBar( min = 0, max = nRep, style = 3 )
  
  # Simulate data and conduct parameter recovery 
  # for each repetition
  for ( nr in 1:nRep ) {
    
    # Simulate responses
    theta = mSDT_prob( dp[IT], crt, lmb[B], Co )
    resp = rbinom( No, 1, theta )
    
    # Carry out parameter recovery
    dat = cbind( resp, IT, B, Co )
    res = MLE( dat, mle_f, st_f, nRep = 20, maxit = 20000 )
    bfp = c( res$param[1:3], logistic( res$param[4:5] ) )
    allPar[nr,] = bfp
    diffsc[nr,] = bfp - gp
    
    # Update the progress bar
    setTxtProgressBar(pb,nr)
  }
  close(pb)
  
  # Exclude abberant d' estimates
  excl = allPar[,1] > 5 | allPar[,2] > 5
  
  # Plot the residuals
  yl = lowerUpper( .2, as.vector( diffsc[!excl,] ) )
  x11()
  blankPlot( c(1,5), yl )
  abline( h = 0, lty = 2 )
  axis( 2, seq(yl[1],yl[2],.2) )
  mtext( 'Estimated - generating', side = 2, line = 2 )
  
  ui = apply( diffsc[ !excl, ], 2, quantile, 
              prob = c(.025,.16,.5,.84,.975) )
  segments( 1:5, ui[1,], 1:5, ui[5,] )
  segments( 1:5, ui[2,], 1:5, ui[4,], lwd = 3 )
  points( 1:5, ui[3,], pch = 21, bg = 'white', cex = 2 )
  
  # In conclusion, with binary data, it is hard to 
  # properly recover unbiased estimates of the model
  
}

###
### Parameter recovery via Bayesian estimation
###
# Lookup - 04

if ( runCode[2] ) {
  
  # Define generating parameters
  gp = c( 1.0, 1.0, -.2, .95, .7 )
  
  plotYes = F
  # plotYes = T
  if ( plotYes ) {
    
    pred = mSDT_prob( 
      gp[ rep( 1:2, 4 ) ],
      gp[ 3 ],
      gp[ rep( rep( 4:5, each = 2 ), 2 ) ], 
      rep( 0:1, each = 4 ) )
    
    x11()
    blankPlot( c(.8,2.2) )
    axis( 2, seq(0,1,.25) )
    mtext( 'P(Right)', side = 2, line = 2.25 )
    axis( 1, 1:2, c('Selective retrieval', 'Baseline' ), tick = F )
    abline( h = .5 )
    
    draw_figure( pred )
    
  }
  
  # Define design
  No = 144
  IT = rep( 2:1, each = 72 )
  B = rep( c( rep(1,54),rep(2,18) ), 2 )
  Co = sample( rep(0:1,each=72) )
  
  # Define probability of picking right
  dp = gp[1:2];
  crt = gp[3]
  lmb = gp[4:5]
  
  # Simulate responses
  theta = mSDT_prob( dp[B], crt, lmb[IT], Co )
  resp = rbinom( No, 1, theta )
  Cond = cbind( B, IT )
  
  stanDat = list(
    No = length( resp ),
    Cond = cbind( B, IT ),
    Y = resp,
    Co = Co
  )
  
  warm = 250 # Warm-up
  niter = 1250 # Number of samples to approximate posterior
  chains = 8 # Number of chains to run
  
  setwd('Stan_scripts')
  
  startTime = Sys.time() # To assess run-time
  sm = stan_model(stanc_ret = stanc_builder("mSDT_model_OS.stan"))
  fit = sampling( sm, data = stanDat, 
                  warmup = warm, iter = warm+niter, 
                  chains = chains )
  
  post = extract(fit)
  
  # Create numerical index of all levels of each variable
  cnd = cbind( rep( rep( 1:2, each = 2 ), 2 ),
               rep( 1:2, 4 ),
               rep( 0:1, each = 4 ) )
  ind = rep(1,No)
  for (i in 2:8) {
    ind[ IT == cnd[i,1] & B == cnd[i,2] & Co == cnd[i,3] ] = i
  }
  
  # Retrodictive checks
  f = function(x) {
    
    tmp = rbinom( No, 1, x )
    out = aggregate( tmp, list( B, IT, Co ), mean )$x
    
    return( out )
  }
  retCheck = apply( post$theta, 1, f )
  
  ui = apply( retCheck, 1, quantile, prob = c(.025,.975) )
  pred = rowMeans( retCheck )
  
  obs = aggregate( resp, list( B, IT, Co ), mean )
  colnames( obs ) = c('B','IT','Co','P')
  x11();
  plot( c(.8,2.2), c(0,1), type = 'n', xaxt = 'n',
        xlab = ' ', ylab = 'P(Right)' )
  abline( h = .5 )
  axis( 1, 1:2, c('Selective retrieval', 'Baseline'), tick =F )
  
  draw_figure( pred, clr = c('blue','white') )
  draw_figure( obs$P )
  
  setwd( orig_dir )
}

###
### Simulation with hierarchical model (w/ item/subject effects)
###
# Lookup - 05

if ( runCode[3] ) {
  
  # Define generating parameters
  beta = c( .18, -.22 )
  sigma_zeta = .25
  Omega = rbind( c(1, .3), c(.3, 1) )
  tau = c(.54,.54)
  phi = c(.9,.7)
  nu = c(20,20)
  
  # Define design
  Ni = 144 # Number of items
  Ns = 10 # Number of subjects
  No = Ni*Ns # Total number of observations
  Co = numeric(No) # Position of correct choice
  for (ns in 1:Ns) Co[ 1:Ni + Ni*(ns-1) ] = sample( rep(0:1,each=72) )
  # Condition indices
  IT = rep( 2:1, each = 72 ); IT = rep(IT,Ns)
  B = rep( c( rep(1,54),rep(2,18) ), 2 ); B = rep(B,Ns)
  Cnd = cbind( B, IT )
  
  subjIndex = rep( 1:Ns, each = Ni )
  itemIndex = numeric( No )
  for (ns in 1:Ns) 
    itemIndex[ 1:Ni + Ni*(ns-1) ] = sample( 1:Ni )
  
  # Simulate data using Stan
  
  setwd('Stan_scripts')
  
  # Create input for Stan
  stanDat = list(
    No = No,
    Ns = Ns,
    Ni = Ni,
    Co = Co,
    subjIndex = subjIndex,
    itemIndex = itemIndex,
    Cond = Cnd,
    beta = beta,
    Omega = Omega,
    tau = tau,
    sigma_zeta = sigma_zeta,
    phi = phi,
    nu = nu
  )
  
  startTime = Sys.time() # To assess run-time
  
  # Compile model
  sm = stan_model(stanc_ret = stanc_builder("mSDT_model_rng.stan"))
  
  # Draw samples
  fit = sampling( sm, data = stanDat, 
                  warmup = 0, iter = 1, 
                  chains = 1, algorithm = 'Fixed_param' )
  
  sim = extract(fit)
  resp = as.vector( sim$Y ) # Extract simulated responses
  
  runTime = Sys.time() - startTime
  print( runTime ); rm( startTime )
  
  # Parameter recovery in Stan
  
  # Create input for Stan
  stanDat = list(
    No = No,
    Ns = Ns,
    Ni = Ni,
    Y = resp, 
    Co = Co, 
    subjIndex = subjIndex,
    itemIndex = itemIndex,
    Cond = Cnd
  )
  
  warm = 250 # Warm-up
  niter = 1250 # Number of samples to approximate posterior
  chains = 8 # Number of chains to run
  
  startTime = Sys.time() # To assess run-time
  
  # Compile model
  sm = stan_model(stanc_ret = stanc_builder("mSDT_model.stan"))
  
  # Draw samples
  fit = sampling( sm, data = stanDat, 
                  warmup = warm, iter = warm+niter, 
                  chains = chains )
  
  post = extract(fit)
  
  runTime = Sys.time() - startTime
  print( runTime ); rm( startTime )
  
  plotYes = F
  # plotYes = T
  if ( plotYes ) {
    
    # Note that the model fails to converge
    x11()
    plot.ts( post$beta )
    
  }
  
}
