#--------------------#
# SDT model check    #
# Kevin Potter       #
# Updated 01/19/2017 #
#--------------------#

# Clear workspace
rm( list = ls() )

# Save current directory
orig_dir = getwd()

# Define code segments to run
runCode = c( F, T )

# Index
# Lookup - 01:  Load in useful packages and data
# Lookup - 02:  Simulation example (One subject)
# Lookup - 03:  Simulation example (Hierarchical)

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
### Simulation example (One subject)
###
# Lookup - 02

if ( runCode[1] ) {
  
  # Define generating parameters
  gp = c( 1.2, .8, .8, 1.2, 0 )
  
  # Define design
  No = 144
  IT = rep( 2:1, each = 72 )
  B = rep( c( rep(1,54),rep(2,18) ), 2 )
  Co = sample( rep(0:1,each=72) )
  Cnd = c( rep(3,54), rep(4,18), rep(1,54), rep(2,18) )
  
  # Define probability of picking right
  dp = gp[1:4]
  crt = gp[5]
  
  # Simulate responses
  theta = SDT_prob( dp[Cnd], crt, Co )
  resp = rbinom( No, 1, theta )
  Cond = cbind( B, IT )
  
  # Carry out parameter recovery using Stan
  
  setwd('Stan_scripts')
  setwd('SDT')
  
  # Create input for Stan
  stanDat = list(
    No = length( resp ),
    Cond = Cnd,
    Y = resp,
    Co = Co
  )
  
  warm = 250 # Warm-up
  niter = 1250 # Number of samples to approximate posterior
  chains = 8 # Number of chains to run
  
  startTime = Sys.time() # To assess run-time
  
  # Compile model
  sm = stan_model(stanc_ret = stanc_builder("SDT_model_OS.stan"))
  
  # Draw samples
  fit = sampling( sm, data = stanDat, 
                  warmup = warm, iter = warm+niter, 
                  chains = chains )
  
  post = extract(fit)
  
  runTime = Sys.time() - startTime
  print( runTime ); rm( startTime )
  
  # Retrodictive checks
  pst = cbind( post$dp, post$crt )
  
  f = function( x ) {
    theta_hat = SDT_prob( x[Cnd], x[5], Co )
    yhat = rbinom( No, 1, theta_hat )
    out = aggregate( yhat, list( B, IT, Co ), mean )$x
    return(out)
  }
  rc = apply( pst, 1, f )
  
  plotYes = F
  # plotYes = T
  if ( plotYes ) {
    
    x11()
    plot( c(.8,2.2), c(0,1), type = 'n', xlab = 'Condition',
          ylab = 'P(Right)', xaxt = 'n', bty = 'n' )
    axis( 1, 1:2, c('Selective retrieval','Baseline'), tick = F )
    abline(h=.5)
    legend('topright','Right correct',bty='n')
    legend('bottomright','Left correct',bty='n')
    legend( 'topleft', c('Target','Competitor'),
            pch = c(19,22), pt.bg = 'white', bty = 'n' )
    
    obs = aggregate( resp, list( B, IT, Co ), mean )$x
    draw_figure( obs )
    
    pred = SDT_prob( c(gp[1:4],gp[1:4]), gp[5], rep(0:1,each=4) )
    draw_figure( pred, clr = c('blue','white') )
    
    est = rowMeans(rc)
    draw_figure( est, clr = c('red','white'), lnTp = 2  )
    
  }
  
  setwd( orig_dir )
}

###
### Simulation example (Hierarchical)
###
# Lookup - 03

if ( runCode[2] ) {
  
  # Define generating parameters
  beta_dp = c( 1.2, .8, .8, 1.2 )
  beta_c = 0.2
  sigma_zeta = .25
  Omega = rbind( c(1, .3), c(.3, 1) )
  tau = c(.71,.71)
  
  # Define design
  Ni = 144 # Number of items
  Ns = 10 # Number of subjects
  No = Ni*Ns # Total number of observations
  Nd = 4; # Total number of conditions for d'
  Nc = 1; # Total number of conditions for criterion
  Co = numeric(No) # Position of correct choice
  for (ns in 1:Ns) Co[ 1:Ni + Ni*(ns-1) ] = sample( rep(0:1,each=72) )
  # Condition indices
  IT = rep( 2:1, each = 72 ); IT = rep(IT,Ns)
  B = rep( c( rep(1,54),rep(2,18) ), 2 ); B = rep(B,Ns)
  Cnd = c( rep(3,54), rep(4,18), rep(1,54), rep(2,18) )
  Cnd = rep(Cnd,Ns)
  
  # Define design matrices
  Xd = matrix( 0, No, 4 );
  for (i in 1:4) Xd[ Cnd == i, i ] = 1
  Xc = matrix( 1, No, 1 );
  
  subjIndex = rep( 1:Ns, each = Ni )
  itemIndex = numeric( No )
  for (ns in 1:Ns) 
    itemIndex[ 1:Ni + Ni*(ns-1) ] = sample( 1:Ni )
  
  # Simulate data using Stan
  
  setwd('Stan_scripts')
  setwd('SDT')
  
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
    sigma_zeta = sigma_zeta
  )
  
  startTime = Sys.time() # To assess run-time
  
  # Compile model
  sm = stan_model(stanc_ret = stanc_builder("SDT_model_rng.stan"))
  
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
    c( 0, 1 ),  # beta_dp
    c( 0, 1 ),
    c( 0, 1 ),
    c( 0, 1 ),
    c( 0, 1 ),  # beta_c
    c( 2, 4 ),  # sigma_zeta
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
  sm = stan_model(stanc_ret = stanc_builder("SDT_model.stan"))
  
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
    
    x11(width=12);
    layout( matrix( 1:6, 2, 3, byrow = T ) )
    
    for ( i in 1:4 ) {
      hist( post$beta_dp[,i], breaks = 40, col = 'grey', 
            border = 'white', freq = F, 
            main = paste( "d'[", i, ']', sep = '' ),
            xlab = ' ' )
      abline( v = beta_dp[i], col = 'blue' )
    }
    hist( as.numeric( post$beta_c[,1] ), breaks = 40, col = 'grey', 
          border = 'white', freq = F, 
          main = 'Bias',
          xlab = ' ' )
    abline( v = beta_c, col = 'blue' )
    
  }
  
  setwd( orig_dir )
}
