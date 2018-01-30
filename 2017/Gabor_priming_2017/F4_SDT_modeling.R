#-----------------------------#
# Script to create R data set #
# Kevin Potter                #
# Updated 12/16/2017          #
#-----------------------------#

# Clear workspace
rm(list = ls())

# Save current working directory
orig_dir = getwd()

# Indicator if data from Experiment 1 or 2 should be loaded in
Exp = 1

# Indicate whether to fit models or load 
# in saved posterior samples
modelFit = F


# Indicate which code segments to run
runCode = c(
  F,  # Basic example for SDT
  F,  # Basic example for hierarchical multivariate normal
  F,  # Basic example for hierarchical SDT
  T
)

# Index
# Lookup - 01:  Load in data, useful packages, and functions
# Lookup - 02:  Parameter recovery for base SDT model
# Lookup - 03:  Example of fitting multivariate normal hierarchy
# Lookup - 04:  Parameter recovery for hierarchical SDT model
# Lookup - 05:  Model fit to experimental data

###
### Load in useful packages, functions, and data
###
# Lookup - 01

### Load in useful packages

# For geting github packages
# install.packages(devtools)
# library(devtools)

# Miscellanous functions for modeling/plotting
# install_github("rettopnivek/utilityf")
library(utilityf)

# Load in package for multivariate normal
# install.packages( 'MASS' )
library(MASS)

# Load in package for C++ interface
# install.packages( 'Rcpp' )
library( Rcpp )

# Load in package for Bayesian modeling
# install.packages( 'rstan' )
library(rstan)

### Define useful functions

source( 'F1_useful_functions.R' )

setwd( 'Stan_scripts/SDT' )
sourceCpp( 'SDT_functions.cpp' )
setwd( orig_dir )

quickUI = function( pos, x, coverage, 
                    flip = T, clr = c('grey','white'), lnSz = 2 ) {
  
  i = quickCoverage( coverage[1] )
  
  lmnts = quantile( x, i )
  errorBars( pos, lmnts, flip = flip, col = clr[1], length = .05,
             lwd = lnSz )
  
  bx = matrix( NA, 2, 4 )
  bx[1,] = pos + c(-.1,-.1,.1,.1)
  i = quickCoverage( coverage[2] )
  bx[2,] = rep( quantile( x, i ), 2 )
  bx[2,3:4] = rev( bx[2,3:4] )
  
  if ( flip ) {
    polygon( bx[2,], bx[1,], border = NA, col = clr[1] )
    segments( mean( x ), bx[1,1], mean( x ), bx[1,3],
            col = clr[2], lwd = lnSz )
  } else {
    polygon( bx[1,], bx[2,], border = NA, col = clr[1] )
    segments( bx[1,1], mean( x ), bx[1,3], mean( x ),
              col = clr[2], lwd = lnSz )
  }
  
}

### Load in data

if ( Exp == 1 ) {
  dName = 'Gabor_Exp_1_11-17.RData'
} else {
  dName = 'Gabor_Exp_2_11-17.RData'
}

setwd( 'Data' )
load( dName )
rm( dName )
setwd( orig_dir )

# For easy variable access
ad = curDat
colnames( ad ) = simplify_variable_names( curDat )

# Location of where posteriors are saved
post_location = 
  'C:/Users/Kevin/Documents/Posteriors from Stan/Gabor_priming_2017'

###
### Parameter recovery for base SDT model
###
# Lookup - 02

if ( runCode[1] ) {
  
  ### Design ###
  
  # Prime durations
  pd = sort( unique( ad$PD[ ad$BT == 2 ] ) )
  # Prime types
  pt = unique( ad$PTL )
  # Correct orientation (0 = left, 1 = right)
  co = 0:1
  # Combinations of levels
  allCond = expand.grid( co, pd, pt )
  colnames( allCond ) = c( 'Co', 'PD', 'PT' )
  # Number of trials per condition
  nTrials = 15
  
  sim = data.frame(
    Co = rep( NA, nTrials * nrow( allCond ) ),
    PD = rep( NA, nTrials * nrow( allCond ) ),
    PT = rep( NA, nTrials * nrow( allCond ) ) )
  for ( i in 1:3 )
    sim[,i] = as.vector( sapply( allCond[,i], rep, times = nTrials ) )
  
  # Number of parameters
  Nd = length( pd ) * length( pt )
  Nc = 1
  
  # Design matrix for criterion
  X_crt = matrix( 1, nrow( sim ), Nc )
  # Design matrix for d'
  X_d = matrix( 0, nrow( sim ), Nd )
  cnd = aggregate( rep(1,nrow(sim)), list( 
    sim$PD, sim$PT ), sum )
  colnames( cnd ) = c( 'PD', 'PT', 'N' )
  for ( i in 1:nrow( cnd ) ) {
    sel = sim$PD == cnd$PD[i] & 
      sim$PT == cnd$PT[i]
    X_d[sel,i] = 1
  }
  
  # Generating parameters
  beta = c(
    # d' values
    rnorm( Nd, 1.684, .5 ), 
    # criterion values
    rnorm( Nc, 0, .5 )
  )
  
  ### Simulate data ###
  
  L = list(
    No = nrow( sim ),
    Nd = Nd,
    Nc = Nc,
    X_d = matrix( X_d, nrow( sim ), Nd ),
    X_crt = matrix( X_crt, nrow( sim ), Nc ),
    Co = sim$Co, nrow( sim )
  )
  
  sim$Y = as.vector( SDT_rng( L, rbind( beta ) ) )
  sim$Ac = as.numeric( sim$Co == sim$Y )
  
  # Clean up workspace
  rm( L )
  
  ### Parameter recovery ###
  
  # Compile model
  setwd( 'Stan_scripts/SDT' )
  sdtmos = stan_model(stanc_ret = 
                    stanc_builder("SDT_model_OS.stan"))
  setwd( orig_dir )
  
  # Define priors
  Priors = rbind(
    matrix( c( 1.349, 1 ), Nd, 2, byrow = T ),
    matrix( c( 0, 1 ), Nc, 2, byrow = T )
  )
  
  # Define input
  stanDat = list(
    No = nrow( sim ),
    Nd = Nd,
    Nc = Nc,
    X_d = X_d,
    X_crt = X_crt,
    Co = sim$Co,
    Y = sim$Y,
    Priors = Priors
  )
  
  # Attempt to sample from posterior
  fit_stan = sampling( sdtmos, data = stanDat, 
                       warmup = 1000, 
                       iter = 2250, 
                       cores = 8, 
                       chains = 8 )
  
  # Extract posterior samples
  post = extract( fit_stan )
  
  # Convergence diagnostics
  cnv = convergenceExtract( fit_stan, 'beta' )
  # Plot convergence
  plotConvergence( cnv, FALSE, 'beta' )
  
  # Clean workspace
  rm( fit_stan, stanDat, sdtmos )
  
  ### Posterior predictive check ###
  
  # Define input
  L = list(
    No = nrow( sim ),
    Nd = Nd,
    Nc = Nc,
    X_d = matrix( X_d, nrow( sim ), Nd ),
    X_crt = matrix( X_crt, nrow( sim ), Nc ),
    Co = sim$Co
  )
  
  # Simulate data based on posterior samples
  ppc = SDT_rng( L, post$beta )
  cnd = createConditionIndex( c('PD','PT'), sim )
  ppc = apply( ppc, 2, function(x) 
    as.numeric( by( x == sim$Co, list(cnd), mean ) ) )
  
  # New plotting window
  x11( width = 12 );
  layout( cbind( 1, 2 ) )
  
  # Plotting characteristics
  lnSz = 2
  ptSz = 1.3
  axSz = 1.2
  txtSz = 1.2
  
  ### Posterior retrodictive checks ###
  
  ui = list(
    outer = apply( ppc, 1, quantile, prob = quickCoverage(.9) ),
    inner = apply( ppc, 1, quantile, prob = quickCoverage(.5) )
  )
  
  # Set margins
  par( mar = c(3,15,.5,.5) )
  
  # Create a blank plot
  yl = c( .5, max( cnd ) + .5 )
  xl = c( 0, 1 )
  blankPlot( xl, yl )
  
  obs = aggregate( sim$Ac, list( cnd ), mean )
  colnames( obs ) = c( 'cnd', 'Ac' )
  
  # Add axes
  vertLines( c( .25, .5, .75 ), yl, col = 'grey80', lwd = lnSz )
  customAxes( xl, yl, lnSz = lnSz )
  axis( 1, seq(0,1,.25),
        tick = F, line = -1.4, cex.axis = axSz )
  allCond = aggregate( rep(1,nrow(sim)), 
                       list( sim$PT, sim$PD, cnd ), sum )
  colnames( allCond ) = c( 'PT', 'PD', 'Cnd', 'N' )
  lbls = paste( allCond$PT, allCond$PD, 'ms' )
  axis( 2, 1:max(cnd), lbls, 
        tick = F, cex.axis = axSz, las = 1 )
  mtext( 'P(Correct)', side = 1, line = 2, cex = txtSz )
  
  # Credible intervals
  for ( i in 1:max( cnd ) ) {
    quickUI( i, ppc[i,], c(.9,.5), 
             clr = c( 'grey50','white'), lnSz = lnSz )
    points( obs$Ac[i], i, pch = 19, cex = ptSz )
  }
  
  ### Marginal posteriors ###
  
  # Set margins
  par( mar = c(3,5,.5,.5) )
  
  # Create a blank plot
  yl = c( .5, ( Nd + Nc ) +.5 )
  xl = lowerUpper( .5, as.vector( as.numeric( post$beta ) ) )
  blankPlot(xl,yl)
  
  # Add axes
  vertLines( 0, yl, col = 'grey', lwd = lnSz )
  customAxes( xl, yl )
  lbls = c( rep( "d'", Nd ), rep( "crt", Nc ) )
  lbls = paste( lbls, c( 1:Nd, 1:Nc ), sep = '-' )
  axis( 2, 1:(Nd+Nc), lbls, 
        tick = F, cex.axis = axSz, las = 1 )
  axis( 1, round( seq( xl[1], xl[2], length = 4 ), 1 ),
        tick = F, line = -1.4, cex.axis = axSz )
  
  # Credible intervals
  for ( i in 1:( Nd + Nc ) ) {
    quickUI( i, post$beta[,i], c( .9, .5 ), 
             clr = c( 'grey50', 'white' ), lnSz = lnSz )
    points( beta[i], i, pch = 19, cex = ptSz )
  }
  
  mtext( 'Posterior samples', side = 1, line = 2, cex = txtSz )
  
}

###
### Example of fitting multivariate normal hierarchy
###
# Lookup - 03

if ( runCode[2] ) {
  
  # Load in package for simulating multivariate normal
  # install.packages( 'MASS' )
  library( MASS )
  
  # Design
  Ns = 50; # Number of subjects
  K = 3; # Number of IVs
  n = K*20; # Number of trials per subject
  
  # Design matrix
  X = matrix( 0, n, K )
  X[,1] = 1;
  for ( k in 2:K ) {
    sel = 1:(n/K)
    shft = (n/K)*(k-1)
    X[ sel, k ] = -1
    X[ sel + shft, k ] = 1
  }
  
  # Define generating parameters
  
  # Group-level
  Mu = matrix( runif( K, -1, 1 ), K, 1 ); # Mean vector
  
  # Random correlation matrix
  ev = runif(K) # Positive eigenvalues
  U = svd( matrix( runif(K*K), nc = K ) )$u # Orthogonal matrix
  pdM = U %*% diag(ev) %*% t(U) # Positive definite matrix
  Omega = cov2cor( pdM ) # Convert to correlation matrix
  # Clean up workspace
  rm( ev, U, pdM )
  
  # Standard deviations
  Tau = runif( K, .5, 1.5 );
  
  # Create covariance matrix
  Sigma = diag( Tau ) %*% Omega %*% diag( Tau )
  
  # Simulate subject-level parameters
  beta = mvrnorm( Ns, Mu, Sigma )
  
  # Error term
  sigma = runif( 1, .5, 1.5 )
  
  # Simulate data
  tmp = apply( beta, 1, function(x) X %*% matrix(x,K,1) )
  mu = as.vector( tmp ); rm( tmp )
  Y = rnorm( length( mu ), mu, sigma )
  
  # Parameter recovery
  
  # Compile model
  setwd( 'Stan_scripts/MVN' )
  sdtmos = stan_model(stanc_ret = 
                        stanc_builder("MVN_model.stan"))
  setwd( orig_dir )
  
  tmp = X
  X = matrix( NA, length(Y), K )
  for ( ns in 1:Ns ) {
    X[ 1:nrow(tmp) + nrow( tmp )*( ns - 1 ), ] = tmp
  }
  
  # Define input
  stanDat = list(
    No = length( Y ),
    Ns = Ns,
    K = K,
    X = X,
    subj_index = rep( 1:Ns, each = n ),
    Y = Y
  )
  
  # Attempt to sample from posterior
  fit_stan = sampling( sdtmos, data = stanDat, 
                       warmup = 1000, 
                       iter = 2250, 
                       cores = 8, 
                       chains = 8 )
  
  # Convergence diagnostics
  cnv = convergenceExtract( fit_stan, 'sigma' )
  plotConvergence( cnv, FALSE, 'sigma' )
  
  # Extract posterior samples
  post = extract( fit_stan )
  
  quickDesc = function(x) {
    out = numeric( 5 )
    out[1] = mean(x)
    out[2] = median(x)
    out[3] = findMode(x)
    out[4:5] = quantile( x, c(.05,.95) )
    
    return( out )
  }
  
  print( t( apply( post$Mu, 2, 
                   function(x) round( quickDesc(x), 3 ) ) ) )
  print( round( Mu, 3 ) )
  
}

###
### Parameter recovery for hierarchical SDT model
###
# Lookup - 04

if ( runCode[3] ) {
  
  ### Design ###
  
  # Prime durations
  pd = sort( unique( ad$PD[ ad$BT == 2 ] ) )
  # Prime types
  pt = unique( ad$PTL )
  # Correct orientation (0 = left, 1 = right)
  co = 0:1
  # Combinations of levels
  allCond = expand.grid( co, pd, pt )
  colnames( allCond ) = c( 'Co', 'PD', 'PT' )
  # Number of trials per condition
  nTrials = 15
  # Number of subjects
  Ns = 48
  
  sim = data.frame(
    Co = rep( NA, Ns * nTrials * nrow( allCond ) ),
    PD = rep( NA, Ns * nTrials * nrow( allCond ) ),
    PT = rep( NA, Ns * nTrials * nrow( allCond ) ) )
  for ( i in 1:3 )
    sim[,i] = rep( 
      as.vector( sapply( allCond[,i], rep, times = nTrials ) ), Ns )
  sim$S = rep( 1:Ns, each = nTrials * nrow( allCond ) )
  
  # Number of parameters
  Nd = length( pd ) * length( pt )
  Nc = 1
  
  # Design matrix for criterion
  X_crt = matrix( 1, nrow( sim ), Nc )
  # Design matrix for d'
  X_d = matrix( 0, nrow( sim ), Nd )
  cnd = aggregate( rep(1,nrow(sim)), list( 
    sim$PD, sim$PT ), sum )
  colnames( cnd ) = c( 'PD', 'PT', 'N' )
  for ( i in 1:nrow( cnd ) ) {
    sel = sim$PD == cnd$PD[i] & 
      sim$PT == cnd$PT[i]
    X_d[sel,i] = 1
  }
  
  # Define generating parameters
  
  # Group-level
  Mu = c(
    # d' values
    rnorm( Nd, 1.684, .5 ), 
    # criterion values
    rnorm( Nc, 0, .5 )
  )
  
  # Random correlation matrix
  K = Nd + Nc
  ev = runif(K) # Positive eigenvalues
  U = svd( matrix( runif(K*K), nc = K ) )$u # Orthogonal matrix
  pdM = U %*% diag(ev) %*% t(U) # Positive definite matrix
  Omega = cov2cor( pdM ) # Convert to correlation matrix
  # Clean up workspace
  rm( ev, U, pdM )
  
  # Standard deviations
  Tau = runif( K, .5, 1.5 );
  
  # Create covariance matrix
  Sigma = diag( Tau ) %*% Omega %*% diag( Tau )
  
  # Simulate subject-level parameters
  beta = mvrnorm( Ns, Mu, Sigma )
  
  ### Simulate data ###
  
  no = nrow(sim)/Ns
  L = list(
    No = no,
    Nd = Nd,
    Nc = Nc,
    X_d = matrix( X_d[ 1:no, ], no, Nd ),
    X_crt = matrix( X_crt[ 1:no, ], no, Nc ),
    Co = sim$Co[1:no]
  )
  
  sim$Y = as.vector( SDT_rng( L, beta ) )
  sim$Ac = as.numeric( sim$Co == sim$Y )
  
  # Clean up workspace
  rm( L )
  
  ### Parameter recovery ###
  
  if ( modelFit ) {
    
    # Track run time
    time_to_fit_model = Sys.time()
    
    # Compile model
    setwd( 'Stan_scripts/SDT' )
    sm = stan_model(stanc_ret = 
                      stanc_builder("SDT_model_MS.stan"))
    setwd( orig_dir )
    
    # Define priors
    Priors = rbind(
      matrix( c( 1.349, 1 ), Nd, 2, byrow = T ),
      matrix( c( 0, 1 ), Nc, 2, byrow = T )
    )
    Priors = cbind( Priors, 4, 4 )
    Priors = rbind( Priors, c(2,0,0,0) )
    
    # Define input
    stanDat = list(
      Ns = Ns, 
      No = nrow( sim ),
      Nd = Nd,
      Nc = Nc,
      X_d = X_d,
      X_crt = X_crt,
      subj_index = sim$S, 
      Co = sim$Co,
      Y = sim$Y,
      Priors = Priors
    )
    
    # Attempt to sample from posterior
    fit_stan = sampling( sm, data = stanDat, 
                         warmup = 1000, 
                         iter = 2250, 
                         cores = 8, 
                         chains = 8 )
    
    # Extract posterior samples
    post = extract( fit_stan )
    
    # Convergence diagnostics
    cnv = convergenceExtract( fit_stan, 'Tau' )
    
    # Compute run time
    time_to_fit_model = Sys.time() - time_to_fit_model
    
    setwd(post_location)
    save( post, cnv, time_to_fit_model,
          file = 'Test_of_hierarchical_SDT.RData' )
    setwd( orig_dir )
    
    # Clean workspace
    rm( fit_stan, stanDat, sm )
    
  } else {
    
    setwd(post_location)
    load( file = 'Test_of_hierarchical_SDT.RData' )
    setwd( orig_dir )
    
  }
  
  ### Plot results ###
  
  # Plot convergence
  plotConvergence( cnv, FALSE, 'Tau' )
  
  # Posterior predictive checks
  
  # Extract conditions
  cnd = createConditionIndex( c('PD','PT'), sim )
  
  # Simulate from group-level means and 
  # covariance matrix
  f = function(Mu,Omega,Tau) {
    out = 
      mvrnorm( 1, matrix( Mu, length(Mu),1), 
               diag( Tau  ) %*% Omega %*% diag( Tau ) )
    return( out )
  }
  beta = matrix( NA, nrow( post$Mu ), ncol( post$Mu ) )
  for ( i in 1:nrow( post$Mu ) ) {
    beta[i,] = f( post$Mu[i,], post$Omega[i,,], post$Tau[i,] )
  }
  
  # Generate observed accuracies
  no = nrow(sim)/Ns
  L = list(
    No = no,
    Nd = Nd,
    Nc = Nc,
    X_d = matrix( X_d[ 1:no, ], no, Nd ),
    X_crt = matrix( X_crt[ 1:no, ], no, Nc ),
    Co = sim$Co[1:no]
  )
  ppc = SDT_rng( L, beta )
  rm( beta )
  ppc = apply( ppc, 2, function(x) 
    as.numeric( by( x == sim$Co[1:no], list(cnd[1:no]), mean ) ) )
  
  # Observed values
  obs = aggregate( sim$Ac, list( sim$PD, sim$PT, cnd ), mean )
  colnames( obs ) = c( 'PD', 'PT', 'Cnd', 'Ac' )
  
  # New plotting window
  x11( width = 12 )
  layout( cbind( 1,1,2,2) )
  
  # Plotting characteristics
  lnSz = 2
  ptSz = 1.5
  axSz = 1.5
  lbSz = 1.5
  
  # Create blank plot
  pd = sort( unique( obs$PD ) )
  pt = unique( obs$PT )
  xl = c( .5, length( pd ) + .5 )
  yl = c( 0, 1 )
  blankPlot( xl, yl )
  
  # Add axes
  horizLines( c( .25, .5, .75 ), xl, col = 'grey80', lwd = lnSz )
  customAxes( xl, yl, lnSz = lnSz )
  axis( 1, 1:length( pd ), 
        tick = F, line = -1.4, cex.axis = axSz )
  axis( 2, seq(0,1,.25), 
        tick = F, line = -1.4, cex.axis = axSz )
  mtext( 'Prime durations', side = 1, line = 2, cex = lbSz )
  mtext( 'P(Correct)', side = 2, line = 2, cex = lbSz )
  
  xa = rep( 1:length( pd ), length( pt ) )
  pts = rep( c(19,15,17), each = length( pd ) )
  shft = rep( 
    seq( -.25, .25, length = length( pt ) ),
    each = length( pd ) )
  
  for ( i in 1:nrow( obs ) ) {
    
    quickUI( xa[i] + shft[i], ppc[i,], c(.9,.5), flip = F )
    
    points( xa[i] + shft[i], 
            obs$Ac[i], pch = pts[i], cex = ptSz, col = 'blue' )
    
  }
  
  # Add legend
  blankPlot()
  legend( 'left', pt, pch = unique( pts ),
          bty = 'n', cex = lbSz )
  legend( 'topleft', c('Observed','Predicted'),
          fill = c( 'blue', 'grey' ),
          bty = 'n', cex = lbSz )
  
  mtext( 'Posterior retrodictive checks',
         side = 3, line = -2, outer = T, cex = lbSz )
  
  
  ### Plot marginal posteriors ###
  
  # Determine labels for parameters
  assignment = list(
    dp = aggregate( X_d %*% cbind( 1:Nd ),
                    list( sim$PD, sim$PT ), unique ),
    crt = aggregate( X_crt %*% cbind( 1:Nc + Nd ),
                     list( sim$PD, sim$PT ), unique )
  )
  colnames( assignment$dp ) = c( 'PD', 'PT', 'dp' )
  colnames( assignment$crt ) = c( 'PD', 'PT', 'crt' )
  
  # New plotting window
  x11( width = 12 )
  
  # Plotting characteristics
  lnSz = 2
  ptSz = 1.5
  axSz = 1.5
  lbSz = 1.5
  
  # Create a blank plot
  K = Nd + Nc
  xl = c( .5, K + .5 )
  yl = lowerUpper( .5, 
                   as.vector( apply( post$Mu, 2, quantile, 
                          prob = quickCoverage(.9) ) ) )
  blankPlot( xl, yl )
  
  # Add axes
  vertLines( Nd + .5, yl, lty = 2, lwd = lnSz )
  horizLines( seq( yl[1], yl[2], .5 ), xl, col = 'grey80', lwd = lnSz )
  customAxes( xl, yl )
  
  lab = sapply( strsplit( assignment$dp$PT, split = ' ' ),
                function(x) unlist(x)[1] )
  lab = paste( lab, ' (', assignment$dp$PD, ')', sep = '' )
  lab = c( lab, 'Bias left' )
  rotatedAxesLabels( 1:K, lab, tick = F, adj = 0 )
  axis( 2, seq( yl[1], yl[2], .5 ), 
        tick = F, line = -1.4, cex.axis = axSz )
  mtext( 'Posterior samples', side = 2, line = 2, cex = lbSz )
  
  # Draw credible intervals and generating parameters
  for ( k in 1:K ) {
    
    quickUI( k, post$Mu[,k], coverage = c(.9,.5), flip = F )
    points( k, Mu[k], pch = 19, cex = ptSz, col = 'blue' )
  }
  
  # Plot correlations and standard deviations for hierarchical 
  # parameters
  
  # Determine position of variables to extract
  corr_pos = matrix( NA, K*(K+1)/2, 2 )
  
  # Intialize variables
  nc = 0
  inc = 1
  # Loop through rows
  for ( r in 1:K ) {
    nc = nc + 1
    
    # Loop through defined columns
    for ( c in 1:nc ) {
      corr_pos[inc,] = c(r,c);
      inc = inc + 1
    }
    
  }
  
  inc = 1
  inc2 = 1
  nc = 0 # Number of columns
  nWindows = ceiling( K / 4 )
  for ( nw in 1:nWindows ) {
    
    x11( width = 12 )
    layout( cbind( 1:4 ) )
    
    for ( i in 1:4 ) {
      
      # Increment number of columns
      nc = nc + 1
      
      # Create blank plot
      par( mar = c( 3, 4, .5, .5 ) )
      xl = c( .5, K + .5 )
      yl = c(-1, 1 )
      blankPlot( xl, yl )
      
      # Add axes
      horizLines( 0, xl, lwd = lnSz, col = 'grey' )
      customAxes( xl, yl )
      
      # Add credible intervals
      for ( c in 1:nc ) {
      
      index = corr_pos[inc,]
      if ( index[1] == index[2] ) {
        x = post$Tau[,inc2]
        tp = Tau[inc2]
        inc2 = inc2 + 1
      } else {
        x = post$Omega[ , index[1], index[2] ]
        tp = Omega[index[1],index[2]]
      }
      inc = inc+1
      
      quickUI( index[2], x, c(.9,.5), flip = F )
      points( index[2], tp, pch = 19, col = 'blue' )
      }
      
    }
    
  }
   
}

###
### Model fit to experimental data
###
# Lookup - 05

if ( runCode[4] ) {
  
  ### Extract data ###
  
  # Isolate experimental trials
  cd = ad[ ad$BT == 2, ]
  
  # Remove time-out responses
  cd = cd[ cd$Ch != 2, ]
  
  # Number of subjects
  Ns = max( cd$S )
  
  createSDTStanInputs = function( cd, model = 1 ) {
    
    # Prime durations
    pd = sort( unique( cd$PD ) )
    # Prime types
    pt = unique( cd$PTL )
    
    # - Examine interaction of prime duration and 
    #   prime type for drift rate
    # - Intercept only for criterion
    if ( model == 1 ) {
      
      # Number of parameters
      Nd = length( pd ) * length( pt )
      Nc = 1
      
      # Design matrix for criterion
      X_crt = matrix( 1, nrow( cd ), Nc )
      
      # Design matrix for d'
      X_d = matrix( 0, nrow( cd ), Nd )
      cnd = aggregate( rep(1,nrow(cd)), list( 
        cd$PD, cd$PT ), sum )
      colnames( cnd ) = c( 'PD', 'PT', 'N' )
      for ( i in 1:nrow( cnd ) ) {
        sel = cd$PD == cnd$PD[i] & 
          cd$PT == cnd$PT[i]
        X_d[sel,i] = 1
      }
      
    }
    
    # - Examine interaction of prime duration and 
    #   prime type for drift rate
    # - Criterion changes based on previously 
    #   primed choice
    if ( model == 2 ) {
      
      # Number of parameters
      Nd = length( pd ) * length( pt )
      Nc = 2
      
      # Design matrix for criterion
      X_crt = matrix( 1, nrow( cd ), Nc )
      sel = cd$PRL == 'Right'
      X_crt[sel,2] = -1
      sel = cd$PRL == 'Vertical'
      X_crt[sel,2] = 0
      
      # Design matrix for d'
      X_d = matrix( 0, nrow( cd ), Nd )
      cnd = aggregate( rep(1,nrow(cd)), list( 
        cd$PD, cd$PT ), sum )
      colnames( cnd ) = c( 'PD', 'PT', 'N' )
      for ( i in 1:nrow( cnd ) ) {
        sel = cd$PD == cnd$PD[i] & 
          cd$PT == cnd$PT[i]
        X_d[sel,i] = 1
      }
      
    }
    
    # Define priors
    Priors = rbind(
      matrix( c( 1.5, 1 ), Nd, 2, byrow = T ),
      matrix( c( 0, 1 ), Nc, 2, byrow = T )
    )
    Priors = cbind( Priors, 4, 4 )
    Priors = rbind( Priors, c(2,0,0,0) )
    
    # Define input
    stanDat = list(
      Ns = Ns, 
      No = nrow( cd ),
      Nd = Nd,
      Nc = Nc,
      X_d = X_d,
      X_crt = X_crt,
      subj_index = cd$S, 
      Co = cd$Co,
      Y = cd$Ch,
      Priors = Priors
    )
    
    return( stanDat )
  }
  
  ### Fit model to data ###
  
  modelType = 1
  fName = paste( 'Hierarchical_SDT_M', modelType, '.RData', sep = '' )
  
  if ( modelFit ) {
    
    # Track run time
    time_to_fit_model = Sys.time()
    
    # Compile model
    setwd( 'Stan_scripts/SDT' )
    sm = stan_model(stanc_ret = 
                      stanc_builder("SDT_model_MS.stan"))
    setwd( orig_dir )
    
    # Define input
    stanDat = createSDTStanInputs( cd, modelType )
    
    # Attempt to sample from posterior
    fit_stan = sampling( sm, data = stanDat, 
                         warmup = 1000, 
                         iter = 2250, 
                         cores = 8, 
                         chains = 8 )
    
    # Extract posterior samples
    post = extract( fit_stan )
    
    # Convergence diagnostics
    cnv = convergenceExtract( fit_stan, 'Tau' )
    
    # Compute run time
    time_to_fit_model = Sys.time() - time_to_fit_model
    
    setwd(post_location)
    save( post, cnv, time_to_fit_model,
          file = fName )
    setwd( orig_dir )
    
    # Clean workspace
    rm( fit_stan, stanDat, sm )
    
  } else {
    
    setwd(post_location)
    load( file = fName )
    setwd( orig_dir )
    
  }
  
  ### Plot results ###
  
  # Plot convergence
  plotConvergence( cnv, FALSE, 'Tau' )
  
  # Simulate from group-level means and 
  # covariance matrix
  f = function(Mu,Omega,Tau) {
    out = 
      mvrnorm( 1, matrix( Mu, length(Mu),1), 
               diag( Tau  ) %*% Omega %*% diag( Tau ) )
    return( out )
  }
  beta = matrix( NA, nrow( post$Mu ), ncol( post$Mu ) )
  for ( i in 1:nrow( post$Mu ) ) {
    beta[i,] = f( post$Mu[i,], post$Omega[i,,], post$Tau[i,] )
  }
  
  # Generate observed accuracies
  stanDat = createSDTStanInputs( cd, modelType )
  sel = stanDat$subj_index == 1
  no = sum(sel)
  L = list(
    No = no,
    Nd = stanDat$Nd,
    Nc = stanDat$Nc,
    X_d = stanDat$X_d[sel,],
    X_crt = stanDat$X_d[sel,],
    Co = stanDat$Co[sel]
  )
  ppc = SDT_rng( L, beta )
  rm( beta )
  ppc = apply( ppc, 2, function(x) 
    as.numeric( by( x == stanDat$Co[sel], 
                    list(cd$DP[sel]), mean ) ) )
  
  # Observed values
  obs = aggregate( cd$Ac, list( cd$PD, cd$PTL, cd$DP ), mean )
  colnames( obs ) = c( 'PD', 'PT', 'Cnd', 'Ac' )
  
  # New plotting window
  x11( width = 12 )
  layout( cbind( 1,1,2,2) )
  
  # Plotting characteristics
  lnSz = 2
  ptSz = 1.5
  axSz = 1.5
  lbSz = 1.5
  
  # Create blank plot
  pd = sort( unique( obs$PD ) )
  pt = unique( obs$PT )
  xl = c( .5, length( pd ) + .5 )
  yl = c( 0, 1 )
  blankPlot( xl, yl )
  
  # Add axes
  horizLines( c( .25, .5, .75 ), xl, col = 'grey80', lwd = lnSz )
  customAxes( xl, yl, lnSz = lnSz )
  axis( 1, 1:length( pd ), pd, 
        tick = F, line = -1.4, cex.axis = axSz )
  axis( 2, seq(0,1,.25), 
        tick = F, line = -1.4, cex.axis = axSz )
  mtext( 'Prime durations', side = 1, line = 2, cex = lbSz )
  mtext( 'P(Correct)', side = 2, line = 2, cex = lbSz )
  
  xa = rep( 1:length( pd ), length( pt ) )
  pts = rep( c(19,15,17), each = length( pd ) )
  shft = rep( 
    seq( -.25, .25, length = length( pt ) ),
    each = length( pd ) )
  
  for ( i in 1:nrow( obs ) ) {
    
    quickUI( xa[i] + shft[i], ppc[i,], c(.9,.5), flip = F )
    
    points( xa[i] + shft[i], 
            obs$Ac[i], pch = pts[i], cex = ptSz, col = 'blue' )
    
  }
  
  # Add legend
  blankPlot()
  legend( 'left', pt, pch = unique( pts ),
          bty = 'n', cex = lbSz )
  legend( 'topleft', c('Observed','Predicted'),
          fill = c( 'blue', 'grey' ),
          bty = 'n', cex = lbSz )
  
  mtext( 'Posterior retrodictive checks',
         side = 3, line = -2, outer = T, cex = lbSz )
  
  ### Plot marginal posteriors ###
  
  Nd = stanDat$Nd
  Nc = stanDat$Nc
  
  # Determine labels for parameters
  assignment = list(
    dp = aggregate( stanDat$X_d %*% cbind( 1:Nd ),
                    list( cd$PD, cd$PTL ), unique ),
    crt = aggregate( stanDat$X_crt %*% 
                       cbind( 1:Nc + Nd),
                     list( cd$PD, cd$PTL ), unique )
  )
  colnames( assignment$dp ) = c( 'PD', 'PT', 'dp' )
  colnames( assignment$crt ) = c( 'PD', 'PT', 'crt' )
  
  # New plotting window
  x11( width = 12 )
  
  # Plotting characteristics
  lnSz = 2
  ptSz = 1.5
  axSz = 1.5
  lbSz = 1.5
  
  # Create a blank plot
  K = Nd + Nc
  xl = c( .5, K + .5 )
  yl = lowerUpper( .5, 
                   as.vector( apply( post$Mu, 2, quantile, 
                                     prob = quickCoverage(.9) ) ) )
  blankPlot( xl, yl )
  
  # Add axes
  vertLines( Nd + .5, yl, lty = 2, lwd = lnSz )
  horizLines( seq( yl[1], yl[2], .5 ), xl, col = 'grey80', lwd = lnSz )
  customAxes( xl, yl )
  
  lab = sapply( strsplit( assignment$dp$PT, split = ' ' ),
                function(x) unlist(x)[1] )
  lab = paste( lab, ' (', assignment$dp$PD, ')', sep = '' )
  lab = c( lab, 'Bias left' )
  rotatedAxesLabels( 1:K, lab, tick = F, adj = 0 )
  axis( 2, seq( yl[1], yl[2], .5 ), 
        tick = F, line = -1.4, cex.axis = axSz )
  mtext( 'Posterior samples', side = 2, line = 2, cex = lbSz )
  
  # Draw credible intervals and generating parameters
  for ( k in 1:K ) {
    
    quickUI( k, post$Mu[,k], coverage = c(.9,.5), flip = F,
             clr = c('grey40','white') )
  }
  
  ### Difference between target and foil conditions ###
  
  x11( width = 12 )
  
  # Plotting characteristics
  lnSz = 2
  ptSz = 1.5
  axSz = 1.5
  lbSz = 1.5
  
  # Create a blank plot
  xl = c( .5, length( pd ) + .5 )
  yl = c( -3, 3 )
  blankPlot( xl, yl )
  
  # Add axes
  vertLines( Nd + .5, yl, lty = 2, lwd = lnSz )
  horizLines( seq( yl[1], yl[2], 1 ), xl, col = 'grey80', lwd = lnSz )
  customAxes( xl, yl )
  axis( 1, 1:length(pd), pd,
        tick = F, line = -1.4, cex.axis = axSz )
  axis( 2, round( seq( yl[1], yl[2], length = 5 ), 2 ),
        tick = F, line = -1.4, cex.axis = axSz )
  mtext( 'Prime duration (ms)', side = 1, line = 2, cex = lbSz )
  mtext( 'Target - Foil', side = 2, line = 2, cex = lbSz )
  
  for ( i in 1:length( pd ) ) {
    
    sel1 = assignment$dp$PT == 'Target primed' & 
      assignment$dp$PD == pd[i]
    sel2 = assignment$dp$PT == 'Foil primed' & 
      assignment$dp$PD == pd[i]
    
    x = post$Mu[,sel1] - post$Mu[,sel2]
    quickUI( i, x, coverage = c(.9,.5), flip = F,
             clr = c('grey40','white') )
    rm( x )
    
  }
  
}