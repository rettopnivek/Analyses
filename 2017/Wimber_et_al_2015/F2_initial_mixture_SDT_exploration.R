#------------------------------------------#
# Initial exploration of mixture SDT model #
# Kevin Potter                             #
# Updated 01/14/2017                       #
#------------------------------------------#

# Clear workspace
rm( list = ls() )

# Save current directory
orig_dir = getwd()

# Define code segments to run
runCode = c( F, F, F, T, T, T )

# Indicate whether to save a pdf file
savePlot = T
if (savePlot) {
  setwd( 'Plots' )
  fName = 'mSDT_model_results.pdf'
  pdf( fName, width = 12, height = 6 )
  rm( fName )
  setwd( orig_dir )
}

# Index
# Lookup - 01:  Load in useful packages and data
# Lookup - 02:  Define useful functions
# Lookup - 03:  Parameter recovery (MLE)
# Lookup - 04:  Parameter recovery via Bayesian estimation
# Lookup - 05:  Simulation with hierarchical model 
#               (w/ item/subject effects)
# Lookup - 06:  Fit to data (Hierarchical, no item effects)
# Lookup - 07:  Plots of model fit results
# Lookup - 08:  Posterior simulations for false-positives

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
# install.packages( 'rstan' )
library(rstan)
# For parallel processing
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# Load in package for truncated normal
# install.packages( 'truncnorm' )
library( truncnorm )

# Load in package for linear modeling
# install.packages( lme4 )
library( lme4 )

# Load in data
setwd( 'Data' )
load( 'Original_all_data.RData' )
setwd( orig_dir )

# Extract observations for final recogniton memory 
# test and pre-process data

d = OriginalAllData[ OriginalAllData$Cond == 6, ]
# For easy manipulation
colnames( d ) = c('S','Tr','C','IN','Co','Ch','RT',
                  'Ac','IT','B','Bl','Cat','CR','fN')

# Use dummy coding for responses
d$Ch = d$Ch - 1; d$Co = d$Co - 1

# Define a variable denoting conditions

# Targets that underwent selective retrieval (T-SR)
d$Cnd = 1;
# Targets in the baseline condition (T-B)
d$Cnd[ d$IT == 1 & d$B == 1 ] = 2
# Competitors that underwent selective retrieval (C-SR)
d$Cnd[ d$IT == 2 & d$B == 0 ] = 3
# Competitors in the baseline condition (C-B)
d$Cnd[ d$IT == 2 & d$B == 1 ] = 4

# Set missing data to be incorrect
d$Ac[ is.na(d$RT) ] = 0
d$Ch[ is.na(d$RT) ] = 1 - d$Co[ is.na(d$RT) ]

###
### Define useful functions
###
# Lookup - 02

mSDT_prob = function( dp, crt, lmb, Co ) {
  # Purpose:
  # A function to calculate the probability of picking 
  # an alternative on the right in a two-alternative 
  # forced-choice task based on an mixture SDT model.
  # Arguments:
  # dp  - The d' value(s), the separation between means for the 
  #       strength distributions for the left and right 
  #       choices.
  # crt - The criterion value(s), denoting the bias towards 
  #       left (positive values) or right (negative values)
  # lmb - The mixture probability (represents properly attending 
  #       to the stimulus)
  # Co  - The position(s) of the correct alternative, where 
  #       0 = left, and 1 = right
  # Returns: 
  # The probability (or if vectors are inputs, a vector of 
  # probabilities) for picking the alternative on the right.
  
  # Calculate the probability of picking right given the 
  # correct answer is on the right.
  P_RR = lmb*( 1.0 - pnorm( crt, dp/2.0, 1.0 ) ) + 
    ( 1.0 - lmb )*( 1.0 - pnorm( crt, 0.0, 1.0 ) )
  
  # Calculate the probability of picking right given the 
  # correct answer is on the left.
  P_RL = lmb*( 1.0 - pnorm( crt, -dp/2.0, 1.0 ) ) + 
    ( 1.0 - lmb )*( 1.0 - pnorm( crt, 0.0, 1.0 ) )
  
  # Select the appropriate probability for the given trial
  theta = Co*P_RR + (1-Co)*P_RL
  
  return( theta )
}

draw_figure = function( prp, lnSz = 1, lnTp = 1, ptSz = 1, 
                        pts = c(19,22,19,22),
                        clr = c('black','white'),
                        xa = NULL ) {
  # Purpose:
  # Adds a set of lines giving the proportion of times the 
  # response on the right was picked.
  # Arguments:
  # prp  - A set of 8 proportions, sorted by the correct position 
  #        of the target, the type of image (target versus competitor),
  #        and the condition (selective retrieval versus baseline)
  # lnSz - The width of the line
  # lnTp - The type of line
  # ptSz - The size of the point
  # pts  - A vector of 4 point types
  # clr  - The foreground and background color for the points
  # xa   - An optional list giving the x-axis positions for each 
  #        condition
  
  if ( length( xa ) == 0 ) 
    xa = list(
      c(.9,2.1),
      c(1.1,1.9),
      c(.9,2.1),
      c(1.1,1.9) )
  
  for ( i in 1:4 ) {
    
    sel = 1:2 + 2*(i-1)
    lines( xa[[i]], prp[sel], lwd = lnSz, lty = lnTp, col = clr[1] )
    points( xa[[i]], prp[sel], pch = pts[i], col = clr[1],
            bg = clr[2], cex = ptSz )
    
  }
  
}

convergence_extract = function( fit, par_name = NULL ) {
  # Purpose:
  # Extract convergence diagnostics from a Stan fit object.
  # Arguments:
  # fit      - A Stan fit object
  # par_name - An optional string giving the final parameter label 
  #            of the subset of the output to include
  # Notes:
  # Extracting the convergence statistics can be slow, especially 
  # when a large number of parameters were stored.
  # Returns:
  # A list with the Gelman-Rubin convergence statistic, the 
  # effective number of samples, and the total number of samples.
  
  Rhat = summary(fit)$summary[,"Rhat"]
  n_eff = summary(fit)$summary[,"n_eff"]
  totSampSize = length(extract(fit, pars = "lp__")[[1]])
  # We're only interested in a subset of parameters
  if ( length( par_name ) == 0 ) 
    par_name = names( Rhat )[ 
      which( names( Rhat ) == "logLik[1]" ) - 1 ]
  sel = 1:which( names(Rhat) == par_name )
  Rhat = Rhat[sel]; n_eff = n_eff[sel];
  
  return( list( Rhat = Rhat, n_eff = n_eff, totSampSize = totSampSize ) )
}

find_dec = function( x, spacing = 10 ) {
  # Purpose:
  # Determines the rounded leading digit and the 
  # number of trailing zeros for a number or the 
  # number of decimal places.
  # Arguments:
  # x       - A vector of values
  # spacing - The value whose exponent should be increased
  # Returns:
  # A vector giving the leading digit, the number of 
  # trailing/leading zeros, the same but in scientific 
  # notation, and 1 if it's trailing zeros, -1 if it's 
  # decimal places.
  
  mx = max( x )
  rnd = mx
  
  if ( round( mx ) > 1 ) {
    
    inc = 0;
    
    while ( rnd > 1 ) {
      inc = inc + 1;
      rnd = round( mx/( spacing^inc ) )
    }
    
    v = round( mx/spacing^(inc-1) )
    f = spacing^(inc-1)
    out = c( v,f,inc-1,1)
  }
  
  if ( round( mx ) == 1 ) {
    
    out = c( 1, 1, 1, 0 )
    
  }
  
  if ( round( mx ) == 0 ) {
    
    inc = 0;
    
    while ( rnd < 1 ) {
      inc = inc + 1;
      rnd = round( mx*( spacing^inc ) )
    }
    
    v = round( mx*spacing^(inc) )
    f = spacing^(inc)
    out = c( v,f,inc,-1)
    
  }
  
  return( out )
}

###
### Parameter recovery (MLE)
###
# Lookup - 03

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

###
### Fit to data (Hierarchical, no item effects)
###
# Lookup - 06

if ( runCode[4] ) {
  
  setwd( 'Stan_scripts' )
  
  # Define input for Stan
  stanDat = list(
    Ns = length( sort( unique( d$S ) ) ), 
    Ni = length( sort( unique( d$IN ) ) ), 
    No = nrow( d ), 
    Y = d$Ch, 
    Co = d$Co, 
    subjIndex = d$S, 
    Cond = cbind( d$B+1, d$IT )
  )
  
  warm = 750 # Warm-up
  niter = 1250 # Number of samples to approximate posterior
  chains = 8 # Number of chains to run
  
  startTime = Sys.time() # To assess run-time
  
  # Compile model
  sm = stan_model(stanc_ret = stanc_builder("mSDT_model_no_item.stan"))
  
  # Draw samples
  fit = sampling( sm, data = stanDat, 
                  warmup = warm, iter = warm+niter, 
                  chains = chains,
                  control = list( adapt_delta = .92 ) )
  
  post = extract(fit)
  
  setwd( orig_dir )
  
  # Save results to a .RData file
  fName = "C:/Users/Kevin/Documents/Posteriors from Stan"
  setwd(fName)
  setwd( "Wimber_et_al" )
  save( stanDat, fit, post, file = 'mSDT_model_no_item_post.RData' )
  setwd( orig_dir )
  
}

###
### Plots of model fit results
###
# Lookup - 07

if ( runCode[5] ) {
  
  # Load in posterior
  fName = "C:/Users/Kevin/Documents/Posteriors from Stan"
  setwd(fName)
  setwd( "Wimber_et_al" )
  load( 'mSDT_model_no_item_post.RData' )
  setwd( orig_dir )
  
  ### Check convergence ###
  conv = convergence_extract(fit,'nu[2]')
  
  if (!savePlot) x11(width=12)
  layout( cbind(1,2) )
  
  # Plot a histogram of the Gelman-Rubin statistics for the 
  # marginal posterior samples of the parameters
  
  tmp = hist( conv$Rhat, plot = F )
  scl = find_dec( tmp$density )
  if ( scl[4] == 1 ) scl = scl[1]*(scl[2]/10) else 
    scl = scl[1]/(scl[2])
  
  yl = lowerUpper( scl, tmp$density )
  yl[1] = 0
  
  xl = lowerUpper( .1, conv$Rhat )
  xl[2] = max( xl[2], 1.12 )
  xl[1] = min( xl[1], .98 )
  
  plot( xl, yl, type = 'n', cex.axis = 1.5, cex.lab = 1.5,
        xlab = expression( hat(R) ), ylab = 'Density',
        bty = 'l', main = 'Gelman-Rubin Statistic' )
  
  segments( tmp$mids, rep(0,length(tmp$mids)),
            tmp$mids, tmp$density, lwd = 3,
            col = 'grey' )
  abline( v = 1.1, lty = 2, lwd = 2 )
  
  # Plot a histogram of the effective sample size for the 
  # set of parameters
  
  tmp = hist( conv$n_eff, plot = F )
  scl = find_dec( tmp$density )
  if ( scl[4] == 1 ) scl = scl[1]*(scl[2]/10) else 
    scl = scl[1]/(scl[2])
  
  yl = lowerUpper( scl, tmp$density )
  yl[1] = 0
  
  xl=c(0,conv$totSampSize)
  plot( xl, yl, type = 'n', cex.axis = 1.5, cex.lab = 1.5,
        xlab = expression(N[sample]), ylab = 'Density',
        bty = 'l', main = 'Effective sample size' )
  
  segments( tmp$mids, rep(0,length(tmp$mids)),
            tmp$mids, tmp$density, lwd = 3,
            col = 'grey' )
  
  ### Posterior retrodictive checks ###
  
  # Extract data
  x = aggregate( d$Ch, list( d$B, d$IT, d$Co, d$S ), mean )
  colnames( x ) = c('B','IT','Co','S','P')
  
  # Group-level performance
  
  # Generate retrodictive checks for group-level performance
  f = function(s) {
    yhat = rbinom( stanDat$No, 1, post$theta[s,] )
    out = aggregate( yhat, list( stanDat$Cond[,1], 
                                 stanDat$Cond[,2], 
                                 stanDat$Co ), mean )$x
    
    return( out )
  }
  rc = sapply( 1:nrow(post$beta), f )
  
  # Create a plot of the proportions
  if (!savePlot) x11(width=12)
  layout( cbind(1,2) )
  plot( c(.8,2.2), c(0,1), type = 'n', xlab = 'Condition',
        ylab = 'P(Right)', xaxt = 'n', bty = 'n' )
  axis( 1, 1:2, c('Selective retrieval','Baseline'), tick = F )
  abline(h=.5)
  legend('topright','Right correct',bty='n')
  legend('bottomright','Left correct',bty='n')
  legend( 'topleft', c('Target','Competitor'),
          pch = c(19,22), pt.bg = 'white', bty = 'n' )
  
  ya = seq( .1, .9, .01 )
  points( rep(1,length(ya)), ya, pch = '.', cex = .5 )
  points( rep(2,length(ya)), ya, pch = '.', cex = .5 )
  
  pred = apply( rc, 1, findMode )
  
  ui = apply( rc, 1, quantile, prob = c( .025,.16,.25,.75,.84,.975 ) )
  ui = ui[ c(1:3,6:4), ]
  
  xa = c(
    c(.9,2.1),
    c(1.1,1.9),
    c(.9,2.1),
    c(1.1,1.9) )
  clr = c( rgb( .5, .5, .5, .4 ),
           rgb( .5, .5, .5, .4 ),
           rgb( .5, .5, .5, .4 ) )
  sz = c( .05, .04, .03 )
  
  for (i in 1:8) {
    
    for (j in 1:3) {
      
      drawEllipse( sz[j], diff( ui[c(j,j+3),i] ),
                   Xc = xa[i], Yc = ui[j,i] + diff( ui[c(j,j+3),i] )/2,
                   border = NA, col = clr[j] )
      
    }
    
  }
  points( xa, pred, pch = 19, col = 'grey40' )
  
  obs = aggregate( x$P, list( x$B, x$IT, x$Co ), 
                   function(x) c( mean(x), sd(x)/sqrt(length(x)) ) )
  colnames( obs ) = c( 'B', 'IT', 'Co', 'DS' )
  draw_figure( obs$DS[,1] )
  
  mtext( 'Group-level performance', side = 3, outer = T, line = -2 )
  
  # Subject-level performance
  
  # Generate retrodictive checks for subject-level performance
  indexS = NULL
  f = function(s) {
    yhat = rbinom( stanDat$No, 1, post$theta[s,] )
    out = aggregate( yhat, list( stanDat$Cond[,1], 
                                 stanDat$Cond[,2], stanDat$Co,
                                 stanDat$subjIndex ), mean )
    if ( length( indexS ) == 0 ) indexS <<- out[,4]
    
    return( out$x )
  }
  rc = sapply( 1:nrow(post$beta), f )
  
  if (!savePlot) x11(width=12)
  layout( matrix( 1:12, 3, 4, byrow = T ) )
  
  for ( s in 1:12 ) {
    blankPlot( c(.8,2.2), c(0,1) )
    title( paste('Subject', s ) )
    if ( s == 1 | s == 5 | s == 9 ) 
      axis( 2, seq(0,1,.25) )
    if ( s == 9 | s == 10 | s == 11 | s == 12 ) 
      axis( 1, 1:2, c('SR','B' ), tick = F )
    
    draw_figure( x$P[ x$S == s ] )
    draw_figure( apply( rc[ indexS == s, ], 1, findMode ),
                 clr = c('blue','white'), lnTp = 2 )
  }
  
  if (!savePlot) x11(width=12)
  layout( matrix( 1:12, 3, 4, byrow = T ) )
  
  for ( s in 13:24 ) {
    blankPlot( c(.8,2.2), c(0,1) )
    title( paste('Subject', s ) )
    if ( s == 1 | s == 5 | s == 9 ) 
      axis( 2, seq(0,1,.25) )
    if ( s == 9 | s == 10 | s == 11 | s == 12 ) 
      axis( 1, 1:2, c('SR','B' ), tick = F )
    
    
    draw_figure( x$P[ x$S == s ] )
    draw_figure( apply( rc[ indexS == s, ], 1, findMode ),
                 clr = c('blue','white'), lnTp = 2 )
  }
  
  # Posterior distributions
  
  if (!savePlot) x11(width=12)
  layout( matrix( 1:6, 2, 3, byrow = T ) )
  
  ttl = c( "d' for SR", "d' for B" )
  for ( i in 1:2 ) {
    tmp = hist( post$mu_dp[,i], col = 'grey', border = 'white', 
                breaks = 40, freq = F, main = ttl[i], xlab = ' ' )
    ui = quantile( post$mu_dp[,i], c(.025,.975) )
    dn = c( max( which( tmp$mids < ui[1] ) ),
            min( which( tmp$mids > ui[2] ) ) )
    dn = tmp$density[ dn ]
    segments( ui, c(0,0), ui, dn, col = 'blue' )
  }
  
  ttl = c( "Mixture proportion for targets", 
           "Mixture proportion for competitors" )
  for ( i in 1:2 ) {
    tmp = hist( post$phi[,i], col = 'grey', border = 'white', 
                breaks = 40, freq = F, main = ttl[i], xlab = ' ' )
    ui = quantile( post$phi[,i], c(.025,.975) )
    dn = c( max( which( tmp$mids < ui[1] ) ),
            min( which( tmp$mids > ui[2] ) ) )
    dn = tmp$density[ dn ]
    segments( ui, c(0,0), ui, dn, col = 'blue' )
  }
  
  # Criterion
  crt = post$mu_crt
  tmp = hist( crt, col = 'grey', 
              border = 'white', breaks = 40, freq = F, 
              main = 'Mean criterion', xlab = ' ' )
  ui = quantile( crt, c(.025,.975) )
  dn = c( max( which( tmp$mids < ui[1] ) ),
          min( which( tmp$mids > ui[2] ) ) )
  dn = tmp$density[ dn ]
  segments( ui, c(0,0), ui, dn, col = 'blue' )
  
  blankPlot()
  
}

###
### Posterior simulations for false-positives
###
# Lookup - 08

if ( runCode[6] ) {
  
  # Load in posterior
  fName = "C:/Users/Kevin/Documents/Posteriors from Stan"
  setwd(fName)
  setwd( "Wimber_et_al" )
  load( 'mSDT_model_no_item_post.RData' )
  setwd( orig_dir )
  
  # Define design
  Ni = 144 # Number of items
  Ns = 24 # Number of subjects
  No = Ni*Ns # Total number of observations
  Co = numeric(No) # Position of correct choice
  for (ns in 1:Ns) Co[ 1:Ni + Ni*(ns-1) ] = sample( rep(0:1,each=72) )
  # Condition indices
  IT = rep( 2:1, each = 72 ); IT = rep(IT,Ns)
  B = rep( c( rep(1,54),rep(2,18) ), 2 ); B = rep(B,Ns)
  Cnd = cbind( B, IT )
  indexS = rep( 1:Ns, each = Ni )
  
  # Vectors for parameter estimates
  dp_all = matrix( NA, Ns, 2 )
  l_all = matrix( NA, Ns, 2 )
  dp = rep( NA, No )
  crt = rep( NA, No )
  lmb = rep( NA, No )
  
  # Effects coded variables
  SR_ef = rep(-1,No); SR_ef[ B == 1 ] = 1
  A_ef = rep(-1,No); A_ef[ IT == 1 ] = 1
  
  # Matrix to store p-values
  res = matrix( NA, 1000, 3 )
  
  # Create a progress bar using a base R function
  pb = txtProgressBar( min = 0, max = 1000, style = 3 )
  
  # Loop over a subset of the samples from the hyper-parameters
  # generating fake subjects and calculates the pair of t-tests
  # the interaction
  inc = 1
  for ( s in sample( 1:nrow(post$beta), 1000 ) ) {
    
    # Update the progress bar
    setTxtProgressBar(pb,inc)
    
    # Simulate data
    dp_all[,1] = rtruncnorm( Ns, 
                             post$mu_dp[s,1], 
                             post$sigma_dp[s,1], a = 0 )
    dp_all[,2] = rtruncnorm( Ns, 
                             post$mu_dp[s,2], 
                             post$sigma_dp[s,2], a = 0 )
    crt_all = rnorm( Ns, post$mu_crt[s], post$sigma_crt[s] )
    l_all[,1] = rbeta( Ns, post$phi[s,1]*post$nu[s,1],
                 (1 - post$phi[s,1])*post$nu[s,1] )
    l_all[,2] = rbeta( Ns, post$phi[s,2]*post$nu[s,2],
                 (1 - post$phi[s,2])*post$nu[s,2] )
    
    dp[ B == 1 ] = dp_all[ indexS[ B == 1 ], 1 ]
    dp[ B == 2 ] = dp_all[ indexS[ B == 2 ], 2 ]
    lmb[ IT == 1 ] = l_all[ indexS[ IT == 1 ], 1 ]
    lmb[ IT == 2 ] = l_all[ indexS[ IT == 2 ], 2 ]
    crt = crt_all[ indexS ]
    
    theta = mSDT_prob( dp, crt, lmb, Co )
    
    Y = rbinom( No, 1, theta )
    Ac = as.numeric( Y == Co )
    
    # Calculate p-values
    obs = aggregate( Ac, list( B, IT, indexS ), mean )
    colnames( obs ) = c('B','IT','S','P')
    
    # T-tests
    res1 = t.test( obs$P[ obs$B == 1 & obs$IT == 1 ] - 
              obs$P[ obs$B == 2 & obs$IT == 1 ] )
    res2 = t.test( obs$P[ obs$B == 1 & obs$IT == 2 ] - 
                     obs$P[ obs$B == 2 & obs$IT == 2 ] )
    # Interaction
    res3 = lmer( Ac ~ SR_ef + A_ef + SR_ef*A_ef + 
                   (1|indexS) )
    res3 = pt(summary(res3)$coefficients[4,3],Ns-1)
    res[inc,] = c(res1$p.value,res2$p.value,res3)
    inc = inc + 1
    
  }
  close(pb)
  
  # Calculate how many times desired pattern occurred
  chk = apply( res, 1, function(x) 
    c( x[1] > .05 & x[2] < .05,
       x[3] < .05 ) )
  
  # Select p-values that met necessary pattern
  sel = chk[1,] == T
  
  # Plot of t-values
  if (!savePlot) x11(width=12)
  layout( cbind(1) )
  blankPlot( c(1,1000), c(0,1) )
  axis( 2, seq(0,1,.25) )
  abline( h = .05, lty = 2 )
  mtext( 'p-value', side = 2, line = 2 )
  axis( 3, sum( sel ), 
        paste( 100*(sum(sel)/1000), '%', sep = '' ),
        tick = F )
  title( 'T-test p-values' )
  
  polygon( c(.5,sum(sel)+.5,sum(sel)+.5,.5),
           c(0,0,1,1), col = 'grey', border = NA )
  
  ord = 1:sum( chk[1,] )
  points( ord, res[sel,1], pch = 19 )
  points( ord, res[sel,2], pch = 19, col = 'red' )
  
  ord = ( sum( chk[1,] ) + 1 ):1000
  points( ord, res[!sel,1], pch = 19 )
  points( ord, res[!sel,2], pch = 19, col = 'red' )
  
  # so turn off clipping:
  par(xpd=TRUE)
  legend(350,-.08,c("Test for targets", 
                    "Test for competitors"), 
         pch = 19, col = c('black','red'), horiz = T, bty = 'n' )
  par(xpd=FALSE)
}

if (savePlot) dev.off()