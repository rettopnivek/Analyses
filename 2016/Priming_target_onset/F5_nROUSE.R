#-----------------------------#
# nROUSE parameter estimation #
# Kevin Potter                #
# Updated 01/03/2017          #
#-----------------------------#

# Initialize script
source('F4_starting_script.R')

# Load in package for simulating nROUSE model
# install_github("rettopnivek/nROUSE")
library(nROUSE)

# Load in package for adaptive Bayesian estimation
# install.packages('MHadaptive')
library( MHadaptive )

# Load in data
setwd( 'Data' )
load( 'Priming_offset.RData' )
setwd( orig_dir )

# Indicate which code segments to run
runCode = c( F, T, T, T )

# Indicate whether to save a pdf file
savePlot = T
if (savePlot) {
  
  setwd( 'Plots' )
  
  if ( runCode[1] ) pdf( file = 'nROUSE_model_sim.pdf', 
                         width = 12, height = 6 )
  
  if ( !runCode[1] ) pdf( file = 'nROUSE_model.pdf', 
                         width = 12, height = 6 )
  setwd( orig_dir )
}

# For easy manipulation
d = allDat
colnames( d ) = c('S','PD','PT','O','Co','RT',
                  'Ch','Ac','OT','DT','CDT','Cnd','TD')

# Index
# Lookup - 01:  Define useful functions
# Lookup - 02:  Simulation example for one subject
# Lookup - 03:  Model estimation at the group level
# Lookup - 04:  Model estimation at the subject level
# Lookup - 05:  Save covariates from nROUSE estimation

###
### Define additional functions
###
# Lookup - 01

pH = function( pst, ttl, alpha = .68, dec = 2 ) {
  # Purpose:
  # Creates a histogram of a vector of posterior samples 
  # and draws line segments denoting a desired credible interval.
  # Arguments:
  # pst   - A vector with the posterior samples
  # ttl   - The title for the histogram
  # alpha - The width of the credible interval
  # dec   - The number of decimal places to report
  
  tmp = hist( pst, breaks = 40, freq = F, col = 'grey',
              border = 'white', xlab = ' ',
              main = ttl )
  bnd = (1-alpha)/2
  bnd = c( bnd, 1-bnd)
  ui = quantile( pst, bnd )
  h = c( max( which( tmp$mids <= ui[1] ) ),
         min( which( tmp$mids >= ui[2] ) ) )
  segments( ui, c(0,0), ui, tmp$density[h], lwd = 2 )
  text( ui, max(tmp$density[h]) + max(tmp$density[h])*.1, 
        as.character( round(ui,dec) ), pos=c(2,4) )
}

pP = function( pst, lbls ) {
  # Purpose:
  # Creates a plot of the bivariate density for a pair of 
  # vectors of posterior samples.
  # Arguments:
  # pst  - A matrix with two columns, one for each vector of 
  #        posterior samples
  # lbls - A vector giving the x-axis label, the y-axis label,
  #        and the title
  
  # Extract posterior estimates
  x = pst[,1]
  y = pst[,2]
  
  # Determine limits for plotting window
  zx = scale( pst[,1] )
  zy = scale( pst[,2] )
  
  rn = c( min( min(zx), min(zy) ),
          max( max(zx), max(zy) ) )
  
  xl = rn*sd(x) + mean(x)
  yl = rn*sd(y) + mean(y)
  
  plot( xl, yl, type = 'n', bty = 'l', xlab = lbls[1],
        ylab = lbls[2], main = lbls[3] )
  points( x, y, pch = 19, col = rgb( 0, 0, 1, .2 ) )
  
}

# Define a function to generate starting values
st_f = function() {
  sv = c( N = .0302, I = .9844, Ta = 1 )
  sv = sv + runif( 3, c(-.015,-.5,-.5),
                   c(.01,.75,.75) )
  sv = log( sv )
  return( sv )
}

# Define a function that calculates the summed log-likelihoods for MLE
mle_f = function( par, dat, priors = NULL ) {
  
  out = nROUSE_logLik( par, dat )
  
  return( out )
}

###
### Simulation example for one subject
###
# Lookup - 02

if ( runCode[1] ) {
  
  # Generating parameters
  mu = c( N = .0283, I = 1.2866, Ta = 0.8931, TD = 65 )
  sig = c( N = .002, I = .35, Ta = .3, TD = 15 )
  # Jitter by some noise
  gp = rnorm( 4, mu, sig ); gp[4] = round( gp[4] )
  
  # Create data set to obtain nROUSE model estimates
  nDat = cbind( TarDur = gp[4], MaskDur = 500 - gp[4],
                PrimeDur = c(50,2000,50,2000),
                Type = c(2,2,-2,-2), Y = rep(0,4),
                N = rep(0,4) )
  nDat = as.data.frame( nDat )
  theta = nROUSE_logLik( log(gp[1:3]), nDat, estimate = F )
  
  # Simulate frequency correct
  nDat$N = rep( 80, 4 )
  nDat$Y = rbinom( 4, nDat$N, theta )
  
  # Define a set of priors
  priors = cbind( c(.0283, 1.2866, .8931 ),
                  c(.01, .5, .4 ) )
  
  # Define a log-likelihood function
  ll_f = function( prm, nDat, priors ) {
    
    ll = nROUSE_logLik( log(prm), nDat )
    ll = ll + sum( dnorm( prm, priors[,1], priors[,2], log = T ) )
    if (is.na(ll)) ll = -Inf
    
    return( ll )
  }
  
  # Generate random starting values
  st_val = runif( 3, priors[,1] - priors[,2],
                  priors[,1] + priors[,2] )
  
  # Carry out Bayesian estimation
  res = Metro_Hastings( ll_f, st_val, nDat = nDat, priors = priors, 
                        #prop_sigma = ps, 
                        par_names = c('N','I','A'),
                        burn_in = 3000, 
                        iterations = 10000 )
  
  # Check convergence
  if ( !savePlot ) x11(width=12)
  plot.ts( res$trace, main = 'Convergence' )
  
  # Plot the posterior distributions
  if ( !savePlot ) x11(width=12)
  layout( cbind( c(1,7,8), c(4,2,9), c(5,6,3) ) )
  
  # Marginal posterior distributions
  pH( res$trace[,1], 'Noise multiplier' )
  pH( res$trace[,2], 'Inibition' )
  pH( res$trace[,3], 'Temporal attention' )
  
  # Bivariate distributions
  pP( res$trace[,2:1], c('I','N',' ') )
  pP( res$trace[,c(3,1)], c('T','N',' ') )
  pP( res$trace[,c(3,2)], c('T','I',' ') )
  
  # Title
  blankPlot()
  blankPlot()
  legend('left','Marginal and bivariate \nposterior distributions',
         bty = 'n', cex = 1.5 )
  blankPlot()
  
}

###
### Model estimation at the group level
###
# Lookup - 03

if ( runCode[2] ) {
  
  # Determine aggregate performance
  nDat = aggregate( d$Ac, list( d$PD, d$PT ), mean )
  
  # Incorporate prime duration and type
  colnames( nDat ) = c( 'PD','PT','P' )
  nDat$PrimeDur = 50; nDat$PrimeDur[ nDat$PD == 1] = 2000
  nDat$Type = -2; nDat$Type[ nDat$PT == 0 ] = 2;
  
  nDat$TarDur = round( mean( d$TD ) )
  nDat$MaskDur = 500 - nDat$TarDur
  
  tmp = aggregate( rep(1,nrow(d)), 
                   list( d$PD, d$PT, d$S ), 
                   sum )
  nDat$N = round( aggregate( tmp$x, list( tmp[,1], tmp[,2] ), mean )$x )
  nDat$Y = round( nDat$P * nDat$N )
  
  # Estimate the best fitting parameters over multiple runs with 
  # dispersed starting points
  startTime = Sys.time()
  model_fit = MLE( nDat, mle_f, st_f, nRep = 10, maxit = 10000 )
  runTime = Sys.time() - startTime
  print( runTime )
  
  # Extract best-fitting parameters
  bfp = exp( model_fit$param )
  
  # Generate predicted latencies
  predLat = matrix(NA,4,2)
  for (i in 1:4) {
    pres = c(nDat$PrimeDur[i],nDat$TarDur[i],
             nDat$MaskDur[i],500)
    if ( nDat$Type[i] == -2 ) pi = c(0,2) else pi = c(2,0)
    prm = c(0.25, 0.0302, 0.15, 0.324, 0.022, 0.9844, 
            0.15, 1, 0.0294, 0.0609, 0.015)
    prm[ c(1,6,8) ] = bfp
    sim = simulate_nROUSE( pres, pi, prm )
    predLat[i,] = sim$Latencies[1:2]
  }
  
  # Generate model predictions
  pred = nROUSE_logLik( model_fit$param, nDat, estimate = F )
  
  # Uncertainty intervals
  ui = rbind(
    qbinom( .025, nDat$N, nDat$Y/nDat$N )/nDat$N,
    qbinom( .975, nDat$N, nDat$Y/nDat$N )/nDat$N
  )
  
  # Plotting dimensions
  if (!savePlot) x11(width=12);
  layout( cbind(1,2) )
  par( mar = c(4, 5, 1, 1 ) )
  blankPlot( xDim = c( -5, 1 ) )
  axis( 2, seq(0,1,.25), cex.axis = 1.5, lwd = 2 )
  x = log( c( .05, .4 ) )
  axis( 1, x,
        c( 50, 400 ), cex.axis = 1.5, lwd = 2 )
  mtext( 'Prime duration (ms)', side = 1, cex = 1.5, 
         line = 2.5 )
  mtext( 'P(Correct)', side = 2, cex = 1.5, 
         line = 2.5 )
  
  # Add observed data
  segments( x, ui[ 1, nDat$PT == 0 ],
            x, ui[ 2, nDat$PT == 0 ], col = 'orange', lwd = 2 )
  segments( x, ui[ 1, nDat$PT == 1 ],
            x, ui[ 2, nDat$PT == 1 ], col = 'red', lwd = 2 )
  lines( x, nDat$P[ nDat$PT == 0 ], col = 'orange', lwd = 2 )
  points( x, nDat$P[ nDat$PT == 0 ], col = 'orange', pch = 19, cex = 1.5 )
  lines( x, nDat$P[ nDat$PT == 1 ], col = 'red', lwd = 2 )
  points( x, nDat$P[ nDat$PT == 1 ], col = 'red', pch = 19, cex = 1.5 )
  
  abline( h = .5, lty = 2, lwd = 2 )
  
  # Add in model predictions
  lines( x, pred[ nDat$PT == 0 ], lwd = 2, lty = 2, col = 'orange' )
  points( x, pred[ nDat$PT == 0 ], pch = 21, bg = 'white', col = 'orange',
          cex = 1.5, lwd = 2 )
  lines( x, pred[ nDat$PT == 1 ], lwd = 2, lty = 2, col = 'red' )
  points( x, pred[ nDat$PT == 1 ], pch = 21, bg = 'white', col = 'red',
          cex = 1.5, lwd = 2 )
  
  # Add in legend
  legend( 'bottomleft', c('Target','Foil'),
          fill = c('orange','red'), cex = 1.5, bty = 'n' )
  legend( 'bottomright', c('Observed','nROUSE'), 
          pch = c(19,21),pt.bg = 'white',
          pt.cex = 1.5, lty = c(1,2), lwd = 2, 
          cex = 1.5, bty = 'n' )
  
  par( mar=c(0,0,0,0) )
  blankPlot()
  legend( 'left', c('Estimation of nROUSE',
                    'parameters at the group level'),
          cex = 2, bty = 'n' )
  
}

###
### Model estimation at the subject level
###
# Lookup - 03

if ( runCode[3] ) {
  
  # Create a progress bar using a base R function
  pb = txtProgressBar( min = 0, max = N, style = 3 )
  
  obs = matrix( NA, N, 12 )
  nROUSE_res = matrix( NA, N, 3 + 8 + 4 )
  
  # Loop over subjects
  for ( n in 1:N ) {
    
    # Determine aggregate performance
    sel = d$S == n
    nDat = aggregate( d$Ac[sel], list( d$PD[sel], d$PT[sel] ), mean )
    
    # Incorporate prime duration and type, as well 
    # as target flash and mask duration
    colnames( nDat ) = c( 'PD','PT','P' )
    nDat$PrimeDur = 50; nDat$PrimeDur[ nDat$PD == 1] = 2000
    nDat$Type = -2; nDat$Type[ nDat$PT == 0 ] = 2;
    
    nDat$TarDur = round( unique(d$TD[sel] ) )
    nDat$MaskDur = 500 - nDat$TarDur
    
    nDat$N = aggregate( rep(1,nrow(d))[sel], 
                        list( d$PD[sel], d$PT[sel] ),
                        sum )$x
    nDat$Y = round( nDat$P * nDat$N )
    
    # Store results for plotting
    obs[n,] = c( nDat$P, nDat$Y, nDat$N )
    
    # Obtain maximum likelihood estimates
    model_fit = MLE( nDat, mle_f, st_f, nRep = 10, maxit = 10000 )
    
    # Extract best-fitting parameters
    bfp = exp( model_fit$param )
    
    # Generate predicted latencies
    predLat = matrix(NA,4,2)
    for (i in 1:4) {
      pres = c(nDat$PrimeDur[i],nDat$TarDur[i],
               nDat$MaskDur[i],500)
      if ( nDat$Type[i] == -2 ) pi = c(0,2) else pi = c(2,0)
      prm = c(0.25, 0.0302, 0.15, 0.324, 0.022, 0.9844, 
              0.15, 1, 0.0294, 0.0609, 0.015)
      prm[ c(1,6,8) ] = bfp
      sim = simulate_nROUSE( pres, pi, prm )
      predLat[i,] = sim$Latencies[1:2]
    }
    
    # Generate model predictions
    pred = nROUSE_logLik( model_fit$param, nDat, estimate = F )
    
    # Store results for plotting
    nROUSE_res[n,] = c(
      bfp, predLat[,1], predLat[,2], pred )
    
    # Update the progress bar
    setTxtProgressBar(pb,n)
    
  }
  close(pb)
  
  # Add meaningful column names
  colnames( nROUSE_res ) = c(
    'NoiseMult', 'Inhibit', 'TemporalAtten',
    'TL_SFP','TL_LFP','TL_STP','TL_LTP',
    'FL_SFP','FL_LFP','FL_STP','FL_LTP',
    'P_SFP','P_LFP','P_STP','P_LTP')
  
  if (!savePlot) x11(width=12)
  layout( cbind( 1, 2, 3 ) )
  
  par( mar = c(5,6,3,1) )
  xa = 1:N
  yl = lowerUpper( .005, nROUSE_res[,1] )
  plot( c(0,N+1), yl, type = 'n', bty = 'l',
        xlab = 'Subjects', xaxt='n', ylab = 'Noise multiplier',
        cex.lab = 2, cex.axis = 2 )
  abline( h = .0302, lty = 2, lwd = 2 )
  points( xa, nROUSE_res[,1], pch = 19, cex = 1.5 )
  
  yl = lowerUpper( .1, nROUSE_res[,2] )
  plot( c(0,N+1), yl, type = 'n', bty = 'l',
        xlab = 'Subjects', xaxt='n', ylab = 'Inhibition',
        cex.lab = 2, cex.axis = 2 )
  abline( h = .9844, lty = 2, lwd = 2 )
  points( xa, nROUSE_res[,2], pch = 19, cex = 1.5 )
  
  yl = lowerUpper( .1, nROUSE_res[,3] )
  plot( c(0,N+1), yl, type = 'n', bty = 'l',
        xlab = 'Subjects', xaxt='n', ylab = 'Temporal attention',
        cex.lab = 2, cex.axis = 2 )
  abline( h = 1, lty = 2, lwd = 2 )
  points( xa, nROUSE_res[,3], pch = 19, cex = 1.5 )
  
  # Plot observed versus predicted performance
  if (!savePlot) x11(width=12)
  
  lyt = matrix( 1:28, 4, 7, byrow = T )
  layout( lyt )
  
  x = log( c(.05,.4) )
  
  # Loop over subjects
  for ( n in 1:N ) {
    
    par( mar = c(3,3,1,1) )
    blankPlot( xDim = c(-3.5,-.5) )
    
    axis(1,x,c('.05','.4'),cex.axis = 1.5,
         lwd = 2 )
    axis(2,c(0,.5,1),cex.axis = 1.5,
         lwd = 2 )
    
    ya = obs[n,1:4]
    segments( rep(x[1],2), ya[ nDat$PD==0],
              rep(x[2],2), ya[ nDat$PD==1], lwd = 2 )
    points( x, ya[nDat$PT == 0], pch = 19, cex = 1.5 )
    points( x, ya[nDat$PT == 1], pch = 17, cex = 1.5 )
    
    ya = nROUSE_res[n,12:15]
    segments( rep(x[1],2), ya[ nDat$PD==0],
              rep(x[2],2), ya[ nDat$PD==1], lwd = 2, lty = 2 )
    points( x, ya[nDat$PT == 0], pch = 21, cex = 1.5, bg = 'white' )
    points( x, ya[nDat$PT == 1], pch = 24, cex = 1.5, bg = 'white' )
    
    legend('bottomleft',as.character(n),bty='n',cex=1.5)
  }
  
  par( mar = c(0,0,0,0) )
  blankPlot()
  legend( 'topleft', c('Target primed','Foil primed'),
          pch = c(19,17), cex = 1.5, bty = 'n' )
  legend( 'bottomleft', c('Observed','nROUSE'),
          pch = c(19,21), cex = 1.5, bty = 'n' )
  
  if (savePlot) dev.off()
}

###
### Save covariates from nROUSE estimation
###
# Lookup - 05

if ( runCode[3] & runCode[4] ) {
  
  # Define set of covariates for the predicted latencies 
  # from the nROUSE model (raw, and the standardized 
  # inverses).
  allDat$TargetLat = 0; allDat$FoilLat = 0;
  allDat$zTargetLat = 0; allDat$zFoilLat = 0;
  
  cnd = cbind( nDat$PD, nDat$PT ) # Priming conditions
  # Loop over subjects
  for ( n in 1:N ) {
    
    # Match latencies to conditions
    for (i in 1:4) {
      sel = allDat$Subject == n & 
        allDat$PrimeDur == cnd[i,1] & allDat$PrimeType == cnd[i,2]
      
      allDat$TargetLat[sel] = nROUSE_res[n,3+i]
      allDat$FoilLat[sel] = nROUSE_res[n,3+4+i]
    }
    
    # Take the inverse and standardize the latencies on a 
    # subject to subject basis
    sel = allDat$Subject == n
    allDat$zTargetLat[sel] = scale( 1/allDat$TargetLat[sel] )
    allDat$zFoilLat[sel] = scale( 1/allDat$FoilLat[sel] )
    
  }
  
  # Include the new covariates and nROUSE estimation results 
  # with the data
  setwd( 'Data' )
  save(rawDat,allDat,N,nROUSE_res,file='Priming_offset.RData')
  setwd( orig_dir )
}

if (savePlot) {
  setwd( orig_dir )
  setwd( 'Plots' )
  dev.off()
}

setwd( orig_dir )