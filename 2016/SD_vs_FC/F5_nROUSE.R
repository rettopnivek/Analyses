#-----------------------------#
# nROUSE parameter estimation #
# Kevin Potter                #
# Updated 11/26/2016          #
#-----------------------------#

# Initialize script
source('F4_starting_script.R')

# Load in package for simulating nROUSE model
# install_github("rettopnivek/nROUSE")
library(nROUSE)

# Load in data
setwd( 'Data' )
load( 'SD_v_FC.RData' )
# Load estimates of nROUSE parameters/latencies from Dave
load( 'nROUSE_est_from_Dave.RData' )
setwd( orig_dir )

# For easy manipulation
d = allDat
colnames( d ) = c('S','Cnd','Ac','RT','E','PD',
                  'Ta','Co','PT','Ch')

# Define code segments to run
runCode = c( F, T, T )

# Indicate if a pdf should be generated
savePlot = T
if (savePlot) {
  setwd( 'Plots' )
  pdf( 'nROUSE_results.pdf', width = 12, height = 6 )
  setwd( orig_dir )
}

# Index
# Lookup - 01:  Define useful functions
# Lookup - 02:  Model estimation at the group level
# Lookup - 03:  Model estimation at the subject level
# Lookup - 04:  Save covariates from nROUSE estimation

###
### Define useful functions
###
# Lookup - 01

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
### Model estimation at the group level
###
# Lookup - 02

if ( runCode[1] ) {
  
  # Determine aggregate performance
  sel = d$Ta == 0
  nDat = aggregate( d$Ac[sel], list( d$PD[sel], d$PT[sel] ), mean )
  
  # Incorporate prime duration and type
  colnames( nDat ) = c( 'PD','PT','P' )
  nDat$PrimeDur = 50; nDat$PrimeDur[ nDat$PD == 1] = 400
  nDat$Type = -2; nDat$Type[ nDat$PT == 1 ] = 2;
  
  nDat$TarDur = round( mean( nROUSE_prev_est$target_duration ) )
  nDat$MaskDur = 500 - nDat$TarDur
  
  tmp = aggregate( rep(1,nrow(d))[sel], 
                   list( d$PD[sel], d$PT[sel], d$S[sel] ), 
                   sum )
  nDat$N = round( aggregate( tmp$x, list( tmp[,1], tmp[,2] ), mean )$x )
  nDat$Y = round( nDat$P * nDat$N )
  
  # Estimate the best fitting parameters over multiple runs with 
  # dispersed starting points
  startTime = Sys.time()
  model_fit = MLE( nDat, mle_f, st_f, nRep = 10 )
  runTime = Sys.time() - startTime
  print( runTime )
  
  # Extract best-fitting parameters
  bfp = exp( model_fit$param )
  
  # Generate predicted latencies
  predLat = matrix(NA,4,2)
  for (i in 1:4) {
    pres = c(nDat$PrimeDur[i],nDat$TarDur[i],
             nDat$MaskDur[i],500)
    pi = c(0,0);
    pi[ nDat$PT[i] + 1 ] = -2*(1-nDat$PT[i]) + 2*nDat$PT[i]
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
          bty = 'l', cex = 2, bty = 'n' )
  
}

###
### Model estimation at the subject level
###
# Lookup - 03

if ( runCode[2] ) {
  
  # Create a progress bar using a base R function
  pb = txtProgressBar( min = 0, max = N, style = 3 )
  
  obs = matrix( NA, N, 12 )
  nROUSE_res = matrix( NA, N, 3 + 8 + 4 )
  
  # Loop over subjects
  for ( n in 1:N ) {
    
    # Determine aggregate performance
    sel = d$Ta == 0 & d$S == n
    nDat = aggregate( d$Ac[sel], list( d$PD[sel], d$PT[sel] ), mean )
    
    # Incorporate prime duration and type, as well 
    # as target flash and mask duration
    colnames( nDat ) = c( 'PD','PT','P' )
    nDat$PrimeDur = 50; nDat$PrimeDur[ nDat$PD == 1] = 400
    nDat$Type = -2; nDat$Type[ nDat$PT == 1 ] = 2;
    
    nDat$TarDur = round( nROUSE_prev_est$target_duration[n] )
    nDat$MaskDur = 500 - nDat$TarDur
    
    nDat$N = aggregate( rep(1,nrow(d))[sel], 
                        list( d$PD[sel], d$PT[sel] ),
                        sum )$x
    nDat$Y = round( nDat$P * nDat$N )
    
    # Store results for plotting
    obs[n,] = c( nDat$P, nDat$Y, nDat$N )
    
    # Obtain maximum likelihood estimates
    model_fit = MLE( nDat, mle_f, st_f, nRep = 20 )
    
    # Extract best-fitting parameters
    bfp = exp( model_fit$param )
    
    # Generate predicted latencies
    predLat = matrix(NA,4,2)
    for (i in 1:4) {
      pres = c(nDat$PrimeDur[i],nDat$TarDur[i],
               nDat$MaskDur[i],500)
      pi = c(0,0);
      pi[ nDat$PT[i] + 1 ] = -2*(1-nDat$PT[i]) + 2*nDat$PT[i]
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
  lyt[4,5:7] = 26
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
  legend( 'topright', c('Observed','nROUSE'),
          pch = c(19,21), cex = 1.5, bty = 'n' )
  
  if (savePlot) dev.off()
}

###
### Save covariates from nROUSE estimation
###
# Lookup - 04

if ( runCode[2] & runCode[3] ) {
  
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
  save(rawDat,allDat,N,nROUSE_res,file='SD_v_FC.RData')
  setwd( orig_dir )
}