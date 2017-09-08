#-----------------------------#
# nROUSE parameter estimation #
# Kevin Potter                #
# Updated 06/03/2017          #
#-----------------------------#

# Initialize script
source('F3_starting_script.R')

# Load in package for simulating nROUSE model
# install_github("rettopnivek/nROUSE")
library(nROUSE)

# Load estimates of nROUSE parameters/latencies from Dave
setwd( 'Data' )
load( 'nROUSE_est_from_Dave.RData' )
setwd( orig_dir )

# Define code segments to run
runCode = c( T, T, T )

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

# Define a function for MLE
mle_f = function( prm, dat, sum = T, priors = NULL ) { 
  
  # Calculate the log-likelihoods 
  mp = c( 2, 6, 8 )
  p_index = NULL
  if ( any( dat$p_index != 0 ) ) {
    p_index = dat$p_index[ dat$p_index != 0 ]
  }
  
  if ( is.null( p_index ) ) { prm = 1; mp = NULL } else 
    mp = mp[ p_index ]
  
  ll = nROUSE_logLik( prm, nDat, mapping = mp, estimate = F ) 
  if ( !sum ) return( ll )
  
  # Sum the log-likelihoods 
  sll = sum( ll ) 
  
  # Check for NA values 
  if ( is.na( sll ) ) sll = -Inf 
  
  return( sll ) 
}

# Define a function to generate starting values
st_f = function() {
  
  sv = c( N = .027, I = 2, Ta = .7 )
  sv = log( sv )
  
  sv = sv + runif( 3, c(-1.5,-1.5,-1),
                   c(1.5,.5,1) )
  if ( is.null( nDat$p_index ) ) out = NULL else 
    out = sv[nDat$p_index]
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
  
  # Incorporate target and mask durations
  nDat$TarDur = round( mean( nROUSE_prev_est$target_duration ) )
  nDat$MaskDur = 500 - nDat$TarDur
  
  # Compute the frequency correct
  tmp = aggregate( rep(1,nrow(d))[sel], 
                   list( d$PD[sel], d$PT[sel], d$S[sel] ), 
                   sum )
  nDat$N = round( aggregate( tmp$x, list( tmp[,1], tmp[,2] ), mean )$x )
  nDat$Y = round( nDat$P * nDat$N )
  
  # Define a scaling adjustment parameter
  pscl = c( 3.571, 0.231, 0.166 )
  # Define index for parameters
  nDat$p_index = 0; nDat$p_index[1:3] = 1:3
  
  # Clean up workspace
  rm( tmp )
  
  ### Estimate parameters of full model ###
  
  # Estimate the best fitting parameters over multiple runs with 
  # dispersed starting points
  model_fit = MLE( nDat, mle_f, st_f, nRep = 100, 
                   control = list( maxit = 10000,
                                   parscale = pscl ),
                   method = 'Nelder-Mead', 
                   hessian = F, parNames = c('Noise-multiplier',
                                             'Inhibition',
                                             'Temporal-attention' ) )
  
  # Density plot of the sum of the log-likelihoods
  sll = model_fit$track_value # Exclude infinite likelihoods
  sll = sll[ sll != -Inf ]
  sll = densityPoints( sll )
  
  # Plot empirical density of sum of the log-likelihoods 
  # over the different estimation attempts
  if ( !savePlot ) x11( width = 12 )
  
  # Plotting characteristics
  layout( cbind( 2, 1, 1, 2 ) )
  lnSz = 2
  txtSz = 2
  axSz = 1.25
  ptSz = 2
  yl = lowerUpper( .01, sll$y )
  xl = lowerUpper( 10, sll$x )
  
  blankPlot( xl, yl ) # Create blank plot
  # Axes
  abline( h = yl[1], lwd = lnSz )
  abline( v = xl[1], lwd = lnSz )
  axis( 1, seq( xl[1], xl[2], 10 ), tick = F, 
        line = -.5, cex.axis = axSz )
  mtext( 'Sum of log-likelihoods for nROUSE', side = 1, 
         line = 3, cex = txtSz )
  mtext( 'Density', side = 2, line = 3, cex = txtSz )
  
  legend( 'topleft', paste( 'Number of estimation attempts:', 
                            length( sll$x ) ),
          bty = 'n', cex = txtSz )
  
  # Add lines/points
  lines( sll$x, sll$y, lwd = lnSz )
  points( sll$x, sll$y, pch = 19, cex = ptSz )
  
  blankPlot()
  
  # Extract best-fitting parameters
  bfp = exp( coef( model_fit ) )
  
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
  pred = nROUSE_logLik( log( bfp ), nDat, predict = T )
  
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
  legend( 'bottomleft', c('Foil','Target'),
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
  
  ### Estimate models with reduced number of parameters ###
  
  mdl_pred = matrix( NA, 9, 4 )
  rownames( mdl_pred ) = c( 'O', 'D', 'N', 'I', 'TA',
                            'N-TA', 'N-I', 'I-TA',
                            'N-I-TA' )
  colnames( mdl_pred ) = c( 'S-FP','L-FP','S-TP','L-TP' )
  mdl_pred[1,] = nDat$P
  
  # 1) All parameters are fixed
  nDat$p_index = 0;
  m1 = MLE( nDat, mle_f, NULL, 
            hessian = F )
  # Generate predictions
  mdl_pred[2,] = nROUSE_logLik( 1, nDat, mapping = NULL, 
                                predict = T )
  
  # 2) Noise multiplier only
  nDat$p_index = 0; nDat$p_index[1] = 1
  m2 = MLE( nDat, mle_f, st_f, 
            control = list( lower = log( .0001 ), upper = log( .5 ) ), 
            hessian = F, parNames = c('Noise-multiplier') )
  mdl_pred[3,] = nROUSE_logLik( coef(m2), nDat, mapping = 2, 
                                predict = T )
  
  # 3) Inhibition only
  nDat$p_index = 0; nDat$p_index[1] = 2
  m3 = MLE( nDat, mle_f, st_f, 
            control = list( lower = log( .01 ), upper = log( 3 ) ), 
            hessian = F, parNames = c('Inhibition') )
  mdl_pred[4,] = nROUSE_logLik( coef(m3), nDat, mapping = 6, 
                                predict = T )
  
  # 4) Temporal attention only
  nDat$p_index = 0; nDat$p_index[1] = 3
  m4 = MLE( nDat, mle_f, st_f, 
            control = list( lower = log( .01 ), upper = log( 3 ) ),
            hessian = F, parNames = c('Temporal-attention') )
  mdl_pred[5,] = nROUSE_logLik( coef(m2), nDat, mapping = 8, 
                                predict = T )
  
  # 5) Noise multiplier and temporal attention
  nDat$p_index = 0; nDat$p_index[1:2] = c(1,3)
  m5 = MLE( nDat, mle_f, st_f, nRep = 100, 
            method = 'Nelder-Mead', 
            hessian = F, parNames = c('Noise-multiplier', 
                                      'Temporal-attention') )
  mdl_pred[6,] = nROUSE_logLik( coef(m5), nDat, mapping = c(2,8), 
                           predict = T )
  
  # 6) Noise multiplier and inhibition
  nDat$p_index = 0; nDat$p_index[1:2] = c(1,2)
  m6 = MLE( nDat, mle_f, st_f, nRep = 100, 
            method = 'Nelder-Mead', 
            hessian = F, parNames = c('Noise-multiplier', 
                                      'Inhibition') )
  mdl_pred[7,] = nROUSE_logLik( coef(m6), nDat, mapping = c(2,6), 
                           predict = T )
  
  # 7) Temporal attention and inhibition
  nDat$p_index = 0; nDat$p_index[1:2] = c(2,3)
  m7 = MLE( nDat, mle_f, st_f, nRep = 100, 
            method = 'Nelder-Mead', 
            hessian = F, parNames = c('Inhibition', 
                                      'Temporal-attention') )
  mdl_pred[8,] = nROUSE_logLik( coef(m7), nDat, mapping = c(6,8), 
                           predict = T )
  mdl_pred[9,] = nROUSE_logLik( coef(model_fit), nDat, 
                                predict = T )
  
  # Carry out model comparison
  mc = anova( m1, m2, m3, m4, m5, m6, m7, model_fit, 
                modelNames = c( 'D', 'N', 'I', 'TA', 'N-TA', 'N-I', 
                                'I-TA', 'N-I-TA' ), finite = F )
  
  # Store results
  overall_fit = list(
    obs = mdl_pred[1,],
    pred = mdl_pred[-1,],
    comp = mc )
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
    
    # Define a scaling adjustment parameter
    pscl = c( 3.571, 0.231, 0.166 )
    # Define index for parameters
    nDat$p_index = 0; nDat$p_index[1:3] = 1:3
    
    # Store results for plotting
    obs[n,] = c( nDat$P, nDat$Y, nDat$N )
    
    # Estimate the best fitting parameters over multiple runs with 
    # dispersed starting points
    model_fit = MLE( nDat, mle_f, st_f, nRep = 100, 
                     control = list( maxit = 10000,
                                     parscale = pscl ),
                     method = 'Nelder-Mead', 
                     hessian = F, parNames = c('Noise-multiplier',
                                               'Inhibition',
                                               'Temporal-attention' ) )
    
    # Extract best-fitting parameters
    bfp = exp( coef( model_fit ) )
    
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
    pred = nROUSE_logLik( log(bfp), nDat, predict = T )
    
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
  
  # Convert nROUSE results to data frame
  nROUSE_res = as.data.frame( nROUSE_res )
  
  # Include the new covariates and nROUSE estimation results 
  # with the data
  setwd( 'Data' )
  save(rawDat,allDat,N,nROUSE_res,overall_fit,file='SD_v_FC.RData')
  setwd( orig_dir )
}

if (savePlot) dev.off()