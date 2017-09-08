#--------------------#
# Model development  #
# Kevin Potter       #
# Updated 06/28/2017 #
#--------------------#

# Initialize script
source('F3_starting_script.R')

# Load in package for Bayes factor tests
library( BayesFactor )

# Indicate which code to run
runCode = c( F, T, F, T, F, F, F, F )

# Indicate whether a PDF file should be generated
savePlot = F

# Create PDF file
if ( savePlot ) {
  setwd( 'Plots' )
  pdf( 'nROUSE_crossover_effect.pdf', width = 12 )
  setwd( orig_dir )
}

# Load in most flexible but psychologically valid 
# diffusion race model
setwd( 'RT_and_choice_MLE' )
setwd( 'Estimation_results' )
load( 'DRM_SD_Full_results.RData' )
setwd( orig_dir )

# Index
# Lookup - 01:  Magnitude of cross-over interaction
# Lookup - 02:  Set up for analyzing minimum response times
# Lookup - 03:  Individual minimum response times
# Lookup - 04:  Group-level minimum response times (Forced-choice)
# Lookup - 05:  Group-level minimum response times (Same-different)

###
### Magnitude of cross-over interaction
###
# Lookup - 01

# Compute accuracy over prime type and duration for 
# 2AFC
sel = d$TaL == 'Forced-choice'
ac = aggregate( d$Ac[sel], 
                list( d$PDL[sel], d$PTL[sel], d$S[sel] ), mean )
colnames( ac ) = c( 'PD', 'PT', 'S', 'P' )
ac$Cnd = rep( 1:4, N )

# Compute interaction (difference of differences)
fp_v_tp_s = (ac$P[ ac$Cnd == 1 ] - ac$P[ ac$Cnd == 3 ])
fp_v_tp_l = (ac$P[ ac$Cnd == 2 ] - ac$P[ ac$Cnd == 4 ])
fp_v_tp = fp_v_tp_s - fp_v_tp_l

# Rank order subjects based on magnitude of interaction
int_rank_ord = order( fp_v_tp )

if ( runCode[1] ) {
  
  # Compute accuracy over prime type and duration for 
  # same-different task
  sel = d$TaL == 'Same-different'
  ac_sd = aggregate( d$Ac[sel], 
                  list( d$PDL[sel], d$PTL[sel], d$S[sel] ), mean )
  colnames( ac_sd ) = c( 'PD', 'PT', 'S', 'P' )
  ac_sd$Cnd = rep( 1:4, N )
  
  # Create panels for each subject
  if ( !savePlot ) x11( width = 12 )
  lyt = matrix(seq(1,4*7,1), 4, 7, byrow = T )
  lyt[4,5:7] = 26
  layout( lyt )
  
  # Loop over subjects
  for ( n in 1:N ) {
    
    # Create blank plot
    par( mar = c( 1, 1, 1, 1 ) )
    xl = c( .5, 2.5 )
    yl = c( 0, 1 )
    blankPlot( xl, yl )
    
    # Add axes and grid-lines
    segments( rep( xl[1], 3 ),
              seq(.25,.75,.25),
              rep( xl[2], 3 ),
              seq(.25,.75,.25),
              lwd = 2, col = 'grey80' )
    customAxes( xl, yl )
    
    # Plot accuracy (Same-different)
    ya = ac_sd$P[ ac_sd$S == int_rank_ord[n] ]
    xa = c( 1, 2, 1, 2 )
    pts = c( 21, 24, 21, 24 )
    clr = rep( c('white','grey50'), each = 2 )
    segments( rep( 1, 2 ), ya[c(1,3)], rep( 2, 2 ), ya[c(2,4)],
              lty = 1:2, col = 'grey50' )
    points( xa, ya, pch = pts, bg = clr, col = 'grey50' )
    
    # Plot accuracy (Forced-choice )
    ya = ac$P[ ac$S == int_rank_ord[n] ]
    xa = c( 1, 2, 1, 2 )
    pts = c( 21, 24, 21, 24 )
    clr = rep( c('white','black'), each = 2 )
    segments( rep( 1, 2 ), ya[c(1,3)], rep( 2, 2 ), ya[c(2,4)],
              lty = 1:2 )
    points( xa, ya, pch = pts, bg = clr )
    
    # Include subject ID and magnitude of interaction
    shft = .08
    text( xl[1], yl[1] + shft, paste( 'S', 
                                      int_rank_ord[n], sep='-' ),
          pos = 4 )
    text( xl[2], yl[1] + shft, round( fp_v_tp[int_rank_ord[n]], 2 ),
          pos = 2 )
    
  }
  
  # Title and legends
  blankPlot()
  legend( 'topleft', '2AFC accuracy by prime duration and type',
          bty = 'n', cex = 2 )
  legend( 'bottomleft', c('Foil-primed', 'Target-primed' ),
          fill = unique( clr ), bty = 'n' )
  legend( 'bottom', c('50 ms prime', '400 ms prime' ),
          pch = unique( pts ), bty = 'n', pt.bg = 'black' )
  legend( 'bottomright', c('2AFC', 'SD' ),
          fill = c('black','grey50'), bty = 'n' )
          
  
  # Create panel for each free parameter for the nROUSE model
  if ( !savePlot ) x11( width = 12 )
  # Reset margins
  par( mar = c( 4, 5, 3, 1 ) )
  layout( cbind( 1, 2, 3 ) )
  
  # Plotting characteristics
  lnSz = 2
  ptSz = 2
  txtSz = 1.5
  
  ttl = c( 'Noise multiplier', 'Inhibition',
           'Temporal attention' )
  
  adj = c( .005, .25, .25 )
  val = c( 'NoiseMult', 'Inhibit', 'TemporalAtten' )
  
  # Loop over free parameters
  for ( i in 1:3 ) {
    
    # Create a blank plot
    xl = c( .5, N + .5 )
    ya = as.numeric( nROUSE_res[,val[i]] )
    yl = lowerUpper( adj[i], ya )
    blankPlot( xl, yl )
    
    # Add axes
    customAxes( xl, yl, lnSz = lnSz )
    axis( 2, seq( yl[1], yl[2], adj[i] ), tick = F, 
          line = -1, cex.axis = txtSz )
    
    # Add points
    points( 1:N, ya[int_rank_ord], pch = 19,
            cex = ptSz )
    title( ttl[i], cex = txtSz )
    
    if ( i == 1 ) 
      mtext( 'Parameter values', side = 2, line = 2, cex = txtSz )
    
  }
  
  mtext( 'Subjects', side = 1, line = -2, outer = T,
         cex = txtSz )
  
  dtbf = data.frame(
    y = fp_v_tp,
    NM = scale( as.numeric( nROUSE_res$NoiseMult ) ),
    I = scale( as.numeric( nROUSE_res$Inhibit ) ),
    TA = scale( as.numeric( nROUSE_res$TemporalAtten ) )
  )
  
  bf = regressionBF( y ~ NM + I + TA, data = dtbf )
  if ( !savePlot ) x11( width = 12 ); plot( bf )
  
  print( bf[ 'NM + I + TA' ] / bf[ 'NM + I' ] )
  
}

###
### Set up for analyzing minimum response times
###
# Lookup - 02

# Extract minimum response times for all conditions
mt = aggregate( d$RT, list( d$PDL, d$PTL, d$CoL, 
                            d$ChL, d$TaL, d$S ),
                min )
colnames( mt ) = c( 'PD', 'PT', 'Co', 'Ch', 'Ta', 'S', 'RT' )
mt$Ac = 0; mt$Ac[ mt$Co == mt$Ch ] = 1
mt_ag = aggregate( mt$RT, list( mt$PD, mt$PT, mt$Co,
                                mt$Ac, mt$Ta ), 
                   function(x) return( c( mean = mean(x), 
                                          sd = sd(x),
                                          sem = sem(x) ) ) )
colnames( mt_ag ) = c('PD','PT','Co','Ac','Ta','T_x')
mt_ag$Cnd = rep( 1:8, 4 )
mt$Cnd = NA
for ( i in 1:8 ) {
  mt$Cnd[ mt$PD == mt_ag$PD[i] & 
            mt$PT == mt_ag$PT[i] & 
            ( mt$Co == mt_ag$Co[i] | 
                mt$Co == mt_ag$Co[i+16] ) ] = i
}

# Extract conditions
Cnd = mt_ag[,c('PD','PT','Co','Ac','Ta','Cnd')]

group_means = function( bf_object ) {
  # Purpose: 
  # A function to extract posterior samples and 
  # the predicted group means from a bayes factor object.
  # Arguments: 
  # bf_object - A bayes factor object.
  # Returns: 
  # A list with the posterior samples and means, as well as 
  # the predicted group means.
  
  # Extract posterior estimates
  post = posterior( bf_object, iterations = 10000, progress = F )
  post = as.matrix( post )
  sel = grep( 'S', colnames( post ) )
  sel = c( sel, grep( 'g', colnames( post ) ) )
  est = post[ , -sel ]
  if ( is.matrix( est ) ) { 
    out = colMeans( est )
  } else {
    out = mean( est ); names( out ) = 'mu'
  }
  vrb = names( out )
  out = out[ vrb != 'sig2' ]
  
  if ( length( vrb ) > 1 ) {
    X = aggregate( dtbf[,vrb[-1]], list( dtbf$PD, dtbf$PT, dtbf$Co,
                                         dtbf$Ac ), unique )
    X = cbind( 1, X[,-(1:4)] )
    X = as.matrix( X )
  } else {
    X = matrix( 1, 16, 1 )
  }
  pred = exp( X %*% cbind( out[ vrb ] ) )
  
  pred = data.frame( mean = pred )
  pred$Cnd = rep( 1:8, 2 )
  pred$Ch = rep( val, each = 8 )
  
  return( list( prm = out, pred = pred, post = post[,vrb] ) )
}

###
### Individual minimum response times
###
# Lookup - 03

if ( runCode[2] ) {
  
  # Forced-choice task
  task = mt$Ta == 'Forced-choice'
  
  # Create panels for each subject
  if ( !savePlot ) x11( width = 12 )
  lyt = matrix(seq(1,4*7,1), 4, 7, byrow = T )
  lyt[4,5:7] = 26
  layout( lyt )
  
  # Loop over subjects
  for ( n in 1:N ) {
    
    # Create blank plot
    par( mar = c( 1, 4, 3, 1 ) )
    xl = c( .5, 8.5 )
    yl = c( .2, .6 )
    blankPlot( xl, yl )
    
    sbj = int_rank_ord[n]
    title( paste( 'S-', sbj, ' (', 
                  round( fp_v_tp[sbj], 2 ), ')',
                  sep = '' ) )
    
    # Add axes and grid-lines
    segments( c( 2.5, 4.5, 6.5 ), rep( yl[1], 3 ),
              c( 2.5, 4.5, 6.5 ), rep( yl[2], 3 ),
              col = 'grey80' )
    customAxes( xl, yl, inc = c( 0, .2 ), axSz = 1 )
    
    pos = c( 'Left', 'Right' )
    clr = c( 'black', 'grey50' )
    for ( j in 1:2 ) {
      xa = mt$Cnd[ task & mt$S == sbj & 
                    mt$Ch == pos[j] ]
      ya = mt$RT[ task & mt$S == sbj & 
                    mt$Ch == pos[j] ]
      lines( xa, ya, type = 'b', pch = 19, col = clr[j] )
    }
    
  }
  
  # Add title and legends
  par( mar = c( 0, 0, 0, 0 ) )
  blankPlot()
  legend( 'topleft', 'Minimum response times for forced-choice',
          bty = 'n', cex = 2 )
  
  # Forced-choice task
  task = mt$Ta == 'Same-different'
  
  # Create panels for each subject
  if ( !savePlot ) x11( width = 12 )
  lyt = matrix(seq(1,4*7,1), 4, 7, byrow = T )
  lyt[4,5:7] = 26
  layout( lyt )
  
  # Loop over subjects
  for ( n in 1:N ) {
    
    # Create blank plot
    par( mar = c( 1, 4, 3, 1 ) )
    xl = c( .5, 8.5 )
    yl = c( .2, .6 )
    blankPlot( xl, yl )
    
    sbj = int_rank_ord[n]
    title( paste( 'S-', sbj, ' (', 
                  round( fp_v_tp[sbj], 2 ), ')',
                  sep = '' ) )
    
    # Add axes and grid-lines
    segments( c( 2.5, 4.5, 6.5 ), rep( yl[1], 3 ),
              c( 2.5, 4.5, 6.5 ), rep( yl[2], 3 ),
              col = 'grey80' )
    customAxes( xl, yl, inc = c( 0, .2 ), axSz = 1 )
    
    pos = c( 'Different', 'Same' )
    clr = c( 'black', 'grey50' )
    for ( j in 1:2 ) {
      xa = mt$Cnd[ task & mt$S == sbj & 
                     mt$Ch == pos[j] ]
      ya = mt$RT[ task & mt$S == sbj & 
                    mt$Ch == pos[j] ]
      lines( xa, ya, type = 'b', pch = 19, col = clr[j] )
    }
    
  }
  
  # Add title and legends
  par( mar = c( 0, 0, 0, 0 ) )
  blankPlot()
  legend( 'topleft', 'Minimum response times for same-different',
          bty = 'n', cex = 2 )
}

###
### Group-level minimum response times (Forced-choice)
###
# Lookup - 04

if ( runCode[3] ) {
  
  ### Plot group means over conditions (Forced-choice) ###
  
  task = 'Forced-choice'
  val = c('Error','Correct')
  
  pts = c( 15, 19 )
  clr = c( 'grey70', 'black' )
  lnSz = 2
  ptSz = 1.5
  txtSz = 1.5
  
  if ( !savePlot ) x11( width = 12 );
  # Create blank plot
  xl = c( .5, 8.5 ); yl = c( .2, .6 );
  blankPlot( xl, yl )
  
  # Add in axes and grid-lines
  segments( c( 2.5, 4.5, 6.5 ), rep( yl[1], 3 ),
            c( 2.5, 4.5, 6.5 ), rep( yl[2], 3 ),
            lty = 2, lwd = lnSz )
  customAxes( xl, yl, lnSz = lnSz )
  axis( 1, 1:8, rep( c('50','400'), 4 ), tick = F, line = -1,
        cex.axis = txtSz )
  axis( 1, seq(1,7,2) + .5, 
        rep( c('Foil','Target'), 2 ), tick = F, line = 0,
        cex.axis = txtSz )
  axis( 1, c(2.5,6.5), c('Left','Right'), tick = F, line = 1,
        cex.axis = txtSz )
  axis( 2, seq( yl[1], yl[2], .2 ), tick = F, line = -1.2, 
        cex.axis = txtSz )
  
  # Add points
  for ( i in 1:2 ) {
    sel = mt_ag$Ta == task & mt_ag$Ac == i - 1
    lines( 1:8, mt_ag$T_x[sel,'mean'], type = 'b', pch = pts[i],
           lwd = lnSz, col = clr[i] )
  }
  
  # Add legend
  legend( 3, yl[2], val, pch = pts, col = clr, bty = 'n', 
          cex = txtSz )
  mtext( 'Average minimum RT', side = 2, line = 2, 
         cex = txtSz )
  
  add_nR = F
  if ( add_nR ) {
    
    Lat = colMeans( nROUSE_res )[1:8 + 3]
    zLat = scale( Lat[ c(5:8,1:4) ] )
    sel = mt_ag$Ta == 'Forced-choice'
    tmp = aggregate( mt_ag$T_x[sel,1], list(
      mt_ag$PD[sel], mt_ag$PT[sel], mt_ag$Ac[sel] ),
      mean )$x
    rLat = zLat * sd( tmp ) + mean( tmp )
    lines( 1:8, rLat[ c(1:4,1:4) ], lty = 2, 
           col = 'red', lwd = lnSz )
    lines( 1:8, rLat[ c(5:8,5:8) ], lty = 1, 
           col = 'red', lwd = lnSz )
    legend( 'topright', c( 'nROUSE (Foil)', 'nROUSE (Target)' ),
            lty = 2:1, col = 'red', lwd = lnSz, cex = txtSz * .9, 
            bty = 'n' )
  }
  
  ### Linear model comparison (Forced-choice) ###
  
  # 1) Null model
  #    Intercept and random effect of subject
  
  # 2) Saturated model (No bias)
  #    Dummy-coded effects for every combination 
  #    of condition excluding left/right
  
  # 3) Saturated model
  #    Dummy-coded effects for every combination of 
  #    condition
  
  # 4) Post-hoc model
  
  # Set up data
  dtbf = mt[ mt$Ta == 'Forced-choice', ]
  # Fit log of the minimum response time
  dtbf$lRT = log( dtbf$RT )
  # Set subject up as a factor
  dtbf$S = as.factor( dtbf$S )
  
  # Design matrix for saturated model (no bias)
  W =  matrix( 0, nrow(dtbf), 8 )
  for ( i in 1:4 ) {
    sel = ( dtbf$Cnd == i | dtbf$Cnd == i + 4 ) & 
      dtbf$Ac == 0
    W[sel,i] = 1;
    sel = ( dtbf$Cnd == i | dtbf$Cnd == i + 4 ) & 
      dtbf$Ac == 1
    W[sel,i+4] = 1;
  }
  colnames( W ) = paste( 'W', 1:8, sep = '' )
  dtbf = cbind( dtbf, W )
  dtbf = as.data.frame( dtbf )
  
  # Design matrix for saturated model
  X = matrix( 0, nrow( dtbf ), 16 )
  for ( i in 1:8 ) {
    sel = dtbf$Cnd == i & dtbf$Ac == 0
    X[sel,i] = 1;
    sel = dtbf$Cnd == i & dtbf$Ac == 1
    X[sel,i+8] = 1
  }
  colnames( X ) = paste( 'X', 1:16, sep = '' )
  dtbf = cbind( dtbf, X )
  dtbf = as.data.frame( dtbf )
  
  # Covariates for post-hoc model
  
  # Target primed for 50 ms (Errors)
  dtbf$V1 = 0 
  sel = dtbf$Cnd %in% c(3,7) & dtbf$Ac == 0
  dtbf$V1[sel] = 1
  
  # Target primed for 400 ms (Errors)
  dtbf$V2 = 0
  sel = dtbf$Cnd %in% c(4,8) & dtbf$Ac == 0
  dtbf$V2[sel] = 1
  
  # Foil primed for 400 ms, Target-primed for 
  # 50/400 ms (Correct)
  dtbf$V3 = 0
  sel = ( dtbf$Cnd %in% c(2:4,6:8) & dtbf$Ac == 1 )
  dtbf$V3[sel] = 1
  
  # For quick checking of design matrices
  tmp = list( dtbf$PD, dtbf$PT, dtbf$Co, dtbf$Ch )
  
  # 1) Null model
  null = generalTestBF( lRT ~ S,
                        whichRandom = 'S', data = dtbf,
                        neverExclude = c( 'S' ) )
  
  # 2) Saturated model (No bias)
  saturated_nb = generalTestBF( lRT ~ W2 + W3 + W4 + 
                                      W5 + W6 + W7 + W8 + 
                                      S,
                         whichRandom = 'S', data = dtbf,
                         neverExclude = c('S',
                                          colnames(W) ) )
  
  # Group-level predictions for full model
  mp = group_means( saturated_nb )
  segments( 1:8 - .1, mp$pred$mean[1:8],
            1:8 + .1, mp$pred$mean[1:8],
            lwd = 2, col = 'blue' )
  segments( 1:8 - .1, mp$pred$mean[9:16],
            1:8 + .1, mp$pred$mean[9:16],
            lwd = 2, col = 'blue', lty = 2 )
  
  # 3) Saturated model
  saturated = generalTestBF( lRT ~ X2 + X3 + X4 + 
                                   X5 + X6 + X7 + X8 + 
                                   X9 + X10 + X11 + X12 + 
                                   X13 + X14 + X15 + X16 + 
                                   S,
                         whichRandom = 'S', data = dtbf,
                         neverExclude = c( 'S', 
                                           colnames( X )[-1] ) )
  # Group-level predictions
  mp = group_means( saturated )
  segments( 1:8 - .1, mp$pred$mean[1:8],
            1:8 + .1, mp$pred$mean[1:8],
            lwd = 2, col = 'red' )
  segments( 1:8 - .1, mp$pred$mean[9:16],
            1:8 + .1, mp$pred$mean[9:16],
            lwd = 2, col = 'red', lty = 2 )
  
  # Post-hoc model
  
  # Fit model using Bayes factor package
  posthoc = generalTestBF( lRT ~ V1 + V2 + V3 + S,
                         whichRandom = 'S', data = dtbf,
                         neverExclude = 
                           c( 'S', paste( 'V', 1:3, sep = '' ) ) )
  # Group-level predictions
  mp = group_means( posthoc )
  segments( 1:8 - .1, mp$pred$mean[1:8],
            1:8 + .1, mp$pred$mean[1:8],
            lwd = 2, col = 'purple' )
  segments( 1:8 - .1, mp$pred$mean[9:16],
            1:8 + .1, mp$pred$mean[9:16],
            lwd = 2, col = 'purple', lty = 2 )
  
  # Model comparison
  allbf = c(
    saturated_nb/null,
    saturated/null,
    posthoc/null
  )
  
}

###
### Group-level minimum response times (Same-different)
###
# Lookup - 05

if ( runCode[4] ) {
  
  ### Plot group means over conditions (Forced-choice) ###
  
  task = 'Same-different'
  val = c('Error','Correct')
  
  pts = c( 15, 19 )
  clr = c( 'grey70', 'black' )
  lnSz = 2
  ptSz = 1.5
  txtSz = 1.5
  
  if ( !savePlot ) x11( width = 12 );
  # Create blank plot
  xl = c( .5, 8.5 ); yl = c( .2, .6 );
  blankPlot( xl, yl )
  
  # Add in axes and grid-lines
  segments( c( 2.5, 4.5, 6.5 ), rep( yl[1], 3 ),
            c( 2.5, 4.5, 6.5 ), rep( yl[2], 3 ),
            lty = 2, lwd = lnSz )
  customAxes( xl, yl, lnSz = lnSz )
  axis( 1, 1:8, rep( c('50','400'), 4 ), tick = F, line = -1,
        cex.axis = txtSz )
  axis( 1, seq(1,7,2) + .5, 
        rep( c('Foil','Target'), 2 ), tick = F, line = 0,
        cex.axis = txtSz )
  axis( 1, c(2.5,6.5), c('Different','Same'), tick = F, line = 1,
        cex.axis = txtSz )
  axis( 2, seq( yl[1], yl[2], .2 ), tick = F, line = -1.2, 
        cex.axis = txtSz )
  
  # Add points
  for ( i in 1:2 ) {
    sel = mt_ag$Ta == task & mt_ag$Ac == i - 1
    lines( 1:8, mt_ag$T_x[sel,'mean'], type = 'b', pch = pts[i],
           lwd = lnSz, col = clr[i] )
  }
  
  # Add legend
  legend( 3, yl[2], val, pch = pts, col = clr, bty = 'n', 
          cex = txtSz )
  mtext( 'Average minimum RT', side = 2, line = 2, 
         cex = txtSz )
  
  add_nR = F
  if ( add_nR ) {
    
    Lat = colMeans( nROUSE_res )[1:8 + 3]
    zLat = scale( Lat[ c(5:8,1:4) ] )
    sel = mt_ag$Ta == 'Same-different' 
    tmp = mt_ag$T_x[c(1:4,13:16),'mean']
    rLat = zLat * sd( tmp ) + mean( tmp )
    lines( 1:4, rLat[ 1:4 ], lty = 2, 
           col = 'red', lwd = lnSz )
    lines( 5:8, rLat[ 5:8 ], lty = 1, 
           col = 'red', lwd = lnSz )
    legend( 'topright', c( 'nROUSE (Foil)', 'nROUSE (Target)' ),
            lty = 2:1, col = 'red', lwd = lnSz, cex = txtSz * .9, 
            bty = 'n' )
  }
  
  ### Linear model comparison (Forced-choice) ###
  
  # 1) Null model
  #    Intercept and random effect of subject
  
  # 2) Saturated model (No bias)
  #    Dummy-coded effects for every combination 
  #    of condition excluding left/right
  
  # 3) Saturated model
  #    Dummy-coded effects for every combination of 
  #    condition
  
  # 4) Post-hoc model
  
  # Set up data
  dtbf = mt[ mt$Ta == 'Same-different', ]
  # Fit log of the minimum response time
  dtbf$lRT = log( dtbf$RT )
  # Set subject up as a factor
  dtbf$S = as.factor( dtbf$S )
  
  # Design matrix for saturated model (no bias)
  W =  matrix( 0, nrow(dtbf), 8 )
  for ( i in 1:4 ) {
    sel = ( dtbf$Cnd == i | dtbf$Cnd == i + 4 ) & 
      dtbf$Ac == 0
    W[sel,i] = 1;
    sel = ( dtbf$Cnd == i | dtbf$Cnd == i + 4 ) & 
      dtbf$Ac == 1
    W[sel,i+4] = 1;
  }
  colnames( W ) = paste( 'W', 1:8, sep = '' )
  dtbf = cbind( dtbf, W )
  dtbf = as.data.frame( dtbf )
  
  # Design matrix for saturated model
  X = matrix( 0, nrow( dtbf ), 16 )
  for ( i in 1:8 ) {
    sel = dtbf$Cnd == i & dtbf$Ac == 0
    X[sel,i] = 1;
    sel = dtbf$Cnd == i & dtbf$Ac == 1
    X[sel,i+8] = 1
  }
  colnames( X ) = paste( 'X', 1:16, sep = '' )
  dtbf = cbind( dtbf, X )
  dtbf = as.data.frame( dtbf )
  
  # Covariates for post-hoc model
  
  # Target primed for 50 ms (Errors)
  dtbf$V1 = 0 
  sel = dtbf$Cnd %in% c(3,7) & dtbf$Ac == 0
  dtbf$V1[sel] = 1
  
  # Target primed for 400 ms (Errors)
  dtbf$V2 = 0
  sel = dtbf$Cnd %in% c(4,8) & dtbf$Ac == 0
  dtbf$V2[sel] = 1
  
  # Foil primed for 400 ms, Target-primed for 
  # 50/400 ms (Correct)
  dtbf$V3 = 0
  sel = ( dtbf$Cnd %in% c(2:4,6:8) & dtbf$Ac == 1 )
  dtbf$V3[sel] = 1
  
  # For quick checking of design matrices
  tmp = list( dtbf$PD, dtbf$PT, dtbf$Co, dtbf$Ch )
  
  # 1) Null model
  null = generalTestBF( lRT ~ S,
                        whichRandom = 'S', data = dtbf,
                        neverExclude = c( 'S' ) )
  
  # 2) Saturated model (No bias)
  saturated_nb = generalTestBF( lRT ~ W2 + W3 + W4 + 
                                  W5 + W6 + W7 + W8 + 
                                  S,
                                whichRandom = 'S', data = dtbf,
                                neverExclude = c('S',
                                                 colnames(W) ) )
  
  # Group-level predictions for full model
  mp = group_means( saturated_nb )
  segments( 1:8 - .1, mp$pred$mean[1:8],
            1:8 + .1, mp$pred$mean[1:8],
            lwd = 2, col = 'blue' )
  segments( 1:8 - .1, mp$pred$mean[9:16],
            1:8 + .1, mp$pred$mean[9:16],
            lwd = 2, col = 'blue', lty = 2 )
  
  # 3) Saturated model
  saturated = generalTestBF( lRT ~ X2 + X3 + X4 + 
                               X5 + X6 + X7 + X8 + 
                               X9 + X10 + X11 + X12 + 
                               X13 + X14 + X15 + X16 + 
                               S,
                             whichRandom = 'S', data = dtbf,
                             neverExclude = c( 'S', 
                                               colnames( X )[-1] ) )
  # Group-level predictions
  mp = group_means( saturated )
  segments( 1:8 - .1, mp$pred$mean[1:8],
            1:8 + .1, mp$pred$mean[1:8],
            lwd = 2, col = 'red' )
  segments( 1:8 - .1, mp$pred$mean[9:16],
            1:8 + .1, mp$pred$mean[9:16],
            lwd = 2, col = 'red', lty = 2 )
  
  # Post-hoc model
  
  # Fit model using Bayes factor package
  posthoc = generalTestBF( lRT ~ V1 + V2 + V3 + S,
                           whichRandom = 'S', data = dtbf,
                           neverExclude = 
                             c( 'S', paste( 'V', 1:3, sep = '' ) ) )
  # Group-level predictions
  mp = group_means( posthoc )
  segments( 1:8 - .1, mp$pred$mean[1:8],
            1:8 + .1, mp$pred$mean[1:8],
            lwd = 2, col = 'purple' )
  segments( 1:8 - .1, mp$pred$mean[9:16],
            1:8 + .1, mp$pred$mean[9:16],
            lwd = 2, col = 'purple', lty = 2 )
  
  # Model comparison
  allbf = c(
    saturated_nb/null,
    saturated/null,
    posthoc/null
  )
  
}

###
### Examination of thresholds (Same-different)
###
# Lookup - 06

if ( runCode[5] ) {
  
  # Plotting characteristics
  lnSz = 2
  txtSz = 1.5
  ptSz = 2
  
  # Extract thresholds
  kappa = raw_prm[,1:4]
  
  ord = order( rowMeans( kappa ) )
  
  # Create blank plot
  if ( !savePlot ) x11( width = 12 )
  layout( cbind( 1, 1, 2, 3 ) )
  
  xl = c( .5, 4.5 )
  yl = c( .5, 2.5 )
  blankPlot( xl, yl )
  
  # Add in axes
  customAxes( xl, yl, lnSz = lnSz )
  axis( 1, 1:4, c('50', '400', '50', '400' ), tick = F,
        line = -1.2, cex.axis = txtSz )
  axis( 1, c(1.5, 3.5), c('Same correct','Different correct'), 
        tick = F, line = .2, cex.axis = txtSz )
  axis( 2, seq(0,4,1), tick = F, line = -1, cex.axis = txtSz )
  mtext( 'Threshold', side = 2, line = 2, cex = txtSz )
  
  # Plot threshold values
  sel = 1:4
  xa = seq( -.2, .2, length = N )
  for ( n in 1:N ) {
    ya = exp( kappa[ord[n],sel] )
    lines( 1:4 + xa[n],
           ya, col = 'grey80', lty = 2, lwd = lnSz )
  }
  for ( i in 1:4 ) {
    ya = kappa[ord,sel[i]]
    points( xa + i, exp( ya ), pch = 19, col = 'grey60',
            cex = ptSz )
    segments( -.2 + i, exp( mean( ya ) ),
              .2 + i, exp( mean( ya ) ), lwd = 2 )
  }
  
  # Parameterize in terms of response caution
  xl = c( .5, 2.5 )
  yl = c( .5, 2.5 )
  blankPlot( xl, yl )
  customAxes( xl, yl )
  
  xa = seq( -.2, .2, length = N )
  for ( n in 1:N ) {
    ya = exp( kappa[ord[n],] )
    ya = c( mean( ya[c(1,3)] ), mean( ya[c(2,4)] ) )
    lines( 1:2 + xa[n], ya, col = 'grey80', 
           lty = 2, lwd = lnSz )
  }
  for ( i in 1:2 ) {
    sel = c(1,3) + 1*(i-1)
    ya = rowMeans( exp( kappa[ord,sel] ) )
    points( xa + i, ya, col = 'grey60', pch = 19, cex = ptSz )
  }
  mtext( 'Response caution', side = 2, line = 2, cex = txtSz )
  axis( 2, seq( .5, 2.5, .5 ), tick = F, line = -1, 
        cex.axis = txtSz )
  axis( 1, 1:2, c( '50', '400' ), tick = F, line = -1,
        cex.axis = txtSz )
  
  # Parameterize in terms of bias
  tmp = exp( kappa )
  tmp = cbind( tmp[,3] / ( tmp[,1] + tmp[,3] ),
               tmp[,4] / ( tmp[,2] + tmp[,4] ) )
  ord = order( rowMeans( tmp ) )
  clr = rep( 'grey', N )
  clr_ln = rep( 'grey80', N )
  clr[ which( tmp[,2] - tmp[,1] > 0 ) ] = 'grey40'
  clr_ln[ which( tmp[,2] - tmp[,1] > 0 ) ] = 'grey60'
  xl = c( .5, 2.5 )
  yl = c( .4, .6 )
  blankPlot( xl, yl )
  segments( xl[1], .5, xl[2], .5, lwd = lnSz, col = 'grey' )
  customAxes( xl, yl )
  
  xa = seq( -.2, .2, length = N )
  for ( n in 1:N ) {
    ya = exp( kappa[ord[n],] )
    ya = ya[3:4]/c( ya[1] + ya[3], ya[2] + ya[4] )
    lines( 1:2 + xa[n], ya, col = clr_ln[ord][n], 
           lty = 2, lwd = lnSz )
  }
  for ( i in 1:2 ) {
    sel = c(1,3) + 1*(i-1)
    ya = rowSums( exp( kappa[ord,sel] ) )
    ya = exp( kappa[ord,i+2] )/ya
    points( xa + i, ya, col = clr[ord], pch = 19, cex = ptSz )
  }
  mtext( 'Bias towards "Same"', side = 2, line = 2, 
         cex = txtSz )
  axis( 2, seq( .4, .6, .05 ), tick = F, line = -1, 
        cex.axis = txtSz )
  axis( 1, 1:2, c( '50', '400' ), tick = F, line = -1,
        cex.axis = txtSz )
  
  dtbf = data.frame( kappa = rep( 0, nrow( kappa ) * 
                                    ncol( kappa ) ),
                     S = as.factor( rep( 1:N, ncol( kappa ) ) ),
                     Dur = as.factor( rep( 
                       c('Short','Long','Short','Long'), each = N ) ),
                     Pos = as.factor( rep( 
                       c('Same', 'Same', 'Diff','Diff'), 
                       each = N ) ),
                     Cnd = rep( 1:4, each = N ) )
  for ( i in 1:4 ) {
    dtbf$kappa[ dtbf$Cnd == i ] = kappa[,i]
  }
  
  # Apply a repeated measures ANOVA 
  bf = anovaBF( kappa ~ Dur*Pos + S, data = dtbf, 
                whichRandom = "S" )
  
  if ( !savePlot ) x11( width = 12 )
  layout( cbind( 2, 1, 1, 2 ) )
  plot( bf )
  
  # Main effect of position over interaction
  tst = capture.output( print( bf[2] / bf[4] ) )
  print( tst[3] )
  # Main effect of position over both main effects
  tst = capture.output( print( bf[2] / bf[3] ) )
  print( tst[3] )
  # Main effect of prime duration over interaction
  tst = capture.output( print( bf[3] / bf[4] ) )
  print( tst[3] )
  
}

###
### Examination of drift rates (Same-different)
###
# Lookup - 07

if ( runCode[6] ) {
  
  # Extract drift rates
  xi = apply( raw_prm[,1:16 + 4], 2, median )
  
  # Average nROUSE latencies
  Lat = aggregate( d[ ,c('TL','FL') ], 
                   list( d$PDL, d$PTL, d$CoL ), mean )
  colnames( Lat ) = c( 'PD', 'PT', 'Co', 'TL', 'FL' )
  # Same-different task only
  Lat = Lat[ Lat$Co != 'Left' & Lat$Co != 'Right', ]
  # For same racer
  Lat$SL = c( Lat$FL[1:4], Lat$TL[5:8] )
  # For different racer
  Lat$DL = c( Lat$TL[1:4], Lat$FL[5:8] )
  
  # Rescale for same drift rates
  rsLat = scale( 1/Lat$SL ) * sd( exp( xi[1:8] ) ) + 
    mean( exp( xi[1:8] ) )
  
  # Mirroring around same drift rates to 
  # estimate different drift rates
  
  # Original
  f1 = function( beta, fit = T ) {
    if ( fit ) {
      out = exp( xi[9:16] ) - 
        exp( c( 2*beta - xi[5:8], 2*beta - xi[1:4] ) )
      out = sum( out^2 )
    } else {
      out = exp( c( 2*beta - xi[5:8], 2*beta - xi[1:4] ) )
    }
    return( out )
  }
  res = optim( c(.7), f1, fit = T, 
               method = 'Brent', 
               lower = 0, upper = 4, 
               control = list( maxit = 5000 ) )
  
  # Revised
  f4 = function( beta, fit = T ) {
    if ( fit ) {
      out = exp( xi[9:16] ) - 
        exp( c( 2*beta[1] - beta[2]*xi[5:8], 2*beta[1] - xi[1:4] ) )
      out = sum( out^2 )
    } else {
      out = exp( c( 2*beta[1] - beta[2]*xi[5:8], 2*beta[1] - xi[1:4] ) )
    }
    return( out )
  }
  res2 = optim( c(.7,.7), f4, control = list( maxit = 5000 ) )

  # Two panels
  x11( width = 12 )
  layout( cbind( 1, 2 ) )
  
  # By correct answer
  
  # Create blank plot
  xl = c( .5, 8.5 )
  yl = c( 0, 4 )
  blankPlot( xl, yl )
  
  # Axes
  segments( c( 2.5, 4.5, 6.5 ), rep( yl[1], 3 ), 
            c( 2.5, 4.5, 6.5 ), rep( yl[2], 3 ), 
            col = 'grey', lwd = 2 )
  customAxes( xl, yl )
  axis( 1, seq( 1, 7, 2 ) + .5, 
        rep( c('Foil\nprimed','Target\nprimed'), 2 ),
        tick = F, line = 0 )
  axis( 3, c( 2.5, 6.5 ), c( 'Different correct', 'Same correct' ),
        tick = F, line = 0 )
  par( xpd = T )
  legend( 'top', c( 'Same racer', 'Different racer' ),
          pch = c(22,24), pt.bg = 'black', 
          lty = 1:2, bty = 'n', inset = c(0,-.05),
          horiz = T )
  legend( 'bottom', c( 'Diffusion race', 'nROUSE (rescaled)' ),
          fill = c( 'black', 'orange' ), 
          bty = 'n', inset = c(0,-0.5),
          horiz = T )
  par ( xpd = F )
  
  # Same
  pts = rep( c(22,22), 4 )
  clr = rep( c('black','white'), 4 )
  ya = exp( xi[ c( 5:8, 1:4 ) ] )
  lines( 1:8, ya, lwd = 2, col = 'grey80' )
  points( 1:8, ya, pch = pts, bg = clr )
  
  # nROUSE predictions
  lines( 1:8, rsLat, col = 'orange', type = 'b',
         pch = 15 )
  
  # Different
  pts = rep( c(24,24), 4 )
  clr = rep( c('black','white'), 4 )
  ya = exp( xi[ c( 1:8 ) + 8 ] )
  lines( 1:8, ya, lwd = 2, lty = 2, col = 'grey80' )
  points( 1:8, ya, pch = pts, bg = clr )
  
  # By target/foil
  
  # Create blank plot
  xl = c( .5, 8.5 )
  yl = c( 0, 4 )
  blankPlot( xl, yl )
  
  # Axes
  segments( c( 2.5, 4.5, 6.5 ), rep( yl[1], 3 ), 
            c( 2.5, 4.5, 6.5 ), rep( yl[2], 3 ), 
            col = 'grey', lwd = 2 )
  customAxes( xl, yl )
  axis( 1, seq( 1, 7, 2 ) + .5, 
        rep( c('Foil\nprimed','Target\nprimed'), 2 ),
        tick = F, line = 0 )
  axis( 3, c( 2.5, 6.5 ), c( 'Target', 'Foil' ),
        tick = F, line = 0 )
  par( xpd = T )
  legend( 'top', c( 'Same racer', 'Different racer' ),
          pch = c(22,24), pt.bg = 'black', 
          lty = 1:2, bty = 'n', inset = c(0,-.05),
          horiz = T )
  par ( xpd = F )
  
  # Same racer
  pts = rep( c(22,22), 4 )
  clr = rep( c('black','white'), 4 )
  ya = exp( xi[ 1:8 ] )
  lines( 1:8, ya, lwd = 2, col = 'grey80' )
  points( 1:8, ya, pch = pts, bg = clr )
  
  # Different racer
  pts = rep( c(24,24), 4 )
  clr = rep( c('black','white'), 4 )
  ya = exp( xi[ 1:8 + 8 ] )
  lines( 1:8, ya, lwd = 2, lty = 2, col = 'grey80' )
  points( 1:8, ya, pch = pts, bg = clr )
  
  # Estimating different racer from same racer
  lines( 1:8, f1( res$par, fit = F ), type = 'b',
         col = 'blue', pch = 17, lty = 2 )
  lines( 1:8, f4( res2$par, fit = F ), type = 'b',
         col = 'red', pch = 17, lty = 2 )
  
}

###
###
###
# Lookup - 08 

if ( runCode[7] ) {
  
  # Extract drift rates
  xi = exp( raw_prm[,1:16 + 4] )
  
  x11( width = 12 )
  lyt = matrix( 1:(4*7), 4, 7, byrow = T )
  lyt[4,5:7] = 26
  layout( lyt )
  
  for ( n in 1:N ) {
    
    par( mar = c(1,1,1,1) )
    xl = c(.5,8.5)
    yl = c(0,5)
    blankPlot( xl, yl )
    customAxes( xl, yl )
    
    sel = int_rank_ord[n]
    
    Lat = as.numeric( nROUSE_res[sel,1:8+3] )
    rLat = scale( 1/Lat )*sd( xi[sel,1:8] ) + 
                                  mean( xi[sel,1:8] )
    
    lines( 1:8, xi[sel,1:8], type = 'b', pch = 19, lwd = 2, 
           cex = 1.2 )
    lines( 1:8, xi[sel,1:8+8], type = 'b', pch = 17, lwd = 2, 
           cex = 1.2, col = 'grey50' )
    lines( 1:8, rLat, col = 'orange', lty = 2, lwd = 2 )
  }
  
  x11( width = 12 )
  lyt = matrix( 1:(4*7), 4, 7, byrow = T )
  lyt[4,5:7] = 26
  layout( lyt )
  
  for ( n in 1:N ) {
    
    par( mar = c(1,1,1,1) )
    xl = c(.5,8.5)
    yl = c(-1.5,1.5)
    blankPlot( xl, yl )
    segments( xl[1], 0, xl[2], 0, lty = 2, col = 'grey', lwd = 2 )
    customAxes( xl, yl )
    
    sel = int_rank_ord[n]
    
    Lat = as.numeric( nROUSE_res[sel,1:8+3] )
    rLat = scale( 1/Lat )*sd( xi[sel,1:8] ) + 
      mean( xi[sel,1:8] )
    
    lines( 1:8, xi[sel,1:8] - rLat, type = 'b', pch = 19, lwd = 2, 
           cex = 1.2 )
  }
  
}

###
###
###
#

if ( runCode[8] ) {
  
  # Extract nROUSE latencies
  Lat = nROUSE_res[,1:8+3]
  
  # Change directory to location of modeling results
  setwd( 'RT_and_choice_MLE' )
  setwd( 'Estimation_results' )
  
  # Extract drift rates (Diffusion race)
  load( 'DRM_SD_Full_results.RData' )
  xi_DRM_SD = exp( raw_prm[ , 
                            grep( 'xi', 
                                colnames( raw_prm ) ) ] )
  # Extract drift rates (Wiener)
  load( 'WP_SD_Comp_v2_results.RData' )
  xi_WP_SD = raw_prm[ , 
                      grep( 'xi', 
                            colnames( raw_prm ) ) ]
  
  setwd( orig_dir )
  
  # Standardize nROUSE latencies
  Lat = as.matrix( Lat )
  iLat = 1/Lat
  Lat = t( apply( Lat, 1, function(x) as.numeric( scale(x) ) ) )
  iLat = t( apply( iLat, 1, function(x) as.numeric( scale(x) ) ) )
  colnames( Lat ) = c(
    'TR-SFP', 'TR-LFP', 'TR-STP', 'TR-LTP',
    'FR-SFP', 'FR-LFP', 'FR-STP', 'FR-LTP'
  )
  colnames( iLat ) = c(
    'TR-SFP', 'TR-LFP', 'TR-STP', 'TR-LTP',
    'FR-SFP', 'FR-LFP', 'FR-STP', 'FR-LTP'
  )
  # Standardize diffusion race drift rates
  xi_DRM_SD = as.matrix( xi_DRM_SD )
  xi_DRM_SD = t( apply( xi_DRM_SD, 1, 
                        function(x) as.numeric( scale(x) ) ) )
  colnames( xi_DRM_SD ) = paste( 
    rep( c( 'TR','FR', 'TR', 'FR' ), each = 4 ),
    rep( c( 'SFP','LFP', 'STP', 'LTP' ), 4 ),
    rep( c( 'S', 'D' ), each = 8 ), sep = '-' )
  # Standardize wiener drift rates
  xi_WP_SD = as.matrix( xi_WP_SD )
  xi_WP_SD = t( apply( xi_WP_SD, 1, 
                        function(x) as.numeric( scale(x) ) ) )
  colnames( xi_WP_SD ) = paste( 
    rep( c( 'TR','FR', 'TR', 'FR' ), each = 2 ),
    rep( c( 'SFP','LFP', 'STP', 'LTP' ), 2 ), sep = '-' )
  
  f = function( s, x, lab ) {
    sel = grep(lab,colnames(x))
    return( x[s,sel] )
  }
  
  ttl = c( 'Foil primed for 50 ms',
           'Foil primed for 400 ms',
           'Target primed for 50 ms',
           'Target primed for 400 ms' )
  inc = 1
  for ( lbl in c('SFP','LFP','STP','LTP') ) {
    if ( !savePlot ) x11( width = 12 )
    lyt = matrix( 1:(4 * 7), 4, 7, byrow = T )
    lyt[4,5:7] = 26
    layout( lyt )
    
    for ( j in 1:N ) {
      
      s = int_rank_ord[j]
      
      par( mar = c( 1, 1, 1, 1 ) )
      xl = c( .5, 10.5 )
      yl = c( -2, 2 )
      blankPlot( xl, yl )
      segments( c(4.5, 8.5), rep( yl[1], 2 ),
                c(4.5, 8.5), rep( yl[2], 2 ), lwd = 2, 
                col = 'grey' )
      customAxes( xl, yl )
      
      x = list(
        Lat,
        iLat,
        xi_DRM_SD,
        xi_DRM_SD,
        xi_WP_SD )
      v = list( 1:2, 1:2, 1:2, 3:4, 1:2 )
      clr = c( 'black', 'grey60', 'red', 'orange', 'purple' )
      
      for ( i in 1:length(x) ) {
        
        xa = 1:2 + 2*(i-1)
        ya = f(s,x[[i]],lbl)
        lines( xa, ya[v[[i]]], col = clr[i], lwd = 2 )
        points( xa, ya[v[[i]]], pch = c(21,22), bg = c(clr[i],'white'), col = clr[i],
                cex = 1.5 )
        
      }
      
      legend( 'topright', paste('S',s,sep='-'), bty = 'n' )
      
    }
    
    blankPlot()
    legend('left', ttl[inc], bty = 'n', cex = 1.5 )
    inc = inc + 1
    
    legend( 'bottomleft', c('Target racer','Foil racer'),
            pch = c(21,22), pt.bg = c('black','white'),
            cex = 1.5, bty = 'n' )
    
    legend( 'topright', c( 'Latencies', 'Inverse latencies', 'DRM-Same', 'DRM-Diff.', 'Wiener'),
            fill = clr,
            cex = 1.5, bty = 'n' )
    
  }

}

if ( savePlot ) dev.off()

setwd( orig_dir )