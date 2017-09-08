#------------------------#
# Plots of model results #
# Kevin Potter           #
# Updated 08/21/2017     #
#------------------------#

# Initialize script
source('F3_starting_script.R')

# Indicate which model results to plot

# Diffusion race model (Same-different)
# model_sel = 'DRM_SD_M1_Null'
# model_sel = 'DRM_SD_M2_Prior'
# model_sel = 'DRM_SD_M3_Orig'
# model_sel = 'DRM_SD_M4_Thresh'
# model_sel = 'DRM_SD_M5_Full'
model_sel = 'DRM_SD_M7_PM'

# Wiener process
# model_sel = 'WP_SD_M1_Null'
# model_sel = 'WP_SD_M2_Prior'
# model_sel = 'WP_SD_M3_Drift'
# model_sel = 'WP_SD_M4_Thresh'
# model_sel = 'WP_SD_M5_Sat'

# Indicate which code to run
runCode = c( T, T, T, F )

# Indicate whether a pdf file should be generated
savePlot = T

if ( savePlot ) {
  setwd( 'Plots' )
  pdf( paste( model_sel, '_figures.pdf', sep = '' ) )
  setwd( orig_dir )
}

# Index
# Lookup - 01:  Initial setup
# Lookup - 02:  Predicted versus observed (joint CDF)
# Lookup - 03:  Observed minus predicted (joint CDF)
# Lookup - 04:  Correlation with nROUSE

###
### Initial setup
###
# Lookup - 01

# Change directory to where model results are saved
setwd( 'RT_and_choice_MLE' )
setwd( 'Estimation_results' )

# Load in reference data
load( 'DRM_SD_M1_Null_results.RData' )
mdl_fit_ref = mdl_fit$LogLik

# Load in data
load( paste( model_sel, 'results.RData', sep = '_' ) )

# Plotting characteristics
lnSz = 2
ptSz_obs = 1.5
ptSz_pred = 1
txtSz = 1.2

residual_boxplot = function( e, all_k, all_cnd, yl = NULL ) {
  # Purpose: 
  # A function that generates a set of boxplots over 
  # conditions for the residuals with median correct 
  # and error RTs and accuracy.
  # Arguments: 
  # e       - An array of residuals
  # all_k   - The specific set of residuals to plot 
  #           (1-5 residuals for 'different/right'
  #            RT quantiles; 6-10 residuals for 
  #            'same/left' RT quantiles; 11 residuals 
  #            accuracy).
  # all_cnd - A vector of the 8 conditions to plot
  
  # Extract current residuals
  cur_resid = matrix( NA, N, 8 )
  for ( i in 1:8 ) {
    cur_resid[,i] = e[ , all_cnd[i], all_k[i] ]
  }
  # Compute a bootstrap sample to estimate uncertainty
  boot = resid_bootstrap( cur_resid, 10000 )
  ui = apply( boot, 2, quantile, prob = c( .025, .975 ),
              na.rm = T )
  
  xl = c( .5, 8.5 )
  
  if ( all( k < 11 ) ) {
    # Response times
    if ( is.null( yl ) ) {
      yl = lowerUpper( .25, 
                       na.omit( as.vector( cur_resid ) ) )
    }
  } else {
    # P( Y == 1 )
    if ( is.null( yl ) ) {
      yl = lowerUpper( .12, 
                       na.omit( as.vector( cur_resid ) ) )
    }
  }
  blankPlot( xl, yl )
  customAxes( xl, yl, pos = 2, lnSz = lnSz )
  
  # Sampling error
  ya = numeric(8)
  for ( i in 1:8 ) ya[i] = sampling_error$ui[1,all_cnd[i],all_k[i]]
  lines( 1:8, ya,
         lty = 2, col = 'grey40', lwd = lnSz )
  ya = numeric(8)
  for ( i in 1:8 ) ya[i] = sampling_error$ui[2,all_cnd[i],all_k[i]]
  lines( 1:8, ya,
         lty = 2, col = 'grey40', lwd = lnSz )
  
  if ( all( k < 11 ) ) {
    # Response times
    axis( 2, seq( yl[1], yl[2], .1 ), 
          cex.axis = txtSz, tick = F, line = -1 )
  } else {
    # P( Y == 1 )
    axis( 2, seq( -.1, .1, .1 ), 
          cex.axis = txtSz, tick = F, line = -1 )
  }
  lbls = rep( 'Foil\nprimed', 8 )
  lbls[ ms$Cnd$PT != 'Foil-primed' ] = 'Target\nprimed'
  axis( 1, c( 1, 3, 5, 7 ) + .5,
        lbls[ c( 1, 3, 5, 7 ) ],
        tick = F, cex.axis = txtSz, 
        line = 0 )
  
  segments( c( 2.5, 4.5, 6.5 ), rep( yl[1], 3 ), 
            c( 2.5, 4.5, 6.5 ), rep( yl[2], 3 ), 
            col = 'grey80', lwd = lnSz )
  
  lbls = unique( ms$Cnd$Co[ all_cnd ] )
  axis( 3, c( 2.5, 6.5 ), lbls, tick = F, line = -1.4,
        cex.axis = txtSz )
  # text( 4.5, yl[2], lbls[1], cex = txtSz, pos = 2 )
  # text( 4.5, yl[2], lbls[2], cex = txtSz, pos = 4 )
  
  # Color-code durations and priming
  clr = rep( c( 'grey60', 'white' ), 4 )
  
  # Interquartile range
  iqr = apply( cur_resid, 2, quantile, prob = c( .25, .75 ),
               na.rm = T )
  # Range
  rng = apply( cur_resid, 2, range, na.rm = T )
  
  segments( 1:8, iqr[1,], 1:8, iqr[2,], lwd = lnSz )
  # points( 1:8, rng[1,], pch = 25, bg = clr )
  # points( 1:8, rng[2,], pch = 24, bg = clr )
  
  width = .15
  for ( i in 1:8 ) {
    
    polygon( i + width * c(-1,-1,1,1),
             ui[c(1,2,2,1),i],
             lwd = lnSz, col = clr[i] )
    
  }
  xa = colMeans( cur_resid, na.rm = T )
  segments( 1:8 - width, xa, 1:8 + width, lwd = 2 )
  
  
  segments( xl[1], 0, xl[2], 0, lty = 3, lwd = lnSz )
  
}

###
### Predicted versus observed (joint CDF)
###
# Lookup - 02

if ( runCode[1] ) {
  
  # Progression of plots
  # Left, Right, Different, Same
  
  # If only the same-different task is being plotted
  if ( length( grep( 'SD', model_sel ) ) > 0 ) {
    ms = extract_cnd( d, task = 'Same-different' )
    
  # If only forced-choice task is being plotted
  } else if ( length( grep( 'FC', model_sel ) ) > 0 ) {
    ms = extract_cnd( d, task = 'Forced-choice' )
    
  # If both tasks are being plotted
  } else {
    ms = extract_cnd( d )
  }
  lbls = unique( ms$Cnd$Co )
  
  # Determine number of plotting windows
  n_win = nrow( ms$Cnd ) / 4
  
  # Loop over conditions, creating figures 
  # with 4 panels
  inc = 1
  for ( i in 1:n_win ) {
    
    if (!savePlot) x11( width = 6 )
    layout( matrix( 1:4, 2, 2, byrow = T ) )
    
    for ( cnd in inc:(inc+3) ) {
      
      # Generate blank plot
      xl = c( .1, 1.7 )
      yl = c( 0, 1 )
      par( mar = c( 4, 1, 2, 4 ) )
      blankPlot( xl, yl )
      # Add in axes
      customAxes( xl, yl, pos = c(1,4), lnSz = lnSz, type = 'extend' )
      axis( 1, seq( .2, 1.4, .4 ), 
            tick = F, cex.axis = txtSz, line = -1 )
      axis( 4, seq( 0, 1, .25 ), cex.axis = txtSz, 
            tick = F, line = -1 )
      
      if ( cnd %in% c( 1, 5, 9, 13 ) ) {
        legend( 'topright', '50 ms', bty = 'n', 
                cex = txtSz, inset = c(.1,0) )
        text( xl[1], diff(yl)/2 + yl[1], 
              'Foil-primed', srt = 90,
              cex = txtSz )
      }
      if ( cnd %in% c( 2, 6, 10, 14 ) ) {
        legend( 'topright', '400 ms', bty = 'n', 
                cex = txtSz, inset = c(.1,0) )
        legend( 'topleft', c( 'Observed', 'Model' ),
                pch = c( 22, 22 ), 
                pt.bg = c( 'white', 'grey' ),
                cex = txtSz, 
                bty = 'n' )
      }
      if ( cnd %in% c( 3, 7, 11, 15 ) ) {
        text( xl[1], diff(yl)/2 + yl[1], 
              'Target-primed', srt = 90,
              cex = txtSz )
      }
      if ( cnd %in% c( 4, 8, 12, 16 ) ) {
        legend( 'topleft', c( 'Same', 'Different' ),
                pch = c( 22, 21 ), 
                pt.bg = c( 'white', 'black' ),
                cex = txtSz, 
                bty = 'n' )
      }
      
      # Define function for central tendency
      T_x = function(x) median( na.omit( x ) )
      
      # Quantiles
      prb = seq( .1, .9, .2 )
      
      # Extract predicted results and determine measure of 
      # central tendency
      avg_pred = apply( pred[ , cnd, ], 2, T_x )
      p1 = avg_pred[11]
      p0 = 1 - p1
      
      # Plot predicted response times
      lines( avg_pred[1:5], prb * p1, lwd = lnSz )
      lines( avg_pred[6:10], prb * p0, lty = 2, lwd = lnSz )
      points( avg_pred[1:5], prb * p1, pch = 21, bg = 'grey',
              cex = ptSz_pred )
      points( avg_pred[6:10], prb * p0, pch = 22, bg = 'grey',
              cex = ptSz_pred )
      
      # Plot predicted choice proportions
      segments( xl[2], p1, xl[2] - .05, p1, lwd = lnSz )
      segments( xl[2], p0, xl[2] - .05, p0, lwd = lnSz )
      
      # Extract observed results and determine measure 
      # of central tendency
      avg_obs = apply( obs[ , cnd, ], 2, T_x )
      
      # Plot observed response times
      p1 = avg_obs[11]
      p0 = 1 - p1
      points( avg_obs[1:5], prb * p1, pch = 21, bg = 'black',
              cex = ptSz_obs )
      points( avg_obs[6:10], prb * p0, pch = 22, bg = 'white', 
              cex = ptSz_obs )
      
      # Plot observed choice proportions
      points( xl[2] - .05, p1, pch = 21, bg = 'black', 
              cex = ptSz_obs )
      points( xl[2] - .05, p0, pch = 22, bg = 'white', 
              cex = ptSz_obs )
    }
    
    mtext( 'Time (s)', side = 1, outer = T, 
           cex = txtSz, line = -1.5 )
    mtext( 'Cumulative probability', side = 4, outer = T, 
           cex = txtSz, line = -1.5 )
    mtext( paste( lbls[i], 'correct' ), side = 3, outer = T, 
           cex = txtSz, line = -1.3 )
    
    inc = inc + 4
  }
  
}

###
### Observed minus predicted (joint CDF)
###
# Lookup - 03

if ( runCode[2] ) {
  
  # Progression of plots
  # 2AFC, Same-different
  
  # Compute all residuals
  e = obs - pred
  
  # If only the same-different task is being plotted
  if ( length( grep( 'SD', model_sel ) ) > 0 ) {
    ms = extract_cnd( d, task = 'Same-different' )
    
    # If only forced-choice task is being plotted
  } else if ( length( grep( 'FC', model_sel ) ) > 0 ) {
    ms = extract_cnd( d, task = 'Forced-choice' )
    
    # If both tasks are being plotted
  } else {
    ms = extract_cnd( d )
  }
  
  # Determine number of plotting windows
  n_win = nrow( ms$Cnd ) / 8
  
  for ( i in 1:n_win ) {
    
    # Create figure with panels for median, 
    # accuracy, 10%, and 90% quantiles for 
    # each task
    if (!savePlot) x11( width = 6 )
    layout( rbind( c(1,1,2,2), c(4,3,3,5) ) )
    
    all_cnd = 1:8 + 8 * ( i - 1 )
    
    par( mar = c( 3, 4, 3, .5 ) )
    
    # Median ( Correct )
    k = 3
    all_k = rep( k, 8 )
    all_k[ ms$Cnd$Co[all_cnd] == 'Same' | 
             ms$Cnd$Co[all_cnd] == 'Left' ] = k + 5
    residual_boxplot( e, all_k, all_cnd )
    mtext( 'Median RT (Correct)', side = 3, 
           line = 1, cex = txtSz * .8 )
    
    # Median ( Error )
    k = 3
    all_k = rep( k, 8 )
    all_k[ ms$Cnd$Co[all_cnd] == 'Different' | 
             ms$Cnd$Co[all_cnd] == 'Right' ] = k + 5
    residual_boxplot( e, all_k, all_cnd )
    mtext( 'Median RT (Error)', side = 3, 
           line = 1, cex = txtSz * .8 )
    
    # Accuracy
    k = 11
    all_k = rep( k, 8 )
    residual_boxplot( e, all_k, all_cnd )
    mtext( 'P(Correct)', side = 3, 
           line = 1, cex = txtSz * .8 )
    
    # Legend
    par( mar = c(0,0,0,0) )
    blankPlot()
    legend( 0,.5, paste( 'Primed for ', 
      unique( ms$Cnd$PD ) * 1000, ' ms', sep = '' ),
            fill = c( 'grey50', 'white' ), 
            bty = 'n', cex = txtSz )
    legend( 0,.4, 'Sampling error',
            lty = 2, col = 'grey40', bty = 'n',
            cex = txtSz, inset = -.4 )

    # Empty regions
    blankPlot()
    
  }
  
}

###
### Correlation with nROUSE
###
# Lookup - 04

if ( runCode[3] ) {
  
  if ( length( grep( 'xi', colnames( raw_prm ) ) ) >= 8 ) {
    # Create data frame
    xi = data.frame( 
      S = rep( 1:N, each = 8 ), 
      R = rep( rep( c('Target','Foil'), each = 4 ), N ), 
      PT = rep( rep( rep( c('Foil','Target'), each = 2 ), 2 ), N ), 
      PD = rep( rep( c('50','400'), 4 ), N ), 
      L = rep( paste( 'xi-', 1:8, sep = '' ), N ) )
    
    # Initialize variables for nROUSE latencies
    xi$Lat = 0;
    
    # Loop over subjects
    for ( n in 1:N ) {
      # Extract drift rates
      xi$X[ xi$S == n ] = all_prm[ n, paste( 'xi-', 1:8, sep = '' ) ]
      
      if ( length( grep( 'SD', model_sel ) ) > 0 ) 
        sel = d$S == n & d$TaL == 'Same-different' else
          sel = d$S == n
        tmp = by( d$TL[sel], list( d$PDL[sel], d$PTL[sel], d$CoL[sel] ), 
                  unique )
        tmp = as.numeric( unlist( tmp ) )
        xi$Lat[ xi$S == n & xi$R == 'Target' ] = tmp[1:4]
        tmp = by( d$FL[sel], list( d$PDL[sel], d$PTL[sel], d$CoL[sel] ), 
                  unique )
        tmp = as.numeric( unlist( tmp ) )
        xi$Lat[ xi$S == n & xi$R == 'Foil' ] = tmp[1:4]
    }
    
    cor_val = data.frame( S = 1:N, 
                          R = rep( NA, N ),
                          Tau = rep( NA, N ),
                          Rho = rep( NA, N ) )
    
    if ( !savePlot ) x11()
    layout( matrix( 1:N, round( sqrt( N ) ), 
                    round( sqrt( N ) ), byrow = T ) )
    
    for ( s in 1:N ) {
      par( mar = c(3,3,.5,.5) )
      xl = c( -3, 3 ); yl = c( -3, 3 )
      blankPlot( xl, yl )
      sel = xi$S == s
      xa = scale( 1/xi$Lat[sel] )
      ya = scale( xi$X[sel] )
      crd = data.frame( ya = ya, xa = xa )
      reg = lm( ya ~ xa, data = crd )
      nd = data.frame( xa = seq( xl[1], xl[2], length = 100 ),
                       ya = 0 )
      ln = predict( reg, nd )
      sel = ln > yl[1] & ln < yl[2]
      lines( nd$xa[sel], ln[sel], lwd = 2, col = 'grey' )
      points( xa, ya, pch = 19 )
      customAxes( xl, yl )
      
      cor_val$R[s] = cor( xa, ya )
      cor_val$Tau[s] = cor( xa, ya, method = 'kendall' )
      cor_val$Rho[s] = cor( xa, ya, method = 'spearman' )
      legend( 'bottomright', paste( 'R = ', round( cor_val$R[s], 2 ),
                                    sep = '' ), bty = 'n' )
      legend( 'topleft', as.character( s ), bty = 'n' )
    }
    mtext( 'Inverse latencies', side = 1, outer = T, line = -1.5 )
    mtext( 'Drift rates', side = 2, outer = T, line = -1.5 )
  }
  
  if ( length( grep( 'kappa', colnames( raw_prm ) ) ) >= 8 ) {
    # Create data frame
    kappa = data.frame( 
      S = rep( 1:N, each = 8 ), 
      R = rep( rep( c('Target','Foil'), each = 4 ), N ), 
      PT = rep( rep( rep( c('Foil','Target'), each = 2 ), 2 ), N ), 
      PD = rep( rep( c('50','400'), 4 ), N ), 
      L = rep( paste( 'kappa-', 1:8, sep = '' ), N ) )
    
    # Initialize variables for nROUSE latencies
    kappa$Lat = 0;
    
    # Loop over subjects
    for ( n in 1:N ) {
      # Extract drift rates
      kappa$X[ kappa$S == n ] = all_prm[ n, 
                                         paste( 
                                           'kappa-', 1:8, sep = '' ) ]
      
      if ( length( grep( 'SD', model_sel ) ) > 0 ) 
        sel = d$S == n & d$TaL == 'Same-different' else
          sel = d$S == n
        tmp = by( d$TL[sel], list( d$PDL[sel], d$PTL[sel], d$CoL[sel] ), 
                  unique )
        tmp = as.numeric( unlist( tmp ) )
        kappa$Lat[ kappa$S == n & kappa$R == 'Target' ] = tmp[1:4]
        tmp = by( d$FL[sel], list( d$PDL[sel], d$PTL[sel], d$CoL[sel] ), 
                  unique )
        tmp = as.numeric( unlist( tmp ) )
        kappa$Lat[ kappa$S == n & kappa$R == 'Foil' ] = tmp[1:4]
    }
    
    cor_val = data.frame( S = 1:N, 
                          R = rep( NA, N ),
                          Tau = rep( NA, N ),
                          Rho = rep( NA, N ) )
    
    if ( !savePlot ) x11()
    layout( matrix( 1:N, round( sqrt( N ) ), 
                    round( sqrt( N ) ), byrow = T ) )
    
    for ( s in 1:N ) {
      par( mar = c(3,3,.5,.5) )
      xl = c( -3, 3 ); yl = c( -3, 3 )
      blankPlot( xl, yl )
      sel = kappa$S == s
      xa = scale( 1/kappa$Lat[sel] )
      ya = scale( kappa$X[sel] )
      crd = data.frame( ya = ya, xa = xa )
      reg = lm( ya ~ xa, data = crd )
      nd = data.frame( xa = seq( xl[1], xl[2], length = 100 ),
                       ya = 0 )
      ln = predict( reg, nd )
      sel = ln > yl[1] & ln < yl[2]
      lines( nd$xa[sel], ln[sel], lwd = 2, col = 'grey' )
      points( xa, ya, pch = 19 )
      customAxes( xl, yl )
      
      cor_val$R[s] = cor( xa, ya )
      cor_val$Tau[s] = cor( xa, ya, method = 'kendall' )
      cor_val$Rho[s] = cor( xa, ya, method = 'spearman' )
      legend( 'bottomright', paste( 'R = ', round( cor_val$R[s], 2 ),
                                    sep = '' ), bty = 'n' )
      legend( 'topleft', as.character( s ), bty = 'n' )
    }
    mtext( 'Inverse latencies', side = 1, outer = T, line = -1.5 )
    mtext( 'Thresholds', side = 2, outer = T, line = -1.5 )
  }
  
  if ( ( length( grep( 'alpha', colnames( raw_prm ) ) ) >= 8 ) & 
       ( length( grep( 'theta', colnames( raw_prm ) ) ) >= 8 ) ) {
    
    # Create data frame
    kappa = data.frame( 
      S = rep( 1:N, each = 8 ), 
      R = rep( rep( c('Target','Foil'), each = 4 ), N ), 
      PT = rep( rep( rep( c('Foil','Target'), each = 2 ), 2 ), N ), 
      PD = rep( rep( c('50','400'), 4 ), N ), 
      L = rep( paste( 'kappa-', 1:8, sep = '' ), N ) )
    
    # Initialize variables for nROUSE latencies
    kappa$Lat = 0;
    
    # Loop over subjects
    for ( n in 1:N ) {
      # Extract drift rates
      tmp = all_prm[ n, paste( 'theta-', 1:8, sep = '' ) ]
      tmp = tmp * all_prm[ n, paste( 'alpha-', 1:8, sep = '' ) ]
      kappa$X[ kappa$S == n ] = tmp
      
      if ( length( grep( 'SD', model_sel ) ) > 0 ) 
        sel = d$S == n & d$TaL == 'Same-different' else
          sel = d$S == n
        tmp = by( d$TL[sel], list( d$PDL[sel], d$PTL[sel], d$CoL[sel] ), 
                  unique )
        tmp = as.numeric( unlist( tmp ) )
        kappa$Lat[ kappa$S == n & kappa$R == 'Target' ] = tmp[1:4]
        tmp = by( d$FL[sel], list( d$PDL[sel], d$PTL[sel], d$CoL[sel] ), 
                  unique )
        tmp = as.numeric( unlist( tmp ) )
        kappa$Lat[ kappa$S == n & kappa$R == 'Foil' ] = tmp[1:4]
    }
    
    cor_val = data.frame( S = 1:N, 
                          R = rep( NA, N ),
                          Tau = rep( NA, N ),
                          Rho = rep( NA, N ) )
    
    if ( !savePlot ) x11()
    layout( matrix( 1:N, round( sqrt( N ) ), 
                    round( sqrt( N ) ), byrow = T ) )
    
    for ( s in 1:N ) {
      par( mar = c(3,3,.5,.5) )
      xl = c( -3, 3 ); yl = c( -3, 3 )
      blankPlot( xl, yl )
      sel = kappa$S == s
      xa = scale( 1/kappa$Lat[sel] )
      ya = scale( kappa$X[sel] )
      crd = data.frame( ya = ya, xa = xa )
      reg = lm( ya ~ xa, data = crd )
      nd = data.frame( xa = seq( xl[1], xl[2], length = 100 ),
                       ya = 0 )
      ln = predict( reg, nd )
      sel = ln > yl[1] & ln < yl[2]
      lines( nd$xa[sel], ln[sel], lwd = 2, col = 'grey' )
      points( xa, ya, pch = 19 )
      customAxes( xl, yl )
      
      cor_val$R[s] = cor( xa, ya )
      cor_val$Tau[s] = cor( xa, ya, method = 'kendall' )
      cor_val$Rho[s] = cor( xa, ya, method = 'spearman' )
      legend( 'bottomright', paste( 'R = ', round( cor_val$R[s], 2 ),
                                    sep = '' ), bty = 'n' )
      legend( 'topleft', as.character( s ), bty = 'n' )
    }
    mtext( 'Inverse latencies', side = 1, outer = T, line = -1.5 )
    mtext( 'Thresholds', side = 2, outer = T, line = -1.5 )
  }
  
}

###
### Model convergence
###
# Lookup - 05

if ( runCode[4] ) {
  
  lkhd_ratio = mdl_fit$LogLik - mdl_fit_ref
  
  if ( !savePlot ) x11()
  layout( cbind(1) )
  par( mar = c( 4, 5, 3, .5 ) )
  xl = c( .5, N + .5 )
  yl = lowerUpper( 10, lkhd_ratio )
  yl[1] = min( -10, yl[1] )
  blankPlot( xl, yl )
  customAxes( xl, yl )
  axis( 2, round( seq( yl[1], yl[2], length = 5 ), 2 ),
        tick = F, line = -1, cex.axis = txtSz )
  mtext( 'Log-likelihood difference', side = 2, line = 2, 
         cex = txtSz )
  mtext( 'Subjects', side = 1, line = 0, 
         cex = txtSz )
  par( xpd = T )
  legend( 'bottom', c('Converged','Did not converge'),
          fill = c( 'black', 'grey' ), bty = 'n',
          horiz = T, cex = txtSz, inset = c( 0, -.1 ) )
  par( xpd = F )
  
  clr = rep( 'black', N )
  clr[ is.na( mdl_fit$Converged ) ] = 'grey'
  points( 1:N, lkhd_ratio, pch = 19, col = clr, 
          cex = ptSz_obs )
  
}

if ( savePlot ) dev.off()
setwd( orig_dir )