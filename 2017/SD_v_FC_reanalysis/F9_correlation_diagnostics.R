#
#
# Kevin Potter
# Updated 07/16/2017
#

# Initialize script
source('F3_starting_script.R')

# Indicate which code to run
runCode = c( T, F, F, F )

# Indicate whether a pdf file should be generated
savePlot = F

if ( savePlot ) {
  setwd( 'Plots' )
  pdf( 'Correlation_diagnostics.pdf', width = 12 )
  setwd( orig_dir )
}

quick_plot = function( lbl, pos ) {
  # Purpose: 
  # ... 
  # Arguments: 
  # ... 
  # Returns: 
  # ...
  
  vrb = paste( lbl, '_xi', sep = '' )
  
  dn = density( cor_val[,vrb] )
  trm = dn$x > 1 | dn$x < -1;
  dn$x = dn$x[ !trm ];
  dn$x = c( dn$x[1], dn$x, dn$x[length(dn$x)] )
  dn$y = dn$y[ !trm ];
  dn$y = c( 0, dn$y, 0 )
  pts = densityPoints( cor_val[,vrb] )
  scl = max( dn$y )
  polygon( .4 * -dn$y/scl + pos, dn$x, col = 'grey50', lwd = 2 )
  points( .4 * -pts$y/scl + pos, pts$x, pch = 19 )
  
  #vrb = paste( lbl, '_k', sep = '' )
  #dn = density( cor_val[,vrb] )
  #trm = dn$x > 1 | dn$x < -1;
  #dn$x = dn$x[ !trm ];
  #dn$x = c( dn$x[1], dn$x, dn$x[length(dn$x)] )
  #dn$y = dn$y[ !trm ];
  #dn$y = c( 0, dn$y, 0 )
  #pts = densityPoints( cor_val[,vrb] )
  #scl = max( dn$y )
  #polygon( .4 * dn$y/scl + pos, dn$x, col = 'grey80', lwd = 2 )
  #points( .4 * pts$y/scl + pos, pts$x, pch = 19 )
  
}

###
### Diffusion race model
###

if ( runCode[1] ) {
  
  setwd( 'RT_and_choice_MLE' )
  setwd( 'Estimation_results' )
  
  # Load in saturated model
  load( 'DRM_SD_M6_Sat_results.RData' )
  
  # Extract drift rates
  xi = all_prm[,grep('xi',colnames(raw_prm))]
  # Extract thresholds
  kappa = all_prm[,grep('kappa',colnames(raw_prm))]
  # Extract latencies
  Lat = nROUSE_res[,1:8+3]
  
  # Create data-frame with data to plot
  dtbp = data.frame( 
    # Subject index
    S = rep( 1:N, each = 8 ), 
    # Correct response
    R = rep( rep( c('Target','Foil'), each = 4 ), N ), 
    # Prime type
    PT = rep( rep( rep( c('Foil','Target'), each = 2 ), 2 ), N ), 
    # Prime duration
    PD = rep( rep( c('50','400'), 4 ), N ), 
    # Label
    L = rep( paste( 'xi-', 1:8, sep = '' ), N ) )
  
  # Initialize variables for nROUSE latencies
  dtbp$lat = 0
  # Initialize variables for drift rates
  dtbp$xi = 0
  # Initialize variables for thresholds
  dtbp$kappa = 0
  
  # Loop over subjects
  for ( n in 1:N ) {
    
    # Extract drift rates
    dtbp$xi[ dtbp$S == n ] = xi[n,1:8]
    # Extract thresholds
    dtbp$kappa[ dtbp$S == n ] = kappa[n,1:8]
    
    sel = d$S == n & d$TaL == 'Same-different'
    tmp = by( d$TL[sel], list( d$PDL[sel], d$PTL[sel], d$CoL[sel] ), 
              unique )
    tmp = as.numeric( unlist( tmp ) )
    dtbp$lat[ dtbp$S == n & dtbp$R == 'Target' ] = tmp[1:4]
    tmp = by( d$FL[sel], list( d$PDL[sel], d$PTL[sel], d$CoL[sel] ), 
              unique )
    tmp = as.numeric( unlist( tmp ) )
    dtbp$lat[ dtbp$S == n & dtbp$R == 'Foil' ] = tmp[1:4]
  }
  
  # Data frame for correlations
  cor_val = data.frame( S = 1:N, 
                        R_xi = rep( NA, N ),
                        Tau_xi = rep( NA, N ),
                        Rho_xi = rep( NA, N ),
                        R_k = rep( NA, N ),
                        Tau_k = rep( NA, N ),
                        Rho_k = rep( NA, N ) )
  
  ### Individual correlations with threshold ###
  if ( !savePlot ) x11( width = 12 )
  lyt = matrix( 1:N, round( sqrt( N ) ), 
                round( sqrt( N ) ), byrow = T )
  lyt = cbind( lyt, lyt + N )
  layout( lyt )
  
  # Drift rates
  for ( s in 1:N ) {
    par( mar = c(3,3,.5,.5) )
    xl = c( -3, 3 ); yl = c( -3, 3 )
    blankPlot( xl, yl )
    sel = dtbp$S == s
    xa = scale( 1/dtbp$lat[sel] )
    ya = scale( dtbp$xi[sel] )
    crd = data.frame( ya = ya, xa = xa )
    reg = lm( ya ~ xa, data = crd )
    nd = data.frame( xa = seq( xl[1], xl[2], length = 100 ),
                     ya = 0 )
    ln = predict( reg, nd )
    sel = ln > yl[1] & ln < yl[2]
    lines( nd$xa[sel], ln[sel], lwd = 2, col = 'grey' )
    points( xa, ya, pch = 19 )
    customAxes( xl, yl )
    
    cor_val$R_xi[s] = cor( xa, ya )
    cor_val$Tau_xi[s] = cor( xa, ya, method = 'kendall' )
    cor_val$Rho_xi[s] = cor( xa, ya, method = 'spearman' )
    legend( 'bottomright', paste( 'R = ', round( cor_val$R_xi[s], 2 ),
                                  sep = '' ), bty = 'n' )
    legend( 'topleft', as.character( s ), bty = 'n' )
    
    if ( s == 11 ) {
      mtext( paste( paste( rep('_',26), collapse = '' ),
                    ' Drift rates ', 
                    paste( rep('_',26), collapse = '' ),
                    sep = '' ),
                    side = 2, line = 1 )
    }
  }
  
  # Threshold
  for ( s in 1:N ) {
    par( mar = c(3,3,.5,.5) )
    xl = c( -3, 3 ); yl = c( -3, 3 )
    blankPlot( xl, yl )
    sel = dtbp$S == s
    xa = scale( 1/dtbp$lat[sel] )
    ya = scale( dtbp$kappa[sel] )
    crd = data.frame( ya = ya, xa = xa )
    reg = lm( ya ~ xa, data = crd )
    nd = data.frame( xa = seq( xl[1], xl[2], length = 100 ),
                     ya = 0 )
    ln = predict( reg, nd )
    sel = ln > yl[1] & ln < yl[2]
    lines( nd$xa[sel], ln[sel], lwd = 2, col = 'grey' )
    points( xa, ya, pch = 19 )
    customAxes( xl, yl )
    
    cor_val$R_k[s] = cor( xa, ya )
    cor_val$Tau_k[s] = cor( xa, ya, method = 'kendall' )
    cor_val$Rho_k[s] = cor( xa, ya, method = 'spearman' )
    legend( 'bottomright', paste( 'R = ', round( cor_val$R_k[s], 2 ),
                                  sep = '' ), bty = 'n' )
    legend( 'topleft', as.character( s ), bty = 'n' )
    
    if ( s == 11 ) {
      mtext( paste( paste( rep('_',26), collapse = '' ),
                    ' Thresholds ', 
                    paste( rep('_',26), collapse = '' ),
                    sep = '' ),
             side = 2, line = 1 )
    }
  }
  
  
  mtext( 'Inverse latencies', side = 1, outer = T, line = -1.5 )
  mtext( '    Diffusion race model', side = 1, outer = T, 
         line = -1.5, adj = 0 )
  
  # Reset margins
  par( mar = c( 4, 5, 3, .5 ) )
  
  if ( !savePlot ) x11( width = 12 )
  layout( cbind( 1 ) )
  
  xl = c( .5, 3.5 ); yl = c( -1, 1 )
  blankPlot( xl, yl )
  
  quick_plot( 'R', 1 )
  quick_plot( 'Tau', 2 )
  quick_plot( 'Rho', 3 )
  
  segments( xl[1], 0, xl[2], 0, lty = 2, lwd = 2 )
  customAxes( xl, yl )
  axis( 2, seq(-1,1,.25), tick = F, line = -1.75 )
  axis( 1, 1:3, c( "Pearson's R", "Kendall's Tau",
        "Spearman's Rho" ), tick = F, line = -1 )
  mtext( 'Density', side = 1, line = 1.5 )
  mtext( 'Correlation with nROUSE latencies', side = 2, line = 1.5 )
  
  par( xpd = NA )
  legend( 'top', c('Drift rates', 'Thresholds'),
          fill = c('grey50','grey80'), bty = 'n',
          horiz = T, inset = c(0,-.1) )
  par( xpd = T )
  
  mtext( '    Diffusion race model', side = 1, outer = T, 
         line = -1.5, adj = 0 )
  
  setwd( orig_dir )
}

###
### Wiener process
###
# Lookup - 02

if ( runCode[2] ) {
  
  setwd( 'RT_and_choice_MLE' )
  setwd( 'Estimation_results' )
  
  # Load in saturated model
  load( 'WP_SD_M5_Sat_results.RData' )
  
  # Extract drift rates
  xi = all_prm[,grep('xi',colnames(raw_prm))]
  # Extract boundary separation
  kappa = all_prm[,grep('alpha',colnames(raw_prm))]
  # Convert to threshold for 'Different' response
  tmp = all_prm[,grep('theta',colnames(raw_prm))]
  for ( n in 1:N )
    kappa[n,] = kappa[n,] * tmp[n,]
  # Extract latencies
  Lat = nROUSE_res[,1:8+3]
  
  # Create data-frame with data to plot
  dtbp = data.frame( 
    # Subject index
    S = rep( 1:N, each = 8 ), 
    # Correct response
    R = rep( rep( c('Target','Foil'), each = 4 ), N ), 
    # Prime type
    PT = rep( rep( rep( c('Foil','Target'), each = 2 ), 2 ), N ), 
    # Prime duration
    PD = rep( rep( c('50','400'), 4 ), N ), 
    # Label
    L = rep( paste( 'xi-', 1:8, sep = '' ), N ) )
  
  # Initialize variables for nROUSE latencies
  dtbp$lat = 0
  # Initialize variables for drift rates
  dtbp$xi = 0
  # Initialize variables for thresholds
  dtbp$kappa = 0
  
  # Loop over subjects
  for ( n in 1:N ) {
    
    # Extract drift rates
    dtbp$xi[ dtbp$S == n ] = xi[n,1:8]
    # Extract thresholds
    dtbp$kappa[ dtbp$S == n ] = kappa[n,1:8]
    
    sel = d$S == n & d$TaL == 'Same-different'
    tmp = by( d$TL[sel], list( d$PDL[sel], d$PTL[sel], d$CoL[sel] ), 
              unique )
    tmp = as.numeric( unlist( tmp ) )
    dtbp$lat[ dtbp$S == n & dtbp$R == 'Target' ] = tmp[1:4]
    tmp = by( d$FL[sel], list( d$PDL[sel], d$PTL[sel], d$CoL[sel] ), 
              unique )
    tmp = as.numeric( unlist( tmp ) )
    dtbp$lat[ dtbp$S == n & dtbp$R == 'Foil' ] = tmp[1:4]
  }
  
  # Data frame for correlations
  cor_val = data.frame( S = 1:N, 
                        R_xi = rep( NA, N ),
                        Tau_xi = rep( NA, N ),
                        Rho_xi = rep( NA, N ),
                        R_k = rep( NA, N ),
                        Tau_k = rep( NA, N ),
                        Rho_k = rep( NA, N ) )
  
  ### Individual correlations with threshold ###
  if ( !savePlot ) x11( width = 12 )
  lyt = matrix( 1:N, round( sqrt( N ) ), 
                round( sqrt( N ) ), byrow = T )
  lyt = cbind( lyt, lyt + N )
  layout( lyt )
  
  # Drift rates
  for ( s in 1:N ) {
    par( mar = c(3,3,.5,.5) )
    xl = c( -3, 3 ); yl = c( -3, 3 )
    blankPlot( xl, yl )
    sel = dtbp$S == s
    xa = scale( 1/dtbp$lat[sel] )
    ya = scale( dtbp$xi[sel] )
    crd = data.frame( ya = ya, xa = xa )
    reg = lm( ya ~ xa, data = crd )
    nd = data.frame( xa = seq( xl[1], xl[2], length = 100 ),
                     ya = 0 )
    ln = predict( reg, nd )
    sel = ln > yl[1] & ln < yl[2]
    lines( nd$xa[sel], ln[sel], lwd = 2, col = 'grey' )
    points( xa, ya, pch = 19 )
    customAxes( xl, yl )
    
    cor_val$R_xi[s] = cor( xa, ya )
    cor_val$Tau_xi[s] = cor( xa, ya, method = 'kendall' )
    cor_val$Rho_xi[s] = cor( xa, ya, method = 'spearman' )
    legend( 'bottomright', paste( 'R = ', round( cor_val$R_xi[s], 2 ),
                                  sep = '' ), bty = 'n' )
    legend( 'topleft', as.character( s ), bty = 'n' )
    
    if ( s == 11 ) {
      mtext( paste( paste( rep('_',26), collapse = '' ),
                    ' Drift rates ', 
                    paste( rep('_',26), collapse = '' ),
                    sep = '' ),
             side = 2, line = 1 )
    }
  }
  
  # Threshold
  for ( s in 1:N ) {
    par( mar = c(3,3,.5,.5) )
    xl = c( -3, 3 ); yl = c( -3, 3 )
    blankPlot( xl, yl )
    sel = dtbp$S == s
    xa = scale( 1/dtbp$lat[sel] )
    ya = scale( dtbp$kappa[sel] )
    crd = data.frame( ya = ya, xa = xa )
    reg = lm( ya ~ xa, data = crd )
    nd = data.frame( xa = seq( xl[1], xl[2], length = 100 ),
                     ya = 0 )
    ln = predict( reg, nd )
    sel = ln > yl[1] & ln < yl[2]
    lines( nd$xa[sel], ln[sel], lwd = 2, col = 'grey' )
    points( xa, ya, pch = 19 )
    customAxes( xl, yl )
    
    cor_val$R_k[s] = cor( xa, ya )
    cor_val$Tau_k[s] = cor( xa, ya, method = 'kendall' )
    cor_val$Rho_k[s] = cor( xa, ya, method = 'spearman' )
    legend( 'bottomright', paste( 'R = ', round( cor_val$R_k[s], 2 ),
                                  sep = '' ), bty = 'n' )
    legend( 'topleft', as.character( s ), bty = 'n' )
    
    if ( s == 11 ) {
      mtext( paste( paste( rep('_',26), collapse = '' ),
                    ' Thresholds ', 
                    paste( rep('_',26), collapse = '' ),
                    sep = '' ),
             side = 2, line = 1 )
    }
  }
  
  
  mtext( 'Inverse latencies', side = 1, outer = T, line = -1.5 )
  mtext( '    Wiener process', side = 1, outer = T, 
         line = -1.5, adj = 0 )
  
  # Reset margins
  par( mar = c( 4, 5, 3, .5 ) )
  
  if ( !savePlot ) x11( width = 12 )
  layout( cbind( 1 ) )
  
  xl = c( .5, 3.5 ); yl = c( -1, 1 )
  blankPlot( xl, yl )
  
  quick_plot( 'R', 1 )
  quick_plot( 'Tau', 2 )
  quick_plot( 'Rho', 3 )
  
  segments( xl[1], 0, xl[2], 0, lty = 2, lwd = 2 )
  customAxes( xl, yl )
  axis( 2, seq(-1,1,.25), tick = F, line = -1.75 )
  axis( 1, 1:3, c( "Pearson's R", "Kendall's Tau",
                   "Spearman's Rho" ), tick = F, line = -1 )
  mtext( 'Density', side = 1, line = 1.5 )
  mtext( 'Correlation with nROUSE latencies', side = 2, line = 1.5 )
  
  par( xpd = NA )
  legend( 'top', c('Drift rates', 'Thresholds'),
          fill = c('grey50','grey80'), bty = 'n',
          horiz = T, inset = c(0,-.1) )
  par( xpd = T )
  
  mtext( '    Wiener process', side = 1, outer = T, 
         line = -1.5, adj = 0 )
  
  setwd( orig_dir )
}

# Save figures
if ( savePlot ) dev.off()

# Return to original directory
setwd( orig_dir )