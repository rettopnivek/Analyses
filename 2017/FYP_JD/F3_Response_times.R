#------------------------------#
# Mixed effects modeling of RT #
# Kevin Potter                 #
# Updated 06/01/2017           #
#------------------------------#

# Clear workspace
rm(list = ls())

# Save current directory
orig_dir = getwd()

# Indicate whether to create figures
plotYes = T
# Indicate whether to save figures
savePlot = T
if ( savePlot & plotYes ) {
  setwd('Plots')
  pdf( 'Mixed_effects_RT.pdf' )
  setwd(orig_dir)
}

# Indicate whether to carry out model-fitting
modelYes = T

# Indicate whether debugging messages should be printed
debugging = T

# Index
# Lookup - 01:  Initial setup
# Lookup - 02:  Plot effects on mean RT
# Lookup - 03:  Mixed effects modeling

###
### Initial setup
###
# Lookup - 01

# Load in data
setwd( 'Data' )
load( 'FYP_JD.RData' )
setwd( orig_dir )

# For easy manipulation
d = trimDat
colnames( d ) = c( 'S', 'Cnd', 'Tr', 'SS', 'CS',
                   'LS', 'PT', 'PC', 'PL', 'LP',
                   'LT', 'RT1', 'Ch', 'RT2', 'Cnf',
                   'AcL', 'Ac' )

# Load in useful packages

# Load in package for mixed effects modeling
# install.packages( 'lme4' )
library( lme4 )

# Define additional useful functions
source( 'F0_Useful_functions.R' )

###
### Plot effects on mean RT
###
# Lookup - 02

if ( plotYes ) {
  
  if ( debugging ) print( "Descriptive figures (Main effects)" )
  Cnd_labels = c( '1a: NL', '1b: NC', '1c: BN', '1d: BO' )
  
  ### Main effects ###
  if (!savePlot) x11()
  layout( rbind( c( 1, 1, 2, 2 ), c( 4, 3, 3, 4 ) ) )
  
  if ( debugging ) print( "Figure 1 (Condition)" )
  # Condition
  ef = aggregate( d$RT1, list( d$Cnd, d$Ac ), mean )
  colnames( ef ) = c( 'Cnd', 'Ac', 'MRT' )
  xl = c(.5, 4.5 ); yl = lowerUpper( .25, ef$MRT )
  blankPlot( xl, yl )
  xa = seq( yl[1], yl[2] - .25, .25 )
  segments( rep( .5, length( xa ) ), xa,
            rep( 4.5, length( xa ) ), xa, col = 'grey80', lwd = 2 )
  customAxes( xl, yl, label = c( 'Condition', 'Mean RT' ),
              inc = c( 0, .25 ) )
  axis( 1, 1:4, Cnd_labels, tick = F, line = -.5 )
  legend( 'top', c('Correct','Error'), 
          pch = c(19,21), pt.bg = c('white'), bty = 'n' )
  
  points( 1:4, ef$MRT[ ef$Ac == 1 ], 
          pch = c(15,19,17,18), cex = 1.5 )
  points( 1:4, ef$MRT[ ef$Ac == 0 ], 
          pch = c(22,21,24,23), bg = 'white', cex = 1.5 )
  
  if ( debugging ) print( "Figure 2 (Set-size)" )
  # Set-size
  ef = aggregate( d$RT1, list( d$SS, d$Ac ), mean )
  colnames( ef ) = c( 'SS', 'Ac', 'MRT' )
  xl = c(.5, 4.5 ); yl = lowerUpper( .25, ef$MRT )
  blankPlot( xl, yl )
  xa = seq( yl[1], yl[2] - .25, .25 )
  segments( rep( .5, length( xa ) ), xa,
            rep( 4.5, length( xa ) ), xa, col = 'grey80', lwd = 2 )
  customAxes( xl, yl, label = c( 'Set-size', 'Mean RT' ),
              inc = c( 0, .25 ) )
  axis( 1, 1:4, sort( unique( ef$SS ) ), tick = F, line = -.5 )
  
  lines( 1:4, ef$MRT[ ef$Ac == 1 ], 
         type = 'b', pch = 19, cex = 1.5 )
  lines( 1:4, ef$MRT[ ef$Ac == 0 ], 
         type = 'b', pch = 21, bg = 'white', lty = 2, cex = 1.5 )
  
  # Lag
  if ( debugging ) print( "Figure 3 (Lag)" )
  ef = aggregate( d$RT1, list( d$LP, d$Ac ), mean )
  colnames( ef ) = c( 'LP', 'Ac', 'MRT' )
  xl = c(.5, 5.5 ); yl = lowerUpper( .25, ef$MRT )
  blankPlot( xl, yl )
  segments( rep( .5, length( xa ) ), xa,
            rep( 5.5, length( xa ) ), xa, col = 'grey80', lwd = 2 )
  customAxes( xl, yl, label = c( 'Set-size', 'Mean RT' ),
              inc = c( 0, .25 ) )
  axis( 1, 1:5, sort( unique( ef$LP ) )+1, tick = F, line = -.5 )
  
  lines( 1:5, ef$MRT[ ef$Ac == 1 ] )
  points( 1:5, ef$MRT[ ef$Ac == 1 ],
          pch = 19, col = c(5, 1:4 ), cex = 1.5 )
  lines( 1:5, ef$MRT[ ef$Ac == 0 ], lty = 2 )
  points( 1:5, ef$MRT[ ef$Ac == 0 ],
          pch = 21, col = c(5, 1:4 ), bg = 'white', cex = 1.5 )
  
  # Reset layout for figures
  layout( cbind( 1 ) )
  
  ### Interactions ###
  
  # Compute mean RT over all possible conditions
  ef = aggregate( d$RT1, list( d$SS, d$LP, d$Cnd, d$Ac ), mean )
  colnames( ef ) = c( 'SS', 'LP', 'Cnd', 'Ac', 'MRT' )
  ef = ef[ ef$LP != 0, ] # Remove lag 0 condition
  
  if (!savePlot) x11();
  # For 4 separate plotting panels, comment out 
  # the following line:
  layout( matrix( 1:4, 2, 2, byrow = T ) )
  
  if ( debugging ) print( "Descriptive figures (Interactions)" )
  for ( i in 1:4 ) {
    
    if ( debugging ) {
      print( paste( "Figure", i ) )
    }
    xl = c(.5, 4.5 ); yl = lowerUpper( .25, ef$MRT )
    blankPlot( xl, yl )
    customAxes( xl, yl, label = c( ' ', 'Mean RT' ),
                inc = c( 0, .25 ) )
    axis( 1, 1:4, sort( unique( ef$LP ) ) + 1, tick = F, line = -.5 )
    mtext( 'Lag', side = 1, line = 2 )
    xyd = par( 'usr' )
    if ( i == 1 )
      legend( 'top',
              paste( 'Set Size:', sort( unique( ef$SS ) )[1:2] ),
              fill = 1:2, bty = 'n', horiz = T,
              cex = .9 )
    if ( i == 2 )
      legend( 'top',
              paste( 'Set Size:', sort( unique( ef$SS ) )[3:4] ),
              fill = 3:4, bty = 'n', horiz = T,
              cex = .9 )
    
    sel = ef$Cnd == i
    inc = 1
    for ( lp in sort( unique( ef$SS[sel] ) ) ) {
      for ( j in 1:0 ) {
        cur = sel & ef$SS == lp & ef$Ac == j
        x = ef$LP[cur]
        y = ef$MRT[cur]
        pts = c( 21, 19 )
        ln = c( 2, 1 )
        lines( x, y, type = 'b', pch = pts[j+1], 
               bg = 'white', col = inc, lty = ln[ j + 1 ] )
      }
      inc = inc + 1
    }
    title( paste( "Experiment", Cnd_labels[i] ) )
    
  }
  
  # Reset layout
  layout( cbind( 1 ) )
  
}

###
### Mixed effects modeling
###
# Lookup - 03

# Define function that plots the predicted versus observed 
# effects
plot_fit = function( fit, dat, new = T ) {
  # Purpose:
  # A function that plots the effects on mean RT for 
  # condition, set-size, and lag, averaged over subjects.
  # Arguments:
  # fit - A lme4 fit object
  # dat - The data that was fitted
  # new - Logical; if true, a new plotting window is generated
  # Returns:
  # A figure with 4 panels showing the predicted and observed 
  # effects on mean RT.
  
  est = predict( fit )
  pred = aggregate( as.vector( est ), list( dat$SS, dat$LP, 
                                            dat$Cnd, dat$Ac ), 
                    function(x) mean( exp(x) ) )
  colnames( pred ) = c( 'SS', 'LP', 'Cnd', 'Ac', 'MRT' )
  
  # Compute average P(Correct) over all possible conditions
  ef = aggregate( dat$RT1, list( dat$SS, dat$LP, 
                                dat$Cnd, dat$Ac ), mean )
  colnames( ef ) = c( 'SS', 'LP', 'Cnd', 'Ac', 'MRT' )
  ef = ef[ ef$LP != 0, ] # Just in case, remove lag == 0 condition
  
  # Extract condition
  cur_cnd = unique( ef$Cnd )
  
  if ( new ) x11();
  
  Cnd_labels = c( '1a: NL', '1b: NC', '1c: BN', '1d: BO' )
  
  xl = c(.5, 4.5 ); yl = lowerUpper( .25, ef$MRT )
  yl[1] = yl[1] - .25;
  blankPlot( xl, yl )
  customAxes( xl, yl, label = c( ' ', 'Mean RT' ),
              inc = c( 0, .25 ) )
  axis( 1, 1:4, sort( unique( ef$LP ) ) + 1, tick = F, line = -.5 )
  mtext( 'Lag', side = 1, line = 2 )
  xyd = par( 'usr' )
  
  legend( 'top',
          paste( 'Set Size:', sort( unique( ef$SS ) ) ),
          fill = 1:4, bty = 'n', horiz = T )
  
  sel = ef$Cnd == cur_cnd
  inc = 1
  if ( debugging ) {
    string = paste( "Model fit figure (", Cnd_labels[cur_cnd], 
                    ")", sep = "" )
    print( string )
  }
  for ( lp in sort( unique( ef$SS ) ) ) {
    
    # Plot observed
    for ( j in 0:1 ) {
      cur = ef$SS == lp & ef$Ac == j
      x = ef$LP[cur]
      y = ef$MRT[cur]
      pts = c( 21, 19 )
      ln = c( 2, 1 )
      lines( x, y, type = 'b', pch = pts[j+1], col = inc,
             bg = 'white', lty = ln[j+1] )
      xp = pred$LP[cur]
      yp = pred$MRT[cur]
      pts = c( 24, 17 )
      ln = c( 2, 1 )
      lines( xp, yp, type = 'b', pch = pts[j+1], col = inc,
             bg = 'white', lty = ln[j+1] )
    }
    inc = inc + 1
  }
  title( paste( "Experiment", Cnd_labels[i] ) )
  
}

if ( modelYes ) {
  cnd_iter = 1:4
  
  for (j in cnd_iter) {
    
    if ( debugging ) {
      print( paste( "Condition", j ) )
    }
    
    # Pull out relevant data
    dtbf = d[ d$PT == 0 & d$Cnd == j & d$LP %in% 1:4, 
              c('S','RT1','Ac','Cnd','SS','LP') ]
    
    # Take log of RT
    dtbf$LRT = log( dtbf$RT1 )
    
    # Shift numeric representation of set-size
    dtbf$SSn = dtbf$SS - 2
    
    # Dummy-coded variables for SS
    dtbf$SS3i = 0
    dtbf$SS3i[ dtbf$SS == 3 ] = 1
    dtbf$SS4i = 0
    dtbf$SS4i[ dtbf$SS == 4 ] = 1
    dtbf$SS5i = 0
    dtbf$SS5i[ dtbf$SS == 5 ] = 1
    
    # Create set-size slopes by lag position (Dummy coded)
    dtbf$SSnxL1 = 0;
    dtbf$SSnxL2 = 0; dtbf$SSnxL3 = 0;
    dtbf$SSnxL4 = 0;
    for ( i in 1:4 ) {
      sel = dtbf$LP == i
      ind = grep( 'SSnxL1', colnames( dtbf ) )
      dtbf[ sel, ind + i - 1 ] = dtbf$SSn[ sel ]
    }
    
    if (debugging) print( "Null model" )
    # Have separate random intercepts for correct and 
    # error response times
    m0 = lmer( LRT ~ 1 + Ac + (Ac|S), data = dtbf )
    if ( plotYes ) {
      if (!savePlot) np = T else np = F
      plot_fit( m0, dtbf, new = np )
      mtext( 'Null Model (No fixed effects)', side = 3, 
             outer = T, line = -1 )
    }
    
    if (debugging) print( "Main effect of set size" )
    m1 = lmer( LRT ~ 1 + SSn + Ac + SSn * Ac + (Ac|S), data = dtbf )
    if ( plotYes ) {
      if (!savePlot) np = T else np = F
      plot_fit( m1, dtbf, new = np )
      mtext( 'Main Effect of Set Size', side = 3, 
             outer = T, line = -1 )
    }
    
    if (debugging) print( "Main effect of lag" )
    m2 = lmer( LRT ~ 1 + LP + Ac + LP*Ac + (Ac|S), data = dtbf )
    if ( plotYes ) {
      if (!savePlot) np = T else np = F
      plot_fit( m2, dtbf, new = np )
      mtext( 'Main Effect of Lag', side = 3, 
             outer = T, line = -1 )
    }
    
    if (debugging) print( "All main effects" )
    m3 = lmer( LRT ~ 1 + SSn + LP + Ac + 
                 SSn * Ac + LP * Ac + 
                 SSn * LP * Ac + (Ac|S), data = dtbf )
    if ( plotYes ) {
      if (!savePlot) np = T else np = F
      plot_fit( m3, dtbf, new = np )
      mtext( 'All Main Effects', side = 3, 
             outer = T, line = -1 )
    }
    
    print( 'Interaction of set size and lag' )
    m4 = lmer( LRT ~ 1 + SS3i + SS4i + SS5i + 
                 SS3i * Ac + SS4i * Ac + SS5i * Ac + 
                 SSnxL2 + SSnxL3 + SSnxL4 + 
                 SSnxL2 * Ac + SSnxL3 * Ac + SSnxL4 * Ac + 
                 (Ac|S), data = dtbf )
    if ( plotYes ) {
      if (!savePlot) np = T else np = F
      plot_fit( m4, dtbf, new = np )
      mtext( 'Interaction of Set Size and Lag', side = 3, 
             outer = T, line = -1 )
    }
    
    # Obtain summary and significance of effects
    if (0) { # Hack equivalent to block quote
      results = list("Null" = m0,
                     "SS" = m1,
                     "Lag" = m2,
                     "SS+Lag" = m3,
                     "SS*Lag" = m4,
                     "SS*Lag(2)" = m5)
      sink(paste0("model_results_", j, ".txt"))
      print(lapply( results , summary ))
      sink(NULL)
      sink(paste0("model_comp_", j, ".txt"))
      print(anova( m0, m1, m2, m3, m4, m5 ))
      sink(NULL)
    }
    
  }
  
}

if ( savePlot ) dev.off()
