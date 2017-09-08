#------------------------------------#
# Mixed effects modeling of accuracy #
# Kevin Potter                       #
# Updated 06/01/2017                 #
#------------------------------------#

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
  pdf( 'Mixed_effects_accuracy.pdf' )
  setwd(orig_dir)
}

# Indicate whether to carry out model-fitting
modelYes = T

# Indicate whether debugging messages should be printed
debugging = T

# Index
# Lookup - 01:  Initial setup
# Lookup - 02:  Plot effects on P(Correct)
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
### Plot effects on P(Correct)
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
  ef = aggregate( d$Ac, list( d$Cnd ), mean )
  colnames( ef ) = c( 'Cnd', 'P' )
  xl = c(.5, 4.5 ); yl = c(0,1)
  blankPlot( xl, yl )
  segments( rep(.5,3), c(.25,.5,.75), rep(4.5,3),
            c(.25,.5,.75), col = 'grey80', lwd = 2 )
  customAxes( xl, yl, label = c( 'Condition', 'P(Correct)' ),
              inc = c( 0, .25 ) )
  axis( 1, 1:4, Cnd_labels, tick = F, line = -.5 )
  
  points( 1:4, ef$P, pch = c(15,19,17,18), cex = 1.5 )
  
  if ( debugging ) print( "Figure 2 (Set-size)" )
  # Set-size
  ef = aggregate( d$Ac, list( d$SS ), mean )
  colnames( ef ) = c( 'SS', 'P' )
  xl = c(.5, 4.5 ); yl = c(0,1)
  blankPlot( xl, yl )
  segments( rep(.5,3), c(.25,.5,.75), rep(4.5,3),
            c(.25,.5,.75), col = 'grey80', lwd = 2 )
  customAxes( xl, yl, label = c( 'Set-size', 'P(Correct)' ),
              inc = c( 0, .25 ) )
  axis( 1, 1:4, sort( unique( ef$SS ) ), tick = F, line = -.5 )
  
  lines( 1:4, ef$P, type = 'b', pch = 19, cex = 1.5 )
  
  # Lag
  if ( debugging ) print( "Figure 3 (Lag)" )
  ef = aggregate( d$Ac, list( d$LP ), mean )
  colnames( ef ) = c( 'LP', 'P' )
  xl = c(.5, 5.5 ); yl = c(0,1)
  blankPlot( xl, yl )
  segments( rep(.5,3), c(.25,.5,.75), rep(5.5,3),
            c(.25,.5,.75), col = 'grey80', lwd = 2 )
  customAxes( xl, yl, label = c( 'Lag', 'P(Correct)' ),
              inc = c( 0, .25 ) )
  axis( 1, 1:5, sort( unique( ef$LP ) )+1, tick = F, line = -.5 )
  
  lines( 1:5, ef$P, type = 'b', pch = 19, cex = 1.5,
         col = c(5,1:4) )
  
  # Reset layout for figures
  layout( cbind( 1 ) )
  
  ### Interactions ###
  
  # Compute average P(Correct) over all possible conditions
  ef = aggregate( d$Ac, list( d$SS, d$LP, d$Cnd ), mean )
  colnames( ef ) = c( 'SS', 'LP', 'Cnd', 'P' )
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
    xl = c(.5, 4.5 ); yl = c(0,1)
    blankPlot( xl, yl )
    customAxes( xl, yl, label = c( ' ', 'P(Correct)' ),
                inc = c( 0, .25 ) )
    axis( 1, 1:4, sort( unique( ef$LP ) ) + 1, tick = F, line = -.5 )
    mtext( 'Lag', side = 1, line = 2 )
    xyd = par( 'usr' )
    if ( i == 3 )
      legend( xyd[1] + (xyd[2]-xyd[1])*.05,
              xyd[3] + (xyd[4]-xyd[3])*.6,
              paste( 'Set Size:', sort( unique( ef$SS ) ) ),
              fill = 1:5, bty = 'n' )
    
    sel = ef$Cnd == i
    inc = 1
    for ( lp in sort( unique( ef$SS[sel] ) ) ) {
      cur = sel & ef$SS == lp
      x = ef$LP[cur]
      y = ef$P[cur]
      lines( x, y, type = 'b', pch = 19, col = inc )
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
  # A function that plots the effects on P(Correct) for 
  # condition, set-size, and lag, averaged over subjects.
  # Arguments:
  # fit - A lme4 fit object
  # dat - The data that was fitted
  # new - Logical; if true, a new plotting window is generated
  # Returns:
  # A figure with 4 panels showing the predicted and observed 
  # effects on P(Correct).
  
  est = predict( fit )
  pred = aggregate( as.vector( est ), list( dat$SS, 
                                            dat$LP, dat$Cnd ), 
                    function(x) mean( logistic(x) ) )
  colnames( pred ) = c( 'SS', 'LP', 'Cnd', 'P' )
  
  # Compute average P(Correct) over all possible conditions
  ef = aggregate( dat$Ac, list( dat$SS, dat$LP, dat$Cnd ), mean )
  colnames( ef ) = c( 'SS', 'LP', 'Cnd', 'P' )
  ef = ef[ ef$LP != 0, ]
  
  # Extract condition
  cur_cnd = unique( ef$Cnd )
  
  if ( new ) x11();
  
  Cnd_labels = c( '1a: NL', '1b: NC', '1c: BN', '1d: BO' )
  
  xl = c(.5, 4.5 ); yl = c(0,1)
  blankPlot( xl, yl )
  customAxes( xl, yl, label = c( ' ', 'P(Correct)' ),
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
    cur = ef$SS == lp
    x = ef$LP[cur]
    y = ef$P[cur]
    lines( x, y, type = 'b', pch = 19, col = inc )
    xp = pred$LP[cur]
    yp = pred$P[cur]
    lines( xp, yp, type = 'b', pch = 21, lty = 2, 
           col = inc, cex = 1.25 )
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
              c('S','Ac','Cnd','SS','LP') ]
    
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
    # Coefficients represent the change in the slope of the 
    # linear trend for set-size based on each type of lag
    # relative to the no-lag condition
    
    if (debugging) print( "Null model" )
    m0 = glmer( Ac ~ 1 + (1|S), 
                family = binomial('logit'), data = dtbf )
    if ( plotYes ) {
      if (!savePlot) np = T else np = F
      plot_fit( m0, dtbf, new = np )
      mtext( 'Null Model (No fixed effects)', side = 3, 
             outer = T, line = -1 )
    }
    
    if (debugging) print( "Main effect of set size" )
    m1 = glmer( Ac ~ 1 + SSn + (1|S), 
                family = binomial('logit'), data = dtbf )
    if ( plotYes ) {
      if (!savePlot) np = T else np = F
      plot_fit( m1, dtbf, new = np )
      mtext( 'Main Effect of Set Size', side = 3, 
             outer = T, line = -1 )
    }
    
    if (debugging) print( "Main effect of lag" )
    m2 = glmer( Ac ~ 1 + LP + (1|S), 
                family = binomial('logit'), data = dtbf )
    if ( plotYes ) {
      if (!savePlot) np = T else np = F
      plot_fit( m2, dtbf, new = np )
      mtext( 'Main Effect of Lag', side = 3, 
             outer = T, line = -1 )
    }
    
    if (debugging) print( "All main effects" )
    m3 = glmer( Ac ~ 1 + SSn + LP + (1|S), 
                family = binomial('logit'), data = dtbf )
    if ( plotYes ) {
      if (!savePlot) np = T else np = F
      plot_fit( m3, dtbf, new = np )
      mtext( 'All Main Effects', side = 3, 
             outer = T, line = -1 )
    }
    
    print( 'Interaction of set size and lag (1)' )
    m4 = glmer( Ac ~ 1 + SSn + SSnxL2 + SSnxL3 + SSnxL4 + (1|S), 
                family = binomial('logit'), data = dtbf )
    if ( plotYes ) {
      if (!savePlot) np = T else np = F
      plot_fit( m4, dtbf, new = np )
      mtext( 'Interaction of Set Size and Lag', side = 3, 
             outer = T, line = -1 )
    }
    
    print( 'Interaction of set size and lag (2)' )
    m5 = glmer( Ac ~ 1 + SS3i + SS4i + SS5i + 
                  SSnxL2 + SSnxL3 + SSnxL4 + (1|S), 
                family = binomial('logit'), data = dtbf )
    if ( plotYes ) {
      if (!savePlot) np = T else np = F
      plot_fit( m5, dtbf, new = np )
      mtext( 'Interaction of Set Size and Lag (2)', side = 3, 
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