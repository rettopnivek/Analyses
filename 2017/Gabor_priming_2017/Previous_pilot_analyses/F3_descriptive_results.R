#---------------------#
# Descriptive results #
# Kevin Potter        #
# Updated 12/15/2016  #
#---------------------#

# Clear workspace
rm( list = ls() )

# Save current directory
orig_dir = getwd()

# Indicator for whether to fit pilot data or current study
Pilot = T
if ( Pilot ) dName = 'Pilot' else dName = 'Current'

# Index
# Lookup - 01:  Load in useful packages and data
# Lookup - 02:  Probability-latency plots for initial practice
# Lookup - 03:  Calibration trials
# Lookup - 04:  Target contrast over blocks
# Lookup - 05:  Prime duration by type

# Indicate which code segments to run
runCode = c( T, T, T, T )

# Indicate whether to save a pdf file
savePlot = T
if (savePlot) {
  setwd( 'Plots' )
  pdf( paste( dName, 'descriptive_results.pdf', sep = '_' ), width = 12, height = 6 )
  setwd( orig_dir )
}

###
### Load in useful packages and data
###
# Lookup - 01

# For geting github packages
# install.packages(devtools)
# library(devtools)

# Integration of C++ and R
# install.packages(Rcpp)
library(Rcpp)

# Miscellanous functions for modeling/plotting
# install_github("rettopnivek/utilityf")
library(utilityf)

# Distribution functions for assorted sequential sampling models
# install_github("rettopnivek/seqmodels")
library(seqmodels)

# Functions for plotting response time and choice data
# install_github("rettopnivek/rtplots")
library(rtplots)

# Functions for plotting response time and choice data
# install_github("rettopnivek/nROUSE")
# library(nROUSE)

# Load in data
setwd( 'Data' )
if ( Pilot ) load( 'Gabor_pilot_9-17.RData' ) else load( 'Gabor_current_9-17.RData' ) 
setwd( orig_dir )

# For easy manipulation
d = allDat
colnames( d ) = c('S','RT','Ch','Ac','Co','PT','PD','FC','TC',
                  'PP','A','BT','ID','PTL','CoL','ChL','PS','PSL','BN','Cnd','DxP')

###
### Probability-latency plots for initial practice
###
# Lookup - 02

if ( runCode[1] ) {
  
  new_win = F
  for ( n in 1:N ) {
    
    if ( n > 1 & ( (n-1) %% 6 != 0 ) ) new_win = F else new_win = T
    if ( new_win ) {
      # Obtain dimensions for figure layout
      if ( !savePlot ) x11( width = 12 )
      layout( matrix( 1:6, 2, 3, byrow = T ) )
    }
    
    # Extract practice trials
    tmp = d[ d$S == n & d$BT == 0, ]
    
    # Define function for desired test statistic
    f = function(x) quantile(x,seq(.1,.9,.2))
    
    # Create a blank plot
    par( mar = c(3,3,.5,.5) )
    xl = c(0,1); yl = c(.2,1)
    blankPlot( xl, yl )
    axis( 1, seq(0,1,.25), cex.axis = 1.25 )
    abline( v = .5, lty = 2 )
    axis(2, seq(.2,1,.2), cex.axis = 1.25 )
    
    # Add points
    obj = rtplots( tmp, c( 'RT', 'Ac', 'S', 'TC' ), type = 'PVT', f = f )
    pch = matrix( 19, 8, 5 );
    clr = c( rgb( 0, 0, 0, 1 ),
             rgb( .33, 0, 0, 1 ),
             rgb( .66, 0, 0, 1 ),
             rgb( 1, 0, 0, 1 ) )
    col = matrix( rep( clr, 2 ), 8, 5 )
    pvt_plot( obj, pch = pch, col = col )
    # pvt_points(rt,ac,cvrt=cnd,T_x=T_x,plt=plt)
    
    legend( 'topright', paste( 'Subject', n ), cex = 1.25,
            bty = 'n' )
  }
  par( mar = c(0,0,0,0) )
  blankPlot()
  legend('left',c( 'Target/foil contrast', 
                   as.character(sort(unique(tmp$TC))/unique(tmp$FC)) ), 
         fill = c( NA, clr ), bty='n',cex=1.25)
  legend('topright','Probability v. Latency',bty='n',cex=2)
  legend('right','Initial practice',bty='n',cex=2)
  
}

###
### Calibration trials
###
# Lookup - 03

if ( runCode[2] ) {
  
  # Extract calibration trials
  sel = d$BT == 1
  cD = d[ sel, ]
  
  yl = lowerUpper( .05, cD$TC )
  pst = rep( 1:8, each = 10 )
  par( mar = c( 4, 5, 3, .5 ) )
  
  new_win = F
  for ( n in 1:N ) {
    
    if ( n > 1 & ( (n-1) %% 6 != 0 ) ) new_win = F else new_win = T
    if ( new_win ) {
      # Obtain dimensions for figure layout
      if ( !savePlot ) x11( width = 12 )
      layout( matrix( 1:6, 2, 3, byrow = T ) )
    }
    
    par( mar = c( 3, 5, .5, .5 ) )
    blankPlot( c(1,8), yl )
    abline( h = seq( .05, yl[2] - .05, .05 ), col = 'grey70', lty = 2 )
    axis( 2, seq( yl[1], yl[2], .05 ) )
    lines( pst, cD$TC[ cD$S == n ], lwd = 2 )
    points( pst, cD$TC[ cD$S == n ], pch = 19, cex = 1.25 )
    
    legend( 'topright', paste( 'Subject', n ), bty = 'n' )
    
  }
  mtext( 'Target contrast', side = 2, outer = T, 
         line = -2, cex = 1.5 )
  mtext( 'Segments of 10 trials', side = 1, outer = T, 
         line = -1.6, cex = 1.5 )
  
  par( mar = c( 0, 0, 0, 0 ) )
  blankPlot()
  legend( 'center', 'Calibration trials', cex = 3, bty = 'n' )
  
}

###
### Target contrast over blocks
###
# Lookup - 04

if ( runCode[3] ) {
  
  yl = lowerUpper( .05, d$TC[ d$BT == 2 ] )
  
  pst = 1:6
  new_win = F
  for ( n in 1:N ) {
    
    if ( n > 1 & ( (n-1) %% 6 != 0 ) ) new_win = F else new_win = T
    if ( new_win ) {
      # Obtain dimensions for figure layout
      if ( !savePlot ) x11( width = 12 )
      layout( matrix( 1:6, 2, 3, byrow = T ) )
    }
  
    sel = d$BT == 2 & d$S == n
    tc = aggregate( d$TC[sel], list( d$BN[sel] ), unique )$x
    
    par( mar = c( 4, 5, 3, .5 ) )
    blankPlot( c(1,6), yl )
    axis( 2, seq( yl[1], yl[2], .1 ) )
    abline( h = c( .05, .1, .15 ), col = 'grey70', lty = 2 )
    lines( pst, tc, lwd = 2 )
    points( pst, tc, pch = 19, cex = 1.25 )
    
    legend( 5, yl[1] + .01, paste( 'Subject', n ), bty = 'n' )
    
  }
  mtext( 'Target contrast', side = 2, outer = T, 
         line = -2, cex = 1.5 )
  mtext( 'Blocks', side = 1, outer = T, 
         line = -1.6, cex = 1.5 )
  
  par( mar = c( 0, 0, 0, 0 ) )
  blankPlot()
  legend( 'center', 'Adaption over blocks', cex = 3, bty = 'n' )
  
}

###
### Prime duration by type
###
# Lookup - 05

if ( runCode[4] ) {
  
  par( mar = c( 4, 5, 3, 1 ) )
  if ( !savePlot ) x11( width = 12 )
  layout( cbind( 1, 2, 3 ) )
  
  ### Accuracy
  
  Y = aggregate( d$Ac, list( d$PD, d$PTL ), mean )
  colnames( Y ) = c( 'PD', 'PT', 'P' )
  
  xl = c( .5, 4.5 ); yl = c( 0, 1 );
  blankPlot( xl, yl )
  
  segments( xl[1], .5, xl[2], .5, lty = 2, col = 'grey80' )
  customAxes( xl, yl )
  axis( 1, 1:4, unique( Y$PD ), tick = F, line = -1.25 )
  axis( 2, seq(0,1,.25), tick = F, line = -1.25 )
  mtext( 'Prime duration (ms)', side = 1, line = 1.5 )
  mtext( 'Average accuracy', side = 2, line = 1.5 )
  
  inc = 1
  for ( cnd in unique( Y$PT ) ) {
    lines( 1:4, Y$P[ Y$PT == cnd ], type = 'b', lwd = 2, pch = 19,
           col = inc )
    inc = inc + 1
  }
  
  legend( 'topright', unique( Y$PT ), fill = 1:3, bty = 'n' )
  
  ### Median correct RT
  
  sel = d$Ac == 1
  Y = aggregate( d$RT[sel], list( d$PD[sel], d$PTL[sel] ), median )
  colnames( Y ) = c( 'PD', 'PT', 'P' )
  
  xl = c( .5, 4.5 )
  yl = lowerUpper( .2, Y$P )
  blankPlot( xl, yl )
  
  segments( xl[1], .5, xl[2], .5, lty = 2, col = 'grey80' )
  customAxes( xl, yl )
  axis( 1, 1:4, unique( Y$PD ), tick = F, line = -1.25 )
  axis( 2, seq(yl[1],yl[2],.1), tick = F, line = -1.25 )
  mtext( 'Prime duration (ms)', side = 1, line = 1.5 )
  mtext( 'Median correct RT (s)', side = 2, line = 1.5 )
  
  inc = 1
  for ( cnd in unique( Y$PT ) ) {
    lines( 1:4, Y$P[ Y$PT == cnd ], type = 'b', lwd = 2, pch = 19,
           col = inc )
    inc = inc + 1
  }
  
  ### Median error RT
  
  sel = d$Ac == 0
  Y = aggregate( d$RT[sel], list( d$PD[sel], d$PTL[sel] ), median )
  colnames( Y ) = c( 'PD', 'PT', 'P' )
  
  xl = c( .5, 4.5 )
  blankPlot( xl, yl )
  
  segments( xl[1], .5, xl[2], .5, lty = 2, col = 'grey80' )
  customAxes( xl, yl )
  axis( 1, 1:4, unique( Y$PD ), tick = F, line = -1.25 )
  axis( 2, seq(yl[1],yl[2],.1), tick = F, line = -1.25 )
  mtext( 'Prime duration (ms)', side = 1, line = 1.5 )
  mtext( 'Median error RT (s)', side = 2, line = 1.5 )
  
  inc = 1
  for ( cnd in unique( Y$PT ) ) {
    lines( 1:4, Y$P[ Y$PT == cnd ], type = 'b', lwd = 2, pch = 19,
           col = inc )
    inc = inc + 1
  }
  
}

# Return to original directory
if ( savePlot ) { setwd( 'Plots' ); dev.off() } # Close plotting window
setwd( orig_dir )