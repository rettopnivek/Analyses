#---------------------#
# Descriptive results #
# Kevin Potter        #
# Updated 12/15/2016  #
#---------------------#

# Clear workspace
rm( list = ls() )

# Save current directory
orig_dir = getwd()

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
  pdf( 'Descriptive_results.pdf', width = 12, height = 6 )
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
library(nROUSE)

# Define useful functions
source('F2_useful_functions.R')

# Load in data
setwd( 'Data' )
load( 'Gabor_pilot_Dec.RData' )
setwd( orig_dir )

# For easy manipulation
d = allDat
colnames( d ) = c('S','RT','Ch','Ac','Co','PT','PD','FC','TC',
                  'PP','A','BT','PS','Cnd','DP','BN')

###
### Probability-latency plots for initial practice
###
# Lookup - 02

if ( runCode[1] ) {
  
  if (!savePlot) x11( width = 12 )
  
  lyt = matrix( 1:16, 4, 4, byrow = T )
  lyt[ lyt == 16 ] = 15
  layout( lyt )
  
  for ( n in 1:N ) {
    
    # Extract practice trials
    tmp = d[ d$S == n & d$BT == 0, ]
    
    # Define function for desired test statistic
    T_x = function(x) quantile(x,seq(.1,.9,.2))
    
    # Create a blank plot
    par( mar = c(3,3,.5,.5) )
    blankPlot( c(0,1), c(.5,2) )
    axis( 1, seq(0,1,.25), cex.axis = 1.25 )
    abline( v = .5, lty = 2 )
    axis(2, seq(.5,2,.5), cex.axis = 1.25 )
    
    # Add points
    clr = c(); adj = seq( 1, .2, length = 4 )
    for (i in 1:4) clr = c( clr, rgb(adj[i],0,0,1) )
    cnd = covCreate( cbind( tmp$TC ) )
    rt = tmp$RT; ac = tmp$Ac
    plt = list( pch = 19, col = clr )
    pvt_points(rt,ac,cvrt=cnd,T_x=T_x,plt=plt)
    
    legend( 'topright', paste( 'Subject', n ), cex = 1.25,
            bty = 'n' )
  }
  par( mar = c(0,0,0,0) )
  blankPlot()
  legend('topleft','Target contrast',bty='n',cex=2)
  legend('left',as.character(sort(unique(tmp$TC))),
         fill=clr,bty='n',cex=1.25)
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
  
  if (!savePlot) x11(width=12)
  
  lyt = matrix( 1:16, 4, 4, byrow = T )
  lyt[ lyt == 16 ] = 15
  layout( lyt )
  
  for ( n in 1:N) {
    
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
  
  if (!savePlot) x11(width=12)
  
  lyt = matrix( 1:16, 4, 4, byrow = T )
  lyt[ lyt == 16 ] = 15
  layout( lyt )
  
  yl = lowerUpper( .05, d$TC[ d$BT == 2 ] )
  
  pst = 1:6
  for ( n in 1:N ) {
    
    sel = d$BT == 2 & d$S == n
    tc = aggregate( d$TC[sel], list( d$BN[sel] ), unique )$x
    
    par( mar = c( 3, 3, .5, .5 ) )
    blankPlot( c(1,6), yl )
    axis( 2, seq( yl[1], yl[2], .1 ) )
    abline( h = c( .05, .1, .15 ), col = 'grey70', lty = 2 )
    lines( pst, tc, lwd = 2 )
    points( pst, tc, pch = 19, cex = 1.25 )
    
    legend( 'bottomleft', paste( 'Subject', n ), bty = 'n' )
    
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
  
  if (!savePlot) x11( width = 12 )
  layout( cbind( 1, 2, 3 ) )
  
  txtSz = 2
  lnSz = 2.5
  ptSz = 2
  
  par( mar = c( 4, 5, 3, 1 ) )
  
  cndSel = list(
    d$DP > 0 & d$DP < 3,
    d$DP > 2 & d$DP < 5,
    d$DP > 4 )
  pchVal = list(
    c(19,19,21,21),
    c(15,15,22,22),
    c(17,17,24,24) )
  
  ttl = c(
    'Neither primed',
    'Foil primed',
    'Target primed'
  )
  
  for ( i in 1:3 ) {
    
    # Create a blank plot
    blankRTplot(tDim=c(.6,.8),ver='PvT',bty='n',
                cex.axis = 1.5, cex.lab = 1.5 )
    title( ttl[i] )
    
    segments( seq(.2,.8,.2), rep(.6,4),
              seq(.2,.8,.2), rep(.77,4),
              col = 'grey70', lty = 2 )
    
    # Calculate PvT values for neither primed condition
    cnd = cndSel[[i]]
    rt = d$RT[ cnd ]; ac = d$Ac[ cnd ]; sb = d$S[ cnd ]; cv = d$DP[ cnd ]
    cv = cv - ( min(cv) - 1 )
    outF = pvt_points( rt, ac, cv, grp = sb, 
                       opt = list( out = T, draw = F ) )
    segments( outF$pv$x[1:2], outF$pv$y[1:2],
              outF$pv$x[3:4], outF$pv$y[3:4],
              col = 'grey', lty = 1:2, lwd = lnSz )
    
    # Add in PvT values
    plt = list( pch = pchVal[[i]], bg = 'white' )
    pvt_points( rt, ac, cv, grp = sb, plt = plt, cex = ptSz )
    
    if ( i == 1 ) {
      legend( 'topleft', c('50 ms prime','400 ms prime'),
              lwd = lnSz, lty = 1:2, col = 'grey', bty = 'n', 
              cex = txtSz )
    }
    
    if ( i == 2 ) {
      legend( 'topleft', c('Error','Correct'),
              pch = c(15,22), pt.bg = 'white',
              pt.cex = ptSz, bty = 'n', 
              cex = txtSz )
    }
    
  }
 
}

# Return to original directory
if ( savePlot ) { setwd( 'Plots' ); dev.off() } # Close plotting window
setwd( orig_dir )