#---------------------------#
# Plots of RT distributions #
# Kevin Potter              #
# Updated 01/30/2017        #
#---------------------------#

# Clear workspace
rm( list = ls() )

# Save current directory
orig_dir = getwd()

# Define code segments to run
runCode = c( T, T, T )

# Indicate whether to save a pdf file
savePlot = T
if (savePlot) {
  setwd( 'Plots' )
  fName = 'RT_distribution_plots.pdf'
  pdf( fName, width = 12, height = 6 )
  rm( fName )
  setwd( orig_dir )
}

# Index
# Lookup - 01:  Load in useful packages and data

###
### Load in useful packages and data
###
# Lookup - 01

# For geting github packages
# install.packages(devtools)
# library(devtools)

# Miscellanous functions for modeling/plotting
# install_github("rettopnivek/utilityf")
library(utilityf)

# Load in package for easy plotting of RT distributions
# install_github("rettopnivek/rtplots")
library(rtplots)

# Load in useful functions
source( 'F1_useful_functions.R' )

# Load in data
setwd( 'Data' )
load( 'Recog_mem.RData' )
setwd( orig_dir )

# Remove missing responses
sel = which( is.na( d$RT ) )
d = d[ -sel, ]
N = length( unique( d$S ) )

###
###  Latency by accuracy by conditions
###
# Lookup - 01

if ( runCode[1] ) {
  
  if (!savePlot) x11(width=12)
  lyt = cbind( 1, 2 )
  layout( lyt )
  
  txtSz = 2
  lnSz = 2.5
  ptSz = 2
  
  ### 2AFC ###
  
  yl = c(1.1,1.5)
  
  par( mar=c(5,6,2,1) )
  blankPlot(xDim=c(0,1),yDim=yl)
  axis( 1, seq(0,1,.25), cex.axis = txtSz, lwd = lnSz )
  axis( 2, seq(yl[1],yl[2],.2), cex.axis = txtSz, lwd = lnSz )
  mtext( 'Accuracy', side = 1, cex = txtSz, line = 3.25 )
  mtext( 'MRT (s)', side = 2, cex = txtSz, line = 2.75 )
  
  # Calculate PvT values
  rt = d$RT; ac = d$Ac; sb = d$S; cv = d$Cnd
  outF = pvt_points( rt, ac, cv, grp = sb, 
                     opt = list( out = T, draw = F ) )
  segments( outF$pv$x[1:4], outF$pv$y[1:4],
            outF$pv$x[5:8], outF$pv$y[5:8],
            col = 'grey', lty = 1:2, lwd = lnSz )
  
  # Add in PvT values
  plt = list( pch = c(19,17,21,24), bg = 'white' )
  pvt_points( rt, ac, cv, grp = sb, plt = plt, cex = ptSz )
  
  # Legends
  par( mar=c(1,1,1,1) )
  blankPlot()
  
  legend( 0, .7, c('Selective retrieval', 'Baseline'),
          pch = c(19,17), cex = txtSz, bty = 'n' )
  legend( 0, .5, c('Targets','Competitors'),
          fill = c('black','white'), cex = txtSz, bty = 'n' )
  
  mtext( 'Accuracy by latency plot', side = 3, outer = T, 
         cex = txtSz, line = -2 )
  
}

###
### Latency by accuracy by condition/subject
###
# Lookup - 03

if ( runCode[2] ) {
  
  if (!savePlot) x11(width=12)
  layout( cbind( 1, 2, 3, 4 ) )
  
  # Define a covariate collapsing over correct position
  
  # Determine plotting dimensions
  MRT = aggregate( d$RT, list( d$Cnd, d$Ac, d$S ), mean )
  yl = lowerUpper( .2, MRT$x )
  # Determine color scheme for separating subjects
  MPC = aggregate( d$Ac, list( d$S ), mean )
  ord = order( MPC$x )
  ttl = c('T-SR',
          'T-B',
          'C-SR',
          'C-B' )
  
  # Determine choice probabilities and mean RT per subject
  # for a specific condition
  for (i in 1:4) {
    cnd = d$Cnd == i
    x = aggregate( d$Ac[cnd], list( d$S[cnd] ), p_f )
    colnames( x ) = c( 'S','P' )
    x = as.data.frame(
      cbind( S = c( x$S[ x$P[,1] > 0 ], x$S[ x$P[,2] > 0 ] ),
             A = c( rep(0,sum(x$P[,1] > 0)), rep(1,sum(x$P[,2] > 0)) ),
             P = c( x$P[ x$P[,1] > 0, 1 ], x$P[ x$P[,2] > 0, 2 ] ) )
    )
    y = aggregate( d$RT[cnd], list( d$S[cnd], d$Ac[cnd] ), mean )
    colnames( y ) = c('S','A','MRT')
    
    # Define current plotting dimensions
    yd = c( -diff(yl), diff(yl) )
    par( mar = c( 5, 5, 1, .5 ) )
    blankPlot( yDim = yd )
    # Add meaningful labels
    segments( .5, -1.5*diff(yl), .5, diff(yl)*.95, col = 'grey80', lwd = 2 )
    abline( h = 0,lwd = 2 )
    axis( 1, seq(0,1,.25), cex.axis = 1.5, lwd = 2 )
    lbl = c( seq( yl[2], yl[1], -.5 )[-4], seq( yl[1], yl[2], .5 ) )
    pos = seq( -diff(yl), diff(yl), .5 )
    axis( 2, pos, lbl, cex.axis = 1.5, lwd = 2 )
    
    # Plot points (darker points < overall accuracy)
    clr = ord[ x$S[ x$A == 1 ] ]
    points( x$P[ x$A == 1 ], y$MRT[ y$A == 1 ]-yl[1], pch = 19,
            col = rgb( clr/N, 0, 0, 1 ), cex = 1.5 )
    clr = ord[ x$S[ x$A == 0 ] ]
    points( x$P[ x$A == 0 ], -(y$MRT[ y$A == 0 ]-yl[1]), pch = 17,
            col = rgb( clr/N, 0, 0, 1 ), cex = 1.5 )
    
    # Add legends
    if ( i == 4 ) legend( 'bottomright', c('Correct','Error'),
                          pch = c(19,17), cex = 1.5, bty = 'n' )
    legend( 'topright', ttl[i], cex = 1.5, bty = 'n' )
    
  }
  mtext('MRT (s)', side = 2, outer = T, line = -2.25, cex = 1.5 )
  mtext('Accuracy', side = 1, outer = T, line = -2.25, 
        cex = 1.5 )
  
}

###
###
###
#

if ( runCode[3] ) {
  
  if (!savePlot) x11(width=12)
  lyt = matrix( 1, 10, 20 )
  lyt[,11:20] = 2
  lyt[10,]=3
  layout( lyt )
  
  txtSz = 2
  lnSz = 2.5
  ptSz = 2
  
  ### Left ###
  
  yl = c(1.1,1.5)
  
  par( mar=c(5,6,2,1) )
  blankPlot(xDim=c(0,1),yDim=yl)
  axis( 1, seq(0,1,.25), cex.axis = txtSz, lwd = lnSz )
  axis( 2, seq(yl[1],yl[2],.2), cex.axis = txtSz, lwd = lnSz )
  mtext( 'Accuracy', side = 1, cex = txtSz, line = 3.25 )
  mtext( 'MRT (s)', side = 2, cex = txtSz, line = 2.75 )
  
  # Calculate PvT values
  sel = d$Co == 0
  rt = d$RT[ sel ]; ac = d$Ac[ sel ]; sb = d$S[ sel ]; cv = d$Cnd[ sel ]
  outF = pvt_points( rt, ac, cv, grp = sb, 
                     opt = list( out = T, draw = F ) )
  segments( outF$pv$x[1:4], outF$pv$y[1:4],
            outF$pv$x[5:8], outF$pv$y[5:8],
            col = 'grey', lty = 1:2, lwd = lnSz )
  
  # Add in PvT values
  plt = list( pch = c(19,17,21,24), bg = 'white' )
  pvt_points( rt, ac, cv, grp = sb, plt = plt, cex = ptSz )
  
  ### Right ###
  
  yl = c(1.1,1.5)
  
  par( mar=c(5,6,2,1) )
  blankPlot(xDim=c(0,1),yDim=yl)
  axis( 1, seq(0,1,.25), cex.axis = txtSz, lwd = lnSz )
  axis( 2, seq(yl[1],yl[2],.2), cex.axis = txtSz, lwd = lnSz )
  mtext( 'Accuracy', side = 1, cex = txtSz, line = 3.25 )
  mtext( 'MRT (s)', side = 2, cex = txtSz, line = 2.75 )
  
  # Calculate PvT values
  sel = d$Co == 1
  rt = d$RT[ sel ]; ac = d$Ac[ sel ]; sb = d$S[ sel ]; cv = d$Cnd[ sel ]
  outF = pvt_points( rt, ac, cv, grp = sb, 
                     opt = list( out = T, draw = F ) )
  segments( outF$pv$x[1:4], outF$pv$y[1:4],
            outF$pv$x[5:8], outF$pv$y[5:8],
            col = 'grey', lty = 1:2, lwd = lnSz )
  
  # Add in PvT values
  plt = list( pch = c(19,17,21,24), bg = 'white' )
  pvt_points( rt, ac, cv, grp = sb, plt = plt, cex = ptSz )
  
  # Legends
  par( mar=c(0,0,0,0) )
  blankPlot()
  
  legend( 'bottomleft', c('Selective retrieval', 'Baseline'),
          pch = c(19,17), cex = txtSz, bty = 'n', horiz = T )
  legend( 'bottomright', c('Targets','Competitors'),
          fill = c('black','white'), cex = txtSz, bty = 'n',
          horiz = T )
  
  mtext( 'Accuracy by latency plot', side = 3, outer = T, 
         cex = txtSz, line = -2 )
  
}

if ( savePlot ) dev.off()