#---------------------#
# Descriptive results #
# Kevin Potter        #
# Updated 11/23/2016  #
#---------------------#

# Clear workspace
rm( list = ls() )

# Save current directory
orig_dir = getwd()

# Index
# Lookup - 01:  Load in useful packages and data
# Lookup - 02:  Priming duration x type x accuracy x task
#                 (runCode[1])
# Lookup - 03:  Priming duration x type effect x subject (2AFC)
#                 (runCode[2])
# Lookup - 04:  Priming duration x type effect x subject (SD)
#                 (runCode[3])
# Lookup - 05:  Priming duration x type x correct side effect (2AFC)
#                 (runCode[4])
# Lookup - 06:  Priming duration x type x correct side effect (SD)
#                 (runCode[5])
# Lookup - 07:  Examination of response biases
#                 (runCode[6])
# Lookup - 08:  Plot nROUSE model predictions against observed data
#                 (runCode[7])

# Indicate which code segments to run
runCode = c( T, T, T, T, T, T, T )

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
load( 'SD_v_FC.RData' )
setwd( orig_dir )

# For easy manipulation
d = allDat
colnames( d ) = c('S','Cnd','Ac','RT','E','PD',
                  'Ta','Co','PT','Ch')

###
### Priming duration x type x accuracy x task
###
# Lookup - 02

if ( runCode[1] ) {
  
  if (!savePlot) x11(width=12)
  lyt = matrix( 1, 10, 20 )
  lyt[,11:20] = 2
  lyt[10,]=3
  layout( lyt )
  
  txtSz = 2
  lnSz = 2.5
  ptSz = 2
  
  ### 2AFC ###
  
  yl = c(.4,1)
  
  par( mar=c(5,6,2,1) )
  blankPlot(xDim=c(0,1),yDim=yl)
  axis( 1, seq(0,1,.25), cex.axis = txtSz, lwd = lnSz )
  axis( 2, seq(yl[1],yl[2],.2), cex.axis = txtSz, lwd = lnSz )
  mtext( 'P(Y=y)', side = 1, cex = txtSz, line = 3.25 )
  mtext( 'MRT (s)', side = 2, cex = txtSz, line = 2.75 )
  legend( 'topleft', '2AFC', bty = 'n', cex = txtSz )
  
  # Define a covariate that collapses over position correct
  cvrt = rep( 1, nrow(d) );
  cvrt[ d$Cnd > 8 ] = NA; # Exclude SD task
  cvrt[ d$Cnd == 7 | d$Cnd == 8 ] = 2
  cvrt[ d$Cnd == 1 | d$Cnd == 2 ] = 3
  cvrt[ d$Cnd == 3 | d$Cnd == 4 ] = 4
  
  # Calculate PvT values for foil primed condition (2AFC)
  cnd = d$Cnd >= 5 & d$Cnd <= 8
  rt = d$RT[ cnd ]; ac = d$Ac[ cnd ]; sb = d$S[ cnd ]; cv = cvrt[ cnd ]
  outF = pvt_points( rt, ac, cv, grp = sb, 
                     opt = list( out = T, draw = F ) )
  segments( outF$pv$x[1:2], outF$pv$y[1:2],
            outF$pv$x[3:4], outF$pv$y[3:4],
            col = 'grey', lty = 1:2, lwd = lnSz )
  
  # Calculate PvT values for target primed condition
  cnd = d$Cnd >= 1 & d$Cnd <= 4
  rt = d$RT[ cnd ]; ac = d$Ac[ cnd ]; sb = d$S[ cnd ]; 
  cv = cvrt[ cnd ] - 2
  outT = pvt_points( rt, ac, cv, grp = sb, opt = list( out = T, draw = F ) )
  segments( outT$pv$x[1:2], outT$pv$y[1:2],
            outT$pv$x[3:4], outT$pv$y[3:4],
            col = 'grey', lty = 1:2, lwd = lnSz )
  
  # Add in PvT values
  cnd = d$Cnd >= 5 & d$Cnd <= 8
  rt = d$RT[ cnd ]; ac = d$Ac[ cnd ]; sb = d$S[ cnd ]; cv = cvrt[ cnd ]
  plt = list( pch = c(19,19,21,21), bg = 'white' )
  pvt_points( rt, ac, cv, grp = sb, plt = plt, cex = ptSz )
  
  cnd = d$Cnd >= 1 & d$Cnd <= 4
  rt = d$RT[ cnd ]; ac = d$Ac[ cnd ]; sb = d$S[ cnd ]; 
  cv = cvrt[ cnd ] - 2
  plt = list( pch = c(17,17,24,24), bg = 'white' )
  pvt_points( rt, ac, cv, grp = sb, plt = plt, cex = ptSz )
  
  ### SD ###
  
  yl = c(.4,1)
  
  par( mar=c(5,6,2,1) )
  blankPlot(xDim=c(0,1),yDim=yl)
  axis( 1, seq(0,1,.25), cex.axis = txtSz, lwd = lnSz )
  axis( 2, seq(yl[1],yl[2],.2), cex.axis = txtSz, lwd = lnSz )
  mtext( 'P(Y=y)', side = 1, cex = txtSz, line = 3.25 )
  mtext( 'MRT (s)', side = 2, cex = txtSz, line = 2.75 )
  legend( 'topleft', 'SD', bty = 'n', cex = txtSz )
  
  # Define a covariate that collapses over position correct
  cvrt = rep( 1, nrow(d) );
  cvrt[ d$Cnd < 9 ] = NA; # Exclude SD task
  cvrt[ d$Cnd == 11 | d$Cnd == 12 ] = 2
  cvrt[ d$Cnd == 13 | d$Cnd == 14 ] = 3
  cvrt[ d$Cnd == 15 | d$Cnd == 16 ] = 4
  
  # Calculate PvT values for foil primed condition (2AFC)
  cnd = d$Cnd >= 13 & d$Cnd <= 16
  rt = d$RT[ cnd ]; ac = d$Ac[ cnd ]; sb = d$S[ cnd ]; cv = cvrt[ cnd ]
  outF = pvt_points( rt, ac, cv, grp = sb, 
                     opt = list( out = T, draw = F ) )
  segments( outF$pv$x[1:2], outF$pv$y[1:2],
            outF$pv$x[3:4], outF$pv$y[3:4],
            col = 'grey', lty = 1:2, lwd = lnSz )
  
  # Calculate PvT values for target primed condition
  cnd = d$Cnd >= 9 & d$Cnd <= 12
  rt = d$RT[ cnd ]; ac = d$Ac[ cnd ]; sb = d$S[ cnd ]; 
  cv = cvrt[ cnd ] - 2
  outT = pvt_points( rt, ac, cv, grp = sb, opt = list( out = T, draw = F ) )
  segments( outT$pv$x[1:2], outT$pv$y[1:2],
            outT$pv$x[3:4], outT$pv$y[3:4],
            col = 'grey', lty = 1:2, lwd = lnSz )
  
  # Add in PvT values
  cnd = d$Cnd >= 13 & d$Cnd <= 16
  rt = d$RT[ cnd ]; ac = d$Ac[ cnd ]; sb = d$S[ cnd ]; cv = cvrt[ cnd ]
  plt = list( pch = c(19,19,21,21), bg = 'white' )
  pvt_points( rt, ac, cv, grp = sb, plt = plt, cex = ptSz )
  
  cnd = d$Cnd >= 9 & d$Cnd <= 12
  rt = d$RT[ cnd ]; ac = d$Ac[ cnd ]; sb = d$S[ cnd ]; 
  cv = cvrt[ cnd ] - 2
  plt = list( pch = c(17,17,24,24), bg = 'white' )
  pvt_points( rt, ac, cv, grp = sb, plt = plt, cex = ptSz )
  
  # Legends
  par( mar=c(0,0,0,0) )
  blankPlot()
  
  legend( 'bottomleft', c('Target primed','Foil primed'),
          pch = c(17,19), cex = txtSz, bty = 'n', horiz = T )
  legend( 'bottom', c('.05 s prime','.4 s prime'),
          lty = 1:2, col = 'grey', lwd = lnSz, cex = txtSz, bty = 'n',
          horiz = T )
  legend( 'bottomright', c('Error','Correct'),
          pch = c(19,21), pt.bg = 'white', cex = txtSz, bty = 'n',
          horiz = T )
  
  mtext( 'Prime duration x type x accuracy', side = 3, outer = T, 
         cex = txtSz, line = -2 )
  
}

###
### Priming duration x type effect x subject (2AFC)
###
# Lookup - 03

if ( runCode[2] ) {
  
  if (!savePlot) x11(width=12)
  layout( cbind( 1, 2, 3, 4 ) )
  
  # Define a covariate collapsing over correct position
  cvrt = rep( 1, nrow(d) )
  tmp = (1:8)*2
  for (i in 2:8) {
    cvrt[ d$Cnd > tmp[i-1] & d$Cnd <= tmp[i] ] = i
  }
  
  # Determine plotting dimensions
  MRT = aggregate( d$RT, list( cvrt, d$Ac, d$S ), mean )
  # yl = lowerUpper( .5, MRT$x )
  yl = c(0.0,1.5)
  # Determine color scheme for separating subjects
  MPC = aggregate( d$Ac, list( d$S ), mean )
  ord = order( MPC$x )
  ttl = c('Target primed .05 s',
          'Target primed 2 s',
          'Foil primed .05 s',
          'Foil primed 2 s' )
  
  # Determine choice probabilities and mean RT per subject
  # for a specific condition
  for (i in 1:4) {
    cnd = cvrt == i
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
  mtext('P( Y = y ) for 2AFC', side = 1, outer = T, line = -2.25, 
        cex = 1.5 )

}

###
### Priming duration x type effect x subject (SD)
###
# Lookup - 04

if ( runCode[3] ) {
  
  if (!savePlot) x11(width=12)
  layout( cbind( 1, 2, 3, 4 ) )
  
  # Define a covariate collapsing over correct position
  cvrt = rep( 1, nrow(d) )
  tmp = (1:8)*2
  for (i in 2:8) {
    cvrt[ d$Cnd > tmp[i-1] & d$Cnd <= tmp[i] ] = i
  }
  
  # Determine plotting dimensions
  MRT = aggregate( d$RT, list( cvrt, d$Ac, d$S ), mean )
  # yl = lowerUpper( .5, MRT$x )
  yl = c(0.0,1.5)
  # Determine color scheme for separating subjects
  MPC = aggregate( d$Ac, list( d$S ), mean )
  ord = order( MPC$x )
  ttl = c('Target primed .05 s',
          'Target primed 2 s',
          'Foil primed .05 s',
          'Foil primed 2 s' )
  
  # Determine choice probabilities and mean RT per subject
  # for a specific condition
  for (i in 1:4) {
    cnd = cvrt == i+4
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
  mtext('P( Y = y ) for SD', side = 1, outer = T, line = -2.25, 
        cex = 1.5 )
  
}


###
### Priming duration x type x correct side effect (2AFC)
###
# Lookup - 05

if ( runCode[4] ) {
  
  if (!savePlot) x11(width=12)
  lyt = matrix( 1, 10, 20 )
  lyt[,11:20] = 2
  lyt[10,]=3
  layout( lyt )
  
  txtSz = 2;
  lnSz = 2.5
  
  for ( sd in 0:1) {
    
    par( mar = c( 5, 6, 3, 1 ) )
    blankPlot(xDim=c(0,1),yDim=c(.4,1))
    axis( 1, seq(0,1,.25), cex.axis = txtSz, lwd = 2 )
    axis( 2, seq(.4,1,.2), cex.axis = txtSz, lwd = 2 )
    mtext( 'P(Y=y)', side = 1, cex= txtSz, line = 3.5 )
    mtext( 'MRT (s)', side = 2, cex= txtSz, line = 2.75 )
    if (sd==0) mtext('Left correct',side=3,cex=txtSz,line = -1.25)
    if (sd==1) mtext('Right correct',side=3,cex=txtSz,line = -1.25)
    
    # Calculate PvT values for foil primed condition (left)
    cnd = d$Cnd >= 5 & d$Cnd <= 8 & d$Co == sd
    rt = d$RT[ cnd ]; ac = d$Ac[ cnd ]; sb = d$S[ cnd ]; 
    cv = d$Cnd[ cnd ] - 4
    outF = pvt_points( rt, ac, cv, grp = sb, 
                       opt = list( out = T, draw = F ) )
    segments( outF$pv$x[1:2], outF$pv$y[1:2],
              outF$pv$x[3:4], outF$pv$y[3:4],
              col = 'grey', lty = 1:2, lwd = 2 )
    
    # Calculate PvT values for target primed condition
    cnd = d$Cnd >= 1 & d$Cnd <= 4 & d$Co == sd
    rt = d$RT[ cnd ]; ac = d$Ac[ cnd ]; sb = d$S[ cnd ]; 
    cv = d$Cnd[ cnd ]
    outT = pvt_points( rt, ac, cv, grp = sb, 
                       opt = list( out = T, draw = F ) )
    segments( outT$pv$x[1:2], outT$pv$y[1:2],
              outT$pv$x[3:4], outT$pv$y[3:4],
              col = 'grey', lty = 1:2, lwd = 2 )
    
    # Add in PvT values
    cnd = d$Cnd >= 5 & d$Cnd <= 8 & d$Co == sd
    rt = d$RT[ cnd ]; ac = d$Ac[ cnd ]; sb = d$S[ cnd ]; 
    cv = d$Cnd[ cnd ] - 4
    plt = list( pch = c(19,19,21,21), bg = 'white' )
    pvt_points( rt, ac, cv, grp = sb, plt = plt, cex= txtSz )
    
    cnd = d$Cnd >= 1 & d$Cnd <= 4 & d$Co == sd
    rt = d$RT[ cnd ]; ac = d$Ac[ cnd ]; sb = d$S[ cnd ]; 
    cv = d$Cnd[ cnd ]
    plt = list( pch = c(17,17,24,24), bg = 'white' )
    pvt_points( rt, ac, cv, grp = sb, plt = plt, cex= txtSz )
    
  }
  
  # Legends
  par( mar = c(0,0,0,0) )
  blankPlot()
  
  legend( 'bottomleft', c('Target primed','Foil primed'),
          pch = c(17,19), cex = txtSz*.95, bty = 'n', horiz = T )
  legend( 'bottom', c('.05 s prime','.4 s prime'), horiz = T, 
          lty = 1:2, col = 'grey', lwd = 2, cex = txtSz*.95, bty = 'n' )
  legend( 'bottomright', c('Error','Correct'), horiz = T, 
          pch = c(19,21), pt.bg = 'white', cex = txtSz*.95, bty = 'n' )
  
  mtext( 'Prime duration x type x side x accuracy', side = 3, outer = T, 
         cex= txtSz, line = -2 )
  
}

###
### Priming duration x type x correct side effect (SD)
###
# Lookup - 06

if ( runCode[5] ) {
  
  if (!savePlot) x11(width=12)
  lyt = matrix( 1, 10, 20 )
  lyt[,11:20] = 2
  lyt[10,]=3
  layout( lyt )
  
  txtSz = 2;
  lnSz = 2.5
  
  for ( sd in 0:1) {
    
    par( mar = c( 5, 6, 3, 1 ) )
    blankPlot(xDim=c(0,1),yDim=c(.4,1))
    axis( 1, seq(0,1,.25), cex.axis = txtSz, lwd = 2 )
    axis( 2, seq(.4,1,.2), cex.axis = txtSz, lwd = 2 )
    mtext( 'P(Y=y)', side = 1, cex= txtSz, line = 3.5 )
    mtext( 'MRT (s)', side = 2, cex= txtSz, line = 2.75 )
    if (sd==0) mtext('Same correct',side=3,cex=txtSz,line = -1.25)
    if (sd==1) mtext('Different correct',side=3,cex=txtSz,line = -1.25)
    
    # Calculate PvT values for foil primed condition (left)
    cnd = d$Cnd >= 13 & d$Cnd <= 16 & d$Co == sd
    rt = d$RT[ cnd ]; ac = d$Ac[ cnd ]; sb = d$S[ cnd ]; 
    cv = d$Cnd[ cnd ] - 12
    outF = pvt_points( rt, ac, cv, grp = sb, 
                       opt = list( out = T, draw = F ) )
    segments( outF$pv$x[1:2], outF$pv$y[1:2],
              outF$pv$x[3:4], outF$pv$y[3:4],
              col = 'grey', lty = 1:2, lwd = 2 )
    
    # Calculate PvT values for target primed condition
    cnd = d$Cnd >= 9 & d$Cnd <= 12 & d$Co == sd
    rt = d$RT[ cnd ]; ac = d$Ac[ cnd ]; sb = d$S[ cnd ]; 
    cv = d$Cnd[ cnd ] - 8
    outT = pvt_points( rt, ac, cv, grp = sb, 
                       opt = list( out = T, draw = F ) )
    segments( outT$pv$x[1:2], outT$pv$y[1:2],
              outT$pv$x[3:4], outT$pv$y[3:4],
              col = 'grey', lty = 1:2, lwd = 2 )
    
    # Add in PvT values
    cnd = d$Cnd >= 13 & d$Cnd <= 16 & d$Co == sd
    rt = d$RT[ cnd ]; ac = d$Ac[ cnd ]; sb = d$S[ cnd ]; 
    cv = d$Cnd[ cnd ] - 12
    plt = list( pch = c(19,19,21,21), bg = 'white' )
    pvt_points( rt, ac, cv, grp = sb, plt = plt, cex= txtSz )
    
    cnd = d$Cnd >= 9 & d$Cnd <= 12 & d$Co == sd
    rt = d$RT[ cnd ]; ac = d$Ac[ cnd ]; sb = d$S[ cnd ]; 
    cv = d$Cnd[ cnd ] - 8
    plt = list( pch = c(17,17,24,24), bg = 'white' )
    pvt_points( rt, ac, cv, grp = sb, plt = plt, cex= txtSz )
    
  }
  
  # Legends
  par( mar = c(0,0,0,0) )
  blankPlot()
  
  legend( 'bottomleft', c('Target primed','Foil primed'),
          pch = c(17,19), cex = txtSz*.95, bty = 'n', horiz = T )
  legend( 'bottom', c('.05 s prime','.4 s prime'), horiz = T, 
          lty = 1:2, col = 'grey', lwd = 2, cex = txtSz*.95, bty = 'n' )
  legend( 'bottomright', c('Error','Correct'), horiz = T, 
          pch = c(19,21), pt.bg = 'white', cex = txtSz*.95, bty = 'n' )
  
  mtext( 'Prime duration x type x side x accuracy', side = 3, outer = T, 
         cex= txtSz, line = -2 )
  
}

###
### Examination of response biases
###
# Lookup - 07

if ( runCode[6] ) {
  
  # Calculate accuracy per subject and condition
  ac = aggregate( d$Ac, 
                  list( d$Cnd, d$S, d$PD, d$PT, d$Ta, d$Co ), p_f )
  colnames( ac ) = c('Cnd','S','PD','PT','Ta','Co','A')
  
  # Reorganize data frame
  ac = cbind(
    Cnd = ac$Cnd[ ac$Co == 0 ],
    S = ac$S[ ac$Co == 0 ],
    PD = ac$PD[ ac$Co == 0 ],
    PT = ac$PT[ ac$Co == 0 ],
    Ta = ac$Ta[ ac$Co == 0 ],
    E0 = ac$A[ ac$Co == 0, 1],
    E1 = ac$A[ ac$Co == 1, 1],
    C0 = ac$A[ ac$Co == 0, 2],
    C1 = ac$A[ ac$Co == 1, 2] )
  ac = as.data.frame( ac )
  
  # Determine rows with missing data
  tmp = colSums( apply( ac[,6:9], 1, function(x) x == 0 ) )
  miss = which( tmp > 0 )
  
  # Calculate mean RT per subject and condition
  mrt = aggregate( d$RT, 
                   list( d$Cnd, d$S, d$PD, d$PT, d$Ta, d$Co, d$Ac ), 
                   mean )
  colnames( mrt ) = c('Cnd','S','PD','PT','Ta','Co','A','MRT')
  
  # Color code points based on overall accuracy
  MPC = aggregate( d$Ac, list( d$S ), mean )
  ord = order( MPC$x )
  allOrd=ord[ mrt$S ]
  
  # Trim rows with missing data
  all_rmv = c()
  
  chk = ac[ miss, 2:5 ]
  
  for ( i in 1:nrow( chk ) ) {
    
    comp = as.numeric( chk[i,] )
    
    res = apply( mrt[,2:5], 1, function(x) sum( comp == x ) )
    rmv = which( res == 4 )
    all_rmv = c( all_rmv, rmv )
  }

  mrt = mrt[ -unique(all_rmv), ]
  ord = allOrd[ -unique(all_rmv) ]
  
  # Reorganize data frame
  mrt = cbind(
    Cnd = mrt$Cnd[ mrt$Co == 0 & mrt$A == 0 ],
    S = mrt$S[ mrt$Co == 0 & mrt$A == 0 ],
    PD = mrt$PD[ mrt$Co == 0 & mrt$A == 0 ],
    PT = mrt$PT[ mrt$Co == 0 & mrt$A == 0 ],
    Ta = mrt$Ta[ mrt$Co == 0 & mrt$A == 0 ],
    E0 = mrt$MRT[ mrt$Co == 0 & mrt$A == 0 ],
    E1 = mrt$MRT[ mrt$Co == 1 & mrt$A == 0 ],
    C0 = mrt$MRT[ mrt$Co == 0 & mrt$A == 1 ],
    C1 = mrt$MRT[ mrt$Co == 1 & mrt$A == 1 ] )
  mrt = as.data.frame( mrt )
  
  # Trim row from data frame for accuracy as well
  ac = ac[ -miss, ]
  
  for (ta in 0:1) {
    
    # Plot difference scores for each condition
    if (!savePlot) x11(width=12)
    layout( rbind( 1:4, 5:8 ) )
    
    ttl = c( 'Target primed .05 s',
             'Target primed .4 s',
             'Foil primed .05 s',
             'Foil primed .4 s' )
    
    # Determine color scheme
    clr = ord[1:nrow(ac)]/length( unique( ord ) )
    
    # Define a covariate for conditions
    cvrt=createIncrement(ac$Cnd)
    
    if (ta == 0 ) ind = c(1:4,1:4)
    if (ta == 1) ind = c(5:8,5:8)
    for ( i in 1:8 ) {
      
      cnd = cvrt == ind[i]
      
      yl = c(-0.8,.8)
      xl = c(-.5,.5)
      par( mar = c(5,6,3,.5) )
      blankPlot( yDim = yl, xDim = xl )
      title( ttl[i], cex = 1.5 )
      abline(h=0,lwd=2,col='grey')
      abline(v=0,lwd=2,col='grey')
      axis( 2, round( seq( yl[1], yl[2], .4 ), 2), cex.axis = 1.5,
            lwd = 2 )
      axis( 1, round( seq( xl[1], xl[2], .25 ), 2), cex.axis = 1.5,
            lwd = 2 )
      if (i==1) mtext('Right - left (Error RT in s)',side=2,
                      cex = 1.25, line = 3 )
      if (i==5) mtext('Right - left (Correct RT in s)',side=2,
                      cex = 1.25, line = 3 )
      
      if (i < 5) {
        dft = mrt$E1[cnd] - mrt$E0[cnd]
        dfa = ac$E1[cnd] - ac$E0[cnd]
      } else {
        dft = mrt$C1[cnd] - mrt$C0[cnd]
        dfa = ac$C1[cnd] - ac$C0[cnd]
      }
      
      points( dfa, dft, pch = 19, col = rgb( clr, 0, 0, 1 ) )
      
    }
    if (ta==0) string = 'Right - left (accuracy; 2AFC)' else 
      string = 'Right - left (accuracy; SD)'
    mtext(string,side=1,outer=T,cex=1.5,
          line = -2 )
    
  }
  
}

###
### Plot nROUSE model predictions against observed data
###
# Lookup - 08

if ( runCode[7] ) {
  
  # Calculate P(Correct) by condition
  sel = d$Ta == 0
  y = aggregate( d$Ac[sel], list( d$PD[sel], d$PT[sel] ), mean )
  tot = aggregate( rep(1,nrow(d[sel,])), 
                   list(d$PD[sel],d$PT[sel]), sum )$x
  colnames( y ) = c('PD','PT','P')
  x = log( c( .05, .4 ) )
  # Uncertainty intervals
  tmp = aggregate( d$Ac[sel], list( d$PD[sel], d$PT[sel], d$S[sel] ), 
                   mean )
  colnames( tmp ) = c('PD','PT','S','P')
  ui = aggregate( tmp$P, list( tmp$PD, tmp$PT ), 
                  function(x) {
                    ct = qt( .975, N - 1 );
                    return( c( mean(x) - sem(x)*ct,
                               mean(x) + sem(x)*ct ) ) } )
  colnames( ui ) = c('PD','PT','UI')
  
  # Plotting dimensions
  if (!savePlot) x11(width=12);
  layout( cbind(1,2) )
  par( mar = c(4, 5, 1, 1 ) )
  blankPlot( xDim = c( -5, 1 ) )
  axis( 2, seq(0,1,.25), cex.axis = 1.5, lwd = 2 )
  axis( 1, log( c( .05, .4 ) ),
        c( 50, 400 ), cex.axis = 1.5, lwd = 2 )
  mtext( 'Prime duration (ms)', side = 1, cex = 1.5, 
         line = 2.5 )
  mtext( 'P(Correct)', side = 2, cex = 1.5, 
         line = 2.5 )
  
  # Add observed data
  segments( x, ui$UI[ ui$PT == 0,1 ],
            x, ui$UI[ ui$PT == 0,2 ], col = 'orange', lwd = 2 )
  segments( x, ui$UI[ ui$PT == 1,1 ],
            x, ui$UI[ ui$PT == 1,2 ], col = 'red', lwd = 2 )
  lines( x, y$P[ y$PT == 0 ], col = 'orange', lwd = 2 )
  points( x, y$P[ y$PT == 0 ], col = 'orange', pch = 19, cex = 1.5 )
  lines( x, y$P[ y$PT == 1 ], col = 'red', lwd = 2 )
  points( x, y$P[ y$PT == 1 ], col = 'red', pch = 19, cex = 1.5 )
  
  abline( h = .5, lty = 2, lwd = 2 )
  
  # Define nROUSE parameters
  np = c( .0302, 0.9844, 1 )
  
  # Create data set to obtain nROUSE model estimates
  nDat = cbind( TarDur = 50, MaskDur = 450, PrimeDur = c(50,400,50,400),
                Type = c(2,2,-2,-2), Y = y$P*tot, N = tot )
  pred = nROUSE_logLik( log(np), nDat, estimate = F )
  
  # Add in model predictions
  lines( x, pred[ y$PT == 0 ], lwd = 2, lty = 2, col = 'orange' )
  points( x, pred[ y$PT == 0 ], pch = 21, bg = 'white', col = 'orange',
          cex = 1.5, lwd = 2 )
  lines( x, pred[ y$PT == 1 ], lwd = 2, lty = 2, col = 'red' )
  points( x, pred[ y$PT == 1 ], pch = 21, bg = 'white', col = 'red',
          cex = 1.5, lwd = 2 )
  
  # Add in legend
  legend( 'bottomleft', c('Target','Foil'),
          fill = c('orange','red'), cex = 1.5, bty = 'n' )
  legend( 'bottomright', c('Observed','nROUSE'), 
          pch = c(19,21),pt.bg = 'white',
          pt.cex = 1.5, lty = c(1,2), lwd = 2, 
          cex = 1.5, bty = 'n' )
  
  blankPlot()
}

# Return to original directory
if ( savePlot ) { setwd( 'Data' ); dev.off() } # Close plotting window
setwd( orig_dir )