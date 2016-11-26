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
# Lookup - 02:  Priming duration x type effect 
#                 (runCode[1])
# Lookup - 03:  Priming duration x type effect x subject
#                 (runCode[2])
# Lookup - 04:  Priming duration x type x correct side effect
#                 (runCode[3])
# Lookup - 05:  Examination of biases for responding left versus right
#                 (runCode[4])
# Lookup - 06:  Response time and choice distributions by offset
#                 (runCode[5])
# Lookup - 07:  Fastest RT by offset and condition
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
load( 'Priming_offset.RData' )
setwd( orig_dir )

# For easy manipulation
d = allDat
colnames( d ) = c('S','PD','PT','O','Co','RT',
                  'Ch','Ac','OT','DT','CDT','Cnd')

###
### Priming duration x type effect
###
# Lookup - 02

if ( runCode[1] ) {
  
  if (!savePlot) x11(width=12)
  layout( cbind(1,2) )
  
  blankPlot(xDim=c(0,1),yDim=c(.6,1.2))
  axis( 1, seq(0,1,.25), cex.axis = 1.5, lwd = 2 )
  axis( 2, seq(.6,1.2,.2), cex.axis = 1.5, lwd = 2 )
  mtext( 'P(Y=y)', side = 1, cex = 1.5, line = 2.5 )
  mtext( 'MRT (s)', side = 2, cex = 1.5, line = 2.5 )
  
  # Calculate PvT values for foil primed condition
  cnd = d$DT < 3
  rt = d$RT[ cnd ]; ac = d$Ac[ cnd ]; sb = d$S[ cnd ]; cv = d$DT[ cnd ]
  outF = pvt_points( rt, ac, cv, grp = sb, opt = list( out = T, draw = F ) )
  segments( outF$pv$x[1:2], outF$pv$y[1:2],
            outF$pv$x[3:4], outF$pv$y[3:4],
            col = 'grey', lty = 1:2, lwd = 2 )
  
  # Calculate PvT values for target primed condition
  cnd = d$DT > 2
  rt = d$RT[ cnd ]; ac = d$Ac[ cnd ]; sb = d$S[ cnd ]; cv = d$DT[ cnd ] - 2
  outT = pvt_points( rt, ac, cv, grp = sb, opt = list( out = T, draw = F ) )
  segments( outT$pv$x[1:2], outT$pv$y[1:2],
            outT$pv$x[3:4], outT$pv$y[3:4],
            col = 'grey', lty = 1:2, lwd = 2 )
  
  # Add in PvT values
  cnd = d$DT < 3
  rt = d$RT[ cnd ]; ac = d$Ac[ cnd ]; sb = d$S[ cnd ]; cv = d$DT[ cnd ]
  plt = list( pch = c(19,19,21,21), bg = 'white' )
  pvt_points( rt, ac, cv, grp = sb, plt = plt, cex = 1.5 )
  
  cnd = d$DT > 2
  rt = d$RT[ cnd ]; ac = d$Ac[ cnd ]; sb = d$S[ cnd ]; cv = d$DT[ cnd ] - 2
  plt = list( pch = c(17,17,24,24), bg = 'white' )
  pvt_points( rt, ac, cv, grp = sb, plt = plt, cex = 1.5 )
  
  # Legends
  blankPlot()
  
  legend( 0,.75, c('Target primed','Foil primed'),
          pch = c(19,17), cex = 1.5, bty = 'n' )
  legend( 0,.5, c('.05 s prime','2 s prime'),
          lty = 1:2, col = 'grey', lwd = 2, cex = 1.5, bty = 'n' )
  legend( 0,.25, c('Error','Correct'),
          pch = c(19,21), pt.bg = 'white', cex = 1.5, bty = 'n' )
  
  mtext( 'Prime duration x type x accuracy', side = 3, outer = T, 
         cex = 1.5, line = -2 )
  
}

###
### Priming duration x type effect x subject
###
# Lookup - 03

if ( runCode[2] ) {
  
  if (!savePlot) x11(width=12)
  layout( cbind( 1, 2, 3, 4 ) )
  
  # Determine plotting dimensions
  MRT = aggregate( d$RT, list( d$DT, d$Ac, d$S ), mean )
  # yl = lowerUpper( .5, MRT$x )
  yl = c(.4,3.4)
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
    cnd = d$DT == i
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
    lbl = c( seq( yl[2], yl[1], -1 )[-1], seq( yl[1], yl[2], 1 ) )
    pos = seq( -diff(yl), diff(yl), 1 )
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
  mtext('P( Y = y )', side = 1, outer = T, line = -2.25, cex = 1.5 )

}

###
### Priming duration x type x correct side effect
###
# Lookup - 04

if ( runCode[3] ) {
  
  if (!savePlot) x11(width=12)
  layout( cbind(1,1,1,2,2,2,3) )
  
  for ( sd in 0:1) {
    
    par( mar = c( 4, 5, 3, 1 ) )
    blankPlot(xDim=c(0,1),yDim=c(.6,1.2))
    axis( 1, seq(0,1,.25), cex.axis = 1.5, lwd = 2 )
    axis( 2, seq(.6,1.2,.2), cex.axis = 1.5, lwd = 2 )
    mtext( 'P(Y=y)', side = 1, cex = 1.5, line = 2.5 )
    mtext( 'MRT (s)', side = 2, cex = 1.5, line = 2.5 )
    if (sd==0) mtext('Left correct',side=3,cex=1.5,line = -.5)
    if (sd==1) mtext('Right correct',side=3,cex=1.5,line = -.5)
    
    # Calculate PvT values for foil primed condition (left)
    cnd = d$CDT < 5 & d$Co == sd
    rt = d$RT[ cnd ]; ac = d$Ac[ cnd ]; sb = d$S[ cnd ]; cv = d$CDT[ cnd ]
    outF = pvt_points( rt, ac, cv, grp = sb, opt = list( out = T, draw = F ) )
    segments( outF$pv$x[1:2], outF$pv$y[1:2],
              outF$pv$x[3:4], outF$pv$y[3:4],
              col = 'grey', lty = 1:2, lwd = 2 )
    
    # Calculate PvT values for target primed condition
    cnd = d$CDT > 4 & d$Co == sd
    rt = d$RT[ cnd ]; ac = d$Ac[ cnd ]; sb = d$S[ cnd ]; cv = d$CDT[ cnd ] - 2
    outT = pvt_points( rt, ac, cv, grp = sb, opt = list( out = T, draw = F ) )
    segments( outT$pv$x[1:2], outT$pv$y[1:2],
              outT$pv$x[3:4], outT$pv$y[3:4],
              col = 'grey', lty = 1:2, lwd = 2 )
    
    # Add in PvT values
    cnd = d$CDT < 5 & d$Co == sd
    rt = d$RT[ cnd ]; ac = d$Ac[ cnd ]; sb = d$S[ cnd ]; cv = d$CDT[ cnd ]
    plt = list( pch = c(19,19,21,21), bg = 'white' )
    pvt_points( rt, ac, cv, grp = sb, plt = plt, cex = 1.5 )
    
    cnd = d$CDT > 4 & d$Co == sd
    rt = d$RT[ cnd ]; ac = d$Ac[ cnd ]; sb = d$S[ cnd ]; cv = d$CDT[ cnd ] - 2
    plt = list( pch = c(17,17,24,24), bg = 'white' )
    pvt_points( rt, ac, cv, grp = sb, plt = plt, cex = 1.5 )
    
  }
  
  # Legends
  par( mar = c(0,0,0,0) )
  blankPlot()
  
  legend( 0,.6, c('Target primed','Foil primed'),
          pch = c(19,17), cex = 1.4, bty = 'n' )
  legend( 0,.5, c('.05 s prime','2 s prime'),
          lty = 1:2, col = 'grey', lwd = 2, cex = 1.4, bty = 'n' )
  legend( 0,.4, c('Error','Correct'),
          pch = c(19,21), pt.bg = 'white', cex = 1.4, bty = 'n' )
  
  mtext( 'Prime duration x type x side x accuracy', side = 3, outer = T, 
         cex = 1.5, line = -2 )
  
}

###
### Examination of biases for responding left versus right
###
# Lookup - 05

if ( runCode[4] ) {
  
  # Calculate accuracy per subject and condition
  ac = aggregate( d$Ac, 
                  list( d$DT, d$S, d$PD, d$PT, d$Co ), p_f )
  colnames( ac ) = c('DT','S','PD','PT','Co','A')
  
  # Reorganize data frame
  ac = cbind(
    DT = ac$DT[ ac$Co == 0 ],
    S = ac$S[ ac$Co == 0 ],
    PD = ac$PD[ ac$Co == 0 ],
    PT = ac$PT[ ac$Co == 0 ],
    E0 = ac$A[ ac$Co == 0, 1],
    E1 = ac$A[ ac$Co == 1, 1],
    C0 = ac$A[ ac$Co == 0, 2],
    C1 = ac$A[ ac$Co == 1, 2] )
  ac = as.data.frame( ac )
  
  # Determine rows with missing data
  tmp = colSums( apply( ac[,5:8], 1, function(x) x == 0 ) )
  miss = which( tmp > 0 )
  
  # Calculate mean RT per subject and condition
  mrt = aggregate( d$RT, 
                   list( d$DT, d$S, d$PD, d$PT, d$Co, d$Ac ), mean )
  colnames( mrt ) = c('DT','S','PD','PT','Co','A','MRT')
  
  # Color code points based on overall accuracy
  MPC = aggregate( d$Ac, list( d$S ), mean )
  ord = order( MPC$x )
  allOrd=ord[ mrt$S ]
  
  # Trim rows with missing data
  all_rmv = c()
  
  chk = ac[ miss, 1:4 ]
  
  for ( i in 1:nrow( chk ) ) {
    
    comp = as.numeric( chk[i,] )
    
    res = apply( mrt[,1:4], 1, function(x) sum( comp == x ) )
    rmv = which( res == 4 )
    all_rmv = c( all_rmv, rmv )
  }

  mrt = mrt[ -unique(all_rmv), ]
  ord = allOrd[ -unique(all_rmv) ]
  
  # Reorganize data frame
  mrt = cbind(
    DT = mrt$DT[ mrt$Co == 0 & mrt$A == 0 ],
    S = mrt$S[ mrt$Co == 0 & mrt$A == 0 ],
    PD = mrt$PD[ mrt$Co == 0 & mrt$A == 0 ],
    PT = mrt$PT[ mrt$Co == 0 & mrt$A == 0 ],
    E0 = mrt$MRT[ mrt$Co == 0 & mrt$A == 0 ],
    E1 = mrt$MRT[ mrt$Co == 1 & mrt$A == 0 ],
    C0 = mrt$MRT[ mrt$Co == 0 & mrt$A == 1 ],
    C1 = mrt$MRT[ mrt$Co == 1 & mrt$A == 1 ] )
  mrt = as.data.frame( mrt )
  
  # Trim row from data frame for accuracy as well
  ac = ac[ -miss, ]
  
  # Plot difference scores for each condition
  if (!savePlot) x11(width=12)
  layout( rbind( 1:4, 5:8 ) )
  
  ttl = c( 'Target primed .05 s',
           'Target primed 2 s',
           'Foil primed .05 s',
           'Foil primed 2 s' )
  
  # Determine color scheme
  clr = ord[1:nrow(ac)]/length( unique( ord ) )
  
  ind = c(1:4,1:4)
  for ( i in 1:8 ) {
    
    cnd = ac$DT == ind[i]
    
    yl = c(-1.2,1.2)
    xl = c(-.6,.6)
    par( mar = c(5,6,3,.5) )
    blankPlot( yDim = yl, xDim = xl )
    title( ttl[i], cex = 1.5 )
    abline(h=0,lwd=2,col='grey')
    abline(v=0,lwd=2,col='grey')
    axis( 2, round( seq( -1.2, 1.2, .4 ), 2), cex.axis = 1.5,
          lwd = 2 )
    axis( 1, round( seq( -.6, .6, .3 ), 2), cex.axis = 1.5,
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
  mtext('Right - left (accuracy)',side=1,outer=T,cex=1.5,
        line = -2 )
  
}

###
### Response time and choice distributions by offset
###
# Lookup - 06

if ( runCode[5] ) {
  
  if (!savePlot) x11(width=12)
  layout( rbind(1:5) )
  
  T_x = function(x) quantile(x,prob=seq(.1,.9,.2))
  
  
  ttl = paste( 'Target offset', c('-33','-17','0','+17','+33') )
  for (i in 1:5) {
    
    # Set up each plot
    par( mar = c( 5, 1, 1, 5 ) )
    blankPlot(xDim=c(.2,2.2),yDim=c(-.5,1))
    axis( 4, seq(-1,1,.5), abs( seq(-1,1,.5) ), cex.axis = 1.5,
          lwd = 2 )
    axis( 1, seq(.2,2.2,.5), cex.axis = 1.5, lwd = 2 )
    
    # Condition on offset
    cnd = d$O == i
    rt = d$RT[cnd]; ac = d$Ac[cnd]; sb = d$S[cnd]
    out = cdf_curve( rt, ac, grp = sb, 
                     opt = list( draw = F, out = T ), lwd = 2 )
    add_segments( out, T_x = T_x, col = 'grey', lwd = 2 )
    out = cdf_curve( rt, ac, grp = sb, sel = 0, 
                     opt = list( out = T, draw = F, flip = T ), 
                     lwd = 2, lty = 2 )
    add_segments( out, T_x = T_x, col = 'grey', lwd = 2, lty = 2 )
    abline( h = 0, lwd = 2 )
    cdf_curve( rt, ac, grp = sb, lwd = 2 )
    cdf_curve( rt, ac, grp = sb, sel = 0, opt = list( flip = T ), 
               lwd = 2, lty = 2 )
    legend( 'topleft', ttl[i], bty = 'n', cex = 1.5 )
    
  }
  mtext('RT (s)', side = 1, outer = T, cex = 1.5, line = -2 )
  mtext('Joint CDF', side = 4, outer = T, cex = 1.5, line = -2 )

}

###
### Fastest RT by offset and condition
###
# Lookup - 07

if ( runCode[6] ) {
  
  fRT = aggregate( d$RT, list( d$O, d$Ac, d$S ), min )
  colnames( fRT ) = c('O','A','S','T')
  yl = lowerUpper( .2, fRT$T )
  
  if (!savePlot) x11(width=12);
  dmn = matrix( 1:32, 4, 8, byrow = T )
  dmn[ dmn == N + 1 | dmn == N + 2 ] = N+1
  dmn[ dmn > N+1 ] = N+2
  layout( dmn )
  
  for (n in 1:N) {
    
    par( mar=c( 2, 2, 2, 2 ) )
    blankPlot( xDim=c(.5,5.5), yDim = yl )
    abline(h=.2,lwd=2)
    abline(v=.5,lwd=2)
    
    if ( n == 1 | n == 9 | n == 17 | n == 25 ) {
      axis( 2, seq( yl[1], yl[2], .5 ), 
            tick = F, cex = 1.5, line = -.5 )
    }
    if ( n > N - 3 ) {
      axis( 1, 1:5, c('-33','-17','0','+17','+33'),
            cex = .8 )
    }
    
    sel = fRT$A == 1 & fRT$S == n
    lines( fRT$O[sel], fRT$T[sel], lwd = 2 )
    points( fRT$O[sel], fRT$T[sel], pch = 19, cex = 1.5 )
    sel = fRT$A == 0 & fRT$S == n
    lines( fRT$O[sel], fRT$T[sel], lty = 2, lwd = 2 )
    points( fRT$O[sel], fRT$T[sel], pch = 24, bg = 'white', cex = 1.5 )
    
    if (n!=20) legend( 'topright', as.character(n), bty = 'n', cex = 1.5 ) else 
      legend( 'bottomright', as.character(n), bty = 'n', cex = 1.5 )
    
    
  }
  par( mar = c(0, 0, 0, 0) )
  blankPlot()
  legend( 'left', c('Correct','Error'),
          pch = c(19,24), bg = 'white', bty = 'n',
          cex = 1.5 )
  
  blankPlot()
  legend( 'left', 'Fastest RT by subject and offset',
          bty = 'n', cex = 2 )
  
  
}

###
### Plot nROUSE model predictions against observed data
###
# Lookup - 08

if ( runCode[7] ) {
  
  # Calculate P(Correct) by condition
  y = aggregate( d$Ac, list( d$PD, d$PT ), mean )
  tot = aggregate( rep(1,nrow(d)), list(d$DT), sum )$x
  colnames( y ) = c('PD','PT','P')
  x = log( c( .05, 2 ) )
  # Uncertainty intervals
  tmp = aggregate( d$Ac, list( d$PD, d$PT, d$S ), mean )
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
  axis( 1, log( c( .05, 2 ) ),
        c( 50, 2000 ), cex.axis = 1.5, lwd = 2 )
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
  nDat = cbind( TarDur = 50, MaskDur = 450, PrimeDur = c(50,2000,50,2000),
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