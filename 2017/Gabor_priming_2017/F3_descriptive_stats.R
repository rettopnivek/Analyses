#--------------------#
# Data preprocessing #
# Kevin Potter       #
# Updated 12/06/2017 #
#--------------------#

# Clear workspace
rm(list = ls())

# Save current working directory
orig_dir = getwd()

# Indicator if data from Experiment 1 or 2 should be loaded in
Exp = 1

# Indicate whether to create a PDF of figures
savePlot = F

# Indicate which code segments to run
runCode = c(
  T,  # Accuracy by prime type and duration (Group level)
  F,  # Accuracy by prime type and duration (Subject level)
  T,  # Choice by target type, prime type and duration (Group level)
  F,  # Choice by target type, prime type and duration (Subject level)
  F,  # ???
  F,  # Response time by prime type and duration (Group level)
  F,  # Response time by prime type and duration (Subject level)
  F,  # Response time by target type, prime type and duration (Group level)
  F   # Response time by target type, prime type and duration (Subject level)
)

# Table of contents for figures
table_of_contents = c(
  'Accuracy by prime type and duration (Group level)',
  'Accuracy by prime type and duration (Subject level)',
  'Choice by target type, prime type and duration (Group level)',
  'Choice by target type, prime type and duration (Subject level)',
  '???', 
  'Response time by prime type and duration (Group level)',
  'Response time by prime type and duration (Subject level)',
  'Response time by target type, prime type and duration (Group level)',
  'Response time by target type, prime type and duration (Subject level)'
)

# Index
# Lookup - 01:  Load in data, useful packages, and functions
# Lookup - 02:  Accuracy by prime type and duration (Group level)
# Lookup - 03:  Accuracy by prime type and duration (Subject level)
# Lookup - 04:  Choice by target type, prime type and duration (Group level)
# Lookup - 05:  Choice by target type, prime type and duration (Subject level)
# Lookup - 06:  Response time by prime type and duration (Group level)
# Lookup - 07:  Response time by prime type and duration (Subject level)
# Lookup - 08:  Response time by target type, prime type and duration 
#               (Group level)
# Lookup - 09:  Response time by target type, prime type and duration 
#               (Subject level)
# Lookup - 10:  

###
### Load in useful packages, functions, and data
###
# Lookup - 01

### Load in useful packages

# For geting github packages
# install.packages(devtools)
# library(devtools)

# Miscellanous functions for modeling/plotting
# install_github("rettopnivek/utilityf")
library(utilityf)

### Define useful functions

source( 'F1_useful_functions.R' )

### Load in data

if ( Exp == 1 ) {
  dName = 'Gabor_Exp_1_11-17.RData'
} else {
  dName = 'Gabor_Exp_2_11-17.RData'
}

setwd( 'Data' )
load( dName )
rm( dName )
setwd( orig_dir )

# For easy variable access
ad = curDat
colnames( ad ) = simplify_variable_names( curDat )

### Initialize plot

if ( savePlot ) {
  setwd( 'Plots' )
  fName = paste( 'Descriptive_results_Exp_', Exp, '.pdf', sep = '' )
  pdf( file = fName, width = 12 )
  rm( fName )
  setwd( orig_dir )
}

# Create table of contents
if ( !savePlot ) x11( width = 12 )
shft = c( 0, ceiling( N/6 ) - 1, 0, ceiling( N/6 ) - 1, 0 )
shft = cumsum( shft[ runCode] )
tableContents( table_of_contents[ runCode ], txtSz = 1.5,
               shift = shft )

###
### Accuracy and median RT by prime type and duration
###
# Lookup - 02

if ( runCode[1] ) {
  
  # Experimental trials only
  cd = ad[ ad$BT == 2, ]
  
  # Define test statistics of interest
  T_x = function( x ) {
    out = c( Y = sum( x == 1 ),
             P = mean( x ),
             N = length( x ),
             SE = sd( x )/sqrt( 5 )
    )
    
    return( out )
  }
  
  # Obtain aggregate statistics
  glp = aggregate( cd$Ac, list( 
    cd$PD, cd$PTL ), T_x )
  colnames( glp ) = c( 'PD', 'PT', 'X' )
  glp$Y = glp$X[,1]; glp$P = glp$X[,2]; glp$N = glp$X[,3]
  glp$SE = glp$X[,4];
  glp$X = NULL
  
  tmp = aggregate( cd$Ac, list( 
    cd$PD, cd$PTL, cd$S ), mean )
  tmp = aggregate( tmp$x, list( tmp[,1], tmp[,2] ),
                   function(x) return( quantile( x, c( .25, .75) ) ) )
  glp$UIL = tmp$x[,1]; glp$UIU = tmp$x[,2]; rm( tmp )
  
  if (!savePlot ) x11( width = 12 )
  layout( cbind( 1, 1, 2, 2 ) )
  
  # Plotting characterstics
  lnSz = 2
  lbSz = 1.4
  axSz = 2
  axPos = -1.4
  ptSz = 1.5
  pt = c( 21, 24, 22 )
  clr = cbbPalette[c(2,3,4)]
  
  # Extract prime durations
  pd = sort( unique( glp$PD ) )
  
  # Create blank plot
  xl = c( .5, length(pd) + .5 )
  yl = c( 0, 1 )
  blankPlot( xl, yl )
  
  # Add axes
  horizLines( c(.25,.5,.75), xl, lwd = lnSz, col = 'grey' )
  customAxes( xl, yl, label = c( 'Prime duration (ms)', 'P(Correct)' ),
              lbSz = lbSz, lnSz = lnSz )
  axis( 2, seq(0,1,.25), tick = F, line = axPos, cex.axis = axSz )
  axis( 1, 1:length(pd), pd, tick = F, 
        line = axPos, cex.axis = axSz )
  
  # Add data
  cnd = unique( glp$PT )
  srt = c( grep( 'Target', cnd ), grep( 'Foil', cnd ), grep( 'Neither', cnd ) )
  cnd = cnd[srt]
  adj = seq( -.15, .15, length = length( cnd ) )
  for ( i in 1:length(cnd) ) {
    sel = glp$PT == cnd[i]
    xa = 1:length(pd); ya = glp$P[sel]
    
    ui = matrix( NA, 2, length(pd) )
    ui[1,] = glp$UIL[sel]
    ui[2,] = glp$UIU[sel]
    
    errorBars( xa + adj[i], ui, lwd = lnSz, length = .025,
               col = clr[i] )
    
    lines( xa + adj[i], ya, lwd = lnSz, col = clr[i] )
    points( xa + adj[i], ya, pch = pt[i], bg = clr[i], cex = ptSz,
            col = clr[i] )
  }
  
  # Add legend 
  blankPlot()
  legend( 'center', as.character( cnd ), pch = pt, pt.bg = clr, col = clr,
          bty = 'n', cex = axSz )
  legend( 'top', 'Error bars: Inter-quartile range', cex = axSz,
          bty = 'n' )
  
  mtext( table_of_contents[1], side = 3, outer = T, line = -2, cex = lbSz )
}

###
### Accuracy by prime type and duration (Subject level)
###
# Lookup - 03

if ( runCode[2] ) {
  
  # Experimental trials only
  cd = ad[ ad$BT == 2, ]
  
  # Limits for uncertainty interval
  interval = quickCoverage( pnorm( 1 ) - pnorm( -1 ) )
  
  # Define test statistics of interest
  T_x = function( x ) {
    out = c( Y = sum( x == 1 ),
             P = mean( x ),
             N = length( x )
    )
    
    return( out )
  }
  
  slp = aggregate( cd$Ac, list( 
    cd$PD, cd$PTL, 
    cd$S ), T_x )
  colnames( slp ) = c( 'PD', 'PT', 'S', 'X' )
  slp$Y = slp$X[,1]; slp$P = slp$X[,2]; slp$N = slp$X[,3]
  slp$X = NULL
  slp$UIL = qbinom( interval[1], slp$N, slp$P )/slp$N
  slp$UIU = qbinom( interval[2], slp$N, slp$P )/slp$N
  
  # Extract prime durations
  pd = sort( unique( slp$PD ) )
  
  # Plotting characterstics
  lnSz = 2
  lbSz = 1.4
  axSz = 1.2
  axPos = -1.4
  ptSz = 1.8
  pt = c( 21, 24, 22 )
  clr = cbbPalette[c(2,3,4)]
  
  # Plot 6 subjects at a time
  nWindows = ceiling( N/6 )
  
  for ( j in 1:nWindows ) {
    
    if ( !savePlot ) x11( width = 12 )
    layout( matrix( 1:6, 2, 3, byrow = T ) )
    
    # Indices
    st = 1 + 6 * ( j - 1 )
    en = min( 6 + 6 * ( j - 1 ), N )
    
    for ( n in st:en ) {
      
      # Create blank plot
      xl = c( .5, length(pd) + .5 )
      yl = c( 0, 1 )
      blankPlot( xl, yl )
      
      # Add axes
      horizLines( c(.25,.5,.75), xl, lwd = lnSz, col = 'grey' )
      customAxes( xl, yl, lnSz = lnSz )
      axis( 2, seq(0,1,.25), tick = F, line = axPos, cex.axis = axSz )
      axis( 1, 1:length(pd), pd, tick = F, 
            line = axPos, cex.axis = axSz )
      title( paste( 'Subject', n ) )
      
      # Add data
      cnd = unique( slp$PT )
      srt = c( grep( 'Target', cnd ), 
               grep( 'Foil', cnd ), 
               grep( 'Neither', cnd ) )
      cnd = cnd[srt]
      adj = seq( -.15, .15, length = length( cnd ) )
      
      for ( i in 1:length(cnd) ) {
        
        sel = slp$PT == cnd[i] & slp$S == n
        xa = 1:length(pd); ya = slp$P[sel]
        
        ui = matrix( NA, 2, length(pd) )
        ui[1,] = slp$UIL[sel]
        ui[2,] = slp$UIU[sel]
        
        errorBars( xa + adj[i], ui, lwd = lnSz, length = .025,
                   col = clr[i] )
        
        lines( xa + adj[i], ya, lwd = lnSz, col = clr[i] )
        points( xa + adj[i], ya, pch = pt[i], bg = clr[i], cex = ptSz,
                col = clr[i] )
      }
      
      # Add legends
      if ( n %% 4 == 0 ) {
        
        # Turn off clipping
        par( xpd = TRUE )
        
        plt_sz = dev.size()
        legend( 0, -.075, 
                as.character( cnd ), pch = pt, pt.bg = clr, col = clr,
                bty = 'n', cex = axSz, horiz = T )
        
        legend( 0, -.15, paste( 'Error bars: ', 
                                round( 
                                  diff( interval ) * 100 ), '% (Binomial)',
                                sep = '' ), cex = axSz,
                bty = 'n', hori = T )
        
        # Turn on clipping
        par( xpd = FALSE )
        
      }
      
    }
    
    # Axis labels
    mtext( 'Prime duration', side = 1, outer = T, line = -1.5, cex = lbSz )
    mtext( 'P(Correct)', side = 2, outer = T, line = -1.7, cex = lbSz )
    
  }
  
}

###
### Choice by target type, prime type and duration (Group level)
###
# Lookup - 04

if ( runCode[3] ) {
  
  # Experimental trials only
  cd = ad[ ad$BT == 2, ]
  
  # Define test statistics of interest
  T_x = function( x ) {
    out = c( Y = sum( x == 1 ),
             P = mean( x ),
             N = length( x ),
             SE = sd( x )/sqrt( 5 )
    )
    
    return( out )
  }
  
  # Obtain aggregate statistics
  glp = aggregate( cd$Ch, list( 
    cd$PD, cd$PTL, cd$CoL ), T_x )
  colnames( glp ) = c( 'PD', 'PT', 'Co', 'X' )
  glp$Y = glp$X[,1]; glp$P = glp$X[,2]; glp$N = glp$X[,3]
  glp$SE = glp$X[,4];
  glp$X = NULL
  
  tmp = aggregate( cd$Ch, list( 
    cd$PD, cd$PTL, cd$CoL, 
    cd$S ), mean )
  tmp = aggregate( tmp$x, list( tmp[,1], tmp[,2], tmp[,3] ),
                   function(x) return( quantile( x, c( .25, .75) ) ) )
  glp$UIL = tmp$x[,1]; glp$UIU = tmp$x[,2]; rm( tmp )
  
  # Extract prime durations
  pd = sort( unique( glp$PD ) )
  # Extract prime types
  ptl = unique( glp$PT )
  # Total number of points to plot
  tp = length( pd ) * 2
  
  if (!savePlot ) x11( width = 12 )
  layout( cbind( 1, 1, 2, 2 ) )
  
  # Plotting characterstics
  lnSz = 2
  lbSz = 1.4
  axSz = 1.5
  axPos = -1.4
  ptSz = 1.5
  pt = c( 21, 24, 22 )
  pt = rep( pt[ 1:length(ptl) ], 2 )
  clr = rep( cbbPalette[c(2,3,4)[1:length(ptl)]], each = 2 )
  lnt = rep( 1, tp )
  
  # Create blank plot
  xl = c( .5, tp + .5 )
  yl = c( 0, 1 )
  blankPlot( xl, yl )
  
  # Add axes
  horizLines( c(.25,.5,.75), xl, lwd = lnSz, col = 'grey' )
  customAxes( xl, yl, label = c( 'Prime duration (ms)', 'P(Y = Right)' ),
              lbSz = lbSz, lnSz = lnSz )
  axis( 2, seq(0,1,.25), tick = F, line = axPos, cex.axis = axSz )
  axis( 1, 1:tp, rep( sort( unique( glp$PD ) ), 2 ), tick = F, 
        line = axPos, cex.axis = axSz )
  
  # Add data
  cnd = aggregate( rep( 1, nrow( cd ) ), list( cd$CoL, cd$PTL ), unique )
  cnd$x = NULL
  colnames( cnd ) = c( 'Co', 'PT' )
  srt = c( grep( 'Target', cnd$PT ), 
           grep( 'Foil', cnd$PT ), 
           grep( 'Neither', cnd$PT ) )
  cnd = cnd[srt,]
  adj = rep( seq( -.15, .15, length = nrow( cnd )/2 ), 2 )
  shft = rep( 0, nrow(cnd) );
  shft[ cnd$Co == 'Right' ] = length(pd)
  lnt[ cnd$Co == 'Right' ] = 2
  
  for ( i in 1:nrow(cnd) ) {
    
    sel = glp$PT == cnd$PT[i] & glp$Co == cnd$Co[i]
    xa = 1:length(pd) + shft[i]; ya = glp$P[sel]
    
    ui = matrix( NA, 2, length(pd) )
    ui[1,] = glp$UIL[sel]
    ui[2,] = glp$UIU[sel]
    
    errorBars( xa + adj[i], ui, lwd = lnSz, length = .025,
               col = clr[i] )
    
    lines( xa + adj[i], ya, lwd = lnSz, col = clr[i], lty = lnt[i] )
    points( xa + adj[i], ya, pch = pt[i], bg = clr[i], cex = ptSz,
            col = clr[i] )
    
  }
  
  # Add legend 
  blankPlot()
  legend( 'center', as.character( unique( cnd$PT ) ), 
          pch = pt, pt.bg =  unique( clr ), col = unique( clr ),
          bty = 'n', cex = axSz )
  legend( 'bottom', paste( unique( cnd$Co ), 'correct' ), 
          lty = unique( lnt ), lwd = lnSz, 
          bty = 'n', cex = axSz )
  legend( 'top', 'Error bars: Inter-quartile range', cex = axSz,
          bty = 'n' )
  
  mtext( table_of_contents[3], side = 3, outer = T, line = -2, cex = lbSz )
}

###
###
###
# Lookup - 05

if ( runCode[4] ) {
  
  # Experimental trials only
  cd = ad[ ad$BT == 2, ]
  
  # Limits for uncertainty interval
  interval = quickCoverage( pnorm( 1 ) - pnorm( -1 ) )
  
  # Define test statistics of interest
  T_x = function( x ) {
    out = c( Y = sum( x == 1 ),
             P = mean( x ),
             N = length( x )
    )
    
    return( out )
  }
  
  slp = aggregate( cd$Ch, list( 
    cd$PD, cd$PTL, cd$CoL, 
    cd$S ), T_x )
  colnames( slp ) = c( 'PD', 'PT', 'Co', 'S', 'X' )
  slp$Y = slp$X[,1]; slp$P = slp$X[,2]; slp$N = slp$X[,3]
  slp$X = NULL
  slp$UIL = qbinom( interval[1], slp$N, slp$P )/slp$N
  slp$UIU = qbinom( interval[2], slp$N, slp$P )/slp$N
  
  # Extract prime durations
  pd = sort( unique( slp$PD ) )
  # Extract prime types
  ptl = unique( slp$PT )
  # Total number of points to plot
  tp = length( pd ) * 2
  
  # Plotting characterstics
  lnSz = 2
  lbSz = 1.4
  axSz = 0.9
  axPos = -1.4
  ptSz = 1.5
  pt = c( 21, 24, 22 )
  pt = rep( pt[ 1:length(ptl) ], 2 )
  clr = rep( cbbPalette[c(2,3,4)[1:length(ptl)]], each = 2 )
  lnt = rep( 1, tp )
  
  # Plot 6 subjects at a time
  nWindows = ceiling( N/6 )
  
  for ( j in 1:nWindows ) {
    
    if ( !savePlot ) x11( width = 12 )
    layout( matrix( 1:6, 2, 3, byrow = T ) )
    
    # Indices
    st = 1 + 6 * ( j - 1 )
    en = min( 6 + 6 * ( j - 1 ), N )
    
    for ( n in st:en ) {
      
      # Create blank plot
      xl = c( .5, tp + .5 )
      yl = c( 0, 1 )
      blankPlot( xl, yl )
      
      # Add axes
      horizLines( c(.25,.5,.75), xl, lwd = lnSz, col = 'grey' )
      customAxes( xl, yl, lnSz = lnSz )
      axis( 2, seq(0,1,.25), tick = F, line = axPos, cex.axis = axSz )
      axis( 1, 1:tp, rep(pd,2), tick = F, 
            line = axPos, cex.axis = axSz )
      title( paste( 'Subject', n ) )
      
      # Add data
      cnd = aggregate( rep( 1, nrow( cd ) ), list( cd$CoL, cd$PTL ), unique )
      cnd$x = NULL
      colnames( cnd ) = c( 'Co', 'PT' )
      srt = c( grep( 'Target', cnd$PT ), 
               grep( 'Foil', cnd$PT ), 
               grep( 'Neither', cnd$PT ) )
      cnd = cnd[srt,]
      adj = rep( seq( -.15, .15, length = nrow( cnd )/2 ), 2 )
      shft = rep( 0, nrow(cnd) );
      shft[ cnd$Co == 'Right' ] = length(pd)
      lnt[ cnd$Co == 'Right' ] = 2
      
      for ( i in 1:nrow(cnd) ) {
        
        sel = slp$PT == cnd$PT[i] & slp$Co == cnd$Co[i] & slp$S == n
        xa = 1:length(pd) + shft[i]; ya = slp$P[sel]
        
        ui = matrix( NA, 2, length(pd) )
        ui[1,] = slp$UIL[sel]
        ui[2,] = slp$UIU[sel]
        
        errorBars( xa + adj[i], ui, lwd = lnSz, length = .025,
                   col = clr[i] )
        
        lines( xa + adj[i], ya, lwd = lnSz, col = clr[i], lty = lnt[i] )
        points( xa + adj[i], ya, pch = pt[i], bg = clr[i], cex = ptSz,
                col = clr[i] )
        
        
      }
      
      # Add legends
      if ( n %% 4 == 0 ) {
        
        # Turn off clipping
        par( xpd = TRUE )
        
        plt_sz = dev.size()
        legend( 0, -.075, 
                as.character( unique( cnd$PT ) ), 
                pch = pt, pt.bg = unique( clr ), col = unique( clr ),
                bty = 'n', cex = axSz, horiz = T )
        legend( 0, -.15, 
                paste( unique( cnd$Co ), 'correct' ), 
                lty = unique( lnt ), 
                bty = 'n', cex = axSz, horiz = T )
        legend( 0, -.225, paste( 'Error bars: ', 
                                round( 
                                  diff( interval ) * 100 ), '% (Binomial)',
                                sep = '' ), cex = axSz,
                bty = 'n', hori = T )
        
        # Turn on clipping
        par( xpd = FALSE )
        
      }
      
    }
    
    # Axis labels
    mtext( 'Prime duration', side = 1, outer = T, line = -1.5, cex = lbSz )
    mtext( 'P(Y = Right)', side = 2, outer = T, line = -1.7, cex = lbSz )
    
  }
  
}

###
###
###
# Lookup - 06

if ( runCode[5] ) {
  
  # Experimental trials only
  cd = ad[ ad$BT == 2, ]
  
  # Define test statistics of interest
  T_x = function( x ) {
    out = c( Y = sum( x == 1 ),
             P = mean( x ),
             N = length( x )
    )
    
    return( out )
  }
  
  # Average accuracy performance over prime duration 
  # and type
  slp = aggregate( cd$Ac, list( cd$PD, cd$PTL, cd$S ), T_x )
  colnames( slp ) = c( 'PD', 'PT', 'S', 'X' )
  slp$Y = slp$X[,1]; slp$P = slp$X[,2]; slp$N = slp$X[,3]
  slp$X = NULL
  
  # Number of prime types
  plt = unique( slp$PT )
  np = length(plt)
  
  # Wide-form matrix
  W = matrix( NA, nrow( slp )/np, np*3 )
  lbl = sapply( plt, 
                function(s) return( strsplit( s, split = '' )[[1]][1] ) )
  colnames( W ) = c(
    paste( lbl, 'Y', sep ='-' ),
    paste( lbl, 'P', sep ='-' ),
    paste( lbl, 'N', sep ='-' )
  )
  
  # Fill in matrix
  for ( i in 1:np ) {
    W[,i] = slp$Y[ slp$PT == plt[i] ]
    W[,i+np] = slp$P[ slp$PT == plt[i] ]
    W[,i+np*2] = slp$N[ slp$PT == plt[i] ]
  }
  
  # Compute difference score between target and foil performance
  DbTaF = W[,'T-P'] - W[,'F-P']
  if ( Exp == 1 ) 
    DbNaA = W[,'N-P'] - rowMeans( W[,4:6] ) else 
      DbNaA = rep( NA, nrow(W) )
  
  sel = slp$PT == plt[1]
  slp = data.frame( S = slp$S[sel],
              PD = slp$PD[sel],
              DbTaF = DbTaF,
              DbNaA = DbNaA
  )
  slp = cbind( slp, W )
  
  # Group-level
  i = quickCoverage( pnorm(2) - pnorm(-2) )
  T_x = function(x) {
    out = numeric(4)
    out[1] = mean(x)
    out[2] = sd(x)/sqrt(length(x))
    out[3] = out[1] - out[2]*qt(i[1],N-1)
    out[4] = out[1] + out[2]*qt(i[1],N-1)
    out[5:6] = quantile(x,c(.25,.75))
    
    names( out ) = c( 'Mean', 'SE', 'LUI', 'UUI', 'LIQ', 'UIQ' )
    
    return( out )
  }
  glp = aggregate( slp$DbTaF, list( slp$PD ), T_x )
  colnames( glp ) = c( 'PD', 'X' )
  glp = cbind( glp, glp$X ); glp$X = NULL
  
  # Plotting characterstics
  lnSz = 2
  lbSz = 1.4
  axSz = 1.5
  axPos = -1.4
  ptSz = 1.5
  pt = 21
  clr = 'black'
  
  if ( !savePlot ) x11( width = 12 )
  layout( cbind( 2, 1, 1, 2 ) )
  
  # Number of prime durations
  pd = length( glp$PD )
  
  # Create a blank plot
  xl = c( .5, pd + .5 )
  yl = lowerUpper( .25, as.numeric( unlist( glp[,-1] ) ) )
  blankPlot( xl, yl )
  
  # Add axes
  horizLines( seq( yl[1], yl[2], .1 ), xl, col = 'grey', lwd = lnSz )
  customAxes( xl, yl )
  axis( 2, seq( yl[1], yl[2], .25 ), 
        tick = F, line = axPos, cex.axis = axSz )
  axis( 1, 1:pd, sort( unique( glp$PD ) ), 
        tick = F, line = axPos, cex.axis = axSz )
  mtext( 'Target - foil primed (accuracy)', side = 2, 
         cex = lbSz, line = 1.5 )
  mtext( 'Prime duration (ms)', side = 1, 
         cex = lbSz, line = 1.5 )
  
  # Add data
  lines( 1:pd, glp$Mean, lwd = lnSz, col = clr )
  segments( 1:pd, glp$LIQ, 1:pd, glp$UIQ, lwd = lnSz, col = 'grey50' )
  errorBars( 1:pd, rbind( glp$LUI, glp$UUI ), lwd = lnSz, length = .05 )
  points( 1:pd, glp$Mean, pch = pt, bg = clr, col = clr, cex = ptSz )
  
  par( xpd = NA )
  legend( 'top', 'Error bars: Inner = 95% CI, Outer = IQR', 
          bty = 'n', cex = axSz )
  par( xpd = FALSE )
  
}

###
###
###
# Lookup - 07

if ( runCode[6] ) {
  
  # Experimental trials only
  cd = ad[ ad$BT == 2, ]
  
  # Define test statistics of interest
  T_x = function( x ) {
    out = c( Y = sum( x == 1 ),
             P = mean( x ),
             N = length( x )
    )
    
    return( out )
  }
  
  # Compute accuracy over ratio of target to foil contrast,
  # prime duration and type, and block
  if ( Exp == 2 ) {
    ToFCR = cd$ToR
  } else {
    
    qnt = quantile( cd$ToR, seq( .1, .9, .2 ) )
    ToFCR = rep( qnt[1], nrow( cd ) )
    
    for ( i in 2:length( qnt) ) {
      sel = cd$ToR > qnt[i-1] & cd$ToR <= qnt[i]
      ToFCR[sel] = qnt[i]
    }
    
  }
  slp = aggregate( cd$Ac, list( ToFCR, cd$PD, cd$PTL ), T_x )
  colnames( slp ) = c( 'R', 'PD', 'PT', 'X' )
  slp$Y = slp$X[,1]; slp$P = slp$X[,2]; slp$N = slp$X[,3]
  slp$X = NULL
  
  # Difference scores
  ds = slp$P[ grep( 'Target', slp$PT ) ] - slp$P[ grep( 'Foil', slp$PT ) ]
  slp$DS = NA
  slp$DS[ grep( 'Target', slp$PT ) ] = ds
  slp$DS[ grep( 'Foil', slp$PT ) ] = -ds
  
  # Prime durations
  pd = sort( unique( slp$PD ) )
  
  # Plotting characterstics
  lnSz = 2
  lbSz = 1.4
  axSz = 1.2
  axPos = -1.3
  ptSz = 1.75
  pt = 21
  clr = cbbPalette[ 1:length( pd ) + 1 ]
  
  
  if( !savePlot ) x11( width = 12 )
  layout( cbind( 1, 2 ) )
  
  inc = 1
  for ( ptl in c( 'Foil', 'Target' ) ) {
    
    # Create blank plot
    xl = lowerUpper( .25, slp$R )
    yl = c( 0, 1 )
    blankPlot( xl, yl )
    
    # Add axes
    horizLines( c(.25,.5,.75), xl, col = 'grey', lwd = lnSz )
    customAxes( xl, yl, lnSz = lnSz )
    axis( 1, round( seq( xl[1], xl[2], .5 ), 1 ),
          tick = F, line = axPos, cex.axis = axSz )
    axis( 2, seq(0,1,.25), 
          tick = F, line = axPos, cex.axis = axSz )
    
    cnd = sort( unique( slp$PD ) )
    for ( i in 1:length(cnd) ) {
      
      
      sel = slp$PD == cnd[i]
      sel[ grep( ptl, slp$PT ) ] = FALSE
      sel[ grep( 'Neither', slp$PT ) ] = FALSE
      
      xa = slp$R[sel]
      ya = slp$P[sel]
      no_na = !is.na( ya )
      xa = xa[ no_na ]
      ya = ya[ no_na ]
      o = order( xa )
      
      lines( xa[o], ya[o], col = clr[i], lwd = lnSz, lty = inc )
      points( xa[o], ya[o], pch = pt, bg = clr[i], col = clr[i] )
      
    }
    
    inc = inc + 1
  }
  
  # Add legend
  # Turn off clipping
  par( xpd = NA )
  
  legend( -3.5, -.175, paste ( pd, 'ms prime' ),
          fill = clr, bty = 'n', cex = lbSz, horiz = T )
  
  # Turn on clipping
  par( xpd = FALSE )
  
}

###
###
###
# Lookup - 08

###
###
###
# Lookup - 09

# Close connection to PDf file
if ( savePlot ) dev.off()
setwd( orig_dir )
