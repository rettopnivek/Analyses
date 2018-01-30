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
Exp = 2

# Indicate whether to create a PDF of figures
savePlot = T

# Indicate whether trimmed data should be saved
saveFile = T

# Indicate which code segments to run
runCode = c(
  T,  # Overall accuracy by subject
  T,  # Distribution of fastest/slowest RTs
  T,  # Accuracy predicts fastest/slowest RTs
  T   # Conditional accuracy functions
)

# Table of contents for figures
table_of_contents = c(
  'Overall accuracy by subject',
  'Distributions of fastest/slowest RTs',
  'Accuracy predicts fastest/slowest RTs',
  'Conditional accuracy functions'
)

# Index
# Lookup - 01:  Load in data, useful packages, and functions
# Lookup - 02:  Overall accuracy by subject
# Lookup - 03:  Distribution of fastest/slowest RTs
# Lookup - 04:  Accuracy predicts fastest/slowest RTs
# Lookup - 05:  Conditional accuracy functions

###
### Load in data, useful packages, and functions
###
# Lookup - 01

### Load in useful packages

# For geting github packages
# install.packages('devtools')
# library(devtools)

# Miscellanous functions for modeling/plotting
# install_github("rettopnivek/utilityf")
library(utilityf)

# Package for generalized additive model
# install.packages('mgcv')
library( mgcv )

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
rd = rawDat
colnames( rd ) = simplify_variable_names( rawDat )

### Initialize plot

if ( savePlot ) {
  setwd( 'Plots' )
  fName = paste( 'Data_trimming_results_Exp_', Exp, '.pdf', sep = '' )
  pdf( file = fName, width = 12 )
  rm( fName )
  setwd( orig_dir )
}

# Create table of contents
if ( !savePlot ) x11( width = 12 )
tableContents( table_of_contents[ runCode ], txtSz = 2 )

###
### Overall accuracy by subject
###
# Lookup - 02

if ( runCode[1] ) {
  
  # Consider main trials only
  sel = rd$ BT == 2
  cd = rd[ sel, ]
  
  # Compute observed accuracy for each subject
  ac = aggregate( cd$Ac, list( cd$S ), function(x)
    return( c( F = sum(x), N = length(x), P = mean(x) ) ) )
  ac = data.frame( S = ac[,1], F = ac$x[,1], N = ac$x[,2], P = ac$x[,3] )
  
  # Uncertainty interval based on binomial distribution
  interval = quickCoverage( .95 )
  ac$LB = qbinom( interval[1], ac$N, ac$P )/ac$N
  ac$UB = qbinom( interval[2], ac$N, ac$P )/ac$N
  
  # Sort subjects from lowest to highest accuracy
  ord = order( ac$P )
  
  # Plotting characteristics
  txtSz = 1.4
  lnSz = 2
  ptSz = 1
  
  if ( !savePlot ) x11( width = 12 )
  
  # Create blank plot
  xl = c( .5, N + .5 ); yl = c( 0, 1 );
  blankPlot( xl, yl );
  
  # Add grid lines, axes, and labels
  horizLines( seq( .1, 1, .1 ), xl, col = 'grey80', lwd = lnSz )
  customAxes( xl, yl, lnSz = lnSz )
  axis( 2, seq( 0, 1, .2 ), 
        tick = F, line = -1.5, cex.axis = txtSz )
  mtext( 'Overall Accuracy', side = 2, line = 2, cex = txtSz )
  mtext( 'Subjects', side = 1, line = 1, cex = txtSz )
  
  # Overall average performance
  horizLines( mean( ac$P ), xl, lwd = lnSz, lty = 2, col = 'blue' )
  
  # Desired level of performance via staircase method
  horizLines( .8, xl, lwd = lnSz, lty = 2, col = 'red' )
  
  # Add points and error bars
  errorBars( 1:N, rbind( ac$LB[ord], ac$UB[ord] ), length = .05, lwd = lnSz )
  points( 1:N, ac$P[ord], pch = 19, cex = ptSz )
  
  # Legend
  legend( N*.25, .1, c( 'Average performance', 'Desired criterion level' ), 
         fill = c( 'blue', 'red' ), bty = 'n', cex = txtSz,
         horiz = T )
  
}

###
### Distribution of fastest/slowest RTs
###
# Lookup - 03

if ( runCode[2] ) {
  
  # Consider main trials only
  sel = rd$ BT == 2
  cd = rd[ sel, ]
  
  # Compute RT quantiles from min to max
  rt = aggregate( cd$RT, list( cd$S ),
                  function(x) return( c( min(x), max(x) ) ) )
  rt = data.frame( S = rt[,1], Fast = rt$x[,1], Slow = rt$x[,2] )
  
  # Plotting characteristics
  txtSz = 1.4
  lnSz = 2
  ptSz = 1
  
  if ( !savePlot ) x11( width = 12 )
  layout( cbind( 1, 2 ) )
  
  ### Fastest RTs
  
  # Extract density for each observation
  plt = densityPoints( rt$Fast )
  
  # Create blank plot
  xl = lowerUpper( .5, sort( rt$Fast ) );
  yl = lowerUpper( .5, plt$y )
  yl[1] = 0;
  blankPlot( xl, yl )
  
  segments( .2, 0, .2, yl[2], lwd = lnSz, col = 'grey' )
  
  # Draw variable-width histogram
  lines( plt$x, plt$y, type = 's', lwd = lnSz )
  segments( plt$x, rep( 0, N ),
            plt$x, plt$y, lwd = lnSz )
  
  # Add axes and labels
  customAxes( xl, yl )
  axis( 1, round( seq( xl[1], xl[2], length = 5 ), 2 ),
        tick = F, line = -1.5, cex.axis = txtSz )
  mtext( 'Time (s)', side = 1, line = 1.5, cex = txtSz )
  mtext( 'Density', side = 2, line = 1.5, cex = txtSz )
  
  title( 'Distribution of Fastest RTs', cex = txtSz )
  
  ### Slowest RTs
  
  # Extract density for each observation
  plt = densityPoints( rt$Slow )
  
  # Create blank plot
  xl = lowerUpper( .5, sort( rt$Slow ) );
  yl = lowerUpper( .5, plt$y )
  yl[1] = 0;
  blankPlot( xl, yl )
  
  # Draw variable-width histogram
  lines( plt$x, plt$y, type = 's', lwd = lnSz )
  segments( plt$x, rep( 0, N ),
            plt$x, plt$y, lwd = lnSz )
  
  # Add axes and labels
  customAxes( xl, yl )
  axis( 1, round( seq( xl[1], xl[2], length = 5 ), 2 ),
        tick = F, line = -1.5, cex.axis = txtSz )
  mtext( 'Time (s)', side = 1, line = 1.5, cex = txtSz )
  mtext( 'Density', side = 2, line = 1.5, cex = txtSz )
  
  title( 'Distribution of Slowest RTs', cex = txtSz )
  
  # Reset layout
  layout( cbind(1) )
}

###
### Accuracy predicts fastest/slowest RTs
###
# Lookup - 04

if ( runCode[3] ) {
  
  # Consider main trials only
  sel = rd$ BT == 2
  cd = rd[ sel, ]
  
  # Fastest and slowest RT and accuracy by subject
  ds = aggregate( cd$RT, list( cd$S ), function(x) 
    return( c( min(x), max(x) ) ) )
  ds = data.frame( S = ds[,1], Fast = ds$x[,1], Slow = ds$x[,2] )
  ds$P = aggregate( cd$Ac, list( cd$S ), mean )$x
  
  # Standardize observations
  ds$zF = NA; ds$zS = NA; ds$zP = NA
  for ( i in 1:3 ) {
    ds[,i+4] = ( ds[,i+1] - mean( ds[,i+1]) )/sd( ds[,i+1] )
  }
  
  # Plotting characteristics
  txtSz = 1.4
  lnSz = 2
  ptSz = 1
  
  if ( !savePlot ) x11( width = 12 )
  layout( cbind( 1, 2 ) )
  
  ### Fastest RT predicted by accuracy
  
  # Create blank plot
  xl = c( -3, 3 )
  yl = lowerUpper( .5, ds$Fast );
  yl[1] = -.01
  blankPlot( xl, yl )
  
  # Scatterplot
  points( ds$zP, ds$Fast, pch = 19, cex = ptSz )
  
  # Regression on log-RT
  mF = lm( I( log( Fast ) ) ~ zP, data = ds )
  # Include correction of bias due to transformation
  bias_correction = sum( exp( residuals( mF ) ) )/N
  # Plot model predictions
  prd = list( x = 
                seq( min( ds$zP ), max( ds$zP ), length = 1000 ) )
  prd$y = exp( coef( mF )[1] ) * exp( coef( mF )[2] * prd$x ) * 
    bias_correction
  lines( prd$x, prd$y, lty = 2, lwd = lnSz )
  
  # Axes and labels
  customAxes( xl, yl )
  axis( 2, round( seq( 0, yl[2], length = 5 ), 2 ),
        tick = F, line = -1.7, cex.axis = txtSz )
  mtext( 'Time (s)', side = 2, line = 2, cex = txtSz )
  
  xa = seq( -2.5, 2.5, length = 5 )
  axis( 1, xa, round( xa * sd( ds$P ) + mean( ds$P ), 2 ),
        tick = F, line = -1.7, cex.axis = txtSz )
  mtext( 'Accuracy', side = 1, line = 2, cex = txtSz )
  
  # Slope of for regression on log-RT
  text( -2.2, yl[2],
       substitute(
         paste( beta, " = ", val ),
         list(val = round( coef(mF)[2], 2 ) )
       ), cex = txtSz )
  
  title( 'Fastest RTs by Accuracy', cex = txtSz )
  
  ### Slowest RT predicted by accuracy
  
  # Create blank plot
  xl = c( -3, 3 )
  yl = lowerUpper( .5, ds$Slow );
  blankPlot( xl, yl )
  
  # Scatterplot
  points( ds$zP, ds$Slow, pch = 19, cex = ptSz )
  
  # Regression on log-RT
  mF = lm( I( log( Slow ) ) ~ zP, data = ds )
  # Include correction of bias due to transformation
  bias_correction = sum( exp( residuals( mF ) ) )/N
  # Plot model predictions
  prd = list( x = 
                seq( min( ds$zP ), max( ds$zP ), length = 1000 ) )
  prd$y = exp( coef( mF )[1] ) * exp( coef( mF )[2] * prd$x ) * 
    bias_correction
  lines( prd$x, prd$y, lty = 2, lwd = lnSz )
  
  # Axes and labels
  customAxes( xl, yl )
  axis( 2, round( seq( yl[1], yl[2], length = 5 ), 2 ),
        tick = F, line = -1.7, cex.axis = txtSz )
  mtext( 'Time (s)', side = 2, line = 2, cex = txtSz )
  
  xa = seq( -2.5, 2.5, length = 5 )
  axis( 1, xa, round( xa * sd( ds$P ) + mean( ds$P ), 2 ),
        tick = F, line = -1.7, cex.axis = txtSz )
  mtext( 'Accuracy', side = 1, line = 2, cex = txtSz )
  
  # Slope of for regression on log-RT
  text( 2.2, yl[2],
        substitute(
          paste( beta, " = ", val ),
          list(val = round( coef(mF)[2], 2 ) )
        ), cex = txtSz )
  
  title( 'Slowest RTs by Accuracy', cex = txtSz )

  # Reset layout
  layout( cbind( 1 ) ) 
}

###
### Conditional accuracy functions 
###
# Lookup - 05

if ( runCode[4] ) {
  
  # Initialize progress bar
  pb = txtProgressBar( min = 1, max = N, style = 3 )
  
  # Determine number of plotting windows
  nWindow = ceiling( N / 8 )
  
  inc = 1
  for ( j in 1:nWindow ) {
    
    if ( !savePlot ) x11( width = 12 )
    layout( matrix( 1:8, 2, 4, byrow = T ) )
    
    st = 1 + 8 * ( j - 1 )
    en = min( 8 + 8 * ( j - 1 ), N )
    
    for ( i in st:en ) {
      
      par( mar = c( 4, 4, 1, .5 ) )
      plt = quickCAF( i )
      title( paste( 'Subject', i ) )
      
      mtext( 'RT (s)', side = 1, outer = T, line = -2, cex = 1.4 )
      mtext( 'P(Correct)', side = 2, outer = T, line = -2, cex = 1.4 )
      
      # Update the progress bar
      setTxtProgressBar(pb,inc)
      inc = inc + 1
      
    }
    
  }
  close( pb )
  
  # Reset margins and layout
  layout( cbind( 1 ) ); par( mar = c( 4, 5, 3, 1 ) )
}

###
###
###
#

if ( saveFile ) {
  
  curDat = rawDat
  
  # Save data set
  setwd( orig_dir ); setwd('Data')
  fName = paste( 'Gabor_Exp_', Exp, '_11-17.RData', sep = '' )
  
  save ( curDat, rawDat, N, demographics, experiment_log, 
         file = fName )
  setwd( orig_dir )
}

# Close connection to PDf file
if ( savePlot ) dev.off()
setwd( orig_dir )
