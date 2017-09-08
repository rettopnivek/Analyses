#-----------------------------------------#
# Checks and plots of parameter estimates #
# Kevin Potter                            #
# Updated 02/28/2017                      #
#-----------------------------------------#

# Clear workspace
rm(list = ls())

# Save current working directory
orig_dir = getwd()

# Indicate which code segments to run
runCode = c( T, T )

# Indicate whether a pdf of the figures should be saved
savePlot = T
if ( savePlot ) {
  setwd( 'Plots' )
  pdf( 'Model_results.pdf', width = 12 )
  setwd( orig_dir )
}

# Index
# Lookup - 01:  Load in useful packages, functions, and data
# Lookup - 02:  Parameter recovery (13 conditions)
# Lookup - 03:  Parameter recovery (14 conditions)

###
### Load in useful packages, functions, and data
###
# Lookup - 01

# For geting github packages
# install.packages(devtools)
# library(devtools)

# Miscellanous functions for modeling/plotting
# install_github("rettopnivek/utilityf")
library(utilityf)

# Load in useful functions
source( 'F1_Useful_functions.R' )

# Load in data
setwd( 'Data' )
load( 'DRMdata.RData' )
setwd( orig_dir )

###
### Results for SDT model of 13 observations
###
# Lookup - 02

if ( runCode[1] ) {
  
  # Load in parameter estimates
  setwd( 'Data' )
  parEst = read.csv( file = 'SDT_estimates_13.csv', header = T )
  setwd( orig_dir )
  
  ### Residuals ###
  
  # Initialize array for residuals
  resid = array( NA, dim = c( N, 13, 2 ) )
  
  # Initialize matrix for data
  curDat = matrix( NA, 13, 3 )
  colnames( curDat ) = c('N','Y','P')
  curDat = as.data.frame( curDat )
  
  # Loop over subjects
  for ( s in 1:N ) {
    
    # Convert 14 conditions into 13 conditions
    tmp = ratesDat[ ratesDat$S == s, ]
    curDat[2:13,] = tmp[3:14,c('N','Y','P')]
    curDat[1,] = c( sum( tmp$N[1:2] ),
                    sum( tmp$Y[1:2] ),
                    sum( tmp$Y[1:2] )/sum( tmp$N[1:2] ) )
    
    theta = SDT_4_dist_sim( as.numeric( parEst[s,1:13+3] ) )
    
    # Residuals using uncorrected data
    resid[ s, , 1 ] = curDat$P - theta
    
    # Residuals using corrected data
    tmp = (curDat$Y + .5)/(curDat$N+1)
    resid[ s, , 2 ] = tmp - theta
  }
  
  # Plot uncorrected residuals
  if (!savePlot) x11( width = 12 )
  plot( c( .5, 13.5 ), c( -.3, .3 ), type = 'n',
        yaxt = 'n', xaxt = 'n', ylab = 'Residuals',
        xlab = 'Condition', bty = 'l',
        main = 'Data with no corrections (13)' )
  axis( 1, 1:13, unique( ratesDat$Cond13 ) )
  axis( 2, seq( -.3, .3, .1 ),
        paste( seq( -.3, .3, .1 )*100, '%', sep = '' ) )
  for ( i in 1:13 ) {
    
    xa = seq( -.4, .4, length = N )
    points( xa + i, resid[,i,1], pch = 19, 
            col = rgb( .5, .5, .5, .5 ) )
    
  }
  
  # Plot corrected residuals
  if (!savePlot) x11( width = 12 )
  plot( c( .5, 13.5 ), c( -.3, .3 ), type = 'n',
        yaxt = 'n', xaxt = 'n', ylab = 'Residuals',
        xlab = 'Condition', bty = 'l',
        main = 'Data with corrections (13)' )
  axis( 1, 1:13, unique( ratesDat$Cond13 ) )
  axis( 2, seq( -.3, .3, .1 ),
        paste( seq( -.3, .3, .1 )*100, '%', sep = '' ) )
  for ( i in 1:13 ) {
    
    xa = seq( -.4, .4, length = N )
    points( xa + i, resid[,i,2], pch = 19, 
            col = rgb( .5, .5, .5, .5 ) )
    
  }
  
  ### Parameter distributions ###
  
  if (!savePlot) x11( width = 12 )
  
  # Determine plotting dimensions
  yl = lowerUpper( 1, as.vector( parEst[,-(1:3)] ) )
  
  plot( c( .5, 13.5 ), yl, type = 'n',
        xaxt = 'n', ylab = 'Parameter values',
        xlab = ' ', bty = 'l',
        main = 'Distributions for estimates (13)' )
  
  axis( 1, 1:13, colnames( parEst )[ -(1:3) ] )
  
  for ( i in 1:13 ) {
    
    p = parEst[ parEst$Age == 0, i + 3 ]
    fancy_violin_plot( p, i, clr = 'blue' )
    
    p = parEst[ parEst$Age == 1, i + 3 ]
    fancy_violin_plot( p, i, side = -1, clr = 'orange' )
    
  }
  
  legend( 'bottomleft', c('Young','Old'), 
          fill = c('blue','orange'), bty = 'n' )
  
}

###
### Residuals for SDT model of 14 observations
###
# Lookup - 03

if ( runCode[2] ) {
  
  # Load in data and parameter estimates
  setwd( 'Data' )
  parEst = read.csv( file = 'SDT_estimates_14.csv', header = T )
  setwd( orig_dir )
  
  ### Residuals ###
  
  # Initialize array for residuals
  resid = array( NA, dim = c( N, 14, 2 ) )
  
  # Loop over subjects
  for ( s in 1:N ) {
    
    curDat = ratesDat[ ratesDat$S == s, ]
    
    theta = SDT_4_dist_sim( as.numeric( parEst[s,1:14+3] ) )
    
    # Residuals using uncorrected data
    resid[ s, , 1 ] = curDat$P - theta
    
    # Residuals using corrected data
    tmp = (curDat$Y + .5)/(curDat$N+1)
    resid[ s, , 2 ] = tmp - theta
  }
  
  # Plot uncorrected residuals
  if (!savePlot) x11( width = 12 )
  plot( c( .5, 14.5 ), c( -.3, .3 ), type = 'n',
        yaxt = 'n', xaxt = 'n', ylab = 'Residuals',
        xlab = 'Condition', bty = 'l',
        main = 'Data with no corrections (14)' )
  axis( 1, 1:14, unique( ratesDat$Cond14 ) )
  axis( 2, seq( -.3, .3, .1 ),
        paste( seq( -.3, .3, .1 )*100, '%', sep = '' ) )
  for ( i in 1:14 ) {
    
    xa = seq( -.4, .4, length = N )
    points( xa + i, resid[,i,1], pch = 19, 
            col = rgb( .5, .5, .5, .5 ) )
    
  }
  
  # Plot corrected residuals
  if (!savePlot) x11( width = 12 )
  plot( c( .5, 14.5 ), c( -.3, .3 ), type = 'n',
        yaxt = 'n', xaxt = 'n', ylab = 'Residuals',
        xlab = 'Condition', bty = 'l',
        main = 'Data with corrections (14)' )
  axis( 1, 1:14, unique( ratesDat$Cond14 ) )
  axis( 2, seq( -.3, .3, .1 ),
        paste( seq( -.3, .3, .1 )*100, '%', sep = '' ) )
  for ( i in 1:14 ) {
    
    xa = seq( -.4, .4, length = N )
    points( xa + i, resid[,i,2], pch = 19, 
            col = rgb( .5, .5, .5, .5 ) )
    
  }
  
  ### Parameter distributions ###
  
  if (!savePlot) x11( width = 12 )
  
  # Determine plotting dimensions
  yl = lowerUpper( 1, as.vector( parEst[,-(1:3)] ) )
  
  plot( c( .5, 14.5 ), yl, type = 'n',
        xaxt = 'n', ylab = 'Parameter values',
        xlab = ' ', bty = 'l',
        main = 'Distributions for estimates (14)' )
  
  axis( 1, 1:14, colnames( parEst )[ -(1:3) ], cex.axis = .7 )
  
  for ( i in 1:14 ) {
    
    p = parEst[ parEst$Age == 0, i + 3 ]
    fancy_violin_plot( p, i, clr = 'blue' )
    
    p = parEst[ parEst$Age == 1, i + 3 ]
    fancy_violin_plot( p, i, side = -1, clr = 'orange' )
    
  }
  
  legend( 'bottomleft', c('Young','Old'), 
          fill = c('blue','orange'), bty = 'n' )
  
}

if ( savePlot ) dev.off()

