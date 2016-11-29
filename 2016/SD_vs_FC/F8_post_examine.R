#---------------------------#
# Examination of posteriors #
# Kevin Potter              #
# Updated 11/23/2016        #
#---------------------------#

# Define model type
Type = 1

# Load in estimation results for hierarchical model
source( 'F7_BE_MS.R' )

# Indicate which code segments to run
runCode = c( T, F, F )

# Indicate whether to save a pdf file
savePlot = F
if (savePlot) {
  setwd( 'Plots' )
  fName = paste('Model', Type, 'fit_results.pdf', sep = '_' )
  pdf( fName, width = 12, height = 6 )
  rm( fName )
  setwd( orig_dir )
}

###
### Convergence check
###
# Lookup - 01

if ( runCode[1] ) {
  
  if (!savePlot) x11(width=12)
  layout( cbind(1,2) )
  
  # Plot a histogram of the Gelman-Rubin statistics for the 
  # marginal posterior samples of the parameters
  
  tmp = hist( conv$Rhat, plot = F )
  scl = find_dec( tmp$density )
  if ( scl[4] == 1 ) scl = scl[1]*(scl[2]/10) else 
    scl = scl[1]/(scl[2])
  
  yl = lowerUpper( scl, tmp$density )
  yl[1] = 0
  
  xl = lowerUpper( .1, conv$Rhat )
  xl[2] = max( xl[2], 1.12 )
  xl[1] = min( xl[1], .98 )
  
  plot( xl, yl, type = 'n', cex.axis = 1.5, cex.lab = 1.5,
        xlab = expression( hat(R) ), ylab = 'Density',
        bty = 'l', main = 'Gelman-Rubin Statistic' )
  
  segments( tmp$mids, rep(0,length(tmp$mids)),
            tmp$mids, tmp$density, lwd = 3,
            col = 'grey' )
  abline( v = 1.1, lty = 2, lwd = 2 )
  
  # Plot a histogram of the effective sample size for the 
  # set of parameters
  
  tmp = hist( conv$n_eff, plot = F )
  scl = find_dec( tmp$density )
  if ( scl[4] == 1 ) scl = scl[1]*(scl[2]/10) else 
    scl = scl[1]/(scl[2])
  
  yl = lowerUpper( scl, tmp$density )
  yl[1] = 0
  
  xl=c(0,conv$totSampSize)
  plot( xl, yl, type = 'n', cex.axis = 1.5, cex.lab = 1.5,
        xlab = expression(N[sample]), ylab = 'Density',
        bty = 'l', main = 'Effective sample size' )
  
  segments( tmp$mids, rep(0,length(tmp$mids)),
            tmp$mids, tmp$density, lwd = 3,
            col = 'grey' )
}

###
### Retrodictive checks
###
# Lookup - 02

