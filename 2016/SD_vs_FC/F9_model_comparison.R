#--------------------#
# Model comparisons  #
# Kevin Potter       #
# Updated 12/01/2016 #
#--------------------#

# Clear workspace
rm( list = ls() )

# Save current directory
orig_dir = getwd()

# Package for calculating WAIC and LOO
# install.packages(loo)
library(loo)

# Miscellanous functions for modeling/plotting
# install_github("rettopnivek/utilityf")
library(utilityf)

# Indicate whether to save a pdf file
savePlot = T
if (savePlot) {
  setwd( 'Plots' )
  pdf( file = 'Model_comparison.pdf', width = 12, height = 6 )
  setwd( orig_dir )
}

###
### Calculate WAIC for each model
###

# Change working directory to where posterior samples are stored
folderName = "C:/Users/Kevin/Documents/Posteriors from Stan/SD_vs_FC"
setwd( folderName )

# Initialize a list for the log-likelihood matrices
allWAIC = list()

fileNames = dir()
for ( l in 1:length( fileNames ) ) {
  
  load( fileNames[l] )
  
  if ( l == 1 ) allWAIC$M1 = waic( post$log_lik )
  if ( l == 2 ) allWAIC$M2 = waic( post$log_lik )
  if ( l == 3 ) allWAIC$M3 = waic( post$log_lik )
  
  # Clean up workspace
  rm( conv, input, post )
}

###
### Calculate Akaike weights for ease of interpretation
###

Akaike_weights = c(
  allWAIC$M1$waic,
  allWAIC$M2$waic,
  allWAIC$M3$waic
)
Akaike_weights = Akaike_weights - min( Akaike_weights )
Akaike_weights = exp( -.5*Akaike_weights )
Akaike_weights = Akaike_weights/sum( Akaike_weights )

if (!savePlot) x11();
xl = c( .5, length(Akaike_weights)+.5 )
yl = c(0,1)
plot( xl, yl, type = 'n', xlab = ' ', ylab = 'Akaike weights',
      cex.axis = 1.5, cex.lab = 1.5, bty = 'l', xaxt = 'n' )
# abline( h=c(.25,.5,.75), col = 'grey80', lwd = 2 )
axis( 1, 1:3, c( expression( M[T] ), expression( M[D] ), expression( M[H] ) ),
      tick = F, cex.axis = 1.5, line = 1 )
points( 1:3, Akaike_weights, pch = 19, cex = 2 )

###
### Pair-wise comparisons of models
###

c_2_v_1 = compare( allWAIC$M1, allWAIC$M2 )
c_3_v_1 = compare( allWAIC$M1, allWAIC$M3 )
c_3_v_2 = compare( allWAIC$M2, allWAIC$M3 )

df = c( c_2_v_1[1], c_3_v_1[1], c_3_v_2[1] )
df_se = c( c_2_v_1[2], c_3_v_1[2], c_3_v_2[2] )

if (!savePlot) x11()
xl = c( .5, 3.5 )
yl = c( -20, 220 )

blankPlot( xDim = xl, yDim = yl )

axis( 3, 1:3, c( expression( M[D] ), expression( M[H] ), expression( M[H] ) ),
      tick = F, cex.axis = 1.5 )
axis( 1, 1:3, c( expression( M[T] ), expression( M[T] ), expression( M[D] ) ),
      tick = F, cex.axis = 1.5 )
axis( 2, seq( -20, 220, 40 ), cex = 1.25 )
mtext( 'Difference in predictive accuracy', side = 2, cex = 1.5,
       line = 2.5 )

abline( h = 0, lwd = 2, col = 'grey' )

segments( 1:3, df - 2*df_se, 1:3, df + 2*df_se, lwd = 2 )
points( 1:3, df, pch = 19, cex = 1.5 )

if (savePlot) {
  setwd( orig_dir )
  setwd( 'Plots' )
  dev.off()
}

setwd( orig_dir )
