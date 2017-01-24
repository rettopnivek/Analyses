#--------------------------------#
# Correlations with nROUSE model #
# Kevin Potter                   #
# Updated 12/02/2016             #
#--------------------------------#

# Define model type (can only be 2 or 3)
type = 3

# Load in estimation results for hierarchical model
source( 'F7_BE_MS.R' )

# Load in Dave's original nROUSE data
setwd( 'Data' )
# Load estimates of nROUSE parameters/latencies from Dave
load( 'nROUSE_est_from_Dave.RData' )
setwd( orig_dir )

# Indicate whether to save a pdf file
savePlot = T
if (savePlot) {
  setwd( 'Plots' )
  pdf( file = 'New_priors.pdf', width = 12, height = 6 )
  setwd( orig_dir )
}

###
### Correlation between nROUSE latencies and drift rates
###
# Lookup - 01

# Restrict to forced-choice data
cD = d[ d$Ta == 0, ]

# Determine the aggregate standardized inverse identification latencies
mIL = aggregate( 1/cD$TL, list( cD$Cnd ), mean )
colnames( mIL ) = c( 'Ord', 'IL' )
mIL = mIL[ c(1,3,5,7), ]
tmp = aggregate( 1/cD$FL, list( cD$Cnd ), mean )
tmp = tmp[ c(1,3,5,7), ]
colnames( tmp ) = c( 'Ord', 'IL' )
mIL = rbind( mIL, tmp ); rm( tmp )
mIL$Ord = 1:8

# Define a function to carry out maximum likelihood estimation for simple regression
mle_f = function( prm, dat, priors = NULL ) {
  
  zx = dat[,1]
  zy = dat[,2]
  
  mu = prm[1]*zx
  
  ll = sum( dnorm( zy, mu, prm[2], log = T ) )
  if (is.na(ll)) ll = -Inf
  
  return( ll )
}

# Define a function to generate random starting values
st_f = function() return( runif( 2, c(-1,.2),c(1,2) ) )

# Define a matrix for the data
dat = matrix( NA, 8, 2 ); dat[,1] = scale( mIL$IL )

# Track progress
total = nrow( post$xi_mu )
pb = txtProgressBar( min = 1, max = total, style = 3 )
inc = 0

# Define a function to carry out estimation and extract the best estimates
f = function( y ) {
  
  inc <<- inc + 1
  dat[,2] <<- scale( y )
  out = MLE( dat, mle_f, st_f )$param
  
  # Update the progress bar
  setTxtProgressBar(pb,inc)
  
  return( out )
}

# Calculate how well the latencies predict the group-level drift rates
est_lm = apply( post$xi_mu, 1, f )
close(pb)

# Extract estimate of correlation
est_R = est_lm[1,]

### Plot resulting distributions ###

if (!savePlot) x11(width = 12)
layout( cbind( 1, 2 ) )

hist( est_R, breaks = 40, freq = F, col = 'grey', border = 'white',
      main = ' ', xlab = 'Correlation between nROUSE/drift rates', 
      cex.axis = 1.5, cex.lab = 1.5 )

hist( est_lm[2,], breaks = 40, freq = F, col = 'grey', border = 'white',
      main = ' ', xlab = 'Residual standard deviation', 
      cex.axis = 1.5, cex.lab = 1.5 )

###
### Generate new prior distributions
###

# Determine aggregate parameters
nR_p = colMeans( nROUSE_res[,1:3] )

# Specify new experimental conditions
newPrimeDur = c( 50, 2000, 50, 2000 )
newTarDur = round( median( nROUSE_prev_est$target_duration ) )
newTarDur = rep( newTarDur, 4 )
newType = c( 2, 2, -2, -2 )

# Generate predicted latencies
predLat = matrix(NA,4,2)
for (i in 1:4) {
  pres = c(newPrimeDur[i],newTarDur[i],
           500 - newTarDur, 500 )
  if ( newType[i] == -2 ) pi = c(0,2) else pi = c(2,0)
  prm = c(0.25, 0.0302, 0.15, 0.324, 0.022, 0.9844, 
          0.15, 1, 0.0294, 0.0609, 0.015)
  prm[ c(1,6,8) ] = nR_p
  sim = simulate_nROUSE( pres, pi, prm )
  predLat[i,] = sim$Latencies[1:2]
}
# Determine the inverse latencies
predLat = 1/predLat
predLat = c( predLat[,1], predLat[,2] )

# Determine the predicted posterior distributions for the means
predDist = matrix( NA, conv$totSampSize, 8 )
for ( i in 1:nrow( predDist ) ) {
  
  y = post$xi_mu[i,]
  
  x = predLat;
  zx = scale( x )
  m = est_lm[1,i];
  s = est_lm[2,i];
  
  zy = zx*m + rnorm(length(y),0,s) # Perturb by estimated error
  predDist[i,] = zy*sd(y) + mean(y)
}

est_mu = apply( predDist, 2, mean )
names( est_mu ) = as.character( 1:8 )
est_sig = apply( predDist, 2, sd )
names( est_sig ) = as.character( 1:8 )

print( 'New prior means for group-level drift rates' )
print( round( est_mu, 2 ) )

print( 'New prior SD for group-level drift rates' )
print( round( est_sig, 2 ) )

### Plot a comparison of the old posteriors and the new priors ###

if (!savePlot) x11(width=12)
layout( cbind(1) )
blankPlot( xDim=c(.5,8.5), yDim = c(0,5) )
axis( 1, 1:8, rep( c('STP','LTP','SFP','LFP'), 2 ), tick = F,
      cex = 1.2 )
axis( 1, c(2.5,6.5), c('Target','Foil') , tick = F,
      cex = 1.2, line = 2.5 )
axis( 2, seq(0,5,1), cex = 1.2, lwd = 2 )
mtext( 'Drift rates', side = 2, cex = 1.2, line = 2 )
abline( h = 0, lwd = 2 )
legend( 'topright', c('Old posteriors','New priors'),
        fill = c('grey','white'), bty='n', cex = 1.2 )

for ( i in 1:8 ) {
  
  violinPlot( post$xi_mu[,i], i, border = NA, col = 'grey',
              scaleH = .3 )
  
  v = seq( -4, 4, length = 1000 )
  v = v*est_sig[i] + est_mu[i]
  v = v[ v > 0 ]
  
  h = dnorm( v, est_mu[i], est_sig[i] )
  h = .3*h/max(h)
  lines( i - h, v, lwd = 2 )
  lines( i + h, v, lwd = 2 )
  
  text( i - .4, max( v[ which(h==max(h)) ] ),
        paste( '(', round( est_mu[i], 2 ),', ',
               round( est_sig[2], 2 ),')', sep = '' ),
        srt = 90, cex = 1.2 )
  
}

if (savePlot) {
  setwd( orig_dir )
  setwd( 'Plots' )
  dev.off()
}

setwd( orig_dir )
