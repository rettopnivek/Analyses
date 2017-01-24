#---------------------#
# Plots of model fits #
# Kevin Potter        #
# Updated 01/16/2017  #
#---------------------#

# Clear workspace
rm( list = ls() )

# Save current directory
orig_dir = getwd()

# Indicate model fit to plot
model_num = 3

# Indicate whether to save a pdf file
savePlot = T
if (savePlot) {
  setwd( 'Plots' )
  fName = paste( 'SDT_model_M', model_num, '.pdf', sep = '' )
  pdf( fName, width = 12, height = 6 )
  rm( fName )
  setwd( orig_dir )
}

# Index
# Lookup - 01:  Load in useful packages and data
# Lookup - 02:  Define useful functions
# Lookup - 03:  Check convergence
# Lookup - 04:  Posterior retrodictive checks
# Lookup - 05:  Plots of marginal posterior distributions

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

# Load in package for Bayesian estimation
library(rstan)
# For parallel processing
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# Load in data
setwd( 'Data' )
load( 'Original_all_data.RData' )
setwd( orig_dir )

# Extract observations for final recogniton memory 
# test and pre-process data

d = OriginalAllData[ OriginalAllData$Cond == 6, ]
# For easy manipulation
colnames( d ) = c('S','Tr','C','IN','Co','Ch','RT',
                  'Ac','IT','B','Bl','Cat','CR','fN')
d$I = createIncrement( d$IN )

# Use dummy coding for responses
d$Ch = d$Ch - 1; d$Co = d$Co - 1

# Define a variable denoting conditions

# Targets that underwent selective retrieval (T-SR)
d$Cnd = 1;
# Targets in the baseline condition (T-B)
d$Cnd[ d$IT == 1 & d$B == 1 ] = 2
# Competitors that underwent selective retrieval (C-SR)
d$Cnd[ d$IT == 2 & d$B == 0 ] = 3
# Competitors in the baseline condition (C-B)
d$Cnd[ d$IT == 2 & d$B == 1 ] = 4

# Set missing data to be incorrect
d$Ac[ is.na(d$RT) ] = 0
d$Ch[ is.na(d$RT) ] = 1 - d$Co[ is.na(d$RT) ]

###
### Define useful functions
###
# Lookup - 02

draw_figure = function( prp, lnSz = 1, lnTp = 1, ptSz = 1, 
                        pts = c(19,22,19,22),
                        clr = c('black','white'),
                        xa = NULL ) {
  # Purpose:
  # Adds a set of lines giving the proportion of times the 
  # response on the right was picked.
  # Arguments:
  # prp  - A set of 8 proportions, sorted by the correct position 
  #        of the target, the type of image (target versus competitor),
  #        and the condition (selective retrieval versus baseline)
  # lnSz - The width of the line
  # lnTp - The type of line
  # ptSz - The size of the point
  # pts  - A vector of 4 point types
  # clr  - The foreground and background color for the points
  # xa   - An optional list giving the x-axis positions for each 
  #        condition
  
  if ( length( xa ) == 0 ) 
    xa = list(
      c(.9,2.1),
      c(1.1,1.9),
      c(.9,2.1),
      c(1.1,1.9) )
  
  for ( i in 1:4 ) {
    
    sel = 1:2 + 2*(i-1)
    lines( xa[[i]], prp[sel], lwd = lnSz, lty = lnTp, col = clr[1] )
    points( xa[[i]], prp[sel], pch = pts[i], col = clr[1],
            bg = clr[2], cex = ptSz )
    
  }
  
}

convergence_extract = function( fit, par_name = NULL ) {
  # Purpose:
  # Extract convergence diagnostics from a Stan fit object.
  # Arguments:
  # fit      - A Stan fit object
  # par_name - An optional string giving the final parameter label 
  #            of the subset of the output to include
  # Notes:
  # Extracting the convergence statistics can be slow, especially 
  # when a large number of parameters were stored.
  # Returns:
  # A list with the Gelman-Rubin convergence statistic, the 
  # effective number of samples, and the total number of samples.
  
  Rhat = summary(fit)$summary[,"Rhat"]
  n_eff = summary(fit)$summary[,"n_eff"]
  totSampSize = length(extract(fit, pars = "lp__")[[1]])
  # We're only interested in a subset of parameters
  if ( length( par_name ) == 0 ) 
    par_name = names( Rhat )[ 
      which( names( Rhat ) == "logLik[1]" ) - 1 ]
  sel = 1:which( names(Rhat) == par_name )
  Rhat = Rhat[sel]; n_eff = n_eff[sel];
  
  return( list( Rhat = Rhat, n_eff = n_eff, totSampSize = totSampSize ) )
}

find_dec = function( x, spacing = 10 ) {
  # Purpose:
  # Determines the rounded leading digit and the 
  # number of trailing zeros for a number or the 
  # number of decimal places.
  # Arguments:
  # x       - A vector of values
  # spacing - The value whose exponent should be increased
  # Returns:
  # A vector giving the leading digit, the number of 
  # trailing/leading zeros, the same but in scientific 
  # notation, and 1 if it's trailing zeros, -1 if it's 
  # decimal places.
  
  mx = max( x )
  rnd = mx
  
  if ( round( mx ) > 1 ) {
    
    inc = 0;
    
    while ( rnd > 1 ) {
      inc = inc + 1;
      rnd = round( mx/( spacing^inc ) )
    }
    
    v = round( mx/spacing^(inc-1) )
    f = spacing^(inc-1)
    out = c( v,f,inc-1,1)
  }
  
  if ( round( mx ) == 1 ) {
    
    out = c( 1, 1, 1, 0 )
    
  }
  
  if ( round( mx ) == 0 ) {
    
    inc = 0;
    
    while ( rnd < 1 ) {
      inc = inc + 1;
      rnd = round( mx*( spacing^inc ) )
    }
    
    v = round( mx*spacing^(inc) )
    f = spacing^(inc)
    out = c( v,f,inc,-1)
    
  }
  
  return( out )
}

###
### Check convergence
###
# Lookup - 03

# Load in posterior
fName = "C:/Users/Kevin/Documents/Posteriors from Stan"
setwd(fName)
setwd( "Wimber_et_al" )
fname = paste( 'SDT_model_M', model_num, '_post.RData', sep = '' )
load( fname )
setwd( orig_dir )

conv = convergence_extract(fit,"sigma_zeta")
# Remove NaN for correlation matrix
sel = which( is.na( conv$Rhat ) )
conv$Rhat = conv$Rhat[ -sel ]
conv$n_eff = conv$n_eff[ -sel ]

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

###
### Posterior retrodictive checks
###
# Lookup - 04

# Extract data
x = aggregate( d$Ch, list( d$B, d$IT, d$Co, d$S ), mean )
colnames( x ) = c('B','IT','Co','S','P')

# Group-level performance

# Generate retrodictive checks for group-level performance
f = function(s) {
  yhat = rbinom( stanDat$No, 1, post$theta[s,] )
  out = aggregate( yhat, list( d$Cnd, d$Co ), mean )$x
  
  return( out )
}
rc = sapply( 1:nrow(post$beta), f )

# Create a plot of the proportions
if (!savePlot) x11(width=12)
layout( cbind(1,2) )
plot( c(.8,2.2), c(0,1), type = 'n', xlab = 'Condition',
      ylab = 'P(Right)', xaxt = 'n', bty = 'n' )
axis( 1, 1:2, c('Selective retrieval','Baseline'), tick = F )
abline(h=.5)
legend('topright','Right correct',bty='n')
legend('bottomright','Left correct',bty='n')
legend( 'topleft', c('Target','Competitor'),
        pch = c(19,22), pt.bg = 'white', bty = 'n' )

ya = seq( .1, .9, .01 )
points( rep(1,length(ya)), ya, pch = '.', cex = .5 )
points( rep(2,length(ya)), ya, pch = '.', cex = .5 )

pred = apply( rc, 1, findMode )

ui = apply( rc, 1, quantile, prob = c( .025,.16,.25,.75,.84,.975 ) )
ui = ui[ c(1:3,6:4), ]

xa = c(
  c(.9,2.1),
  c(1.1,1.9),
  c(.9,2.1),
  c(1.1,1.9) )
clr = c( rgb( .5, .5, .5, .4 ),
         rgb( .5, .5, .5, .4 ),
         rgb( .5, .5, .5, .4 ) )
sz = c( .05, .04, .03 )

for (i in 1:8) {
  
  for (j in 1:3) {
    
    drawEllipse( sz[j], diff( ui[c(j,j+3),i] ),
                 Xc = xa[i], Yc = ui[j,i] + diff( ui[c(j,j+3),i] )/2,
                 border = NA, col = clr[j] )
    
  }
  
}
points( xa, pred, pch = 19, col = 'grey40' )

obs = aggregate( x$P, list( x$B, x$IT, x$Co ), 
                 function(x) c( mean(x), sd(x)/sqrt(length(x)) ) )
colnames( obs ) = c( 'B', 'IT', 'Co', 'DS' )
draw_figure( obs$DS[,1] )

mtext( 'Group-level performance', side = 3, outer = T, line = -2 )

# Subject-level performance

# Generate retrodictive checks for subject-level performance
indexS = NULL
f = function(s) {
  yhat = rbinom( stanDat$No, 1, post$theta[s,] )
  out = aggregate( yhat, list( d$Cnd, d$Co, 
                               d$S ), mean )
  if ( length( indexS ) == 0 ) indexS <<- out[,3]
  
  return( out$x )
}
rc = sapply( 1:nrow(post$beta), f )

if (!savePlot) x11(width=12)
layout( matrix( 1:12, 3, 4, byrow = T ) )

for ( s in 1:12 ) {
  blankPlot( c(.8,2.2), c(0,1) )
  title( paste('Subject', s ) )
  if ( s == 1 | s == 5 | s == 9 ) 
    axis( 2, seq(0,1,.25) )
  if ( s == 9 | s == 10 | s == 11 | s == 12 ) 
    axis( 1, 1:2, c('SR','B' ), tick = F )
  
  
  draw_figure( x$P[ x$S == s ] )
  draw_figure( apply( rc[ indexS == s, ], 1, findMode ),
               clr = c('blue','white'), lnTp = 2 )
}

if (!savePlot) x11(width=12)
layout( matrix( 1:12, 3, 4, byrow = T ) )

for ( s in 13:24 ) {
  blankPlot( c(.8,2.2), c(0,1) )
  title( paste('Subject', s ) )
  if ( s == 1 | s == 5 | s == 9 ) 
    axis( 2, seq(0,1,.25) )
  if ( s == 9 | s == 10 | s == 11 | s == 12 ) 
    axis( 1, 1:2, c('SR','B' ), tick = F )
  
  
  draw_figure( x$P[ x$S == s ] )
  draw_figure( apply( rc[ indexS == s, ], 1, findMode ),
               clr = c('blue','white'), lnTp = 2 )
}

###
### Plots of marginal posterior distributions
###
# Lookup - 05

### Model 1

if ( model_num == 1 ) {
  
  if (!savePlot) x11(width=12)
  layout( matrix( 1:8, 2, 4, byrow = T ) )
  
  ttl = c( "d' for T-SR", "d' for T-B", "d' for C-SR", "d' for C-B" )
  for ( i in 1:4 ) {
    tmp = hist( post$beta[,i], col = 'grey', border = 'white', 
                breaks = 40, freq = F, main = ttl[i], xlab = ' ' )
    ui = quantile( post$beta[,i], c(.025,.975) )
    dn = c( max( which( tmp$mids < ui[1] ) ),
            min( which( tmp$mids > ui[2] ) ) )
    dn = tmp$density[ dn ]
    segments( ui, c(0,0), ui, dn, col = 'blue' )
  }
  
  # Mean criterion value
  tmp = hist( post$bias, col = 'grey', 
              border = 'white', breaks = 40, freq = F, 
              main = 'Bias', xlab = ' ' )
  ui = quantile( post$bias, c(.025,.975) )
  dn = c( max( which( tmp$mids < ui[1] ) ),
          min( which( tmp$mids > ui[2] ) ) )
  dn = tmp$density[ dn ]
  segments( ui, c(0,0), ui, dn, col = 'blue' )
  
  # Difference scores for first associates (targets)
  dT = post$beta[,2] - post$beta[,1]
  tmp = hist( dT, col = 'grey', 
              border = 'white', breaks = 40, freq = F, 
              main = 'B - SR for targets', xlab = ' ' )
  ui = quantile( dT, c(.025,.975) )
  dn = c( max( which( tmp$mids < ui[1] ) ),
          min( which( tmp$mids > ui[2] ) ) )
  dn = tmp$density[ dn ]
  segments( ui, c(0,0), ui, dn, col = 'blue' )
  
  # Difference scores for second associates (competitors)
  dC = post$beta[,4] - post$beta[,3]
  tmp = hist( dC, col = 'grey', 
              border = 'white', breaks = 40, freq = F, 
              main = 'B - SR for competitors', xlab = ' ' )
  ui = quantile( dC, c(.025,.975) )
  dn = c( max( which( tmp$mids < ui[1] ) ),
          min( which( tmp$mids > ui[2] ) ) )
  dn = tmp$density[ dn ]
  segments( ui, c(0,0), ui, dn, col = 'blue' )
  
  # Difference of differences
  tmp = hist( dC - dT, col = 'grey', 
              border = 'white', breaks = 40, freq = F, 
              main = 'Difference of differences', xlab = ' ' )
  ui = quantile( dC - dT, c(.025,.975) )
  dn = c( max( which( tmp$mids < ui[1] ) ),
          min( which( tmp$mids > ui[2] ) ) )
  dn = tmp$density[ dn ]
  segments( ui, c(0,0), ui, dn, col = 'blue' )
  
}

### Model 2

if ( model_num == 2 ) {
  
  if (!savePlot) x11(width=12)
  layout( matrix( 1:4, 1, 4, byrow = T ) )
  
  ttl = c( "Intercept", "Image type", "Condition" )
  for ( i in 1:3 ) {
    tmp = hist( post$beta[,i], col = 'grey', border = 'white', 
                breaks = 40, freq = F, main = ttl[i], xlab = ' ' )
    ui = quantile( post$beta[,i], c(.025,.975) )
    dn = c( max( which( tmp$mids < ui[1] ) ),
            min( which( tmp$mids > ui[2] ) ) )
    dn = tmp$density[ dn ]
    segments( ui, c(0,0), ui, dn, col = 'blue' )
  }
  
  # Mean criterion value
  tmp = hist( post$bias, col = 'grey', 
              border = 'white', breaks = 40, freq = F, 
              main = 'Bias', xlab = ' ' )
  ui = quantile( post$bias, c(.025,.975) )
  dn = c( max( which( tmp$mids < ui[1] ) ),
          min( which( tmp$mids > ui[2] ) ) )
  dn = tmp$density[ dn ]
  segments( ui, c(0,0), ui, dn, col = 'blue' )
  
}

### Model 3

if ( model_num == 3 ) {
  
  if (!savePlot) x11(width=12)
  layout( matrix( 1:4, 1, 4, byrow = T ) )
  
  ttl = c( "Intercept", "SR-C" )
  for ( i in 1:2 ) {
    tmp = hist( post$beta[,i], col = 'grey', border = 'white', 
                breaks = 40, freq = F, main = ttl[i], xlab = ' ' )
    ui = quantile( post$beta[,i], c(.025,.975) )
    dn = c( max( which( tmp$mids < ui[1] ) ),
            min( which( tmp$mids > ui[2] ) ) )
    dn = tmp$density[ dn ]
    segments( ui, c(0,0), ui, dn, col = 'blue' )
  }
  
  # Mean criterion value
  tmp = hist( post$bias, col = 'grey', 
              border = 'white', breaks = 40, freq = F, 
              main = 'Bias', xlab = ' ' )
  ui = quantile( post$bias, c(.025,.975) )
  dn = c( max( which( tmp$mids < ui[1] ) ),
          min( which( tmp$mids > ui[2] ) ) )
  dn = tmp$density[ dn ]
  segments( ui, c(0,0), ui, dn, col = 'blue' )
  
  blankPlot()
  
}

if (savePlot) dev.off()