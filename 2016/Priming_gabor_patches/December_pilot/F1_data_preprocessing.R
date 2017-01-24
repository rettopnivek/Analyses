#--------------------#
# Data preprocessing #
# Kevin Potter       #
# Updated 12/15/16   #
#--------------------#

# Clear workspace
rm(list = ls())

# Save current working directory
orig_dir = getwd()

# Index
# Lookup - 01: Load in useful packages
# Lookup - 02: Load in data and set script options
# Lookup - 03: Density plot of minimum and maximum RTs
# Lookup - 04: Remove excessively slow responses
# Lookup - 05: Trim data using mixture model
# Lookup - 06: Present results on trimmed data
# Lookup - 07: Save trimmed data

###
### Load in useful packages
###
# Lookup - 01

# For geting github packages
# install.packages(devtools)
# library(devtools)

# Miscellanous functions for modeling/plotting
# install_github("rettopnivek/utilityf")
library(utilityf)

# Functions for preprocessing response time and choice data
# install_github("rettopnivek/rtclean")
library(rtclean)

# Functions for plotting response time and choice data
# install_github("rettopnivek/rtplots")
library(rtplots)

###
### Load in data and set script options
###
# Lookup - 02

# Load in data
setwd( 'Data' )
load( 'Gabor_pilot_Dec.RData' )
setwd( orig_dir )

# Cut-off rules:
# RTs above 4 s are too slow
# Use a mixture model, a shifted inverse gaussian 
#   and a uniform distribution, to identify remaining 
#   outliers.

# Create PDF with report of results, etc.
savePlot = T
if (savePlot) {
  setwd( 'Plots' )
  pdf('Trimming_results.pdf',width=12,height=6)
}

# Save output
saveOutput = T

# Select which code segments to run
runCode = c(
  T,
  T,
  F,
  T
)

# Pre-process data only for the main study
rawDat$TrialIndex = 1:nrow(rawDat) # For later trimming
allDat = rawDat[ rawDat$BlockType == 2, ]

# Create variable to keep track of original number of trials
nt_old = nrow( allDat )

###
### Density plot of minimum and maximum RTs
###
# Lookup - 03

if (runCode[1]) {
  
  # For easy manipulation
  d = allDat
  
  if (!savePlot) x11(width=12);
  lyt = matrix( 1:16, 4, 4, byrow = T )
  lyt[ lyt == 16 ] = 15
  lyt = cbind( lyt, matrix(16,4,4) )
  layout( lyt )
  
  for ( n in 1:N ) {
    
    par( mar = c(3.5,3,.5,.5) )
    blankPlot( c(0,1.5), c(0,5) )
    abline(v=.5)
    
    if ( n == 13 | n == 14 | n == 11 | n == 12 ) {
      
      axis( 1, seq(0,1.5,.5) )
      mtext('RT (s)', side = 1, line = 2.5 )
      
    }
    
    rt = d$RT[ d$Subject == n ]
    ch = rep(1,length(rt))
    
    pdf_curve( rt, ch, lwd = 2 )
    
  }
  par(mar=c(0,0,0,0))
  blankPlot()
  legend('left','RT densities by subject',bty='n',cex=2)
  
  par( mar=c(4,5,3,1) )
  ttl = c('Fastest RTs by subject')
  
  minRT = aggregate( d$RT, list( d$Subject ), min )
  colnames( minRT ) = c( 'S', 'RT' )
  
  ed = density( minRT$RT )
  af = approxfun( ed )
  ya = af( minRT$RT )
  yl = lowerUpper( .2, ya )
  xl = lowerUpper( .2, minRT$RT )
  
  plot( xl, c(0,yl[2]), type = 'n', xlab = 'Fastest RT (s)',
        ylab = 'Density', bty = 'l',
        main = ttl )
  ord = order( minRT$RT )
  lines( minRT$RT[ ord ], ya[ ord ] )
  points( minRT$RT, ya, pch = 19 )
  
  too_slow = 2.5
  if (!savePlot) x11(width=12);
  layout( cbind(1,2) )
  ttl = c('Slowest RTs by subject',
          'Slowest RTs after trimming')
  
  for (i in 1:2) {
    
    if ( i == 1 ) {
      d = allDat
    } else {
      sel = allDat$RT < too_slow
      d = allDat[ sel, ]
    }
    
    maxRT = aggregate( d$RT, list( d$Subject ), max )
    colnames( maxRT ) = c( 'S', 'RT' )
    
    ed = density( maxRT$RT )
    af = approxfun( ed )
    ya = af( maxRT$RT )
    yl = lowerUpper( .2, ya )
    xl = lowerUpper( .2, maxRT$RT )
    
    plot( xl, c(0,yl[2]), type = 'n', xlab = 'Slowest RT (s)',
          ylab = 'Density', bty = 'l',
          main = ttl[i] )
    ord = order( maxRT$RT )
    lines( maxRT$RT[ ord ], ya[ ord ] )
    points( maxRT$RT, ya, pch = 19 )
    
  }
  
}

###
### Remove excessively slow responses
###
# Lookup - 04

if (runCode[2]) {
  
  # Somewhat arbitrarily, we'll define slow responses to 
  # be equal to or over 2.5 s.
  too_slow = 2.5
  # By subject
  S_too_slow = aggregate( allDat$RT >= too_slow,
                          list( allDat$Subject ), 
                          sum )$x
  # Total number of fast responses
  N_too_slow = sum( allDat$RT >= too_slow )
  keep = allDat$RT < too_slow
  allDat = allDat[ keep, ]
  
}

###
### Trim data using mixture model
###
# Lookup - 05

if (runCode[3]) {
  
  newDat = c() # Initialize empty data set
  unlikely = numeric( nrow( allDat ) )
  
  # Create a progress bar using a base R function
  pb = txtProgressBar( min = 1, max = N, style = 3 )
  
  # Loop over subjects
  for (n in 1:N) {
    
    # Extract data for current subject
    sel = allDat$Subject == n
    cD = allDat[ sel, ]
    
    # Fit a mixture of shifted inverse gaussian and 
    # a uniform distribution. Observations that have 
    # a relative probability of less than .5 under the 
    # inverse gaussian distribution are trimmed.
    keep = rep( F, nrow(cD) )
    inc = 1
    while( sum(keep) == 0 & inc <= 20 ) {
      tst = rtclean( cD$RT, exclude = 1, nRep = 20 )
      keep = !tst$exclude_1
      inc = inc + 1
    }
    
    # Keep a record of trimmed responses
    unlikely[ sel ] = !keep
    
    if (sum(keep) == 0) break();
    
    # Create new, trimmed data set
    newDat = rbind( newDat, 
                    cD[ keep, ] )
    
    # Update the progress bar
    setTxtProgressBar(pb,n)
  }
  close(pb)
  
  # For record keeping
  N_unlikely = sum( unlikely )
  S_unlikely = aggregate( unlikely,
                          list( allDat$Subject ),
                          sum )$x
}

###
### Present results on trimmed data
###
# Lookup - 07

if ( sum( runCode[2:3] ) == 2 ) {
  
  Total_trimmed = 100*( 1 - nrow(newDat)/nt_old )
  
  Totals = c( 100*N_too_slow/nt_old,
              100*N_unlikely/nt_old )
  trimmed_S = cbind( S_too_slow/N_too_slow,
                     S_unlikely/N_unlikely )
  
  if (!savePlot) x11();
  
  plot( c(0,.6), c(1,N+1), type = 'n',
        bty = 'n', xaxt = 'n', yaxt = 'n',
        xlab = ' ',ylab = ' ' )
  axis( 2, 1:N-.5, 1:N, tick = F, cex.axis = .5,
        line = -.5 )
  axis( 1, c(.1,.3),
        c( 'Too slow', 'Unlikely' ),
        tick = F, line = -.5,
        cex.axis = 1.2 )
  axis( 3, c(.1,.3),
        paste( round( Totals, 2 ), '%', sep = '' ),
        tick = F, line = -2,
        cex.axis = 1.2 )
  mtext('Subject',side=2,cex=1.5,line=2)
  mtext( paste( 'Breakdown of ', round( Total_trimmed, 2 ),
                '% data trimmed', sep = '' ),
         side = 3, cex = 1.5, line = 1 )
  
  pst = c(0,.3)
  for (j in 1:2) {
    for (n in 1:N) {
      wght = trimmed_S[n,j]/max( trimmed_S[,j] )
      if (j==1) clr = rgb( 1, 0, 0, wght )
      if (j==2) clr = rgb( 0, 1, 0, wght )
      if (j==3) clr = rgb( 0, 0, 1, wght )
      polygon( c(0,.2,.2,0) + pst[j],
               c(n-1,n-1,n,n),
               col = clr )
    }
  }
  
  if (savePlot) blankPlot()
  
  # Look at distribution of minimum and maximum
  # RTs over subjects after trimming
  
  if (!savePlot) x11(width=12);
  layout( cbind(1,2) )
  
  d = newDat
  
  ttl = c('Distribution of fastest RTs',
          'Distribution of slowest RTs' )
  
  minRT = aggregate( d$RT, list( d$Subject ), min )
  colnames( minRT ) = c( 'S', 'RT' )
  
  ed = density( minRT$RT )
  af = approxfun( ed )
  ya = af( minRT$RT )
  yl = lowerUpper( .2, ya )
  xl = lowerUpper( .2, minRT$RT )
  
  plot( xl, c(0,yl[2]), type = 'n', xlab = 'Fastest RT (s)',
        ylab = 'Density', bty = 'l',
        main = ttl[1] )
  ord = order( minRT$RT )
  lines( minRT$RT[ ord ], ya[ ord ] )
  points( minRT$RT, ya, pch = 19 )
  
  maxRT = aggregate( d$RT, list( d$Subject ), max )
  colnames( maxRT ) = c( 'S', 'RT' )
  
  ed = density( maxRT$RT )
  af = approxfun( ed )
  ya = af( maxRT$RT )
  yl = lowerUpper( .2, ya )
  xl = lowerUpper( .2, maxRT$RT )
  
  plot( xl, c(0,yl[2]), type = 'n', xlab = 'Slowest RT (s)',
        ylab = 'Density', bty = 'l',
        main = ttl[2] )
  ord = order( maxRT$RT )
  lines( maxRT$RT[ ord ], ya[ ord ] )
  points( maxRT$RT, ya, pch = 19 )
  
}

if (savePlot) {
  dev.off()
  setwd( orig_dir )
}

###
### Save trimmed data
###
# Lookup - 08

if ( runCode[2] & !runCode[3] ) newDat = allDat
if ( runCode[4] ) {
  allDat = rawDat[ rawDat$BlockType != 2, ]
  allDat = rbind( allDat, newDat )
  allDat = allDat[ order( allDat$TrialIndex ), ]
  allDat = allDat[,-17]
}

# Save original and new dataset
if (saveOutput) {
  setwd( 'Data' )
  save( rawDat, allDat, N, file = 'Gabor_pilot_Dec.RData' )
}

setwd(orig_dir)