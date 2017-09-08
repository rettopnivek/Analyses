#--------------------#
# Data preprocessing #
# Kevin Potter       #
# Updated 05/13/17   #
#--------------------#

# Clear workspace
rm(list = ls())

# Save current working directory
orig_dir = getwd()

# Index
# Lookup - 01: Load in useful packages
# Lookup - 02: Load in data and set script options
# Lookup - 03: Density plot of minimum and maximum RTs
# Lookup - 04: Remove excessively fast responses
# Lookup - 05: Remove excessively slow responses
# Lookup - 06: Trim data using mixture model
# Lookup - 07: Present results on trimmed data
# Lookup - 08: Save trimmed data

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

###
### Load in data and set script options
###
# Lookup - 02

# Load in data
setwd( 'Data' )
load( 'SD_v_FC.RData' )
setwd( orig_dir )

# Cut-off rules:
# Trim all subjects who failed to meet Dave's initial 
#   inclusion rules (returned for 2 sessions, overall accuracy 
#   was between 60 and 80%, average response time was less than 
#   1 second, and had no computer errors).
# RTs below 200 ms are too fast
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
  T,
  T
)

allDat = rawDat # Initialize data set to be cleaned

# Remove subjects who failed Dave's initial inclusion rules
allDat = allDat[ allDat$excl == '00000', ]
# Adjust number of subjects
N = length( unique( allDat$Subject ) )
# Create variable to keep track of original number of trials
nt_old = nrow( allDat )

###
### Density plot of minimum and maximum RTs
###
# Lookup - 03

if (runCode[1]) {
  
  too_fast = .2
  if (!savePlot) x11(width=12);
  layout( cbind(1,2) )
  ttl = c('Fastest RTs by subject',
          'Fastest RTs after trimming')
  for (i in 1:2) {
    
    if ( i == 1 ) d = allDat else d = d[ d$RT > too_fast, ]
    
    minRT = aggregate( d$RT, list( d$Subject ), min )
    colnames( minRT ) = c( 'S', 'RT' )
    
    ed = density( minRT$RT )
    af = approxfun( ed )
    ya = af( minRT$RT )
    yl = lowerUpper( .2, ya )
    xl = lowerUpper( .2, minRT$RT )
    
    plot( xl, c(0,yl[2]), type = 'n', xlab = 'Fastest RT (s)',
          ylab = 'Density', bty = 'l',
          main = ttl[i] )
    ord = order( minRT$RT )
    lines( minRT$RT[ ord ], ya[ ord ] )
    points( minRT$RT, ya, pch = 19 )
    
  }
  
  too_slow = 4
  if (!savePlot) x11(width=12);
  layout( cbind(1,2) )
  ttl = c('Slowest RTs by subject',
          'Slowest RTs after trimming')
  for (i in 1:2) {
    
    if ( i == 1 ) d = allDat else d = d[ d$RT < too_slow, ]
    
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
### Remove excessively fast responses
###
# Lookup - 04

if (runCode[2]) {
  
  # Somewhat arbitrarily, we'll define fast responses to 
  # be equal to or under 200 ms (Though we definitely 
  # want a cut-off slower than the N170 marker for 
  # word processing ).
  too_fast = .2
  # By subject
  S_too_fast = aggregate( allDat$RT <= too_fast,
                          list( allDat$Subject ), 
                          sum )$x
  # Total number of fast responses
  N_too_fast = sum( allDat$RT <= too_fast )
  keep = allDat$RT > too_fast
  allDat = allDat[ keep, ]
  
}

###
### Remove excessively slow responses
###
# Lookup - 05

if (runCode[3]) {
  
  # Somewhat arbitrarily, we'll define slow responses to 
  # be equal to or over 4 s.
  too_slow = 4
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
# Lookup - 06

if (runCode[4]) {
  
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

if ( sum( runCode[2:4] ) == 3 ) {
  
  Total_trimmed = 100*( 1 - nrow(newDat)/nt_old )
  
  Totals = c( 100*N_too_fast/nt_old,
              100*N_too_slow/nt_old,
              100*N_unlikely/nt_old )
  trimmed_S = cbind( S_too_fast/N_too_fast,
                     S_too_slow/N_too_slow,
                     S_unlikely/N_unlikely )
  
  if (!savePlot) x11();
  
  plot( c(0,.6), c(1,N+1), type = 'n',
        bty = 'n', xaxt = 'n', yaxt = 'n',
        xlab = ' ',ylab = ' ' )
  axis( 2, 1:N-.5, 1:N, tick = F, cex.axis = .5,
        line = -.5 )
  axis( 1, c(.1,.3,.5),
        c( 'Too fast', 'Too slow', 'Unlikely' ),
        tick = F, line = -.5,
        cex.axis = 1.2 )
  axis( 3, c(.1,.3,.5),
        paste( round( Totals, 2 ), '%', sep = '' ),
        tick = F, line = -2,
        cex.axis = 1.2 )
  mtext('Subject',side=2,cex=1.5,line=2)
  mtext( paste( 'Breakdown of ', round( Total_trimmed, 2 ),
                '% data trimmed', sep = '' ),
         side = 3, cex = 1.5, line = 1 )
  
  pst = c(0,.2,.4)
  for (j in 1:3) {
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

if ( runCode[4] ) allDat = newDat;

# Save original and new dataset
if (saveOutput) {
  setwd( 'Data' )
  save( rawDat, allDat, N, file = 'SD_v_FC.RData' )
}