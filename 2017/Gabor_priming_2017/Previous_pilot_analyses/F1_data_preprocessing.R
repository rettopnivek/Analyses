#--------------------#
# Data preprocessing #
# Kevin Potter       #
# Updated 09/29/2017 #
#--------------------#

# Clear workspace
rm(list = ls())

# Save current working directory
orig_dir = getwd()

# Indicator for whether to fit pilot data or current study
Pilot = T

if ( Pilot ) dName = 'Pilot' else dName = 'Current'

# Create PDF with report of results, etc.
savePlot = T
if (savePlot) {
  setwd( 'Plots' )
  pdf( paste( dName, 'trimming_results.pdf', sep = '_' ), width=12, height=6 )
  setwd( orig_dir )
}

# Save output
saveOutput = T

# Select which code segments to run
runCode = c(
  T,
  T,
  T,
  T,
  T,
  F
)


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
if ( Pilot ) load( 'Gabor_pilot_9-17.RData' ) else load( 'Gabor_current_9-17.RData' ) 
setwd( orig_dir )

# Pre-process data only for the main study
rawDat$TrialIndex = 1:nrow(rawDat) # For later trimming
allDat = rawDat[ rawDat$BlockType == 2, ]

# Create variable to keep track of original number of trials
nt_old = nrow( allDat )

# Cut-off rules:
# RTs below .2 s are too fast
# RTs above 4 s are too slow
# Remove subjects whose data consist of over 20% fast/slow responses
# Use a mixture model, a shifted inverse gaussian 
#   and a uniform distribution, to identify remaining 
#   outliers.

too_fast = .2; too_slow = 4;

###
### Density plot of minimum and maximum RTs
###
# Lookup - 03

if (runCode[1]) {
  
  # For easy manipulation
  d = allDat
  
  ### Conditional accuracy figure for fastest RTs

  new_win = F
  for ( n in 1:N ) {
    
    if ( n > 1 & ( (n-1) %% 6 != 0 ) ) new_win = F else new_win = T
    if ( new_win ) {
      # Obtain dimensions for figure layout
      if ( !savePlot ) x11( width = 12 )
      layout( matrix( 1:6, 2, 3, byrow = T ) )
    }
    
    # Extract data for each subject
    s = d$Subject == n
    rt = d$RT[ s ]; ac = d$Accuracy[ s ]
    
    # Determine RTs <= to 2.5% to 20% quantiles
    q = quantile( rt, prob = seq( .025, .2, .025 ) )
    p = numeric( length(q) ) # Determine accuracy at each point
    tot = numeric( length(q) ) # Determine number of trials at each point
    for ( i in 1:length(q) ) {
      sel = rt <= q[i]
      tot[i] = sum( sel )
      p[i] = sum( ac[ sel ] )/sum( sel )
    }
    # Compute 68% uncertainty interval around chance performance
    ui = rbind( qbinom( .16, tot, .5 )/tot,
                qbinom( .84, tot, .5 )/tot )
    
    # Create blank plot
    xl = c(-.1,.75); yl = c(0,1);
    blankPlot(xl,yl);
    
    # Add axes
    customAxes(xl,yl)
    axis( 1, seq( 0, .75, .25 ), tick = F, line = -1.25, cex.axis = 1 )
    axis( 2, seq(0,1,.25), tick = F, line = -1.25, cex.axis = 1)
    
    # Uncertainty interval
    lines( q, ui[1,], col = 'grey80', lwd = 2 )
    lines( q, ui[2,], col = 'grey80', lwd = 2 )
    
    # CAF
    lines( q, p, type = 'b', pch = 19, lwd = 2, cex = 1.2 )
    
    # Gridlines
    segments( min(rt), yl[1], min(rt), yl[2], lty = 2, lwd = 2 )
    segments( xl[1], .5, xl[2], .5, lty = 2, lwd = 2 )
    
    # Subject
    title( paste( 'Subject', n ), cex = .8 )
    
    # Legend
    if ( n == 1 ) legend( -.05, yl[2], 'At chance',
                          fill = 'grey80', bty = 'n' )
    
    if ( n == 1 | n %% 6 == 5 ) {
      mtext( 'RT (s)', side = 1, outer = T, line = -2 )
      mtext( 'Accuracy', side = 2, outer = T, line = -2 )
      mtext( 'Conditional accuracy functions', side = 3, outer = T, line = -1.5 )
    }
    
  }
  
  ### Overall accuracy / Proportion of fast responses
  
  # Overall accuracy
  Y = aggregate( d$Accuracy, list( d$Subject ), 
                 function(x) c( F = sum(x), N = length(x), P = mean(x) ) )
  colnames( Y ) = c( 'S', 'Y' )
  
  if ( !savePlot ) x11( width = 12 )
  layout( cbind( 1, 2 ) )
  par( mar=c(4,5,3,1) )
  
  xl = c( .5, N + .5 ); yl = c( 0, 1 )
  blankPlot( xl, yl ); customAxes( xl, yl )
  
  ui_50 = c( qbinom( .025, min( Y$Y[,'N'] ), .5 )/min( Y$Y[,'N'] ),
             qbinom( .975, min( Y$Y[,'N'] ), .5 )/min( Y$Y[,'N'] ) )
  ui_75 = c( qbinom( .025, min( Y$Y[,'N'] ), .75 )/min( Y$Y[,'N'] ),
             qbinom( .975, min( Y$Y[,'N'] ), .75 )/min( Y$Y[,'N'] ) )
  
  segments( c( xl[1], xl[1] ), ui_50[c(1,2)], 
            c( xl[2], xl[2] ), ui_50[c(1,2)], lwd = 2, lty = 2, col = 'grey80' )
  
  segments( c( xl[1], xl[1] ), ui_75[c(1,2)], 
            c( xl[2], xl[2] ), ui_75[c(1,2)], lwd = 2, lty = 3, col = 'grey80' )
  
  points( 1:N, Y$Y[,'P'], pch = 19, cex = 1.2 )
  
  axis( 1, 1:N, tick = F, line = -1.25 )
  axis( 2, seq( 0, 1, .25 ), tick = F, line = -1.25 )
  mtext( 'Subject', side = 1, line = 1.25 )
  mtext( 'Accuracy', side = 2, line = 1.25 )
  
  # Proportion of fast/slow responses
  Y = aggregate( d$RT, list( d$Subject ), 
                 function(x) 
                   c( F = sum( x < too_fast | x > too_slow ), 
                      P = sum( x < too_fast | x > too_slow )/ length(x) ) )
  colnames( Y ) = c( 'S', 'Y' )
  
  xl = c( .5, N + .5 ); yl = c( 0, 1 )
  blankPlot( xl, yl ); 
  segments( xl[1], .2, xl[2], .2, lwd = 2, lty = 2, col = 'grey80' )
  
  points( 1:N, Y$Y[,'P'], pch = 19, cex = 1.2 )
  
  customAxes( xl, yl )
  axis( 2, seq( 0, 1, .25 ), tick = F, line = -1.25 )
  mtext( 'Subject', side = 1, line = 1.25 )
  mtext( 'Proportion of fast/slow', side = 2, line = 1.25 )
  
  # Define function to create density plot of fastest/slowest RTs
  fast_slow_density = function( d, type ) {
    
    if ( type == 'Fastest' ) f = min
    if ( type == 'Slowest' ) f = max
    
    if (!savePlot) x11(width=12);
    layout( cbind(1,2) )
    par( mar=c(4,5,3,1) )
    ttl = c( paste( type, 'RTs by subject' ),
             paste( type, 'RTs after trimming' ) )
    
    for (i in 1:2) {
      
      if ( i == 1 ) {
        d = allDat
      } else {
        if ( type == 'Fastest' ) {
          sel = allDat$RT >= too_fast
        } else {
          sel = allDat$RT < too_slow
        }
        d = allDat[ sel, ]
      }
      
      val = aggregate( d$RT, list( d$Subject ), f )
      colnames( val ) = c( 'S', 'RT' )
      
      ed = density( val$RT )
      af = approxfun( ed )
      ya = af( val$RT )
      yl = lowerUpper( .2, ya )
      xl = lowerUpper( .2, val$RT )
      
      plot( xl, c(0,yl[2]), type = 'n', xlab = paste( type, 'RT (s)' ),
            ylab = 'Density', bty = 'l',
            main = ttl[i] )
      ord = order( val$RT )
      lines( val$RT[ ord ], ya[ ord ] )
      points( val$RT, ya, pch = 19 )
      
    }
    
  }
  
  ### Fastest RTs by subject
  fast_slow_density( d, 'Fastest' )
  
  ### Slowest RTs by subject
  fast_slow_density( d, 'Slowest' )
}

###
### Remove subjects with too many fast/slow responses
###
# Lookup - 04

if ( runCode[2] ) {
  
  # Compute proportion of fast/slow responses
  Y = aggregate( allDat$RT, list( allDat$Subject ), 
                 function(x) 
                   c( F = sum( x < too_fast | x > too_slow ), 
                      P = sum( x < too_fast | x > too_slow )/ length(x) ) )
  colnames( Y ) = c( 'S', 'Y' )
  
  
  S_too_many_fast_slow = Y$S[ which( Y$Y[,'P'] > .2 ) ]
  keep = !( allDat$Subject %in% S_too_many_fast_slow )
  
}

###
### Remove excessively fast responses
###
# Lookup - 05

if (runCode[3]) {
  
  # Somewhat arbitrarily, we'll define fast responses to 
  # be less than .2 s.
  # By subject
  S_too_fast = aggregate( allDat$RT < too_fast,
                          list( allDat$Subject ), 
                          sum )$x
  # Total number of fast responses
  N_too_fast = sum( allDat$RT < too_fast )
  keep = keep & allDat$RT >= too_fast
  
}

###
### Remove excessively slow responses
###
# Lookup - 06

if (runCode[4]) {
  
  # Somewhat arbitrarily, we'll define slow responses to 
  # be equal to or over 4 s.
  # By subject
  S_too_slow = aggregate( allDat$RT >= too_slow,
                          list( allDat$Subject ), 
                          sum )$x
  # Total number of fast responses
  N_too_slow = sum( allDat$RT >= too_slow )
  keep = keep & allDat$RT < too_slow
  
}

# Trim responses
if ( any( runCode[2:4] ) ) 
  allDat = allDat[ keep, ]

# Update N and subject indices
N = length( unique( allDat$Subject ) )
allDat$Subject = createIncrement( allDat$Subject )

###
### Trim data using mixture model
###
# Lookup - 07

if (runCode[5]) {
  
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

if ( all( runCode[2:5] ) ) {
  
  Total_trimmed = 100*( 1 - nrow(newDat)/nt_old )
  
  # Adjust totals based on subjects who were removed
  if ( length( S_too_many_fast_slow ) > 0 ) {
    N_too_fast = N_too_fast - sum( S_too_fast[ S_too_many_fast_slow ] )
    N_too_slow = N_too_slow - sum( S_too_slow[ S_too_many_fast_slow ] )
    S_too_fast = S_too_fast[ -S_too_many_fast_slow ]
    S_too_slow = S_too_slow[ -S_too_many_fast_slow ]
  }
  
  Totals = c(
    100 * N_too_fast/nt_old,
    100 * N_too_slow/nt_old,
    100 * N_unlikely/nt_old
  )
  
  trimmed_S = cbind( 
    S_too_fast/N_too_fast,
    S_too_slow/N_too_slow,
    S_unlikely/N_unlikely )
  
  if (!savePlot) x11( width = 12 );
  layout( cbind( 2, 1, 1, 2 ) )
  
  plot( c(0,.6), c(0,N+1), type = 'n',
        bty = 'n', xaxt = 'n', yaxt = 'n',
        xlab = ' ',ylab = ' ' )
  axis( 2, 1:N-.5, 1:N, tick = F, cex.axis = 1,
        line = -.5 )
  axis( 1, c(.1,.3,.5),
        c( 'Too fast', 'Too slow', 'Unlikely' ),
        tick = F, line = -.5,
        cex.axis = 1.2 )
  axis( 3, c(.1,.3,.5 ),
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
  
  mtext( paste( '% from trimmed subjects: ',
                round( 100 * ( length( S_too_many_fast_slow ) * 360 )/nt_old, 2 ),
                sep = ' ' ), side = 1, outer = T, line = -1.5 )
  
  if (savePlot) blankPlot()
  
  # Look at distribution of minimum and maximum
  # RTs over subjects after trimming
  
  if (!savePlot) x11(width=12);
  layout( cbind(1,2) )
  
  ### Conditional accuracy figure for fastest RTs
  
  d = newDat
  
  new_win = F
  for ( n in 1:N ) {
    
    if ( n > 1 & ( (n-1) %% 6 != 0 ) ) new_win = F else new_win = T
    if ( new_win ) {
      # Obtain dimensions for figure layout
      if ( !savePlot ) x11( width = 12 )
      layout( matrix( 1:6, 2, 3, byrow = T ) )
    }
    
    # Extract data for each subject
    s = d$Subject == n
    rt = d$RT[ s ]; ac = d$Accuracy[ s ]
    
    # Determine RTs <= to 2.5% to 20% quantiles
    q = quantile( rt, prob = seq( .025, .2, .025 ) )
    p = numeric( length(q) ) # Determine accuracy at each point
    tot = numeric( length(q) ) # Determine number of trials at each point
    for ( i in 1:length(q) ) {
      sel = rt <= q[i]
      tot[i] = sum( sel )
      p[i] = sum( ac[ sel ] )/sum( sel )
    }
    # Compute 68% uncertainty interval around chance performance
    ui = rbind( qbinom( .16, tot, .5 )/tot,
                qbinom( .84, tot, .5 )/tot )
    
    # Create blank plot
    xl = c(-.1,.75); yl = c(0,1);
    blankPlot(xl,yl);
    
    # Add axes
    customAxes(xl,yl)
    axis( 1, seq( 0, .75, .25 ), tick = F, line = -1.25, cex.axis = 1 )
    axis( 2, seq(0,1,.25), tick = F, line = -1.25, cex.axis = 1)
    
    # Uncertainty interval
    lines( q, ui[1,], col = 'grey80', lwd = 2 )
    lines( q, ui[2,], col = 'grey80', lwd = 2 )
    
    # CAF
    lines( q, p, type = 'b', pch = 19, lwd = 2, cex = 1.2 )
    
    # Gridlines
    segments( min(rt), yl[1], min(rt), yl[2], lty = 2, lwd = 2 )
    segments( xl[1], .5, xl[2], .5, lty = 2, lwd = 2 )
    
    # Subject
    title( paste( 'Subject', n ), cex = .8 )
    
    # Legend
    if ( n == 1 ) legend( -.05, yl[2], 'At chance',
                          fill = 'grey80', bty = 'n' )
    
    if ( n == 1 | n %% 6 == 5 ) {
      mtext( 'RT (s)', side = 1, outer = T, line = -2 )
      mtext( 'Accuracy', side = 2, outer = T, line = -2 )
      mtext( 'Conditional accuracy functions', side = 3, outer = T, line = -1.5 )
    }
    
  }
  
  # Define function to create density plot of fastest/slowest RTs
  fast_slow_density = function( d, type ) {
    
    if ( type == 'Fastest' ) f = min
    if ( type == 'Slowest' ) f = max
    
    par( mar=c(4,5,3,1) )
    ttl = c( paste( type, 'RTs by subject' ),
             paste( type, 'RTs after trimming' ) )
    
    for (i in 2:2) {
      
      if ( i == 1 ) {
        d = allDat
      } else {
        if ( type == 'Fastest' ) {
          sel = allDat$RT >= too_fast
        } else {
          sel = allDat$RT < too_slow
        }
        d = allDat[ sel, ]
      }
      
      val = aggregate( d$RT, list( d$Subject ), f )
      colnames( val ) = c( 'S', 'RT' )
      
      ed = density( val$RT )
      af = approxfun( ed )
      ya = af( val$RT )
      yl = lowerUpper( .2, ya )
      xl = lowerUpper( .2, val$RT )
      
      plot( xl, c(0,yl[2]), type = 'n', xlab = paste( type, 'RT (s)' ),
            ylab = 'Density', bty = 'l',
            main = ttl[i] )
      ord = order( val$RT )
      lines( val$RT[ ord ], ya[ ord ] )
      points( val$RT, ya, pch = 19 )
      
    }
    
  }
  
  if ( !savePlot ) x11( width = 12 )
  layout( cbind( 1, 2 ) )
  
  ### Fastest RTs by subject
  fast_slow_density( d, 'Fastest' )
  
  ### Slowest RTs by subject
  fast_slow_density( d, 'Slowest' )
  
}

if (savePlot) {
  dev.off()
  setwd( orig_dir )
}

###
### Save trimmed data
###
# Lookup - 08

if ( all( runCode[2:5] ) ) {
  allDat = rawDat[ rawDat$BlockType != 2, ]
  allDat = rbind( allDat, newDat )
  allDat = allDat[ order( allDat$TrialIndex ), ]
  sel = which( colnames( allDat ) == 'TrialIndex' )
  allDat = allDat[,-sel]
}

# Save original and new dataset
if (saveOutput) {
  setwd( 'Data' )
  if ( Pilot ) save( rawDat, allDat, N, demographics, experiment_log, 
                     file = 'Gabor_pilot_9-17.RData' )
  if ( !Pilot ) save( rawDat, allDat, N, demographics, experiment_log, 
                      file = 'Gabor_current_9-17.RData' )
}

setwd(orig_dir)