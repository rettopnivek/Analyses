# Checks for correlations between self-report and physio
# Written by Kevin Potter (based on work by Adrian Alvarez)
# email: kevin.w.potter@gmail.com
# Please email me directly if you 
# have any questions or comments
# Last updated 2018-05-18

# Table of contents
# 1) Initial setup
# 2) Double-check correlations
# 3) Scatterplots

###
### 1) Initial setup
###

# Indicate whether to save plots as a PDF
savePlot = T

# Save current working directory
orig_dir = getwd()

# Navigate to project directory
setwd( '..' )
proj_dir = getwd()

# Load in useful packages

# Collection of useful functions
# devtools::install_github("rettopnivek/utilityf")
library( utilityf )

# Load in data
setwd( 'Data/Original_files' )
dat_dir = getwd()
dtba = read.csv( file = 'Checks_for_Adrian.csv', 
                 header = T,
                 stringsAsFactors = FALSE )
dtba2 = read.csv( file = 'Checks_for_Adrian_above_0.csv', 
                 header = T,
                 stringsAsFactors = FALSE )
setwd( proj_dir )

###
### 2) Double-check correlations
###

message( 'All observations' )
message( 'Self report vs. Heart rate' )
message( round( cor( dtba$SR, dtba$HR ), 3 ) )
message( 'Self report vs. Systolic' )
message( round( cor( dtba$SR, dtba$Sys ), 3 ) )
message( 'Self report vs. Diastolic' )
message( round( cor( dtba$SR, dtba$Dia ), 3 ) )

message(' ')

message( 'Above 0' )
message( 'Self report vs. Heart rate' )
message( round( cor( dtba2$SR, dtba2$HR ), 3 ) )
message( 'Self report vs. Systolic' )
message( round( cor( dtba2$SR, dtba2$Sys ), 3 ) )
message( 'Self report vs. Diastolic' )
message( round( cor( dtba2$SR, dtba2$Dia ), 3 ) )

###
### 3) Scatterplots
###

# 3.1)
quick_scatter_plot = function( x, y,
                               lbls, 
                               pts = 21,
                               bc = 'black',
                               fc = 'black',
                               ptSz = 1,
                               lnSz = 2,
                               axPos = -1,
                               axSz = 1.25 ) {
  
  # Means and standard deviations
  xm = mean(x)
  xs = sd(x)
  ym = mean(y)
  ys = sd(y)
  
  # Standardized scores
  zx = ( x - xm )/xs
  zy = ( y - ym )/ys
  
  # Create a blank plot
  xl = lowerUpper( 1, zx )
  yl = lowerUpper( 1, zy )
  blankPlot( xl, yl )
  
  # Grid lines
  horizLines( seq( yl[1], yl[2], .5 ), xl,
              lwd = lnSz, col = 'grey' )
  
  # Add observations
  points( zx, zy, pch = pts, 
          bg = bc, col = fc, 
          cex = ptSz )
  
  # Regression line
  lr = lm( zx ~ zy )
  abline( lr, lwd = lnSz, lty = 2, col = 'red' )
  
  # Add axes and labels
  customAxes( xl, yl,
              label = lbls, 
              lnSz = lnSz )
  axis( 1, seq( xl[1], xl[2], 1 ),
        tick = F, line = axPos, cex.axis = axSz )
  axis( 2, seq( yl[1], yl[2], 1 ),
        tick = F, line = axPos, cex.axis = axSz )
  
}

if ( !savePlot ) x11( width = 12 ) else {
  setwd( 'Figures' )
  pdf( 'Checks for Adrian.pdf', width = 12 )
  setwd( proj_dir )
}
layout( matrix( 1:6, 2, 3, byrow = F ) )

sel = c( 'HR', 'Sys', 'Dia' )
lbls = c( 'Heart rate z-scores',
          'Systolic z-scores',
          'Diastolic z-scores' )
for ( i in 1:3 ) {
  
  quick_scatter_plot( dtba$SR, dtba[,sel[i]],
                      c( 'Self-report z-scores',
                         lbls[i] ) )
  if ( i == 2 ) title( 'All' )
  
  quick_scatter_plot( dtba2$SR, dtba2[,sel[i]],
                      c( 'Self-report z-scores',
                         lbls[i] ) )
  if ( i == 2 ) title( 'Above 0' )
  
}

if ( savePlot ) dev.off()

###
### 4) Additional analyses
###

library( dplyr )
library( BayesFactor )

setwd( proj_dir )
setwd( 'Data' )
setwd( 'Original_files' )

# Read in data
dtba = read.csv( file = 'fnirsAdrianDRE.csv',
                 header = T,
                 stringsAsFactors = F )
# Determine subjects without DRE data
dtba$Missing_data = 
  apply( dtba, 1, function(x) any( is.na( x ) ) )

dtba = dtba %>% arrange( Subject )

dtbp = data.frame(
  ID = dtba$Subject[ dtba$Category == 'Placebo' ],
  stringsAsFactors = F
)

f = function( val ) {
  out = val[ dtba$Category == 'Placebo' ] - 
    val[ dtba$Category != 'Placebo' ]
  return( out )
}
x = f( dtba$Room_Light.Left - dtba$Direct_Light.Left )
bf = ttestBF( na.omit( x ) ); 1/bf
t.test( x )
x = f( dtba$Room_Light.Right - dtba$Direct_Light.Right )
bf = ttestBF( na.omit( x ) ); 1/bf
t.test( x )

x = f( dtba$Room_Light.Left - dtba$Near_Total_Darkness.Left )
bf = ttestBF( na.omit( x ) ); 1/bf
t.test( x )
x = f( dtba$Room_Light.Right - dtba$Near_Total_Darkness.Right )
bf = ttestBF( na.omit( x ) ); 1/bf
t.test( x )

setwd( orig_dir )
