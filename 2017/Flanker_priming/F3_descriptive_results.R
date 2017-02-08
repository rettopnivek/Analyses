#---------------------#
# Descriptive results #
# Kevin Potter        #
# Updated 02/05/2017  #
#---------------------#

# Clear workspace
rm( list = ls() )

# Save current directory
orig_dir = getwd()

# Define code segments to run
runCode = c( T, T, T )

# Indicate whether to save a pdf file
savePlot = T
if (savePlot) {
  setwd( 'Plots' )
  fName = 'Descriptive_results.pdf'
  pdf( fName, width = 12, height = 6 )
  rm( fName )
  setwd( orig_dir )
}

# Index
# Lookup - 01:  Load in useful packages and data

###
### Load in useful packages and data
###
# Lookup - 01

# For geting github packages
# install.packages(devtools)
# library(devtools)

# Integration of C++ and R
# install.packages(Rcpp)
library(Rcpp)

# Miscellanous functions for modeling/plotting
# install_github("rettopnivek/utilityf")
library(utilityf)

# Distribution functions for assorted sequential sampling models
# install_github("rettopnivek/seqmodels")
library(seqmodels)

# Functions for plotting response time and choice data
# install_github("rettopnivek/rtplots")
library(rtplots)

# Functions for plotting response time and choice data
# install_github("rettopnivek/nROUSE")
library(nROUSE)

# Define useful functions
source('F2_useful_functions.R')

# Load in data
setwd( 'Data' )
load( 'Flanker_priming.RData' )
setwd( orig_dir )

# For easy manipulation
d = allDat
colnames( d ) = c( 'S', 'Cnd', 'Co', 'Ac', 'Ch', 'RT', 
                   'T', 'P1', 'P3', 'ID', 'F', 'PT', 'D' )

###
###
###
# Lookup - 02

if ( runCode[1] ) {
  
  # Select no priming conditions
  sel = d$D == 0
  rt = d$RT[ sel ]; ch = d$Ac[ sel ];
  grp = d$S[ sel ];
  cvrt = d$Cnd[ sel ]
  
  # Create accuracy by latency plot
  if (!savePlot) x11(width=12)
  
  # Create 3 figure layout
  lyt = matrix( 1, 10, 10 )
  lyt[,6:10] = 2
  lyt[10,] = 3
  layout( lyt )
  
  # Plot MRT by accuracy per condition
  blankPlot( xDim = c(0,1), yDim = c(.4,.6) )
  axis( 1, seq(0,1,.2) )
  axis( 2, seq(.4,.6,.05) )
  mtext( 'MRT(s)', side = 2, line = 2.5 )
  mtext( 'Accuracy', side = 1, line = 2.5 )
  
  out = pvt_points( rt, ch, cvrt, grp = grp, T_x = mean,
              opt = list( draw = F, out = T ) )
  segments( out$pv$x[1:4], out$pv$y[1:4],
            out$pv$x[5:8], out$pv$y[5:8],
            col = 'grey' )
  
  plt = list( pch = c(21,24,25,22), 
              bg = c('white',colors()[26],colors()[616],'orange') )
  pvt_points( rt, ch, cvrt, grp = grp, T_x = mean,
              plt = plt, cex = 2 )
  
  # Plot median RT by accuracy per condition
  blankPlot( xDim = c(0,1), yDim = c(.4,.6) )
  axis( 1, seq(0,1,.2) )
  axis( 2, seq(.4,.6,.05) )
  mtext( 'Median RT(s)', side = 2, line = 2.5 )
  mtext( 'Accuracy', side = 1, line = 2.5 )
  
  out = pvt_points( rt, ch, cvrt, grp = grp, T_x = median,
                    opt = list( draw = F, out = T ) )
  segments( out$pv$x[1:4], out$pv$y[1:4],
            out$pv$x[5:8], out$pv$y[5:8],
            col = 'grey' )
  
  plt = list( pch = c(21,24,25,22), 
              bg = c('white',colors()[26],colors()[616],'orange') )
  pvt_points( rt, ch, cvrt, grp = grp, T_x = median,
              plt = plt, cex = 2 )
  
  par( mar = c(0,0,0,0) )
  blankPlot()
  legend( 'left', c( 'Flanker type: ', 'None', 'Identical',
                    'Same response', 'Incompatible' ),
          pch = c(19,21,24,25,22),
          col = c('white','black','black','black','black'),
          pt.bg = c('white','white',colors()[26],
                 colors()[616],'orange'),
  horiz = T, cex = 2, bty = 'n' )
  
  mtext( 'Prime duration: 0 s', side = 3, outer = T, line = -2,
         cex = 1.5 )
  
}

###
###
###
# Lookup - 03

if ( runCode[2] ) {
  
  # Select no priming conditions
  sel = d$D == 0.1
  rt = d$RT[ sel ]; ch = d$Ac[ sel ];
  grp = d$S[ sel ];
  cvrt = createIncrement( d$Cnd[ sel ] )
  
  # Create accuracy by latency plot
  if (!savePlot) x11(width=12)
  par( mar = c(5, 4, 4, 2) )
  
  # Create 3 figure layout
  lyt = matrix( 1, 10, 10 )
  lyt[,6:10] = 2
  lyt[10,] = 3
  layout( lyt )
  
  # Plot MRT by accuracy per condition
  blankPlot( xDim = c(0,1), yDim = c(.4,.6) )
  axis( 1, seq(0,1,.2) )
  axis( 2, seq(.4,.6,.05) )
  mtext( 'MRT(s)', side = 2, line = 2.5 )
  mtext( 'Accuracy', side = 1, line = 2.5 )
  
  out = pvt_points( rt, ch, cvrt, grp = grp, T_x = mean,
                    opt = list( draw = F, out = T ) )
  segments( out$pv$x[1:6], out$pv$y[1:6],
            out$pv$x[7:12], out$pv$y[7:12],
            col = 'grey' )
  
  plt = list( pch = c(21,21,21,24,25,22), 
              bg = c('white','white','white',
                     colors()[26],colors()[616],'orange') )
  pvt_points( rt, ch, cvrt, grp = grp, T_x = mean,
              plt = plt, cex = 3 )
  
  plt = list( pch = c(24,25,22,24,25,22), 
              bg = c(colors()[26],colors()[616],'orange',
                     colors()[26],colors()[616],'orange') )
  pvt_points( rt, ch, cvrt, grp = grp, T_x = mean,
              plt = plt, cex = 1.5 )
  
  # Plot median RT by accuracy per condition
  blankPlot( xDim = c(0,1), yDim = c(.4,.6) )
  axis( 1, seq(0,1,.2) )
  axis( 2, seq(.4,.6,.05) )
  mtext( 'Median RT(s)', side = 2, line = 2.5 )
  mtext( 'Accuracy', side = 1, line = 2.5 )
  
  out = pvt_points( rt, ch, cvrt, grp = grp, T_x = median,
                    opt = list( draw = F, out = T ) )
  segments( out$pv$x[1:6], out$pv$y[1:6],
            out$pv$x[7:12], out$pv$y[7:12],
            col = 'grey' )
  
  plt = list( pch = c(21,21,21,24,25,22), 
              bg = c('white','white','white',
                     colors()[26],colors()[616],'orange') )
  pvt_points( rt, ch, cvrt, grp = grp, T_x = median,
              plt = plt, cex = 3 )
  
  plt = list( pch = c(24,25,22,24,25,22), 
              bg = c(colors()[26],colors()[616],'orange',
                     colors()[26],colors()[616],'orange') )
  pvt_points( rt, ch, cvrt, grp = grp, T_x = median,
              plt = plt, cex = 1.5 )
  
  par( mar = c(0,0,0,0) )
  blankPlot()
  legend( 'topleft', c( 'Type: ', 'None', 'Identical',
                     'Same response', 'Incompatible' ),
          pch = c(19,21,24,25,22),
          col = c('white','black','black','black','black'),
          pt.bg = c('white','white',colors()[26],
                    colors()[616],'orange'),
          horiz = T, cex = 1.5, bty = 'n' )
  legend( 'bottomright', c('Inner = prime','Outer = flanker'),
          horiz = T, cex = 1.5, bty = 'n' )
  
  mtext( 'Prime duration: .1 s', side = 3, outer = T, line = -2,
         cex = 1.5 )
  
}

###
###
###
# Lookup - 04

if ( runCode[3] ) {
  
  # Select no priming conditions
  sel = d$D == 0.8
  rt = d$RT[ sel ]; ch = d$Ac[ sel ];
  grp = d$S[ sel ];
  cvrt = createIncrement( d$Cnd[ sel ] )
  
  # Create accuracy by latency plot
  if (!savePlot) x11(width=12)
  par( mar = c(5, 4, 4, 2) )
  
  # Create 3 figure layout
  lyt = matrix( 1, 10, 10 )
  lyt[,6:10] = 2
  lyt[10,] = 3
  layout( lyt )
  
  # Plot MRT by accuracy per condition
  blankPlot( xDim = c(0,1), yDim = c(.4,.6) )
  axis( 1, seq(0,1,.2) )
  axis( 2, seq(.4,.6,.05) )
  mtext( 'MRT(s)', side = 2, line = 2.5 )
  mtext( 'Accuracy', side = 1, line = 2.5 )
  
  out = pvt_points( rt, ch, cvrt, grp = grp, T_x = mean,
                    opt = list( draw = F, out = T ) )
  segments( out$pv$x[1:6], out$pv$y[1:6],
            out$pv$x[7:12], out$pv$y[7:12],
            col = 'grey' )
  
  plt = list( pch = c(21,21,21,24,25,22), 
              bg = c('white','white','white',
                     colors()[26],colors()[616],'orange') )
  pvt_points( rt, ch, cvrt, grp = grp, T_x = mean,
              plt = plt, cex = 3 )
  
  plt = list( pch = c(24,25,22,24,25,22), 
              bg = c(colors()[26],colors()[616],'orange',
                     colors()[26],colors()[616],'orange') )
  pvt_points( rt, ch, cvrt, grp = grp, T_x = mean,
              plt = plt, cex = 1.5 )
  
  # Plot median RT by accuracy per condition
  blankPlot( xDim = c(0,1), yDim = c(.4,.6) )
  axis( 1, seq(0,1,.2) )
  axis( 2, seq(.4,.6,.05) )
  mtext( 'Median RT(s)', side = 2, line = 2.5 )
  mtext( 'Accuracy', side = 1, line = 2.5 )
  
  out = pvt_points( rt, ch, cvrt, grp = grp, T_x = median,
                    opt = list( draw = F, out = T ) )
  segments( out$pv$x[1:6], out$pv$y[1:6],
            out$pv$x[7:12], out$pv$y[7:12],
            col = 'grey' )
  
  plt = list( pch = c(21,21,21,24,25,22), 
              bg = c('white','white','white',
                     colors()[26],colors()[616],'orange') )
  pvt_points( rt, ch, cvrt, grp = grp, T_x = median,
              plt = plt, cex = 3 )
  
  plt = list( pch = c(24,25,22,24,25,22), 
              bg = c(colors()[26],colors()[616],'orange',
                     colors()[26],colors()[616],'orange') )
  pvt_points( rt, ch, cvrt, grp = grp, T_x = median,
              plt = plt, cex = 1.5 )
  
  par( mar = c(0,0,0,0) )
  blankPlot()
  legend( 'topleft', c( 'Type: ', 'None', 'Identical',
                        'Same response', 'Incompatible' ),
          pch = c(19,21,24,25,22),
          col = c('white','black','black','black','black'),
          pt.bg = c('white','white',colors()[26],
                    colors()[616],'orange'),
          horiz = T, cex = 1.5, bty = 'n' )
  legend( 'bottomright', c('Inner = prime','Outer = flanker'),
          horiz = T, cex = 1.5, bty = 'n' )
  
  mtext( 'Prime duration: .8 s', side = 3, outer = T, line = -2,
         cex = 1.5 )
  
}

if (savePlot) dev.off()