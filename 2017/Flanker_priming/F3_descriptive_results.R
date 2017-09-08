#---------------------#
# Descriptive results #
# Kevin Potter        #
# Updated 05/01/2017  #
#---------------------#

# Clear workspace
rm( list = ls() )

# Save current directory
orig_dir = getwd()

# Define code segments to run
runCode = c( T, T, T, T, T, T, T )

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
# Lookup - 02:  Accuracy-latency figures (No prime)
# Lookup - 03:  Conditional accuracy functions (No prime)
# Lookup - 04:  Accuracy-latency figures (100 ms prime)
# Lookup - 05:  Conditional accuracy functions (100 prime)
# Lookup - 06:  Accuracy-latency figures (800 ms prime)
# Lookup - 07:  Conditional accuracy functions (800 prime)

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

# Define useful functions
source('F2_useful_functions.R')

# Load in data
setwd( 'Data' )
load( 'Flanker_priming.RData' )
setwd( orig_dir )

# For easy manipulation
d = allDat
colnames( d ) = c( 'S', 'Cnd', 'Co', 'Ac', 'Ch', 'RT', 
                   'T', 'P1', 'P3', 'ID', 'FL', 'F', 
                   'PT', 'PTL', 'PD', 'CL', 'RL' )

# Sort experimental conditions
exp_str = aggregate( d$Cnd, list( d$PD, d$PTL, d$FL, d$Cnd ), unique )
exp_str = exp_str[,-4]
colnames( exp_str ) = c( 'PD', 'PT', 'F','Cnd' )

###
### Accuracy-latency figures (No prime)
###
# Lookup - 02

if ( runCode[1] ) {
  
  # Select no priming conditions
  sel = d$PD == 0
  rt = d$RT[ sel ]; ch = d$Ac[ sel ];
  grp = d$S[ sel ];
  cvrt = d$Cnd[ sel ]
  
  # Create accuracy by latency plot
  if (!savePlot) x11(width=12)
  
  # Create 3 figure layout
  lyt = matrix( 2, 10, 10 )
  lyt[,6:10] = 3
  lyt[1:2,] = 1
  layout( lyt )
  
  # Plot experimental conditions
  display = list(
    None = c( "+++++", "+++++","++A++"),
    Identical = c( "+++++", "+++++","AAAAA"),
    Same = c( "+++++", "+++++","EEAEE"),
    Incompatiable = c( "+++++", "+++++","RRARR") )
  cnd_clr = c('black',colors()[26],colors()[616],'orange')
  
  par( mar = c(0,0,0,0) )
  blankPlot()
  design( 0, display, cnd_color )
  
  # Plot MRT by accuracy per condition
  par( mar = c( 4, 5, 3, 1 ) )
  xl = c(0,1); yl = c(.4,.6)
  blankPlot( xDim = xl, yDim = yl )
  
  axis( 1, seq(0,1,.2), cex.axis = 1.5 )
  axis( 2, seq(.4,.6,.05), cex.axis = 1.5 )
  mtext( 'MRT(s)', side = 2, line = 2.5 )
  mtext( 'Accuracy', side = 1, line = 2.5 )
  
  create_grid( xl, yl, .1, .025, col = 'grey80' )
  
  out = pvt_points( rt, ch, cvrt, grp = grp, T_x = mean,
              opt = list( draw = F, out = T ) )
  segments( out$pv$x[1:4], out$pv$y[1:4],
            out$pv$x[5:8], out$pv$y[5:8],
            col = 'grey' )
  
  plt = list( pch = c(21,24,25,22), 
              bg = cnd_clr )
  pvt_points( rt, ch, cvrt, grp = grp, T_x = mean,
              plt = plt, cex = 2 )
  
  # Plot median RT by accuracy per condition
  blankPlot( xDim = c(0,1), yDim = c(.4,.6) )
  
  axis( 1, seq(0,1,.2), cex.axis = 1.5 )
  axis( 2, seq(.4,.6,.05), cex.axis = 1.5 )
  mtext( 'Median RT(s)', side = 2, line = 2.5 )
  mtext( 'Accuracy', side = 1, line = 2.5 )
  
  create_grid( xl, yl, .1, .025, col = 'grey80' )
  
  out = pvt_points( rt, ch, cvrt, grp = grp, T_x = median,
                    opt = list( draw = F, out = T ) )
  segments( out$pv$x[1:4], out$pv$y[1:4],
            out$pv$x[5:8], out$pv$y[5:8],
            col = 'grey' )
  
  plt = list( pch = c(21,24,25,22), 
              bg = cnd_clr )
  pvt_points( rt, ch, cvrt, grp = grp, T_x = median,
              plt = plt, cex = 2 )
  
}

###
### Conditional accuracy functions (No prime)
###
# Lookup - 03

if ( runCode[2] ) {
  
  # Select no priming conditions
  sel = d$PD == 0
  rt = d$RT[ sel ]; ch = d$Ac[ sel ];
  grp = d$S[ sel ];
  cvrt = d$Cnd[ sel ]
  
  # Create plot of conditional accuracy functions
  if (!savePlot) x11(width=12)
  
  # Create 3 figure layout
  lyt = matrix( 3, 10, 10 )
  lyt[3:10,3:7] = 2
  lyt[1:2,] = 1
  layout( lyt )
  
  # Plot experimental conditions
  display = list(
    None = c( "+++++", "+++++","++A++"),
    Identical = c( "+++++", "+++++","AAAAA"),
    Same = c( "+++++", "+++++","EEAEE"),
    Incompatiable = c( "+++++", "+++++","RRARR") )
  cnd_clr = c('black',colors()[26],colors()[616],'orange')
  
  par( mar = c(0,0,0,0) )
  blankPlot()
  design( 0, display, cnd_color )
  
  # CAF
  xl = c( .3, .8 )
  yl = c( .5, 1 )
  par( mar = c( 4, 5, 3, 1 ) )
  blankPlot( xl, yl )
  
  axis( 1, seq(xl[1],xl[2],.1), cex.axis = 1.5 )
  axis( 2, seq(yl[1],yl[2],.1), cex.axis = 1.5 )
  mtext( 'Conditional accuracy function', side = 2, line = 2.5 )
  mtext( 'RT (s)', side = 1, line = 2.5 )
  
  create_grid( xl, yl, .05, .05, col = 'grey80' )
  
  pts = c(21,24,25,22)
  
  for ( i in unique( cvrt ) ) {
    sel = cvrt == i
    caf_points( rt[sel], ch[sel], grp = grp[sel],
                lwd = 2, col = cnd_clr[i],
                opt = list( pts = F ) )
    caf_points( rt[sel], ch[sel], grp = grp[sel],
                bg = cnd_clr[i], cex = 2, pch = pts[i],
                opt = list( jnt = T, pts = T ) )
  }
}


###
### Accuracy-latency figures (100 ms prime)
###
# Lookup - 04

if ( runCode[3] ) {
  
  # Select no priming conditions
  sel = d$PD == 0.1
  rt = d$RT[ sel ]; ch = d$Ac[ sel ];
  grp = d$S[ sel ];
  cvrt = createIncrement( d$Cnd[ sel ] )
  
  # Create accuracy by latency plot
  if (!savePlot) x11(width=12)
  
  # Create 4 figure layout
  lyt = matrix( 4, 10, 10 )
  lyt[,2:5] = 2
  lyt[,6:9] = 3
  lyt[1:3,] = 1
  layout( lyt )
  
  # Plot experimental conditions
  display = list(
    Identical_None = c( "+++++", "++a++","++A++"),
    Same_None = c( "+++++", "++e++","++A++"),
    Incompatible_None = c( "+++++", "++r++","++A++"),
    Identical = c( "+++++", "++A++","AAAAA"),
    Same = c( "+++++", "++E++","EEAEE"),
    Incompatiable = c( "+++++", "++R++","RRARR") )
  cnd_clr = c(colors()[26],colors()[616],'orange',
              colors()[26],colors()[616],'orange')
  
  par( mar = c(0,0,0,0) )
  blankPlot()
  design( 100, display, cnd_color, yax = c(.075,.1),
          shft = c(.15,.1),txtSz=1.2 )
  
  legend( 'topleft', c( 'Inner = prime', 'Outer = flanker' ),
          bty = 'n', cex =  1.2 )
  
  # Plot MRT by accuracy per condition
  par( mar = c(5, 4, 4, 2) )
  xl = c( 0, 1 ); yl = c(.4,.6)
  blankPlot( xl, yl )
  
  axis( 1, seq(0,1,.2) )
  axis( 2, seq(.4,.6,.05) )
  mtext( 'MRT(s)', side = 2, line = 2.5 )
  mtext( 'Accuracy', side = 1, line = 2.5 )
  
  create_grid( xl, yl, .1, .025, col = 'grey80' )
  
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
  blankPlot( xl, yl )
  
  axis( 1, seq(0,1,.2) )
  axis( 2, seq(.4,.6,.05) )
  mtext( 'Median RT(s)', side = 2, line = 2.5 )
  mtext( 'Accuracy', side = 1, line = 2.5 )
  
  create_grid( xl, yl, .1, .025, col = 'grey80' )
  
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
}

###
### Conditional accuracy functions (100 prime)
###
# Lookup - 05

if ( runCode[4] ) {
  
  # Select no priming conditions
  sel = d$PD == .1
  rt = d$RT[ sel ]; ch = d$Ac[ sel ];
  grp = d$S[ sel ];
  cvrt = createIncrement( d$Cnd[ sel ] )
  
  # Create plot of conditional accuracy functions
  if (!savePlot) x11(width=12)
  
  # Create 3 figure layout
  lyt = matrix( 3, 10, 10 )
  lyt[4:10,3:7] = 2
  lyt[1:3,] = 1
  layout( lyt )
  
  # Plot experimental conditions
  display = list(
    Identical_None = c( "+++++", "++a++","++A++"),
    Same_None = c( "+++++", "++e++","++A++"),
    Incompatible_None = c( "+++++", "++r++","++A++"),
    Identical = c( "+++++", "++A++","AAAAA"),
    Same = c( "+++++", "++E++","EEAEE"),
    Incompatiable = c( "+++++", "++R++","RRARR") )
  cnd_clr = c(colors()[26],colors()[616],'orange',
              colors()[26],colors()[616],'orange')
  
  par( mar = c(0,0,0,0) )
  blankPlot()
  design( 100, display, cnd_color, yax = c(.075,.1),
          shft = c(.15,.1),txtSz=1.2 )
  
  legend( 'topleft', c( 'Inner = prime', 'Outer = flanker' ),
          bty = 'n', cex =  1.2 )
  
  # CAF
  xl = c( .3, .8 )
  yl = c( .5, 1 )
  par( mar = c( 4, 5, 3, 1 ) )
  blankPlot( xl, yl )
  
  axis( 1, seq(xl[1],xl[2],.1), cex.axis = 1.5 )
  axis( 2, seq(yl[1],yl[2],.1), cex.axis = 1.5 )
  mtext( 'Conditional accuracy function', side = 2, line = 2.5 )
  mtext( 'RT (s)', side = 1, line = 2.5 )
  
  create_grid( xl, yl, .05, .05, col = 'grey80' )
  
  cnd_clr1 = cnd_clr; cnd_clr1[1:3] = 'white'
  pts1 = c(21,21,21,24,25,22)
  pts2 = c(24,25,22,24,25,22)
  
  for ( i in unique( cvrt ) ) {
    sel = cvrt == i
    caf_points( rt[sel], ch[sel], grp = grp[sel],
                lwd = 2, col = cnd_clr[i],
                opt = list( pts = F ) )
    caf_points( rt[sel], ch[sel], grp = grp[sel],
                bg = cnd_clr1[i], cex = 3, pch = pts1[i],
                opt = list( pts = T ) )
    caf_points( rt[sel], ch[sel], grp = grp[sel],
                bg = cnd_clr[i], cex = 1.5, pch = pts2[i],
                opt = list( pts = T ) )
  }
}


###
### Accuracy-latency figures (800 ms prime)
###
# Lookup - 06

if ( runCode[5] ) {
  
  # Select no priming conditions
  sel = d$PD == 0.8
  rt = d$RT[ sel ]; ch = d$Ac[ sel ];
  grp = d$S[ sel ];
  cvrt = createIncrement( d$Cnd[ sel ] )
  
  # Create accuracy by latency plot
  if (!savePlot) x11(width=12)
  
  # Create 4 figure layout
  lyt = matrix( 4, 10, 10 )
  lyt[,2:5] = 2
  lyt[,6:9] = 3
  lyt[1:3,] = 1
  layout( lyt )
  
  # Plot experimental conditions
  display = list(
    Identical_None = c( "+++++", "++a++","++A++"),
    Same_None = c( "+++++", "++e++","++A++"),
    Incompatible_None = c( "+++++", "++r++","++A++"),
    Identical = c( "+++++", "++A++","AAAAA"),
    Same = c( "+++++", "++E++","EEAEE"),
    Incompatiable = c( "+++++", "++R++","RRARR") )
  cnd_clr = c(colors()[26],colors()[616],'orange',
              colors()[26],colors()[616],'orange')
  
  par( mar = c(0,0,0,0) )
  blankPlot()
  design( 800, display, cnd_color, yax = c(.075,.1),
          shft = c(.15,.1),txtSz=1.2 )
  
  legend( 'topleft', c( 'Inner = prime', 'Outer = flanker' ),
          bty = 'n', cex =  1.2 )
  
  # Plot MRT by accuracy per condition
  par( mar = c(5, 4, 4, 2) )
  xl = c( 0, 1 ); yl = c(.4,.6)
  blankPlot( xl, yl )
  
  axis( 1, seq(0,1,.2) )
  axis( 2, seq(.4,.6,.05) )
  mtext( 'MRT(s)', side = 2, line = 2.5 )
  mtext( 'Accuracy', side = 1, line = 2.5 )
  
  create_grid( xl, yl, .1, .025, col = 'grey80' )
  
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
  blankPlot( xl, yl )
  
  axis( 1, seq(0,1,.2) )
  axis( 2, seq(.4,.6,.05) )
  mtext( 'Median RT(s)', side = 2, line = 2.5 )
  mtext( 'Accuracy', side = 1, line = 2.5 )
  
  create_grid( xl, yl, .1, .025, col = 'grey80' )
  
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
}

###
### Conditional accuracy functions (800 prime)
###
# Lookup - 07

if ( runCode[6] ) {
  
  # Select no priming conditions
  sel = d$PD == .8
  rt = d$RT[ sel ]; ch = d$Ac[ sel ];
  grp = d$S[ sel ];
  cvrt = createIncrement( d$Cnd[ sel ] )
  
  # Create plot of conditional accuracy functions
  if (!savePlot) x11(width=12)
  
  # Create 3 figure layout
  lyt = matrix( 3, 10, 10 )
  lyt[4:10,3:7] = 2
  lyt[1:3,] = 1
  layout( lyt )
  
  # Plot experimental conditions
  display = list(
    Identical_None = c( "+++++", "++a++","++A++"),
    Same_None = c( "+++++", "++e++","++A++"),
    Incompatible_None = c( "+++++", "++r++","++A++"),
    Identical = c( "+++++", "++A++","AAAAA"),
    Same = c( "+++++", "++E++","EEAEE"),
    Incompatiable = c( "+++++", "++R++","RRARR") )
  cnd_clr = c(colors()[26],colors()[616],'orange',
              colors()[26],colors()[616],'orange')
  
  par( mar = c(0,0,0,0) )
  blankPlot()
  design( 800, display, cnd_color, yax = c(.075,.1),
          shft = c(.15,.1),txtSz=1.2 )
  
  legend( 'topleft', c( 'Inner = prime', 'Outer = flanker' ),
          bty = 'n', cex =  1.2 )
  
  # CAF
  xl = c( .3, .8 )
  yl = c( .5, 1 )
  par( mar = c( 4, 5, 3, 1 ) )
  blankPlot( xl, yl )
  
  axis( 1, seq(xl[1],xl[2],.1), cex.axis = 1.5 )
  axis( 2, seq(yl[1],yl[2],.1), cex.axis = 1.5 )
  mtext( 'Conditional accuracy function', side = 2, line = 2.5 )
  mtext( 'RT (s)', side = 1, line = 2.5 )
  
  create_grid( xl, yl, .05, .05, col = 'grey80' )
  
  cnd_clr1 = cnd_clr; cnd_clr1[1:3] = 'white'
  pts1 = c(21,21,21,24,25,22)
  pts2 = c(24,25,22,24,25,22)
  
  for ( i in unique( cvrt ) ) {
    sel = cvrt == i
    caf_points( rt[sel], ch[sel], grp = grp[sel],
                lwd = 2, col = cnd_clr[i],
                opt = list( pts = F ) )
    caf_points( rt[sel], ch[sel], grp = grp[sel],
                bg = cnd_clr1[i], cex = 3, pch = pts1[i],
                opt = list( pts = T ) )
    caf_points( rt[sel], ch[sel], grp = grp[sel],
                bg = cnd_clr[i], cex = 1.5, pch = pts2[i],
                opt = list( pts = T ) )
  }
}

###
### Joint CDF plots
###
# Lookup - 08

if ( runCode[7] ) {
  
  # Function for uncertainty intervals
  f = function(x,alpha) {
    xbar = mean(x)
    se = sem(x)
    n = length(x)
    crt = abs( qt( alpha, n - 1 ) )
    return( c( xbar - crt * se, xbar + crt * se ) )
  }
  
  # Prime duration = 0 ms
  if ( !savePlot ) x11( width = 12 )
  
  lyt = matrix( 5, 2, 12 )
  lyt[,3:10] = matrix( rep( 1:4, each = 2 ), 2, 8, byrow = T )
  layout( lyt )
  
  v = list(
    PD = rep( 0, 4 ),
    PT = rep( 'None', 4 ),
    FL = c( 'None', 'Identical','Same','Incompatible')
  )
  
  for ( j in 1:4 ) {
    cnd = list(
      c1 = list( PD = v$PD[j], PT = v$PT[j], FL = v$FL[j], Co = 1, flip = F ),
      c2 = list( PD = v$PD[j], PT = v$PT[j], FL = v$FL[j], Co = 0, flip = T ) )
    
    plot_jcdf( d, cnd, f = f, opt = list( yl = c(-1,1), xl = c(.2,1.2) ) )
    title( paste( 'Flanker:', v$FL[j] ) )
  }
  
  # Prime duration = 100 ms
  if ( !savePlot ) x11( width = 12 )
  
  lyt = matrix( rep( 1:6, each = 2 ), 2, 12, byrow = T )
  layout( lyt )
  
  v = list(
    PD = rep( .1, 6 ),
    PT = rep( c( 'Identical','Same','Incompatible' ), 3 ),
    FL = c( rep( 'None', 3 ), 'Identical','Same','Incompatible' )
  )
  
  for ( j in 1:6 ) {
    cnd = list(
      c1 = list( PD = v$PD[j], PT = v$PT[j], FL = v$FL[j], Co = 1, flip = F ),
      c2 = list( PD = v$PD[j], PT = v$PT[j], FL = v$FL[j], Co = 0, flip = T ) )
    
    plot_jcdf( d, cnd, f = f, opt = list( yl = c(-1,1), xl = c(.2,1.2) ) )
    title( paste( 'Flanker:', v$FL[j] ) )
  }
  
  # Prime duration = 800 ms
  if ( !savePlot ) x11( width = 12 )
  
  lyt = matrix( rep( 1:6, each = 2 ), 2, 12, byrow = T )
  layout( lyt )
  
  v = list(
    PD = rep( .8, 6 ),
    PT = rep( c( 'Identical','Same','Incompatible' ), 3 ),
    FL = c( rep( 'None', 3 ), 'Identical','Same','Incompatible' )
  )
  
  for ( j in 1:6 ) {
    cnd = list(
      c1 = list( PD = v$PD[j], PT = v$PT[j], FL = v$FL[j], Co = 1, flip = F ),
      c2 = list( PD = v$PD[j], PT = v$PT[j], FL = v$FL[j], Co = 0, flip = T ) )
    
    plot_jcdf( d, cnd, f = f, opt = list( yl = c(-1,1), xl = c(.2,1.2) ) )
    title( paste( 'Flanker:', v$FL[j] ) )
  }
  
}

if (savePlot) dev.off()