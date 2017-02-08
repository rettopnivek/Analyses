#---------------------#
# Model comparisons   #
# Kevin Potter        #
# Updated 01/29/2017  #
#---------------------#

# Clear workspace
rm( list = ls() )

# Save current directory
orig_dir = getwd()

# Define code segments to run
runCode = c( T, T )

# Indicate whether to save a pdf file
savePlot = T
if (savePlot) {
  setwd( 'Plots' )
  fName = 'Model_comparisons.pdf'
  pdf( fName, width = 12, height = 6 )
  rm( fName )
  setwd( orig_dir )
}

# Index
# Lookup - 01:  Load in useful packages and data
# Lookup - 02:  Comparison of standard SDT models
# Lookup - 03:  Comparison across model types

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

# Load in package for calculating WAIC
library(loo)

# Load in useful functions
source( 'F1_useful_functions.R' )

# Load in data
setwd( 'Data' )
load( 'Recog_mem.RData' )
setwd( orig_dir )

###
### Comparison of standard SDT models
###
# Lookup - 02

if ( runCode[1] ) {
  
  ### Load in posterior samples ###
  
  mNames = c(
    'SDT_model_M1_post.RData',
    'SDT_model_M2_post.RData',
    'SDT_model_M3_post.RData' )
  
  waic_all = c()
  loo_all = c()
  
  ### Calculate WAIC ###
  
  for (i in 1:3 ) {
    
    # Load in posterior
    fName = "C:/Users/Kevin/Documents/Posteriors from Stan"
    setwd(fName)
    setwd( "Wimber_et_al" )
    load( mNames[i] )
    setwd( orig_dir )
    
    waic_all = c( waic_all, list( waic( post$logLik ) ) )
    loo_all = c( loo_all, list( loo( post$logLik ) ) )
    
  }
  names( waic_all ) = c('M1','M2','M3')
  
  ### Aikake weights ###
  
  ll = c(
    waic_all$M1$waic,
    waic_all$M2$waic,
    waic_all$M3$waic
  )
  
  lbls = c(
    'Full model',
    'Main effects only',
    'SR-C' )
  
  ll = ll - min( ll )
  w = exp( -.5*ll )
  w = w/sum(w)
  
  if (!savePlot) x11()
  plot( c(1,length(w)), c(0,1), type = 'n', xlab = ' ',
        ylab = 'Akaike weights', xaxt = 'n', bty = 'n',
        main = 'Standard SDT' )
  abline( h = 0 )
  axis( 1, 1:3, lbls, tick = F )
  points( 1:length(w), w, pch = 19 )
  
  ### Pair-wise comparisons ###
  
  m1_v_m2 = compare( waic_all$M1, waic_all$M2 )
  m1_v_m3 = compare( waic_all$M1, waic_all$M3 )
  m2_v_m3 = compare( waic_all$M2, waic_all$M3 )
  
  ya = c(m1_v_m2[1],
         m1_v_m3[1],
         m2_v_m3[1])
  u = c(m1_v_m2[2],
         m1_v_m3[2],
         m2_v_m3[2])
  ui = cbind(
    ya - qnorm( .975 )*u,
    ya + qnorm( .975 )*u )
  
  if (!savePlot) x11()
  par( mar = c( 5, 4, 4, .5 ) )
  plot( c(0,4), c(-10,10), type = 'n', bty = 'n',
        xlab = ' ', xaxt = 'n', ylab = 'WAIC difference' )
  abline( h = 0 )
  axis( 1, 1:3, lbls[c(1,1,2)], tick = F )
  axis( 3, 1:3, lbls[c(2,3,3)], tick = F )
  
  segments( 1:3, ui[,1], 1:3, ui[,2] )
  points( 1:3, ya, pch = 19 )
  
}

###
### Comparison across model types
###
# Lookup - 03

if ( runCode[1] ) {
  
  ### Load in posterior samples ###
  
  mNames = c(
    'SDT_model_M3_post.RData',
    'rSDT_model_M1_post.RData',
    'mSDT_model_M1_post.RData',
    'RvSDT_model_post.RData' )
  
  waic_all = c()
  loo_all = c()
  
  ### Calculate WAIC ###
  
  for (i in 1:4 ) {
    
    # Load in posterior
    fName = "C:/Users/Kevin/Documents/Posteriors from Stan"
    setwd(fName)
    setwd( "Wimber_et_al" )
    load( mNames[i] )
    setwd( orig_dir )
    
    waic_all = c( waic_all, list( waic( post$logLik ) ) )
    loo_all = c( loo_all, list( loo( post$logLik ) ) )
    
  }
  names( waic_all ) = c('M1','M2','M3','M4')
  
  ### Aikake weights ###
  
  ll = c(
    waic_all$M1$waic,
    waic_all$M2$waic,
    waic_all$M3$waic,
    waic_all$M4$waic
  )
  
  lbls = c(
    'SDT',
    'rSDT',
    'mSDT',
    'RvSDT' )
  
  ll = ll - min( ll )
  w = exp( -.5*ll )
  w = w/sum(w)
  
  if (!savePlot) x11()
  plot( c(1,length(w)), c(0,1), type = 'n', xlab = ' ',
        ylab = 'Akaike weights', xaxt = 'n', bty = 'n',
        main = 'Comparisons across model types' )
  abline( h = 0 )
  axis( 1, 1:length(w), lbls, tick = F )
  points( 1:length(w), w, pch = 19 )
  
  ### Pair-wise comparisons ###
  
  m1_v_m2 = compare( waic_all$M1, waic_all$M2 )
  m1_v_m3 = compare( waic_all$M1, waic_all$M3 )
  m1_v_m4 = compare( waic_all$M1, waic_all$M4 )
  m2_v_m3 = compare( waic_all$M2, waic_all$M3 )
  m2_v_m4 = compare( waic_all$M2, waic_all$M4 )
  m3_v_m4 = compare( waic_all$M3, waic_all$M4 )
  
  ya = c(m1_v_m2[1],
         m1_v_m3[1],
         m1_v_m4[1],
         m2_v_m3[1],
         m2_v_m4[1],
         m3_v_m4[1] )
  u  = c(m1_v_m2[2],
         m1_v_m3[2],
         m1_v_m4[2],
         m2_v_m3[2],
         m2_v_m4[2],
         m3_v_m4[2] )
  ui = cbind(
    ya - qnorm( .975 )*u,
    ya + qnorm( .975 )*u )
  
  if (!savePlot) x11()
  par( mar = c( 5, 4, 4, .5 ) )
  yl = lowerUpper( 10, ui )
  plot( c(0,length(ya)+1), yl, type = 'n', bty = 'n',
        xlab = ' ', xaxt = 'n', ylab = 'WAIC difference' )
  abline( h = 0 )
  axis( 1, 1:length(ya), 
        lbls[c(1,1,1,2,2,3)], tick = F )
  axis( 3, 1:length(ya), 
        lbls[c(2,3,4,3,4,4)], tick = F )

  segments( 1:length(ya), 
            ui[,1], 
            1:length(ya), ui[,2] )
  points( 1:length(ya), ya, pch = 19 )
  
}

if (savePlot) dev.off()