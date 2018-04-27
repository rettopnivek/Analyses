# Model comparisons
# Written by Kevin Potter
# email: kevin.w.potter@gmail.com
# Please email me directory if you 
# have any questions or comments
# Last updated 2018-04-20

# Table of contents
# 1) Initial setup
# 2) Extract model estimation results
# 3) Compute information criterions

###
### 1) Initial setup
###

# Save current working directory
orig_dir = getwd()

# Navigate to project directory
setwd( '..' )
proj_dir = getwd()

# Indicate whether to save figures in a PDF file
savePlot = T
if ( savePlot ) {
  setwd( 'Figures' )
  pdf( file = 'Model_comparisons.pdf' )
  setwd( proj_dir )
}

# Load in useful packages

# Collection of useful functions
# devtools::install_github("rettopnivek/utilityf")
library( utilityf )

# Package for easier manipulation of data frames
# install.packages( 'dplyr' )
library( dplyr )

# Package for Bayesian estimation
# install.packages( 'rstan' )
library( rstan )
# For parallel processing
options(mc.cores = parallel::detectCores())
# Avoid recompilation
rstan_options(auto_write = TRUE)

# Load in data
setwd( 'Data' )
load( "NBack_3_19_2018.RData" )
setwd( proj_dir )

# Load in useful functions
setwd( 'R' )
source( 'Useful_functions.R' )
setwd( proj_dir )

# Load in estimation functions 
# and compile SDT model script
setwd( 'R' )
source( 'Estimation_functions.R' )

# Load in package for 'Leave one out' 
# information criterion
# install.packages( 'loo' )
library( loo )

# Exclude combined data
dtbf = dat %>% 
  filter( Task != 'Combined' )

# Drop rows with missing 
# self-report data for subject 
# FN_041 (66)
# Subject FN_041 did not have any 
# self report data on how high
# for the T1 drug condition
dtbf = dtbf %>% 
  filter( Self_report_on_high != '-' )

###
### 2) Extract model estimation results
###

# Navigate to folder with posterior estimates
# and posterior predictive checks
setwd( 'Data/Posterior_estimates' )

for ( m in 1:3 ) {
  model = paste( 'm', m, sep = '' )
  # Load in posterior estimates
  fname = paste( 'Posterior_', model, '.RData', sep = '' )
  load( fname )
}
setwd( proj_dir )

###
### 3) Compute information criterions
###

# Compute the Leave-one-out information criterion
loo_m1 = loo( m1$post$logLik )
loo_m2 = loo( m2$post$logLik )
loo_m3 = loo( m3$post$logLik )

loo_ic = numeric( 3 )
loo_ic[1] = loo_m1$looic
loo_ic[2] = loo_m2$looic
loo_ic[3] = loo_m3$looic

# Define function to calculate the Akaike weights
aw = function( ic ) {
  dic = ic - min( ic )
  rl = exp( -.5*dic )
  w = rl / sum( rl )
  
  return( w )
}

# Convert information criterions into Akaike weights
aw_sdt = aw( loo_ic )

# Indicate whether new plotting window should 
# be generated or if figure should be saved as PDF
if (!savePlot) x11()

xl = c( .5, 3.5 ); yl = c( -.05, 1.05 )
blankPlot( xl, yl )
customAxes( xl, yl )
par( mar = c( 7, 3, .5, .5 ) )

for ( i in 1:3 ) vertLines( i, c( 0, aw_sdt[i] ), lwd = 2 )
points( 1:3, aw_sdt, pch = 19, cex = 1.5 )

axis( 2, seq( 0, 1, .25 ), 
      tick = F, line = -1.75, cex.axis = 1.5 )
mtext( 'Relative likelihood to predict new data',
       side = 2, line = 1.5, cex = 1.5 )

xa = 1:3
lbl = c( 'Model 1 (Null)', 
         'Model 2 (Drug)', 
         'Model 3 \n(Drug + \nself-report)' )
axis( 1, at = xa, tick = F, labels = FALSE )
text( x = xa, 
      y = par()$usr[3]+0.0*(par()$usr[4]-par()$usr[3]),
      labels = lbl, srt = 45, adj = 1, xpd = TRUE )


if ( savePlot ) dev.off()
setwd( orig_dir )