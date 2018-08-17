# Examination of fNIRS data
# Written by Kevin Potter
# email: kevin.w.potter@gmail.com
# Please email me directory if you 
# have any questions or comments
# Last updated 2018-05-16

# Table of contents
# 1) Initial setup
# 2) Convert neural data to tidy form
# 3) Plot of drug/placebo HbO per subject for PRE
# 4) 


###
### 1) Initial setup
###

# Save current working directory
orig_dir = getwd()

# Navigate to project directory
setwd( '..' )
proj_dir = getwd()

# Indicate which code segments to run
run_code = c(
  T, # 
  F, # 
  F, # 
  F  # 
)

# Load in useful packages

# Collection of useful functions
# devtools::install_github("rettopnivek/utilityf")
library( utilityf )

# Packages for easier manipulation of data frames
# install.packages( 'tidyr' )
# install.packages( 'dplyr' )
library(dplyr)

# Package for mixed effects models
# install.packages( 'lme4' )
library( lme4 )

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

# Custom functions
setwd( 'R' );
source( 'S02_Useful_functions.R' )
setwd( proj_dir )

# Indicate whether to save figures in a PDF file
savePlot = F
if ( savePlot ) {
  setwd( 'Figures' )
  pdf( 'fNIRS_examination.pdf', width = 12 )
  setwd( proj_dir )
}

###
### 2) Convert neural data to tidy form
###

# Convert tidy version of data frame with HbO estimates
HbO = neurDat %>% 
  tidyr::gather( key = Channel, value = HbO,
          Channel_1, Channel_2, Channel_3,
          Channel_4, Channel_5, Channel_6,
          Channel_7, Channel_8, Channel_9,
          Channel_10, 
          Channel_11, Channel_12, Channel_13,
          Channel_14, Channel_15, Channel_16,
          Channel_17, Channel_18, Channel_19,
          Channel_20 )
# Convert channel variable to numeric
HbO$Channel = sapply( HbO$Channel, function(x) 
  as.numeric( strsplit( x, split = 'Channel_' )[[1]][2] ) )
# Add tidy version of labels for ROI
HbO$R_DLPFC = NULL
HbO$L_DLPFC = NULL
HbO$MPFC = NULL
HbO$R_VLPFC = NULL
HbO$L_VLPFC = NULL
HbO$X = NULL

# Channels 10, 15, 17, 18
HbO$ROI = 'R_DLPFC'
# Channels 1, 2, 5, 8
sel = HbO$Channel %in% c( 1, 2, 5, 8 )
HbO$ROI[ sel ] = 'L_DLPFC'
# Channels 7, 9, 12, 14
sel = HbO$Channel %in% c( 7, 9, 12, 14 )
HbO$ROI[ sel ] = 'MPFC'
# Channels 13, 16, 19, 20
sel = HbO$Channel %in% c( 13, 16, 19, 20 )
HbO$ROI[ sel ] = 'R_VLPFC'
# Channels 3, 4, 6, 11
sel = HbO$Channel %in% c( 3, 4, 6, 11 )
HbO$ROI[ sel ] = 'L_VLPFC'

# Create standarized estimates
HbO$zHbO = ( HbO$HbO - 
               mean( HbO$HbO, na.rm = T ) ) / 
  sd( HbO$HbO, na.rm = T )

# Compute new incremental subject values
HbO$Subject_old = HbO$Subject
HbO$Subject = createIncrement( HbO$Subject )

###
### 3) Plot of drug/placebo HbO per subject for PRE
###

if ( run_code[1] ) {
  
  # Indicate whether new plotting window should 
  # be generated or if figure should be saved as PDF
  if (!savePlot) x11(width = 12 )
  
  # Create a blank plot
  xl = c( 0, max( HbO$Subject ) + 1 )
  yl = c( -10, 10 )
  blankPlot( xl, yl )
  
  ds = list(
    m = mean( HbO$HbO[ HbO$Timepoints == 'T1_Pre_drug' ], 
              na.rm = T ),
    s = sd( HbO$HbO[ HbO$Timepoints == 'T1_Pre_drug' ], 
            na.rm = T ) )
  
  # Convenience function to select subjects 
  # who either had first placebo, then drug 
  # conditions or vice versa, and plot 
  # their data separately
  quick_draw = function( first, shft = .25, adj = 0,
                         clr = c( 'blue', 'red' ),
                         pts = c( 19, 19 ) ) {
    
    visit_val = 1:2
    if ( first == 'Drug' ) visit_val = 2:1
    
    shft = c( -shft, shft )
    cond = c( 'Placebo', 'Drug' )
    
    for ( i in 1:2 ) {
      
      cond_sel = HbO$Condition == cond[i] & 
        HbO$Visit == visit_val[i]
      
      sel = HbO$Condition == cond[i] & 
        HbO$Timepoints == 'T1_Pre_drug'
      
      y = ( HbO$HbO[ sel & cond_sel ] - ds$m )/ds$s
      x = createIncrement( HbO$Subject[ sel & cond_sel ] )
      
      id = unique( HbO$ID[ sel & cond_sel ] )
      id = sapply( id, function(x)
        strsplit( x, split = 'FN_0' )[[1]][2] )
      
      points( x + shft[i] + adj, y, pch = pts[i], col = clr[i] )
      
      if ( i == 1 ) {
        xy = aggregate( cbind( x, y ), 
                        list( HbO$ID[ cond_sel & sel ] ), max )
      } else {
        xy2 = aggregate( cbind( x, y ), 
                        list( HbO$ID[ cond_sel & sel ] ), max )
        
        xy[,3] = apply( cbind( xy[,3], xy2[,3] ), 1, max )
        text( xy[,2] + adj, xy[,3], id, pos = 3,
              cex = .5 )
      }
      
    }
    
    return( length( unique( x ) ) )
  }
  
  horizLines( qnorm( c(.005,.995) ), xl, 
              col = 'grey', lwd = 2 )
  horizLines( 0, xl, lwd = 2 )
  vertLines( 24.5, yl, lwd = 2 )
  adj = quick_draw( 'Placebo' )
  tmp = quick_draw( 'Drug', adj = adj )
  rm( tmp )
  
  customAxes( xl, yl )
  axis( 2, seq( -10, 10, 2 ), 
        tick = F, line = -1.75, cex.axis = 1.25 )
  mtext( 'z-scores for HbO', side = 2, cex = 1.5, line = 2 )
  
  axis( 1, c( 24/2, 24 + (54 - 24)/2 ),
        c( 'Placebo first', 'Drug first' ),
        tick = F, line = -1.25, cex.axis = 1.5 )
  mtext( 'Subjects', side = 1, cex = 1.5, line = 2 )
  
  legend( 'topright', c( 'Placebo', 'Drug' ),
          fill = c( 'blue', 'red' ),
          bty = 'n', cex = 1.5 )
  
}

###
### 4) 
###

if ( run_code[2] ) {
  
  
  
}

if ( savePlot ) dev.off()
setwd( orig_dir )

