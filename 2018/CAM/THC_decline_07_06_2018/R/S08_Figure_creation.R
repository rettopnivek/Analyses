# Figure creation
# Written by Kevin Potter
# email: kevin.w.potter@gmail.com
# Please email me directly if you 
# have any questions or comments
# Last updated 2018-12-05

# Table of contents
# 1) Initial setup
# 2) Define functions and variables
#   2.1) ln_m
#   2.2) ln_sd
#   2.3) pv
# 3) Figure 1
# 4) Figure 2
# 5) Figure S.1
# 6) Figure S.2

source( 'S01_Folder_paths.R' )

# Indicate code segment to run
run_code = c(
  T, # Figure 1
  T, # Figure 2
  T, # Figure S.1
  T  # Figure S.2
)

# Load in useful packages

# Collection of useful functions
my_package_load( 'utilityf', github = T )

# Package for working with data frames
my_package_load( 'dplyr' )

# Package for Bayesian estimation
Sys.setenv(USE_CXX14 = 1) # For Rstan to work
my_package_load( 'brms' )

# Load in packages for exponential decay model
source( 'S03_Exponential_decay_functions.R' )

# Load in data
setwd( dat_dir )
load( 'THC_decay.RData' )
# Load in data with imputted values
load( 'Imputted_data.RData' )

# If present, load in 
# best-fitting model 
# results
if ( 'Best_fitting_model.RData' %in% dir() ) {
  load( 'Best_fitting_model.RData' )
}

###
### 2) Define functions and variables
###

# 2.1) 
ln_m = function( mu, sigma ) {
  # Purpose:
  # Computes the mean for a log-normal distribution.
  # Arguments:
  # mu    - Mean parameter for log-normal distribution
  # sigma - Standard deviation parameter for log-normal 
  #         distribution
  # Returns:
  # The mean for a log-normal distribution.
  
  s2 = pow( sigma, 2 )
  out = exp( mu + s2/2 )
  return( out )
}

# 2.2) 
ln_sd = function( mu, sigma ) {
  # Purpose:
  # Computes the standard deviation for a 
  # log-normal distribution.
  # Arguments:
  # mu    - Mean parameter for log-normal distribution
  # sigma - Standard deviation parameter for log-normal 
  #         distribution
  # Returns:
  # The standard deviation for a log-normal distribution.
  
  s2 = pow( sigma, 2 )
  out = ( exp( s2 ) - 1 ) * exp( 2*mu + s2 )
  return( sqrt( out ) )
}

# 2.3) 
pv = function( x ) {
  # Purpose:
  # Function to compute posterior p-values from 
  # a MCMC sample.
  # Arguments:
  # x - A vector of values
  # Returns:
  # A posterior p-value.
  
  if ( mean(x) > 0 ) out = sum( x < 0 )/length(x)
  if ( mean(x) <= 0 ) out = sum( x > 0 )/length(x)
  
  return( out )
}

setwd( dat_dir )
load( file = 'Best_fitting_model.RData' )
setwd( R_dir )

# Extract posterior estimates for population-level 
# parameters and subject-level variability
post = as.matrix( bfm$model_fit )

# Extract columns for subject-level posteriors
pn = colnames( post )
sp_subj = grep( 'Start_point', pn )[ -(1:3) ]
er_subj = grep( 'Elimination_rate', pn )[ -(1:3) ]

# Subject-level elimination rates
er = colMeans( post[,'b_Elimination_rate'] + post[,er_subj] )
# Subject-level start points
sp = colMeans( post[,'b_Start_point'] + post[,sp_subj] )

###
### 3) Create figure 1
###

setwd( fig_dir )
if ( run_code[1] ) {
  
  png( 'Figure_1.png', width = 480*2, height = 480 )
  
  tmp = names( sp )
  tmp = strsplit( tmp, split = ',' )
  tmp = unlist( lapply( tmp, function(x) x[1] ) )
  tmp = as.character( 
    sapply( tmp, function(x) strsplit( x, 
                                       split = 'r_ID\\[' )[[1]][2] ) )
  check = data.frame(
    ID = tmp,
    log_SP = sp,
    SP = exp( sp ), 
    ER = er,
    stringsAsFactors = F
  )
  rownames( check ) = 1:nrow( check )
  check = check %>% arrange( SP )
  check$N = NA
  check$U = NA
  for ( s in 1:nrow( check ) ) {
    sel = dtbf$ID == check$ID[s]
    check$N[s] = sum(sel)
    check$U[s] = sum( dtbf$THCCOOH[sel] <= 1 )/check$N[s]
  }
  
  # Save model estimates of starting level
  tbl = data.frame(
    ID = check$ID,
    Starting_level = round( 
      ln_m( check$log_SP, mean( post[,'sd_ID__Start_point'] ) ) ),
    stringsAsFactors = F
  )
  setwd( proj_dir )
  setwd( 'Documents' )
  write.table( tbl, 
               row.names = F,
               quote = F,
               sep = ',',
               file = 'Starting_levels.csv'
  )
  setwd( fig_dir )
  
  # Select example subjects
  ex_subj = c(
    # Low starting levels
    '10055',
    '10022',
    # High starting levels
    '10066',
    '10020'
  )
  
  # x11( width = 12 )
  layout( cbind( 1, 2 ) )
  
  xl = c( 0, 35 )
  yl = lowerUpper( 50, dtbf$THCCOOH[ dtbf$ID %in% ex_subj ] )
  yl[1] = -5
  blankPlot( xl, yl )
  
  horizLines( seq( 0, yl[2], 50 ), xl, col = 'grey' )
  
  pts = c( 24, 22, 23, 25 )
  clr = rep( c( 'black', 'white' ), each = 2 )
  
  for ( i in 1:length(ex_subj) ) {
    
    sel = dtbf$ID == ex_subj[i]
    xa = dtbf$Time[ sel ]
    ya = dtbf$THCCOOH[ sel ]
    
    lines( xa, ya, lwd = 2 )
    points( xa, ya, pch = pts[i], bg = clr[i] )
    
  }
  
  # Axes and labels
  customAxes( xl, yl )
  
  # x-axis
  axis( 1, seq( 0, 32, 4 ),
        tick = F, line = -1.5, cex.axis = 1.25 )
  
  # y-axis
  axis( 2, round( seq( 0, yl[2], 50 ) ),
        tick = F, line = -1.5, cex.axis = 1.25 )
  mtext( 'CN-THCCOOH - ng/mg',
         side = 2, line = 2, cex = 1.25 )
  
  yl = lowerUpper( 2.5, dtbf$log_THCCOOH[ dtbf$ID %in% ex_subj ] )
  blankPlot( xl, yl )
  
  horizLines( seq( yl[1], yl[2], 2.5 ), xl, col = 'grey' )
  
  for ( i in 1:length(ex_subj) ) {
    
    sel = dtbf$ID == ex_subj[i]
    xa = dtbf$Time[ sel ]
    ya = dtbf$log_THCCOOH[ sel ]
    
    b("
      lmf = lm( ya ~ xa )
      cf = coef( lmf )
      segments( 0, cf[1],
      max( xa ), cf[1] + cf[2]*max( xa ),
      lwd = 2, col = 'grey' )
      ")
    
    lines( xa, ya, lwd = 2 )
    points( xa, ya, pch = pts[i], bg = clr[i] )
    
  }
  customAxes( xl, yl )
  
  # x-axis
  axis( 1, seq( 0, 32, 4 ),
        tick = F, line = -1.5, cex.axis = 1.25 )
  
  # y-axis
  axis( 2, seq( -7.5, yl[2], 2.5 ),
        tick = F, line = -1.5, cex.axis = 1.25 )
  mtext( 'CN-THCCOOH - log(ng/mg)',
         side = 2, line = 2, cex = 1.25 )
  
  mtext( 'Days Since Last Cannabis Exposure',
         side = 1, line = -2, cex = 1.25, outer = T )
  
  par( xpd = NA )
  legend( -30, 10,
          paste( 'ID =', ex_subj ), 
          pch = pts,
          pt.bg = clr,
          horiz = T,
          bty = 'n',
          cex = 1.25
  )
  par( xpd = F )
  
  dev.off()
  
  
}
setwd( R_dir )


###
### 4) Figure 2
###

setwd( fig_dir )
if ( run_code[2] ) {
  
  png( 'Figure_2.png' )
  
  sel = dtbf$THCCOOH < 1.2
  dtbf$Orig_THCCOOH = dtbf$THCCOOH
  dtbf$Orig_THCCOOH[sel] = 0
  pts = rep( 19, nrow(dtbf) )
  pts[sel] = 1
  
  # Mean over subjects
  plt = dtbf %>% 
    group_by( Time ) %>% 
    summarise(
      R = mean( THCCOOH )
    )
  
  # Obtain model predictions
  sm = summary( bfm$model_fit )
  xa = seq( 0, max( dtbf$Time ) )
  pred = exp( sm$fixed[1,1] - xa*sm$fixed[3,1] )
  
  ### Raw data
  
  # x11()
  
  x = dtbf$Time
  y = dtbf$THCCOOH
  
  # Create a blank plot
  xl = c( -1, 35 )
  yl = c( -7.5, 750 )
  blankPlot( xl, yl )
  
  horizLines( seq( 150, 750, 150 ), xl, col = 'grey80' )
  
  # Add observations
  points( x, y, pch = 19, col = 'grey' )
  
  # Add mean
  lines( plt$Time, plt$R, lwd = 2, col = 'grey40' )
  points( plt$Time, plt$R, pch = 17 )
  # Add estimates
  lines( xa, pred, lwd = 2 )
  
  # Add axes and labels
  customAxes( xl, yl )
  
  # x-axis
  axis( 1, seq( 0, 32, 4 ),
        tick = F, line = -1.5, cex.axis = 1.25 )
  mtext( 'Days Since Last Cannabis Exposure',
         side = 1, line = 2, cex = 1.25 )
  
  # y-axis
  axis( 2, round( seq( 0, yl[2], 150 ) ),
        tick = F, line = -1.5, cex.axis = 1.25 )
  mtext( 'CN-THCCOOH - ng/mg',
         side = 2, line = 2, cex = 1.25 )
  
  dev.off()
  
  
}
setwd( R_dir )

###
### 5) Figure S.1
###

setwd( fig_dir )
if ( run_code[3] ) {
  
  png( 'Figure_S1.png' )
  
  # x11()
  
  # Generating parameters
  gp = c(
    alpha = 100,
    kappa = .25,
    sigma = 1
  )
  
  # Simulate data
  sim_dat = data.frame(
    Time = seq( 0, 34, length = 100 ),
    epsilon = NA, 
    log_THCCOOH = NA,
    THCCOOH.obs = NA,
    log_THCCOOH.obs = NA
  )
  
  set.seed( 393 )
  sim_dat$epsilon = rnorm( nrow(sim_dat), 0, gp[3] )
  sim_dat$log_THCCOOH = ed_model( gp, sim_dat$Time, log = T )
  sim_dat$log_THCCOOH.obs = sim_dat$log_THCCOOH + sim_dat$epsilon
  sim_dat$THCCOOH.obs = exp( sim_dat$log_THCCOOH.obs )
  sim_dat$THCCOOH.obs.censored = sim_dat$THCCOOH.obs
  
  cut_off = 1
  
  sel = sim_dat$THCCOOH.obs.censored < cut_off
  sim_dat$THCCOOH.obs.censored[sel] = 0
  tr = estimate_ed_lm( sim_dat$Time, sim_dat$THCCOOH.obs.censored )
  
  xl = c( -1, 35 )
  yl = lowerUpper( 1, sim_dat$log_THCCOOH.obs )
  blankPlot( xl, yl )
  
  lines( sim_dat$Time, sim_dat$log_THCCOOH,
         lwd = 2, col = 'blue' )
  
  horizLines( log( cut_off ), xl, lwd = 2, lty = 2 )
  
  clr = rep( 'black', nrow( sim_dat ) )
  clr[ sim_dat$log_THCCOOH.obs < log( cut_off ) ] = 'grey'
  points( sim_dat$Time, sim_dat$log_THCCOOH.obs, 
          pch = 19, col = clr )
  
  sel = sim_dat$log_THCCOOH.obs > log( cut_off )
  lmf = lm( log_THCCOOH.obs ~ Time, data = sim_dat[sel,] )
  est = coef( lmf )
  lines( sim_dat$Time, est[1] + est[2]*sim_dat$Time, 
         col = 'red', lwd = 2 )
  lines( tr$x, log( tr$pred ), col = 'purple', lwd = 2 )
  
  customAxes( xl, yl )
  
  axis( 1, seq( 0, 32, 4 ),
        tick = F, line = -1.5, cex.axis = 1.25 )
  mtext( 'Days Since Last Cannabis Exposure', side = 1, 
         line = 2, cex = 1.25 )
  
  axis( 2, round( seq( yl[1], yl[2], length = 5 ), 1 ), 
        tick = F, line = -1.5, cex.axis = 1.25 )
  mtext( 'CN-THCCOOH - ln( ng/mg)', side = 2, 
         line = 2, cex = 1.25 )
  
  legend( 'topright',
          c( 'Observed', 'Censored', 'Generating', 'Estimated', 'Tobit' ),
          fill = c( 'black', 'grey', 'blue', 'red', 'purple' ),
          cex = 1.25, bty = 'n' )
  
  dev.off()
}
setwd( R_dir )

###
### Figure S.2
###

setwd( fig_dir )
if ( run_code[4] ) {
  
  png( 'Figure_S2.png' )
  
  # x11()
  
  x = dtbf$Time
  y = dtbf$THCCOOH
  
  # Create a blank plot
  xl = c( -1, 35 )
  yl = c( -26, 10 )
  blankPlot( xl, yl )
  
  horizLines( seq( -25, 5, 5 ), xl, col = 'grey80' )
  horizLines( log(1), xl, lwd = 2, lty = 2 )
  
  # Add observations
  clr = rep( 'black', length( y ) )
  clr[ y <= 1 ] = 'grey'
  points( x, log( y ), pch = 19, col = clr )
  
  # 
  tmp = data.frame(
    x = x,
    y = log( y )
  )
  tmp$y[ tmp$y <= 0 ] = NA
  plt_2 = tmp %>% 
    group_by( x ) %>% 
    summarize(
      M = mean( y, na.rm = T )
    )
  
  # Add mean
  lines( plt$Time, log( plt$R ), lwd = 2 )
  lines( plt$Time, plt_2$M, lwd = 2, col = 'red' )
  # Add estimates
  lines( xa, log( pred ), col = 'blue', lwd = 2 )
  
  # Add axes and labels
  customAxes( xl, yl )
  
  # x-axis
  axis( 1, seq( 0, 32, 4 ),
        tick = F, line = -1.5, cex.axis = 1.25 )
  mtext( 'Days Since Last Cannabis Exposure',
         side = 1, line = 2, cex = 1.25 )
  
  # y-axis
  axis( 2, round( seq( -25, 10, 5 ) ),
        tick = F, line = -1.5, cex.axis = 1.25 )
  mtext( 'CN-THCCOOH - ln(ng/mg)',
         side = 2, line = 2, cex = 1.25 )
  
  # Legend
  legend( 20, 13,  
          c( 'Mean', 'Mean - censored', 'PK model' ),
          fill = c( 'black', 'red', 'blue' ),
          cex = 1.25, bty = 'n', xpd = T )
  
  dev.off()
}
setwd( R_dir )


