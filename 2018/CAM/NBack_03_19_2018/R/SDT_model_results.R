# Script to display model estimation results
# Written by Kevin Potter
# email: kevin.w.potter@gmail.com
# Please email me directory if you 
# have any questions or comments
# Last updated 2018-04-20

# Table of contents
# 1) Initial setup
# 2) Extract model estimation results
# 3) Convergence check
# 4) Posterior predictive checks for mean hits/correct omissions
# 5) Marginal posteriors (Effect sizes for d')
#   5.1) quick_ci
#   5.2) plt_ci
# 6) Marginal posteriors (Effect sizes for criterion)
# 7) Marginal posterior (Response bias per condition)

###
### 1) Initial setup
###

# Select model results to extract
model = 'm3'

# Save current working directory
orig_dir = getwd()

# Navigate to project directory
setwd( '..' )
proj_dir = getwd()

# Indicate which code segments to run
run_code = c(
  F, # Convergence checks
  T, # Posterior predictive checks
  T, # Marginal posteriors (d')
  T, # Marginal posteriors (criterion)
  F, # Marginal posteriors (response bias)
  F
)

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

# Indicate whether to save figures in a PDF file
savePlot = T
if ( savePlot ) {
  setwd( 'Figures' )
  pdf( 'MR_script_output.pdf', width = 12 )
  setwd( proj_dir )
  
  # Create a table of contents for 
  # figures
  
  tc_index = run_code
  
  # Initialize table of contents
  tbl_cnt = c()
  
  # run_code[1]
  tbl_cnt = c( tbl_cnt,
               'Estimation convergence diagnostics' )
  
  # run_code[2]
  tbl_cnt = c( tbl_cnt,
               'Posterior retrodictive check' )
  
  # run_code[3]
  tbl_cnt = c( tbl_cnt,
               "Marginal posteriors (d')" )
  
  # run_code[4]
  tbl_cnt = c( tbl_cnt,
               "Marginal posteriors (criterion)" )
  
  # run_code[5]
  tbl_cnt = c( tbl_cnt,
               "Marginal posteriors (response bias)" )
  
  # Generate table of contents
  tableContents( tbl_cnt[ tc_index ] )
  
}

# Load in estimation functions 
# and compile SDT model script
setwd( 'R' )
source( 'Estimation_functions.R' )

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

# Create a new data frame that
# a) Excludes the combined trials
# b) Has separate variables for hits and false alarms
all_ptbp = dat %>% 
  filter( Task != 'Combined' & Response_type == 'Hits' ) %>% 
  arrange( Task, Condition, Timepoints, Subject ) %>% 
  mutate( H = Counts / Trials,
          fH = Counts, 
          Npos = Trials )
tmp = dat %>% 
  filter( Task != 'Combined' & Response_type == 'False_alarms' ) %>% 
  arrange( Task, Condition, Timepoints, Subject ) %>% 
  mutate( FA = Counts / Trials,
          fFA = Counts,
          Nneg = Trials )
all_ptbp$FA = tmp$FA
all_ptbp$fFA = tmp$fFA
all_ptbp$Nneg = tmp$Nneg

# Clean up workspace
rm( tmp )

# Convert the hit/false alarm rates into 
# estimates of d' and bias
all_ptbp = all_ptbp %>% 
  mutate( dp = sdt_calc_binary( fH, Npos, fFA, Nneg, 1, T ),
          crt = sdt_calc_binary( fH, Npos, fFA, Nneg, 1, F ) )

# Standardize BOLD estimates

# Variable names for ROI
roi = c( 'R_DLPFC',
         'L_DLPFC',
         'MPFC',
         'R_VLPFC',
         'L_VLPFC' )

# Flip sign of BOLD estimates
all_ptbp[,roi] = -all_ptbp[,roi]

# Mean and standard deviation over all estimates
BOLD_ds = c(
  m = mean( as.vector( unlist( all_ptbp[,roi] ) ), na.rm = T ),
  sd = sd( as.vector( unlist( all_ptbp[,roi] ) ), na.rm = T )
)

# Standardize estimates
all_ptbp = all_ptbp %>% 
  mutate(
    R_DLPFC = ( R_DLPFC - BOLD_ds[1] )/BOLD_ds[2],
    L_DLPFC = ( L_DLPFC - BOLD_ds[1] )/BOLD_ds[2],
    MPFC = ( MPFC - BOLD_ds[1] )/BOLD_ds[2],
    R_VLPFC = ( R_VLPFC - BOLD_ds[1] )/BOLD_ds[2],
    L_VLPFC = ( L_VLPFC - BOLD_ds[1] )/BOLD_ds[2]
  )

# Add in correct rejections and 
# misses
all_ptbp$CR = 1 - all_ptbp$FA
all_ptbp$fCR = all_ptbp$Nneg - all_ptbp$fFA
all_ptbp$M = 1 - all_ptbp$H
all_ptbp$fM = all_ptbp$Npos - all_ptbp$fH

###
### 2) Extract model estimation results
###

# Navigate to folder with posterior estimates
# and posterior predictive checks
setwd( 'Data/Posterior_estimates' )

# Load in posterior estimates
fname = paste( 'Posterior_', model, '.RData', sep = '' )
load( fname )

# Standardize variable names
if ( model == 'm1' ) {
  mdl = m1
  X = X1
}

if ( model == 'm2' ) {
  mdl = m2
  X = X2
}

if ( model == 'm3' ) {
  mdl = m3
  X = X3
}

# Load in posterior predictive checks
fname = paste( 'PPC_', model, '.RData', sep = '' )
load( fname )

setwd( proj_dir )

###
### 3) Convergence check
###

if ( run_code[1] )
  plotConvergence( mdl$conv, savePlot )

###
### 4) Posterior predictive checks for mean hits/correct omissions
###

if ( run_code[2] ) {
  
  # Statistics for boxplots
  ct = all_ptbp %>% 
    group_by( Task, Condition, Timepoints ) %>% 
    summarize( mean_H = mean( H ),
               mean_CR = mean( CR ),
               median_H = median( H ),
               median_CR = median( CR ) )
  
  # Indicate whether new plotting window should 
  # be generated or if figure should be saved as PDF
  if (!savePlot) x11(width = 12 )
  
  # Plotting characteristics
  lnSz = 2
  lnSz2 = 3
  txtSz = 1.5
  ptSz = 1.25
  axPos = -1.75
  pts = rep( rep( c(22,24), each = 3 ), 4 )
  clr = rep( rep( c( 'black', 'white' ), each = 3 ), 4 )
  
  # Extract elements to plot
  mean_val = numeric(24)
  median_val = numeric(24)
  error_bars = matrix( NA, 4, 24 )
  lbls = matrix( ' ', 4, 2 )
  for ( i in 1:4 ) {
    v1 = 1:3 + 6*(i-1)
    v2 = 4:6 + 6*(i-1)
    v3 = 1:3 + 3*(i-1)
    mean_val[v1] = ct$mean_H[v3]
    mean_val[v2] = ct$mean_CR[v3]
    median_val[v1] = ct$median_H[v3]
    median_val[v2] = ct$median_CR[v3]
    error_bars[1,v1] = pred[v3,'CI_2.5','H']
    error_bars[1,v2] = 1 - pred[v3,'CI_2.5','FA']
    error_bars[2,v1] = pred[v3,'CI_97.5','H']
    error_bars[2,v2] = 1 - pred[v3,'CI_97.5','FA']
    error_bars[3,v1] = pred[v3,'CI_16','H']
    error_bars[3,v2] = 1 - pred[v3,'CI_16','FA']
    error_bars[4,v1] = pred[v3,'CI_84','H']
    error_bars[4,v2] = 1 - pred[v3,'CI_84','FA']
    lbls[i,1] = unique( ct$Condition[v3] )
    lbls[i,2] = unique( ct$Task[v3] )
  }
  lbls[ lbls[,2] == 'NBack_0', 2 ] = '0-back'
  lbls[ lbls[,2] == 'NBack_2', 2 ] = '2-back'
  
  # Create a blank plot
  xl = c( .5, 24.5 ); yl = c( .7, 1 )
  blankPlot( xl, yl )
  
  # Add grid lines
  horizLines( seq( .7, 1, .05 ), xl, col = 'grey80', 
              lwd = lnSz )
  customAxes( xl, yl )
  # Separate task and conditions
  vertLines( c( 6.5, 12.5, 18.5 ), yl, lwd = lnSz )
  
  # Add axis labels
  axis( 3, 24/4 * 1:4 - 2.5, 
        lbls[,1],
        tick = F, line = axPos, cex.axis = txtSz )
  axis( 1, 24/2 * 1:2 - 5.5, 
        unique( lbls[,2] ),
        tick = F, line = axPos, cex.axis = txtSz )
  axis( 2, seq( .7, 1, .1 ),
        seq( .7, 1, .1 ) * 100, 
        tick = F, line = axPos, cex.axis = txtSz )
  mtext( 'Average percentage', side = 2, line = 1.5, cex = txtSz )
  
  # Add plotting elements
  
  segments( 1:24, error_bars[1,],
            1:24, error_bars[2,], 
            lwd = lnSz, col = 'grey50' )
  errorBars( 1:24, error_bars[3:4,], 
             lwd = lnSz, length = .05, col = 'grey50' )
  
  # Lines by condition
  for ( i in 1:8 ) {
    v1 = 1:3 + 3*(i-1)
    lines( v1, mean_val[v1], lwd = lnSz )
  }
  
  # Add average percentages
  points( 1:24, mean_val, pch = pts, bg = clr )
  # points( 1:24, median_val, pch = 19, col = 'grey' )
  
  # Add labels for timepoints
  text( 1:3, rep( .69, 3 ), c( expression(T[1]),
                                   expression(T[2]),
                                   expression(T[3]) ),
        pos = 1, cex = txtSz )
  
  legend( -.5, yl[1] - diff(yl)*.1, 
          c( 'Correct Omissions', 'Hits' ),
          pch = rev( unique( pts ) ), 
          pt.bg = c( 'white', 'black' ),
          bty = 'n', horiz = T, xpd = T, cex = txtSz )
  
  legend( 7.5, yl[1] - diff(yl)*.1, 
          c( 'T1: pre-drug', 'T2: post-drug', 'T3: post-drug (2)' ),
          bty = 'n', horiz = T, xpd = T, cex = txtSz )
  
  legend( 7.5, yl[2] + diff(yl)*.2,
          c( 'Observed', 'Predicted' ),
          fill = c( 'black', 'grey50' ), bty = 'n',
          horiz = T, xpd = T, cex = txtSz )
  
  # Indicate which model results are being plotted
  # if ( model == 'm1' ) 
  #   ttl = c( 'Model 1 (Null)' )
  # if ( model == 'm2' ) 
  #   ttl = c( 'Model 2 (Effect of drug)' )
  # if ( model == 'm3' ) 
  #   ttl = c( 'Model 3 (Drug x intoxication)' )
  
  # mtext( ttl, side = 3, line = 1, cex = txtSz )
  
}

###
### 5) Marginal posteriors (Effect sizes for d')
###

# 5.1)
quick_ci = function(x) {
  # Purpose:
  # A function to compute metrics of choice 
  # over posterior samples.
  # Arguments:
  # x - A vector of posterior samples
  # Returns:
  # A vector with the 68% and 95% credible intervals 
  # and the posterior mean and median.
  
  p = c( .025, .16, .5, .84, .975 )
  q = quantile( x, p )
  names( q ) = c(
    'CI_2.5',
    'CI_16',
    'Median',
    'CI_84',
    'CI_97.5' )
  out = c(
    Mean = mean( x ),
    q )
  return( out )
}

# 5.2)
plt_ci = function( pos, ci, clr = 'black' ) {
  # Purpose:
  # A function to quickly plot credible intervals 
  # for posterior samples.
  # Arguments:
  # pos - The x-axis position
  # ci  - The vector with the information for the 
  #       marginal posterior
  # clr - The line color

  segments( pos, ci['CI_2.5'],
            pos, ci['CI_97.5'],
            lwd = lnSz, col = clr )
  errorBars( pos, ci[c('CI_16','CI_84')],
             lwd = lnSz, length = .05, col = clr )
  points( pos, ci['Mean'], pch = 19, cex = ptSz, col = clr )
  points( pos, ci['Median'], pch = '-', cex = ptSz, col = clr )
  
}


if ( run_code[3] ) {
  
  # Details regarding posterior samples
  dps = 1:ncol( X[[1]] )
  cs = 1:ncol( X[[2]] ) + ncol( X[[1]] )
  
  # Parameter labels
  dp_n = colnames( X[[1]] )
  c_n = colnames( X[[2]] )
  
  ptbp = mdl$post$group_param[,dps]
  colnames( ptbp ) = dp_n
  ptbp = ptbp[, -(grep('Baseline', dp_n))]
  cip = apply( ptbp, 2, quick_ci )
  
  # Indicate whether new plotting window should 
  # be generated or if figure should be saved as PDF
  if (!savePlot) x11(width = 12 )
  par( mar = c( 7, 4, 3, 2 ) )
  
  # Plotting characteristics
  lnSz = 2
  txtSz = 1.5
  ptSz = 1.25
  axPos = -1.75
  
  yl = lowerUpper( .5, as.vector( cip ) )
  xl = c( .5, .5 + ncol(cip) )
  blankPlot( xl, yl )
  horizLines( seq( yl[1], yl[2], .5 ), 
              xl, lwd = lnSz, col = 'grey80' )
  horizLines( 0, xl, lwd = lnSz )
  customAxes( xl, yl, lnSz = lnSz )
  
  axis( 2, seq( yl[1], yl[2], .5 ), 
        tick = F, line = -1.75, cex.axis = txtSz )
  mtext( "Estimated effect size (d')", 
         side = 2, line = 1.5, cex = txtSz )
  
  # Identify statistically significant 
  # effects
  ssp = cip[ 'CI_2.5', ] > 0 | cip[ 'CI_97.5', ] < 0
  clr = rep( 'black', ncol( cip ) )
  clr[ ssp ] = 'blue'
  b("
  if ( model %in% c( 'm2', 'm3' ) ) {
    clr[ grep( 'Visit', colnames( cip ) ) ] = 'blue'
  }
  if ( model == 'm3' ) {
    clr[ grep( 'Self', colnames( cip ) ) ] = 'blue'
  }")
  
  for ( i in 1:ncol(cip) ) plt_ci( i, cip[,i], clr = clr[i] )
  
  if ( model == 'm1' ) {
    lbl = c(
      'POST (0-back)',
      'POST 2 (0-back)',
      'POST (2-back)',
      'POST 2 (2-back)',
      'Visit order effect' )
  }
  
  if ( model == 'm2' ) {
    lbl = c(
      'PRE (Placebo, 0-back)',
      'POST (Placebo, 0-back)',
      'POST 2 (Placebo, 0-back)',
      'POST (Drug, 0-back)',
      'POST 2 (Drug, 0-back)',
      'PRE (Placebo, 2-back)',
      'POST (Placebo, 2-back)',
      'POST 2 (Placebo, 2-back)',
      'POST (Drug, 2-back)',
      'POST 2 (Drug, 2-back)',
      'Visit order effect' )
  }
  
  if ( model == 'm3' ) {
    lbl = c(
      'PRE (Placebo, 0-back)',
      'POST (Placebo, 0-back)',
      'POST 2 (Placebo, 0-back)',
      'POST (Drug, 0-back)',
      'POST 2 (Drug, 0-back)',
      'PRE (Placebo, 2-back)',
      'POST (Placebo, 2-back)',
      'POST 2 (Placebo, 2-back)',
      'POST (Drug, 2-back)',
      'POST 2 (Drug, 2-back)',
      'Visit order effect',
      'Intoxication effect')
  }
  
  xa = 1:ncol( cip )
  axis( 1, at = xa, tick = F, labels = FALSE )
  text( x = xa, 
        y = par()$usr[3]+0.0*(par()$usr[4]-par()$usr[3]),
       labels = lbl, srt = 45, adj = 1, xpd = TRUE )
  
  dl = max( cumsum( sapply( lbl, 
                            function(x) 
                              length( grep( '0-back', x ) ) ) ) )
  dl = c( dl,
          which( colnames( cip ) == 'Visit_2' ) - 1 )
  vertLines( dl + .5, yl, lwd = lnSz, col = c( 'grey', 'black' ) )
  axis( 3, dl[1], 'Difference from reference: PRE (Drug)',
        line = 0, tick = F, cex.axis = txtSz )
  
  axis( 4, c( -.5, .5 ),
        c( 'Worse performance',
           ' Better performance' ),
        tick = F, line = -1 )
  
  # Reset plotting margins
  par( mar = c( 5, 4, 3, .5 ) )
  
}

###
### 6) Marginal posteriors (Effect sizes for criterion)
###

if ( run_code[4] ) {
  
  # Details regarding posterior samples
  dps = 1:ncol( X[[1]] )
  cs = 1:ncol( X[[2]] ) + ncol( X[[1]] )
  
  # Parameter labels
  dp_n = colnames( X[[1]] )
  c_n = colnames( X[[2]] )
  
  ptbp = mdl$post$group_param[,cs]
  colnames( ptbp ) = c_n
  cip = apply( ptbp, 2, quick_ci )
  
  # Indicate whether new plotting window should 
  # be generated or if figure should be saved as PDF
  if (!savePlot) x11(width = 12 )
  par( mar = c( 7, 4, 3, 2 ) )
  
  # Plotting characteristics
  lnSz = 2
  txtSz = 1.5
  ptSz = 1.25
  axPos = -1.75
  
  yl = lowerUpper( .5, as.vector( cip ) )
  xl = c( .5, .5 + ncol(cip) )
  blankPlot( xl, yl )
  horizLines( seq( yl[1], yl[2], .5 ), 
              xl, lwd = lnSz, col = 'grey80' )
  horizLines( 0, xl, lwd = lnSz )
  customAxes( xl, yl, lnSz = lnSz )
  
  axis( 2, seq( yl[1], yl[2], .5 ), 
        tick = F, line = -1.75, cex.axis = txtSz )
  mtext( "Estimated effect size (bias)", 
         side = 2, line = 1.5, cex = txtSz )
  
  # Identify significant effects
  ssp = cip[ 'CI_2.5', ] > 0 | cip[ 'CI_97.5', ] < 0
  clr = rep( 'black', ncol( cip ) )
  clr[ ssp ] = 'blue'
  b("
  if ( model %in% c( 'm2', 'm3' ) ) {
    clr[ grep( 'Visit', colnames( cip ) ) ] = 'blue'
  }
  if ( model == 'm3' ) {
    clr[ grep( 'Self', colnames( cip ) ) ] = 'blue'
  }")
  
  for ( i in 1:ncol(cip) ) plt_ci( i, cip[,i], clr = clr[i] )
  
  if ( model == 'm1' ) {
    lbl = c(
      'Bias (0-back)',
      'POST (0-back)',
      'POST 2 (0-back)',
      'Bias (2-back)',
      'POST (2-back)',
      'POST 2 (2-back)',
      'Visit order effect' )
  }
  
  if ( model == 'm2' ) {
    lbl = c(
      'Bias (0-back)',
      'PRE (Placebo, 0-back)',
      'POST (Placebo, 0-back)',
      'POST 2 (Placebo, 0-back)',
      'POST (Drug, 0-back)',
      'POST 2 (Drug, 0-back)',
      'Bias (2-back)',
      'PRE (Placebo, 2-back)',
      'POST (Placebo, 2-back)',
      'POST 2 (Placebo, 2-back)',
      'POST (Drug, 2-back)',
      'POST 2 (Drug, 2-back)',
      'Visit order effect' )
  }
  
  if ( model == 'm3' ) {
    lbl = c(
      'Bias (0-back)',
      'PRE (Placebo, 0-back)',
      'POST (Placebo, 0-back)',
      'POST 2 (Placebo, 0-back)',
      'POST (Drug, 0-back)',
      'POST 2 (Drug, 0-back)',
      'Bias (2-back)',
      'PRE (Placebo, 2-back)',
      'POST (Placebo, 2-back)',
      'POST 2 (Placebo, 2-back)',
      'POST (Drug, 2-back)',
      'POST 2 (Drug, 2-back)',
      'Visit order effect',
      'Intoxication effect')
  }
  
  xa = 1:ncol( cip )
  axis( 1, at = xa, tick = F, labels = FALSE )
  text( x = xa, 
        y = par()$usr[3]+0.0*(par()$usr[4]-par()$usr[3]),
        labels = lbl, srt = 45, adj = 1, xpd = TRUE )
  
  dl = max( cumsum( sapply( lbl, 
                            function(x) 
                              length( grep( '0-back', x ) ) ) ) )
  dl = c( dl,
          which( colnames( cip ) == 'Visit_2' ) - 1 )
  vertLines( dl + .5, yl, lwd = lnSz, col = c( 'grey', 'black' ) )
  
  axis( 4, c( -.25, .25 ),
        c( ' Biased toward',
           'Biased against' ),
        tick = F, line = -1 )
  
}

###
### 7) Marginal posterior (Response bias per condition)
###

if ( run_code[5] ) {
  
  # Extract small version of design matrix
  dtbf = dat %>% 
    filter( Task != 'Combined' )
  dtbf = dtbf %>% 
    filter( Self_report_on_high != '-' )
  check = create_design_mat( dtbf, X[[2]] )
  
  # Details regarding posterior samples
  dps = 1:ncol( X[[1]] )
  cs = 1:ncol( X[[2]] ) + ncol( X[[1]] )
  
  # Parameter labels
  dp_n = colnames( X[[1]] )
  c_n = colnames( X[[2]] )
  
  # Extract posterior samples
  ptbp = mdl$post$group_param[,cs]
  # Convert into estimates of response bias 
  # for the 12 conditions
  ptbp = ptbp %*% t( as.matrix( check$X ) )
  cip = apply( ptbp, 2, quick_ci )
  
  # Indicate whether new plotting window should 
  # be generated or if figure should be saved as PDF
  if (!savePlot) x11(width = 12 )
  par( mar = c( 7, 5, 3, .5 ) )
  
  # Plotting characteristics
  lnSz = 2
  txtSz = 1.5
  ptSz = 1.25
  axPos = -1.75
  
  yl = c( -1, 1 )
  xl = c( .5, .5 + ncol(cip) )
  blankPlot( xl, yl )
  horizLines( seq( yl[1], yl[2], .5 ), 
              xl, lwd = lnSz, col = 'grey80' )
  horizLines( 0, xl, lwd = lnSz )
  customAxes( xl, yl, lnSz = lnSz )
  
  axis( 2, seq( yl[1], yl[2], .5 ), 
        tick = F, line = -1.75, cex.axis = txtSz )
  axis( 2, c( -.5,.5 ), 
        c('Biased toward', 'Biased against' ),
        tick = F, line = 1, cex.axis = txtSz )
  mtext( 'Response bias', side = 2, line = 3.5, 
         cex = txtSz )
  
  for ( i in 1:ncol(cip) ) plt_ci( i, cip[,i], clr = clr[i] )
  
  # Add condition labels
  xa = 1:12
  axis( 1, at = xa, tick = F, labels = FALSE )
  lbl = paste(
    rep( c( 'PRE ', 'POST ', 'POST 2 ' ), 4 ),
    '(', 
    check$Condition, ', ',
    rep( c( '0-back)', '2-back)' ), each = 6 ), sep = '' )
  text( x = xa, 
        y = par()$usr[3]+0.0*(par()$usr[4]-par()$usr[3]),
        labels = lbl, srt = 45, adj = 1, xpd = TRUE )
  
}

if ( savePlot ) dev.off()
setwd( orig_dir )