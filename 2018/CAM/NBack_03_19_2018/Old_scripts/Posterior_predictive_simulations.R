# Posterior predictive simulations
# Written by Kevin Potter
# email: kevin.w.potter@gmail.com
# Please email me directory if you 
# have any questions or comments
# Last updated 2018-04-19

# Table of contents
# 1) Initial setup
# 2) Extract model estimation results
# 3) Group-level predictions
# 4) Predictions over self-reported highs
#   4.1) create_SR_PPC_input
# 5) Save results

# Clear workspace
rm(list = ls())

# Select model results to extract
model = 'm1'

# Indicate whether to append path analysis tagline
PA = T

# Save current working directory
orig_dir = getwd()

# Navigate to project directory
setwd( '..' )
proj_dir = getwd()

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

if ( !PA ) {
  
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
  
} else {
  
  # Focus on 2-back task
  dtbf = dat %>% 
    filter( Task == 'NBack_2' )
  
  # Exclude rows with NA values for predictors
  init_X = create_design_mat( dtbf )
  dtbf = dtbf[ !is.na( init_X$R_DLPFC ) & !is.na( init_X$SR ), ]
  
  # Subject FN_095 has no neural data
  # Therefore, it is necessary to 
  # redefine subject values to be in 
  # sequential order
  dtbf$Subject_old = dtbf$Subject
  dtbf$Subject = createIncrement( dtbf$Subject )
  # The new assignments are
  # FN_096 = 63 (Originally 64)
  # FN_097 = 64 (Originally 65)
  # FN_041 = 65 (Originally 66)
  
}

# Create a new data frame that
# a) Excludes the combined trials
# b) Has separate variables for hits and false alarms
all_ptbp = dtbf %>% 
  filter( Task != 'Combined' & Response_type == 'Hits' ) %>% 
  arrange( Task, Condition, Timepoints, Subject ) %>% 
  mutate( H = Counts / Trials,
          fH = Counts, 
          Npos = Trials )
tmp = dtbf %>% 
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

###
### 2) Extract model estimation results
###

# Navigate to folder with posterior estimates
setwd( 'Data/Posterior_estimates' )
if ( !PA ) {
  fname = paste( 'Posterior_', model, '.RData', sep = '' )
  load( fname )
} else {
  fname = paste( 'Posterior_', model, '_PA.RData', sep = '' )
  load( fname )
}
setwd( proj_dir )

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

###
### 3) Group-level predictions
###

# Generate posterior simulations for hit/false alarm rates
sim = sdt_mix_rng( mdl$post, 
                   nrow( dtbf ),
                   max( dtbf$Subject ),
                   dtbf$Trials,
                   dtbf$Target,
                   dtbf$Subject,
                   X,
                   mdl$random )

# Initialize array to store statistics of interest
if ( !PA ) {
  ppc = array( NA, dim = c( 12, 4, ncol( sim ) ) )
} else {
  ppc = array( NA, dim = c( 6, 4, ncol( sim ) ) )
}

# Create a temporary data frame to store aggregated 
# data
tmp_df = all_ptbp %>% 
  group_by( Task, Condition, Timepoints ) %>% 
  summarize( H = mean( H ), 
             FA = mean( FA ), 
             Npos = round( mean( Npos ) ), 
             Nneg = round( mean( Nneg ) ) )
tmp_df = as.data.frame( tmp_df )
# Create list of variables to aggregate over
vrb = list( Timepoints = dtbf$Timepoints,
            Condition = dtbf$Condition,
            Task = dtbf$Task,
            Response_type = dtbf$Response_type )
sel_H = rep( c(F,T), each = nrow( tmp_df ) )
out = matrix( NA, nrow( tmp_df ), 4 )

# Function to aggregate over conditions 
# and extract hits, false alarms, d',
# and criterion values
fta = function( x ) {
  
  tmp = aggregate( x, vrb, mean )
  tmp_df$H = tmp$x[sel_H]
  tmp_df$FA = tmp$x[!sel_H]
  out[,1] = tmp_df$H
  out[,2] = tmp_df$FA
  out[,3] = sdt_calc_binary( round( tmp_df$H * 
                                      tmp_df$Npos ),
                             tmp_df$Npos,
                             round( tmp_df$FA * 
                                      tmp_df$Nneg ),
                             tmp_df$Nneg,
                             correct = 1,
                             dp = T )
  out[,4] = sdt_calc_binary( round( tmp_df$H * 
                                      tmp_df$Npos ),
                             tmp_df$Npos,
                             round( tmp_df$FA * 
                                      tmp_df$Nneg ),
                             tmp_df$Nneg,
                             correct = 1,
                             dp = F )
  
  return( out )
}

# Create a progress bar using a base R function
pb = txtProgressBar( min = 1, max = ncol( sim ), style = 3 )

# Generate posterior predictive checks 
# over conditions and samples
for ( i in 1:ncol( sim ) ) {
  ppc[,,i] = fta( sim[,i] )
  # Update the progress bar
  setTxtProgressBar(pb,i)
}
close(pb)

pred = array( NA, dim = c( nrow( tmp_df ), 5, 4 ) )
for ( i in 1:4 ) {
  pred[,,i] = t( apply( ppc[,i,], 1, function(x) {
    q = quantile( x, c( .025, .16, .84, .975 ) );
    names( q ) = c(
      'CI_2.5',
      'CI_16',
      'CI_84',
      'CI_97.5' );
    out = c(Mean = mean(x), q );
    return( out )
  } ) )
}
dimnames(pred)[[2]] = c('Mean','CI_2.5','CI_16',
                        'CI_84','CI_97.5' )
dimnames(pred)[[3]] = c( 'H', 'FA', 'dp', 'crt' )

###
### 4) Predictions over self-reported highs
###

# 4.1)
create_SR_PPC_input = function( Cond, mdl, X ) {
  # Purpose:
  # A convenience function to generate the 
  # relevant posterior simulations over 
  # the self-report metrics.
  # Arguments:
  # Cond - A list with the desired response 
  #        type, task, condition, and post-drug 
  #        time point to include.
  # mdl  - A list with the posterior output
  # X    - A list with the two design matrices
  # Returns:
  # An array with the average of the mean difference
  # and the lower and upper 95% prediction intervals 
  # over posterior simulations.
  
  # Cond = list( 'NBack_0', 'Placebo', 'T2_Post_drug' )
  
  # Timepoints
  pre = 'T1_Pre_drug'
  po = Cond[[3]]
  
  # Isolate conditions of interest
  sel = dtbf$Task == Cond[[1]] & 
    dtbf$Condition == Cond[[2]] & 
    ( dtbf$Timepoints == pre | 
        dtbf$ Timepoints == po )
  
  # Pull out design matrices
  Xd = X[[1]][sel,]
  Xc = X[[2]][sel,]
  
  # Extract data to be plotted
  dtbp = dtbf %>% 
    filter( Task == Cond[[1]] & 
              Condition == Cond[[2]] & 
              ( Timepoints == pre | 
                  Timepoints == po )
    )
  
  # Number of subjects
  Ns = max( dtbf$Subject )
  # Number of points for self-reported high
  SR = seq( 0, 1, .1 )
  # Number of group parameters
  Nd = ncol( X[[1]] )
  Nc = ncol( X[[2]] )
  # Parameter labels
  dp_n = colnames( X[[1]] )
  c_n = colnames( X[[2]] )
  
  # Create simplified design matrices
  Xdsp = aggregate( Xd, list(
    dtbp$Timepoints ), mean )
  Xcsp = aggregate( Xc, list(
    dtbp$Timepoints ), mean )
  
  if ( 'Self_report_high' %in% dp_n ) 
    Xdsp[,'Self_report_high'] = c(0,1)
  if ( 'Self_report_high' %in% c_n ) 
    Xcsp[,'Self_report_high'] = c(0,1)
  
  cnd = data.frame(
    rep( Xdsp[,1], each = Ns * length(SR) ) )
  
  Xds = matrix( as.numeric( Xdsp[1,-1] ), 
                Ns * length(SR), 
                ncol(X[[1]]), byrow = T )
  Xds = rbind( Xds, 
               matrix( as.numeric( Xdsp[2,-1] ), 
                       Ns * length(SR), 
                       ncol(X[[1]]), byrow = T ) )
  Xcs = matrix( as.numeric( Xcsp[1,-1] ), 
                Ns * length(SR), 
                ncol(X[[2]]), byrow = T )
  Xcs = rbind( Xcs, 
               matrix( as.numeric( Xcsp[2,-1] ), 
                       Ns * length(SR), 
                       ncol(X[[2]]), byrow = T ) )
  colnames( Xds ) = colnames( X[[1]] )
  colnames( Xcs ) = colnames( X[[2]] )
  
  if ( 'Self_report_high' %in% dp_n ) {
    Xds[,'Self_report_high'] = 
      Xds[,'Self_report_high'] * rep( SR, Ns * 2 )
  }
  
  if ( 'Self_report_high' %in% c_n ) {
    Xcs[,'Self_report_high'] = 
      Xcs[,'Self_report_high'] * rep( SR, Ns * 2 )
  }
  
  X_sr = list( rbind( Xds, Xds ), 
               rbind( Xcs, Xcs ) )
  
  Nt = dtbp %>% 
    group_by( Timepoints ) %>% 
    summarize( Nt = round( mean( Trials ) ) )
  Nt = rep( Nt$Nt, each = Ns * length( SR ) )
  
  # Create a data frame 
  df = data.frame( 
    Subject = rep( rep( 1:Ns, each = length( SR ) ), 4 ),
    Timepoints = rep(
      rep( c(pre,po), each = Ns * length( SR ) ), 2 ), 
    Condition = Cond[[3]], 
    Task = Cond[[2]], 
    Response_type = rep( c( 'Hits', 'False_alarms' ), each = 
                           Ns * length( SR ) * 2 ), 
    Target = rep( c( 1, 0 ), each = Ns * length( SR ) * 2 ),
    Trials = rep( Nt, 2 ), 
    Self_report_on_high = rep( SR * 100, 2 ), 
    stringsAsFactors = FALSE )
  
  # Posterior simulations
  sim = sdt_mix_rng( mdl$post,
                       nrow( df ),
                       Ns,
                       df$Trials,
                       df$Target,
                       df$Subject,
                       X_sr,
                       mdl$random )
  
  # Convert to d' and criterion estimates
  sel_H = df$Target == 1
  sim_dp = matrix( NA, nrow( sim )/2, ncol( sim ) )
  sim_crt = matrix( NA, nrow( sim )/2, ncol( sim ) )
  
  for ( i in 1:ncol( sim ) ) {
    
    sim_dp[,i] = sdt_calc_binary(
      sim[sel_H,i] * df$Trials[sel_H],
      df$Trials[sel_H],
      sim[!sel_H,i] * df$Trials[!sel_H],
      df$Trials[!sel_H],
      correct = 1,
      dp = T )
    
    sim_crt[,i] = sdt_calc_binary(
      sim[sel_H,i] * df$Trials[sel_H],
      df$Trials[sel_H],
      sim[!sel_H,i] * df$Trials[!sel_H],
      df$Trials[!sel_H],
      correct = 1,
      dp = F )
      
  }
  
  # Compute difference scores
  sel_pre = df$Timepoints == pre
  diff_H = sim[ sel_H & sel_pre, ] - 
    sim[ sel_H & !sel_pre, ]
  diff_dp = sim_dp[ sel_pre[ sel_H ], ] - 
    sim_dp[ !sel_pre[ sel_H ], ]
  diff_crt = sim_dp[ sel_pre[ sel_H ], ] - 
    sim_dp[ !sel_pre[ sel_H ], ]
  
  # Compute statistics over posterior simulations
  ppc = array( NA, dim = c( length( SR ), 3, 3, ncol( sim ) ) )
  # Function to compute statistics
  fta = function( x ) {
    aggregate( x, list( df$Self_report_on_high[ 1:nrow( diff_H ) ] ),
               function(x) c( 
                 mean(x), quantile(x,c(.025,.975)) ) )$x
  }
  
  for ( i in 1:ncol( sim ) ) {
    
    # Difference between pre-post hits
    ppc[,,1,i] = fta( diff_H[,i] )
    # Difference between pre-post d'
    ppc[,,2,i] = fta( diff_dp[,i] )
    # Difference between pre-post criterion
    ppc[,,3,i] = fta( diff_crt[,i] )
    
  }
  
  out = array( NA, dim = c( length( SR ), 3, 3 ) )
  dimnames( out )[[1]] = 
    paste( 'SRH_', SR * 100, '%', sep = '' )
  dimnames( out )[[2]] = c( 'Mean', 'PI_2.5', 'PI_97.5' )
  dimnames( out )[[3]] = c( 'H', 'dp', 'crt' )
  
  for ( i in 1:3 )
    out[,,i] = apply( ppc[,,i,], c(1,2), mean )
  
  return(out)
}

if ( !PA ) {
  # Conditions to examine
  pred_2_cnd = rbind(
    c( 'NBack_0', 'Placebo', 'T2_Post_drug' ),
    c( 'NBack_0', 'Drug', 'T2_Post_drug' ),
    c( 'NBack_0', 'Placebo', 'T3_Post_drug' ),
    c( 'NBack_0', 'Drug', 'T3_Post_drug' ),
    c( 'NBack_2', 'Placebo', 'T2_Post_drug' ),
    c( 'NBack_2', 'Drug', 'T2_Post_drug' ),
    c( 'NBack_2', 'Placebo', 'T3_Post_drug' ),
    c( 'NBack_2', 'Drug', 'T3_Post_drug' ) )
} else {
  # Conditions to examine
  pred_2_cnd = rbind(
    c( 'NBack_2', 'Placebo', 'T2_Post_drug' ),
    c( 'NBack_2', 'Drug', 'T2_Post_drug' ),
    c( 'NBack_2', 'Placebo', 'T3_Post_drug' ),
    c( 'NBack_2', 'Drug', 'T3_Post_drug' ) )
}

# Initialize output
pred_2 = array( NA,
                dim = c( nrow( pred_2_cnd ),
                         11, 3, 3 ) )

# Create a progress bar using a base R function
pb = txtProgressBar( min = 0, max = nrow( pred_2_cnd), style = 3 )

# Loop over conditions
for ( i in 1:nrow( pred_2_cnd ) ) {
  
  cond = list(
    pred_2_cnd[i,1],
    pred_2_cnd[i,2],
    pred_2_cnd[i,3] )
  
  # Generat posterior simulations
  pred_2[i,,,] = create_SR_PPC_input( cond, mdl, X )
  
  # Update the progress bar
  setTxtProgressBar(pb,i)
  
}
close(pb)

# Add meaningful labels
dimnames( pred_2 )[[2]] = 
  paste( 'SRH_', seq(0,1,.1) * 100, '%', sep = '' )
dimnames( pred_2 )[[3]] = c( 'Mean', 'PI_2.5', 'PI_97.5' )
dimnames( pred_2 )[[4]] = c( 'H', 'dp', 'crt' )

###
### 5) # Save results
###

setwd( proj_dir )
setwd( 'Data/Posterior_estimates' )
if ( !PA ) {
  save( pred, pred_2_cnd, pred_2, 
        file = paste( 'PPC_', model, '.RData', sep = '' ) )
} else {
  save( pred, pred_2_cnd, pred_2, 
        file = paste( 'PPC_', model, '_PA.RData', sep = '' ) )
}
setwd( orig_dir )
