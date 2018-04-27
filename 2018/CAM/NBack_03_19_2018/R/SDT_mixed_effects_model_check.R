# Simulation checks for SDT stan scripts
# Written by Kevin Potter
# email: kevin.w.potter@gmail.com
# Please email me directory if you 
# have any questions or comments
# Last updated 2018-04-20

# Table of contents
# 1) Initial setup
# 2) Parameter recovery with simulated data (Ver. 1)
# 3) Parameter recovery with simulated data (Ver. 2)

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
  T  # 
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

# Package for random draws from a multivariate normal
# install.packages( 'MASS' )
library( MASS )

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

###
### 2) Parameter recovery with simulated data (Ver. 1)
###

if ( run_code[1] ) {
  
  ### Simulate data ###
  
  # Design
  design = data.frame(
    Subject = 1,
    Condition = rep( c('Placebo', 'Drug'), each = 4 ),
    Timepoints = rep( rep( 
      c('T1_Pre_drug','T2_Post_drug'), each = 2 ), 2 ),
    Task = 'NBack_2',
    Response_type = rep( c('Hits','False_alarms'), 4 ),
    Target = rep( c(1,0), 4 ),
    Trials = rep( c(22,19), 4 ), 
    stringsAsFactors = FALSE
  )
  
  Ns = 66 # Number of subjects
  
  # Create matrix to store data
  dtbf = data.frame(
    # Subject index
    Subject = rep( 1:Ns, each = nrow( design ) ),
    # Placebo versus drug
    Condition = rep( design$Condition, Ns ),
    # Pre-post design
    Timepoints = rep( design$Timepoints, Ns ),
    # One task
    Task = rep( design$Task, Ns ),
    # Hits versus false alarms
    Response_type = rep( design$Response_type, Ns ),
    # 1 = positive trials, 0 = negative trials
    Target = rep( design$Target, Ns ),
    # Number of trials per observation
    Trials = rep( design$Trials, Ns ), 
    stringsAsFactors = FALSE )
  
  # Design matrices
  Xd = matrix( 0, nrow( dtbf ), 4 )
  sel = dtbf$Condition == 'Placebo'; Xd[sel,1] = 1
  sel = dtbf$Condition == 'Drug'; Xd[sel,2] = 1
  sel = dtbf$Condition == 'Placebo' & 
    dtbf$Timepoints == 'T2_Post_drug'; Xd[sel,3] = 1
  sel = dtbf$Condition == 'Drug' & 
    dtbf$Timepoints == 'T2_Post_drug'; Xd[sel,4] = 1
  
  Xc = matrix( 0, nrow( dtbf ), 2 )
  sel = dtbf$Condition == 'Placebo'; Xc[sel,1] = 1
  sel = dtbf$Condition == 'Drug'; Xc[sel,2] = 1
  
  # Generating parameters
  
  # Group-level means
  set.seed( 89784 ) # For reproducibility
  beta_dp = runif( 4, 
                   c( 1.7, 1.7, -.5, -.5 ),
                   c( 3.3, 3.3, .5, .5 ) )
  set.seed( 233 ) # For reproducibility
  beta_crt = runif( 2, -.5, .5 )
  
  # Define number of random effects
  Nef = 4
  # Specify position of random effects
  ref_pos = c( 1:2, 5:6 )
  
  # Generate a covariance matrix
  set.seed( 677 ) # For reproducibility
  A = matrix( runif(Nef^2)*2-1, ncol = Nef )
  Sigma = t(A) %*% A
  
  # Extract the correlation matrix
  Omega = cov2cor( Sigma )
  # Extract the variances
  Tau = diag( Sigma )
  
  # Scale down variances
  set.seed( 397 ) # For reproducibility
  Tau = runif( Nef, .1, .75 )
  Sigma = diag( Tau ) %*% Omega %*% diag( Tau )
  
  # Subject level effects
  eta = mvrnorm( Ns, rep( 0, Nef ), Sigma )
  
  # Create d' and criterion values
  dp = numeric( nrow( dtbf ) )
  crt = numeric( nrow( dtbf ) )
  beta_dp = cbind( beta_dp )
  beta_crt = cbind( beta_crt )
  eta_dp = matrix( 0, nrow( beta_dp ), 1 )
  eta_crt = matrix( 0, nrow( beta_crt ), 1 )
  xd = matrix( 0, 1, nrow( beta_dp ) )
  xc = matrix( 0, 1, nrow( beta_crt ) )
  for ( no in 1:nrow( dtbf ) ) {
    xd = Xd[no,]
    eta_dp[,1] = eta[ dtbf$Subject[no], ref_pos[1:2] ]
    dp[no] = xd %*% ( beta_dp + eta_dp )
    xc = Xc[no,]
    eta_crt[,1] = eta[ dtbf$Subject[no], ref_pos[1:2] ]
    crt[no] = xc %*% ( beta_crt + eta_crt )
  }
  theta = sdt_prob( dp, crt, dtbf$Target )
  # Simulate data
  dtbf$Counts = rbinom( nrow( dtbf ), dtbf$Trials, theta )
  
  dtbf %>% 
    group_by( Response_type, Condition, Timepoints ) %>% 
    summarize( P = mean( Counts / Trials ) )
  
  ### Parameter recovery ###
  
  # Specify priors for model
  Priors = rbind(
    c( 2.5, 1 ),
    c( 2.5, 1 ),
    c( 0, 1 ),
    c( 0, 1 ),
    c( 0, 1 ),
    c( 0, 1 ),
    c( 2, 0 ),
    c( 0, 1 )
  )
  
  # Obtain model estimates
  res = estimate_sdt( dtbf, list( Xd, Xc ),
                      ref_pos, Priors, debug = T )
  
  # Compare posterior estimates 
  # against generating parameters
  x11( width = 12 )
  
  xl = c( .5, 6.5 )
  yl = lowerUpper( .5, range( res$post$group_param ) )
  
  blankPlot( xl, yl )
  horizLines( seq( yl[1], yl[2], .5 ), xl, col = 'grey80',
              lwd = 2 )
  customAxes( xl, yl )
  
  ds_post = apply( res$post$group_param, 2,
                   function(x) {
                     p = c( .025, .16, .5, .84, .975 );
                     q = quantile(x,p);
                     names( q ) = c(
                       'CI_2.5','CI_16',
                       'Median',
                       'CI_84','CI_97.5');
                     out = c(
                       Mode = findMode( x ),
                       Mean = mean( x ),
                       q );
                     return( out )
                   } )
  for ( i in 1:6 ) {
    segments( i, ds_post['CI_2.5',i],
              i, ds_post['CI_97.5',i],
              lwd = 2 )
    errorBars( i, ds_post[c('CI_16','CI_84'),i],
               length = .05, lwd = 2 )
    points( i, ds_post['Median'], pch = 15, cex = 1.25 )
    points( i, ds_post['Mean'], pch = 19, cex = 1.25 )
    points( i, ds_post['Mean'], pch = 17, cex = 1.25 )
    
    if ( i < 5 )
      points( i, beta_dp[i,1], pch = 21, bg = 'white' )
    if ( i > 4 ) 
      points( i, beta_crt[i-4,1], pch = 21, bg = 'white' )
    
  }
  
}

###
### 3) Parameter recovery with simulated data (Ver. 2)
###

if ( run_code[2] ) {
  
  ### Simulate data ###
  
  # Design
  design = data.frame(
    Subject = 1,
    Condition = rep( c('Placebo', 'Drug'), each = 4 ),
    Timepoints = rep( rep( 
      c('T1_Pre_drug','T2_Post_drug'), each = 2 ), 2 ),
    Task = 'NBack_2',
    Response_type = rep( c('Hits','False_alarms'), 4 ),
    Target = rep( c(1,0), 4 ),
    Trials = rep( c(22,19), 4 ), 
    stringsAsFactors = FALSE
  )
  
  Ns = 66 # Number of subjects
  
  # Create matrix to store data
  dtbf = data.frame(
    # Subject index
    Subject = rep( 1:Ns, each = nrow( design ) ),
    # Placebo versus drug
    Condition = rep( design$Condition, Ns ),
    # Pre-post design
    Timepoints = rep( design$Timepoints, Ns ),
    # One task
    Task = rep( design$Task, Ns ),
    # Hits versus false alarms
    Response_type = rep( design$Response_type, Ns ),
    # 1 = positive trials, 0 = negative trials
    Target = rep( design$Target, Ns ),
    # Number of trials per observation
    Trials = rep( design$Trials, Ns ), 
    stringsAsFactors = FALSE )
  
  # Design matrices
  Xd = matrix( 0, nrow( dtbf ), 4 )
  sel = dtbf$Condition == 'Placebo'; Xd[sel,1] = 1
  sel = dtbf$Condition == 'Drug'; Xd[sel,2] = 1
  sel = dtbf$Condition == 'Placebo' & 
    dtbf$Timepoints == 'T2_Post_drug'; Xd[sel,3] = 1
  sel = dtbf$Condition == 'Drug' & 
    dtbf$Timepoints == 'T2_Post_drug'; Xd[sel,4] = 1
  
  Xc = matrix( 0, nrow( dtbf ), 2 )
  sel = dtbf$Condition == 'Placebo'; Xc[sel,1] = 1
  sel = dtbf$Condition == 'Drug'; Xc[sel,2] = 1
  
  # Generating parameters
  
  # Group-level means
  set.seed( 89784 ) # For reproducibility
  beta_dp = runif( 4, 
                   c( 1.7, 1.7, -.5, -.5 ),
                   c( 3.3, 3.3, .5, .5 ) )
  set.seed( 233 ) # For reproducibility
  beta_crt = runif( 2, -.5, .5 )
  
  # Define number of random effects
  Nef = 4
  # Specify position of random effects
  ref_pos = c( 1:2, 5:6 )
  
  # Generate a covariance matrix
  set.seed( 677 ) # For reproducibility
  A = matrix( runif(Nef^2)*2-1, ncol = Nef )
  Sigma = t(A) %*% A
  
  # Extract the correlation matrix
  Omega = cov2cor( Sigma )
  # Extract the variances
  Tau = diag( Sigma )
  
  # Scale down variances
  set.seed( 397 ) # For reproducibility
  Tau = runif( Nef, .1, .75 )
  Sigma = diag( Tau ) %*% Omega %*% diag( Tau )
  
  post_input = list(
    skip = 'Skip',
    group_param = matrix( c( beta_dp, beta_crt ), 
                          2, 6, byrow = T ),
    Omega = array( NA, dim = c( 2, Nef, Nef ) ),
    Tau = matrix( Tau, 2, Nef, byrow = T )
  )
  for ( i in 1:2 ) post_input$Omega[i,,] = Omega
  
  sim = sdt_mix_rng( post_input,
                     nrow( dtbf ),
                     max( dtbf$Subject ),
                     dtbf$Trials,
                     dtbf$Target, 
                     dtbf$Subject,
                     list( Xd, Xc ),
                     ref_pos )
  dtbf$Counts = round( sim[,1] * dtbf$Trials )
  
  dtbf %>% 
    group_by( Response_type, Condition, Timepoints ) %>% 
    summarize( P = mean( Counts / Trials ) )
  
  # Obtain model estimates
  res = estimate_sdt( dtbf, list( Xd, Xc ),
                      ref_pos, Priors, debug = T )
  
  # Compare posterior estimates 
  # against generating parameters
  x11( width = 12 )
  
  xl = c( .5, 6.5 )
  yl = lowerUpper( .5, range( res$post$group_param ) )
  
  blankPlot( xl, yl )
  horizLines( seq( yl[1], yl[2], .5 ), xl, col = 'grey80',
              lwd = 2 )
  customAxes( xl, yl )
  
  ds_post = apply( res$post$group_param, 2,
                   function(x) {
                     p = c( .025, .16, .5, .84, .975 );
                     q = quantile(x,p);
                     names( q ) = c(
                       'CI_2.5','CI_16',
                       'Median',
                       'CI_84','CI_97.5');
                     out = c(
                       Mode = findMode( x ),
                       Mean = mean( x ),
                       q );
                     return( out )
                   } )
  for ( i in 1:6 ) {
    segments( i, ds_post['CI_2.5',i],
              i, ds_post['CI_97.5',i],
              lwd = 2 )
    errorBars( i, ds_post[c('CI_16','CI_84'),i],
               length = .05, lwd = 2 )
    points( i, ds_post['Median'], pch = 15, cex = 1.25 )
    points( i, ds_post['Mean'], pch = 19, cex = 1.25 )
    points( i, ds_post['Mean'], pch = 17, cex = 1.25 )
    
    if ( i < 5 )
      points( i, beta_dp[i], pch = 21, bg = 'white' )
    if ( i > 4 ) 
      points( i, beta_crt[i-4], pch = 21, bg = 'white' )
    
  }
  
  
  
}

setwd( orig_dir )