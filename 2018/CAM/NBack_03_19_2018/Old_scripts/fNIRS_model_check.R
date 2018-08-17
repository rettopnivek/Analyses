# Simulation checks for fNIRS stan scripts
# Written by Kevin Potter
# email: kevin.w.potter@gmail.com
# Please email me directory if you 
# have any questions or comments
# Last updated 2018-04-23

# Table of contents
# 1) Initial setup
# 2) Parameter recovery with simulated data

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
  T  
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

# Labels for ROI
roi = c(
  'R_DLPFC',
  'L_DLPFC',
  'MPFC',
  'R_VLPFC',
  'L_VLPFC' )

# Mappings of ROI to channels
roi_mapping = 
  c( 2, #  1 = L. DLPFC -  2
     2, #  2 = L. DLPFC -  2
     5, #  3 = L. VLPFC -     5
     5, #  4 = L. VLPFC -     5
     2, #  5 = L. DLPFC -  2
     5, #  6 = L. VLPFC -     5
     3, #  7 = MPFC     -   3
     2, #  8 = L. DLPFC -  2
     3, #  9 = MPFC     -   3
     1, # 10 = R. DLPFC - 1
     5, # 11 = L. VLPFC -     5
     3, # 12 = MPFC     -   3
     4, # 13 = R. VLPFC -    4
     3, # 14 = MPFC     -   3
     1, # 15 = R. DLPFC - 1
     4, # 16 = R. VLPFC -    4
     1, # 17 = R. DLPFC - 1
     1, # 18 = R. DLPFC - 1
     4, # 19 = R. VLPFC -    4
     4  # 20 = R. VLPFC -    4
  )

# 2.1)
rstudentt = function( n, nu, mu = 0, sigma = 1 ) {
  # Purpose:
  # A function to generate random draws from a student-t 
  # distribution given the degrees of freedom, mean, and 
  # standard deviation.
  # Arguments:
  # n     - The number of draws
  # nu    - The degrees of freedom (can be a positive non-integer)
  # mu    - The distribution mean
  # sigma - The distribution standard deviation
  # Returns:
  # Random draws from the student-t distribution.
  
  s2 = ( sigma^2 )/rgamma( n, nu/2.0, nu/2.0 )
  out = rnorm( n, mu, sqrt(s2) )
  
  return( out )
}

if ( run_code[1] ) {
  
  ### Simulate data ###
  
  # Design
  design = data.frame(
    Subject = 1,
    Condition = rep( c('Placebo', 'Drug'), each = 20 * 2 ),
    Timepoints = rep( rep( 
      c('T1_Pre_drug','T2_Post_drug'), each = 20 ), 2 ),
    Task = 'NBack_2',
    Channel = rep( 1:20, 4 ),
    ROI = rep( roi[ roi_mapping ], 4 ), 
    ROI_num = rep( roi_mapping, 4 ), 
    zHbO = NA, 
    nu = 10, 
    stringsAsFactors = FALSE
  )
  
  Ns = 54 # Number of subjects
  
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
    # Channel
    Channel = rep( design$Channel, Ns ),
    # ROI
    ROI = rep( design$ROI, Ns ),
    # Numeric coding for ROI
    ROI_num = rep( design$ROI_num, Ns ), 
    # Observed z-scores for HbO
    zHbO = NA, 
    nu = rep( design$nu, Ns ), # Degrees of freedom
    stringsAsFactors = FALSE )
  
  # 5 ROI per each of the 4 conditions
  Nroi = 5 * 4
  Xc = rbind(
    c( 1, 0, 0, 0 ),
    c( 1, 1, 0, 0 ),
    c( 1, 0, 1, 0 ),
    c( 1, 0, 0, 1 )
  )
  X = matrix( 0, Nroi, 5 * 4 )
  for ( i in 1:5 ) {
    xa = 1:4 + 4*(i-1)
    X[xa,xa] = Xc
  }
  
  X_mapping = dtbf %>% 
    group_by( ROI_num, Condition, Timepoints ) %>% 
    summarize( V = unique( ROI_num ) )
  X_mapping$V = 1:20
  X_mapping = as.data.frame( X_mapping )
  
  dtbf$indexROI = NA
  for ( i in 1:nrow( X_mapping ) ) {
    sel = 
      dtbf$ROI_num == X_mapping$ROI_num[i] & 
      dtbf$Condition == X_mapping$Condition[i] & 
      dtbf$Timepoints == X_mapping$Timepoints[i]
    dtbf$indexROI[sel] = X_mapping$V[i]
  }
  
  # Generating parameters
  
  # Group-level means
  set.seed( 20130 ) # For reproducibility
  beta = cbind( runif( 4 * 5, -.5, .5 ) )
  Mu = X %*% beta
  
  # Generate a covariance matrix
  set.seed( 321 ) # For reproducibility
  A = matrix( runif(Nroi^2)*2-1, ncol = Nroi )
  Sigma = t(A) %*% A
  
  # Extract the correlation matrix
  Omega = cov2cor( Sigma )
  # Extract the variances
  Tau = diag( Sigma )
  
  # Scale down variances
  set.seed( 895 ) # For reproducibility
  Tau = runif( Nroi, .1, .75 )
  Sigma = diag( Tau ) %*% Omega %*% diag( Tau )
  
  # Subject-level ROI means
  ROI_subjects = MASS::mvrnorm( Ns, Mu, Sigma )
  # Subject standard deviations
  sigma_channel = runif( Ns, .5, 2 )
  
  # Simulate zHbO values
  for ( i in 1:nrow( dtbf ) ) {
    dtbf$zHbO[i] = rstudentt( 1, 
                              dtbf$nu[i],
                              ROI_subjects[ dtbf$Subject[i],
                                            dtbf$indexROI[i] ],
                              sigma_channel[ dtbf$Subject[i] ] )
  }
  
  ### Parameter recovery ###
  
  # Specify priors for model
  Priors = list(
    priors_beta = matrix( c(0,1), 2, ncol(X), byrow = F ),
    priors_Omega = 4,
    priors_Tau = matrix( c(0,1), 2, 
                         max( X_mapping$V ), byrow = F )
  )
  
  # Obtain model estimates
  res = estimate_fNIRS( dtbf, X, Priors, debug = T )
  
}


setwd( orig_dir )