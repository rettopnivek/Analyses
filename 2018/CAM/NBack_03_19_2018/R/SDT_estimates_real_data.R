# SDT model estimation with real data
# Written by Kevin Potter
# email: kevin.w.potter@gmail.com
# Please email me directory if you 
# have any questions or comments
# Last updated 2018-04-17

# Table of contents
# 1) Initial setup
# 2) Computations for priors
# 3) Model 1 (Null model, Effect of task / Timepoints only)
# 4) Model 2 (Effect of drug, no self report)
# 5) Model 3 (Effect of drug, self-reported high)

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
  T, # Model 1
  F, # Model 2
  F  # Model 3
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

###
### 2) Computations for priors
###

# Determine average performance over tasks

avg = dtbf %>% 
  group_by( Response_type, Task ) %>% 
  summarize( P = mean( Counts / Trials ) ) %>% 
  arrange( desc( Response_type ) )

Priors_task = data.frame(
  Task = c( 'NBack_0', 'NBack_2' ),
  dp = rep( NA, 2 ),
  crt = rep( NA, 2 ) )

Priors_task$crt = -.5 * ( 
  qnorm( avg$P[ avg$Response_type == 'Hits' ] ) + 
    qnorm( avg$P[ avg$Response_type == 'False_alarms' ] ) )
Priors_task$dp = 2 * ( 
  qnorm( avg$P[ avg$Response_type == 'Hits' ] ) + 
    Priors_task$crt )

###
### 3) Model 1 (Null model, Effect of task / Timepoints only)
###

if ( run_code[1] ) {
  
  # Create design matrices for d' and bias
  init_X = create_design_mat( dtbf )
  
  # Specify separate d' values for 
  # each task and timepoint, 
  # and an order effect
  Xd = cbind(
    Baseline_0 = init_X$NBack_0,
    Time_0_2 = init_X$T2_0,
    Time_0_3 =  init_X$T3_0,
    Baseline_2 = init_X$NBack_2,
    Time_2_2 = init_X$T2_2,
    Time_2_3 =  init_X$T3_2,
    Visit_2 = init_X$O
  )
  
  check = create_design_mat( dtbf, Xd )
  
  # Specify separate bias values for 
  # each task and timepoint
  Xc = Xd
  
  # Save design matrices
  X1 = list( Xd, Xc )
  
  # Define priors for model
  Priors = rbind(
    # d'
    c( 4.1, 1.0 ), # Baseline (0) *R - 1*
    c( 0.0, 1.0 ), # Time (0 - T2)
    c( 0.0, 1.0 ), # Time (0 - T3)
    c( 2.9, 1.0 ), # Baseline (2) *R - 4*
    c( 0.0, 1.0 ), # Time (0 - T2)
    c( 0.0, 1.0 ), # Time (0 - T3)
    c( 0.0, 1.0 ), # Visit 2
    # Criterion
    c( 0.1, 1.0 ), # Baseline (0) *R - 8*
    c( 0.0, 1.0 ), # Time (0 - T2)
    c( 0.0, 1.0 ), # Time (0 - T3)
    c( 0.4, 1.0 ), # Baseline (2) *R - 11*
    c( 0.0, 1.0 ), # Time (0 - T2)
    c( 0.0, 1.0 ), # Time (0 - T3)
    c( 0.0, 1.0 ), # Visit 2
    # Random effects
    c( 2.0, 0.0 ), # Omega
    c( 0.0, 1.0 )  # Tau
  )
  
  # Estimate model using Stan
  m1 = estimate_sdt( dtbf, 
                     X1, 
                     c(1,4,8,11), 
                     Priors )
  
  setwd( 'Data/Posterior_estimates' )
  save( m1, X1, file = 'Posterior_m1.RData' )
  setwd( proj_dir )
  
  # Clean up workspace
  rm( m1, X1 )
  
}

###
### 4) Model 2 (Effect of drug, no self report)
###

if ( run_code[2] ) {
  
  # Create design matrices for d' and bias
  init_X = create_design_mat( dtbf )
  
  # Specify separate d' values 
  # for each of the 12 conditions
  Xd = cbind(
    Baseline_0 = init_X$NBack_0,
    Placebo_0_T1 = init_X$T1_0_P,
    Placebo_0_T2 =  init_X$T2_0_P,
    Placebo_0_T3 = init_X$T3_0_P,
    Drug_0_T2 =  init_X$T2_0_D,
    Drug_0_T3 = init_X$T3_0_D,
    Baseline_2 = init_X$NBack_2,
    Placebo_2_T1 = init_X$T1_2_P,
    Placebo_2_T2 =  init_X$T2_2_P,
    Placebo_2_T3 = init_X$T3_2_P,
    Drug_2_T2 =  init_X$T2_2_D,
    Drug_2_T3 = init_X$T3_2_D,
    Visit_2 = init_X$O
  )
  
  check = create_design_mat( dtbf, Xd )
  
  # Position of random effects (d')
  ref_pos = c( 1,7 )
  
  # Specify separate bias values 
  # for each of the 12 conditions
  Xc = Xd
  
  # Position of random effects (criterion)
  ref_pos = c( ref_pos, ncol( Xd ) + ref_pos )
  
  # Save design matrices
  X2 = list( Xd, Xc )
  
  # Define priors for model
  # Specify priors for model
  Priors = rbind(
    # d'
    c( 4.1, 1.0 ), # Baseline (0) *R - 1*
    c( 0.0, 1.0 ), # Placebo (0 - T1 )
    c( 0.0, 1.0 ), # Placebo (0 - T2 )
    c( 0.0, 1.0 ), # Placebo (0 - T3 )
    c( 0.0, 1.0 ), # Drug (0 - T2)
    c( 0.0, 1.0 ), # Drug (0 - T3)
    c( 2.9, 1.0 ), # Baseline (2) *R - 7*
    c( 0.0, 1.0 ), # Placebo (2 - T1 )
    c( 0.0, 1.0 ), # Placebo (2 - T2 )
    c( 0.0, 1.0 ), # Placebo (2 - T3 )
    c( 0.0, 1.0 ), # Drug (2 - T2)
    c( 0.0, 1.0 ), # Drug (2 - T3)
    c( 0.0, 1.0 ), # Visit 2
    # Criterion
    c( 0.1, 1.0 ), # Baseline (0) *R - 14*
    c( 0.0, 1.0 ), # Placebo (0 - T1 )
    c( 0.0, 1.0 ), # Placebo (0 - T2 )
    c( 0.0, 1.0 ), # Placebo (0 - T3 )
    c( 0.0, 1.0 ), # Drug (0 - T2)
    c( 0.0, 1.0 ), # Drug (0 - T3)
    c( 0.4, 1.0 ), # Baseline (2) *R - 20*
    c( 0.0, 1.0 ), # Placebo (2 - T1 )
    c( 0.0, 1.0 ), # Placebo (2 - T2 )
    c( 0.0, 1.0 ), # Placebo (2 - T3 )
    c( 0.0, 1.0 ), # Drug (2 - T2)
    c( 0.0, 1.0 ), # Drug (2 - T3)
    c( 0.0, 1.0 ), # Visit 2
    # Random effects
    c( 2.0, 0.0 ), # Omega
    c( 0.0, 1.0 )  # Tau
  )
  
  m2 = estimate_sdt( dtbf, 
                     X2, 
                     ref_pos, 
                     Priors )
  
  setwd( 'Data/Posterior_estimates' )
  save( m2, X2, file = 'Posterior_m2.RData' )
  setwd( proj_dir )
  
  # Clean up workspace
  rm( m2, X2 )
}

###
### 5) Model 3 (Saturated)
###

if ( run_code[3] ) {
  
  # Create design matrices for d' and bias
  init_X = create_design_mat( dtbf )
  
  # Specify separate d' values 
  # for each of the 12 conditions
  # and an added effect of self-report
  Xd = cbind(
    Baseline_0 = init_X$NBack_0,
    Placebo_0_T1 = init_X$T1_0_P,
    Placebo_0_T2 =  init_X$T2_0_P,
    Placebo_0_T3 = init_X$T3_0_P,
    Drug_0_T2 =  init_X$T2_0_D,
    Drug_0_T3 = init_X$T3_0_D,
    Baseline_2 = init_X$NBack_2,
    Placebo_2_T1 = init_X$T1_2_P,
    Placebo_2_T2 =  init_X$T2_2_P,
    Placebo_2_T3 = init_X$T3_2_P,
    Drug_2_T2 =  init_X$T2_2_D,
    Drug_2_T3 = init_X$T3_2_D,
    Visit_2 = init_X$O, 
    Self_report_high = init_X$SR
  )
  
  check = create_design_mat( dtbf, Xd )
  
  # Position of random effects (d')
  ref_pos = c( 1,7 )
  
  # Specify separate bias values 
  # for each of the 12 conditions
  # and an added effect of self-report
  Xc = Xd
  
  # Position of random effects (criterion)
  ref_pos = c( ref_pos, ncol( Xd ) + ref_pos )
  
  # Save design matrices
  X3 = list( Xd, Xc )
  
  # Specify priors for model
  Priors = rbind(
    # d'
    c( 4.1, 1.0 ), # Baseline (0) *R - 1*
    c( 0.0, 1.0 ), # Placebo (0 - T1 )
    c( 0.0, 1.0 ), # Placebo (0 - T2 )
    c( 0.0, 1.0 ), # Placebo (0 - T3 )
    c( 0.0, 1.0 ), # Drug (0 - T2)
    c( 0.0, 1.0 ), # Drug (0 - T3)
    c( 2.9, 1.0 ), # Baseline (2) *R - 7*
    c( 0.0, 1.0 ), # Placebo (2 - T1 )
    c( 0.0, 1.0 ), # Placebo (2 - T2 )
    c( 0.0, 1.0 ), # Placebo (2 - T3 )
    c( 0.0, 1.0 ), # Drug (2 - T2)
    c( 0.0, 1.0 ), # Drug (2 - T3)
    c( 0.0, 1.0 ), # Visit 2
    c( 0.0, 1.0 ), # Self-reported high
    # Criterion
    c( 0.1, 1.0 ), # Baseline (0) *R - 15*
    c( 0.0, 1.0 ), # Placebo (0 - T1 )
    c( 0.0, 1.0 ), # Placebo (0 - T2 )
    c( 0.0, 1.0 ), # Placebo (0 - T3 )
    c( 0.0, 1.0 ), # Drug (0 - T2)
    c( 0.0, 1.0 ), # Drug (0 - T3)
    c( 0.4, 1.0 ), # Baseline (2) *R - 21*
    c( 0.0, 1.0 ), # Placebo (2 - T1 )
    c( 0.0, 1.0 ), # Placebo (2 - T2 )
    c( 0.0, 1.0 ), # Placebo (2 - T3 )
    c( 0.0, 1.0 ), # Drug (2 - T2)
    c( 0.0, 1.0 ), # Drug (2 - T3)
    c( 0.0, 1.0 ), # Visit 2
    c( 0.0, 1.0 ), # Self-reported high
    # Random effects
    c( 2.0, 0.0 ), # Omega
    c( 0.0, 1.0 )  # Tau
  )
  
  m3 = estimate_sdt( dtbf, 
                     X3, 
                     ref_pos, 
                     Priors )
  
  setwd( 'Data/Posterior_estimates' )
  save( m3, X3, file = 'Posterior_m3.RData' )
  setwd( proj_dir )
  
  # Clean up workspace
  rm( m3, X3 )
  
}

setwd( orig_dir )

