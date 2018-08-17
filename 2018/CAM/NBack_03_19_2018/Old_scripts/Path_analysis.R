# Title
# Written by Kevin Potter
# email: kevin.w.potter@gmail.com
# Please email me directory if you 
# have any questions or comments
# Last updated 2018-04-10

# Table of contents
# 1) Initial setup
# 2) Model 1
# 3) Model 2
# 4) Model 3
# 5) Path from fNIRS to self-report

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
  F, # Model 1
  F, # Model 2
  T, # Model 3
  F  # Path from fNIRS to self-report
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

# Load in package for Bayesian estimation 
# of mixed effects models
# install.packages( 'rstanarm' )
library( rstanarm )

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

###
### 2) Model 1
###

if ( run_code[1] ) {
  
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
  X1 = list( Xd, Xc )
  
  # Specify priors for model
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
  
  m1 = estimate_sdt( dtbf, 
                     X1, 
                     ref_pos, 
                     Priors )
  
  setwd( 'Data/Posterior_estimates' )
  save( m1, X1, file = 'Posterior_m1_PA.RData' )
  setwd( proj_dir )
  
  # Clean up workspace
  rm( m1, X1 )
  
}

###
### 3) Model 2
###

if ( run_code[2] ) {
  
  # Create design matrices for d' and bias
  init_X = create_design_mat( dtbf )
  
  # Specify separate d' values 
  # for each of the 12 conditions
  # and an added effect of 
  # the 5 ROI from the fNIRS data
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
    R_DLPFC = init_X$R_DLPFC,
    L_DLPFC = init_X$L_DLPFC,
    MPFC = init_X$MPFC,
    R_VLPFC = init_X$R_VLPFC,
    L_VLPFC = init_X$L_VLPFC
  )
  
  check = create_design_mat( dtbf, Xd )
  
  # Position of random effects (d')
  ref_pos = c( 1,7 )
  
  # Specify separate bias values 
  # for each of the 12 conditions
  # and an added effect of 
  # the 5 ROI from the fNIRS data
  Xc = Xd
  
  # Position of random effects (criterion)
  ref_pos = c( ref_pos, ncol( Xd ) + ref_pos )
  
  # Save design matrices
  X2 = list( Xd, Xc )
  
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
    c( 0.0, 1.0 ), # R. DLPFC
    c( 0.0, 1.0 ), # L. DLPFC
    c( 0.0, 1.0 ), # MPFC
    c( 0.0, 1.0 ), # R. VLPFC
    c( 0.0, 1.0 ), # L. VLPFC
    # Criterion
    c( 0.1, 1.0 ), # Baseline (0) *R - 19*
    c( 0.0, 1.0 ), # Placebo (0 - T1 )
    c( 0.0, 1.0 ), # Placebo (0 - T2 )
    c( 0.0, 1.0 ), # Placebo (0 - T3 )
    c( 0.0, 1.0 ), # Drug (0 - T2)
    c( 0.0, 1.0 ), # Drug (0 - T3)
    c( 0.4, 1.0 ), # Baseline (2) *R - 25*
    c( 0.0, 1.0 ), # Placebo (2 - T1 )
    c( 0.0, 1.0 ), # Placebo (2 - T2 )
    c( 0.0, 1.0 ), # Placebo (2 - T3 )
    c( 0.0, 1.0 ), # Drug (2 - T2)
    c( 0.0, 1.0 ), # Drug (2 - T3)
    c( 0.0, 1.0 ), # Visit 2
    c( 0.0, 1.0 ), # R. DLPFC
    c( 0.0, 1.0 ), # L. DLPFC
    c( 0.0, 1.0 ), # MPFC
    c( 0.0, 1.0 ), # R. VLPFC
    c( 0.0, 1.0 ), # L. VLPFC
    # Random effects
    c( 2.0, 0.0 ), # Omega
    c( 0.0, 1.0 )  # Tau
  )
  
  m2 = estimate_sdt( dtbf, 
                     X2, 
                     ref_pos, 
                     Priors,
                     adapt_delta = .999,
                     max_treedepth = 15 )
  
  setwd( 'Data/Posterior_estimates' )
  save( m2, X2, file = 'Posterior_m2_PA.RData' )
  setwd( proj_dir )
  
  # Clean up workspace
  rm( m2, X2 )
  
}

###
### 4) Model 3
###

if ( run_code[3] ) {
  
  # Create design matrices for d' and bias
  init_X = create_design_mat( dtbf )
  
  # Specify separate d' values 
  # for each of the 12 conditions
  # and added effects of 
  # self-reported high and 
  # the 5 ROI from the fNIRS data
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
    Self_report_high = init_X$SR, 
    R_DLPFC = init_X$R_DLPFC,
    L_DLPFC = init_X$L_DLPFC,
    MPFC = init_X$MPFC,
    R_VLPFC = init_X$R_VLPFC,
    L_VLPFC = init_X$L_VLPFC
  )
  
  check = create_design_mat( dtbf, Xd )
  
  # Position of random effects (d')
  ref_pos = c( 1,7 )
  
  # Specify separate bias values 
  # for each of the 12 conditions
  # and an added effect of 
  # the 5 ROI from the fNIRS data
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
    c( 0.0, 1.0 ), # R. DLPFC
    c( 0.0, 1.0 ), # L. DLPFC
    c( 0.0, 1.0 ), # MPFC
    c( 0.0, 1.0 ), # R. VLPFC
    c( 0.0, 1.0 ), # L. VLPFC
    # Criterion
    c( 0.1, 1.0 ), # Baseline (0) *R - 20*
    c( 0.0, 1.0 ), # Placebo (0 - T1 )
    c( 0.0, 1.0 ), # Placebo (0 - T2 )
    c( 0.0, 1.0 ), # Placebo (0 - T3 )
    c( 0.0, 1.0 ), # Drug (0 - T2)
    c( 0.0, 1.0 ), # Drug (0 - T3)
    c( 0.4, 1.0 ), # Baseline (2) *R - 26*
    c( 0.0, 1.0 ), # Placebo (2 - T1 )
    c( 0.0, 1.0 ), # Placebo (2 - T2 )
    c( 0.0, 1.0 ), # Placebo (2 - T3 )
    c( 0.0, 1.0 ), # Drug (2 - T2)
    c( 0.0, 1.0 ), # Drug (2 - T3)
    c( 0.0, 1.0 ), # Visit 2
    c( 0.0, 1.0 ), # Self-reported high
    c( 0.0, 1.0 ), # R. DLPFC
    c( 0.0, 1.0 ), # L. DLPFC
    c( 0.0, 1.0 ), # MPFC
    c( 0.0, 1.0 ), # R. VLPFC
    c( 0.0, 1.0 ), # L. VLPFC
    # Random effects
    c( 2.0, 0.0 ), # Omega
    c( 0.0, 1.0 )  # Tau
  )
  
  m3 = estimate_sdt( dtbf, 
                     X3, 
                     ref_pos, 
                     Priors,
                     adapt_delta = .999,
                     max_treedepth = 15 )
  
  setwd( 'Data/Posterior_estimates' )
  save( m3, X3, file = 'Posterior_m3_PA.RData' )
  setwd( proj_dir )
  
  # Clean up workspace
  rm( m3, X3 )
  
}

###
### 5) Path from fNIRS to self-report
###

if ( run_code[4] ) {
  
  init_X = create_design_mat( dtbf )
  
  dtbf_2 = cbind( Subject = dtbf$Subject, 
                  Task = dtbf$Task,
                  Condition = dtbf$Condition,
                  Timepoints = dtbf$Timepoints, 
                  init_X )
  
  res = lmer( SR ~ -1 + Condition + 
                R_DLPFC + 
                L_DLPFC + 
                MPFC + 
                R_VLPFC + 
                L_VLPFC + 
                (1|Subject), data = dtbf_2 )
  
  
  dtbp = dtbf_2
  dtbp$Pred = predict( res )
  dtbp %>% 
    group_by( Task, Condition, Timepoints ) %>% 
    summarize( Obs = mean( SR ),
               Pred = mean( Pred ) )
  
}

setwd( orig_dir )