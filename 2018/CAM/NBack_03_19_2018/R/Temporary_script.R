# 
# Written by Kevin Potter
# email: kevin.w.potter@gmail.com
# Please email me directory if you 
# have any questions or comments
# Last updated 2018-04-24

# Table of contents
# 1) Initial setup
# 2) Computations for priors

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

# Package for estimating mixed effects models
# install.packages( 'rstanarm' )
library( rstanarm )
# For parallel processing
options(mc.cores = parallel::detectCores())

# Load in data
setwd( 'Data' )
load( "NBack_3_19_2018.RData" )
setwd( proj_dir )

# Load in useful functions
setwd( 'R' )
source( 'Useful_functions.R' )
setwd( proj_dir )

# Exclude combined data
cd = dat %>% 
  filter( Task != 'Combined' )

# Drop rows with missing 
# self-report data for subject 
# FN_041 (66)
# Subject FN_041 did not have any 
# self report data on how high
# for the T1 drug condition
cd = cd %>% 
  filter( Self_report_on_high != '-' )

# Determine count data for misses/correct rejections
cd$Z = cd$Trials - cd$Counts
cd$Y = cd$Counts 

###
### 2) Computations for priors
###

# Determine average performance over tasks
avg = cd %>% 
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
###
###

if ( run_code[1] ) {
  
  # Create design matrices for d' and bias
  init_X = create_design_mat( cd )
  
  # Specify separate dummy-coded variables 
  # for each task and timepoint, and an 
  # order effect
  Xc = cbind(
    Baseline_0 = init_X$NBack_0,
    Time_0_2 = init_X$T2_0,
    Time_0_3 =  init_X$T3_0,
    Baseline_2 = init_X$NBack_2,
    Time_2_2 = init_X$T2_2,
    Time_2_3 =  init_X$T3_2,
    Visit_2 = init_X$O
  )
  
  check = create_design_mat( dtbf, Xc )
  
  # For d' variables, weight dummy coded 
  # values by .5 (False alarms) or -.5 (Hits)
  Xd = Xc; Xd = Xd * .5
  sel = dtbf$Response_type == 'Hits'
  Xd[sel] = -1 * Xd[sel]
  
}



