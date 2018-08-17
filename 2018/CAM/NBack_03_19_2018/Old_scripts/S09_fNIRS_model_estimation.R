# 
# Written by Kevin Potter
# email: kevin.w.potter@gmail.com
# Please email me directory if you 
# have any questions or comments
# Last updated 2018-04-24

# Table of contents
# 1) Initial setup



###
### 1) Initial setup
###

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

# Package for Bayesian multilevel modeling
# install.packages( 'brms' )
library( brms )

# Load in data
setwd( 'Data' )
load( "NBack_3_19_2018.RData" )
setwd( proj_dir )

# Load in useful functions
setwd( 'R' )
source( 'Useful_functions.R' )
setwd( proj_dir )

run_code = c(
  F,
  F,
  F
)

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

if ( run_code[1] ) {
  
  ### Design
  
  Ns = 54 # Number of subjects
  Nc = 5 # Number of channels per ROI
  Nroi = length( roi ) # Number of ROI
  # Data frame with design structure
  design = data.frame(
    Subject = 1,
    Channel = rep( 1:Nc, Nroi ),
    ROI = rep( roi, each = Nc ),
    ROI_num = rep( 1:Nroi, each = Nc ),
    zHbO = NA,
    stringsAsFactors = FALSE
  )
  # Data frame with observations
  sim = data.frame(
    Subject = rep( 1:Ns, each = nrow( design ) ),
    Channel = rep( design$Channel, Ns ),
    ROI = rep( design$ROI, Ns ),
    ROI_num = rep( design$ROI_num, Ns ),
    zHbO = NA, 
    stringsAsFactors = FALSE
  )
  No = nrow( sim )
  
  ### Generating parameters
  
  # Group-level
  beta = runif( Nroi, -1, 1 )
  sigma = runif( 2, -.5, 1.5 )
  nu = 10
  # Subject-level
  eta = rnorm( Ns, 0, sigma[2] )
  mu = beta[ sim$ROI_num ] + eta[ sim$Subject ]
  # Simulate data
  sim$zHbO = rt( No, df = nu )*sqrt( 
    (sigma[1]^2) * (nu-2)/nu) + mu
  
  ### Parameter recovery
  
  # Robust linear mixed effects model
  est = brm( formula = 
               zHbO ~ -1 + ROI_num + (1|Subject),
             family = student( link = "identity" ),
             data = sim,
             prior = c(
               set_prior("normal(0,1)", class = "b"),
               set_prior("cauchy(0,1)", class = "sd") ),
             warmup = 1000,
             iter = 2500, 
             chains = 4, 
             cores = 4, 
             control = list(adapt_delta = 0.95) )
  
}


###
###
###


setwd( orig_dir )