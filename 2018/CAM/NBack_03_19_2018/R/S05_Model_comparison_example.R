# Model comparison examples
# Written by Kevin Potter
# email: kevin.w.potter@gmail.com
# Please email me directly if you 
# have any questions or comments
# Last updated 2018-05-10

# Table of contents
# 1) Initial setup
# 2) Bayes factor/LOO-CV example with 'rstanarm'

###
### 1) Initial setup
###

# Save current working directory
orig_dir = getwd()

# Navigate to project directory
setwd( '..' )
proj_dir = getwd()

# Specify directory for posterior estimates
setwd( 'Data/Posterior_estimates' )
post_dir = getwd()
setwd( proj_dir )

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

# Package for estimating mixed effects models
# install.packages( 'lme4' )
library( lme4 )

# Package for estimating mixed effects models (Bayesian)
# install.packages( 'rstanarm' )
library( rstanarm )
# For parallel processing
options(mc.cores = parallel::detectCores())

# Package to compute Bayes Factors from MCMC samples
# install.packages( 'bridgesampling' )
library( bridgesampling )

###
### 2) Bayes factor/LOO-CV example with 'rstanarm'
###

if ( run_code[1] ) {
  
  data( 'turtles' )
  
  # Example of probit regression with with 
  # deviations based on clutch number
  b('
  m0 = glmer( y ~ 1 + (1|clutch),
              family = binomial( link = "probit" ),
              data = turtles )
  m1 = glmer( y ~ 1 + x + (1|clutch),
              family = binomial( link = "probit" ),
              data = turtles )
  ')
  
  # Estimation using rstnarm
  m0 = stan_glmer( y ~ 1 + (1|clutch),
                   prior_intercept = normal( 0, 10 ),
                   family = binomial( link = "probit" ),
                   data = turtles,
                   warm = 1000, 
                   iter = 13500,
                   chains = 4,
                   cores = 4,
                   diagnostic_file = file.path(tempdir(), "df.csv") )
  # Bridge sampling
  bridge_0 = bridge_sampler( m0 )
  
  # Estimation using rstnarm
  m1 = stan_glmer( y ~ 1 + x + (1|clutch),
                   prior_intercept = normal( 0, 10 ),
                   prior = normal( 0, 10 ), 
                   family = binomial( link = "probit" ),
                   data = turtles,
                   warm = 1000, 
                   iter = 13500,
                   chains = 4,
                   cores = 4,
                   diagnostic_file = file.path(tempdir(), "df.csv") )
  # Bridge sampling
  bridge_1 = bridge_sampler( m1 )
  
  # LOO
  l0 = loo( m0 )
  l1 = loo( m1 )
  
  # Tabulate model comparison results
  mc = data.frame(
    Model = c( 0, 1 ),
    LML = c( bridge_0$logml,
             bridge_1$logml ),
    BF = c( NA, NA ),
    LOOIC = c( l0$estimates[3,1],
               l1$estimates[3,1] ), 
    AW = c( NA, NA ),
    LOOMW = c( NA, NA )
  )
  # Bayes factor
  mc$BF = c( bf( bridge_0, bridge_0 )$bf, 
             bf( bridge_1, bridge_0 )$bf )
  # Akaike weights
  mc$AW = icWeights( mc$LOOIC )
  # Stacking method for combining multiple 
  # predictive distributions
  mc$LOOMW = as.numeric( loo_model_weights(
    list( l0, l1 ) ) )
  
}

setwd( orig_dir )