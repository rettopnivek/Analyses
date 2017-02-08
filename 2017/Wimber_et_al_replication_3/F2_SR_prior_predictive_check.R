#----------------------------------------#
# Prior predictive check for selective   #
# retrieval task analysis                #
# Kevin Potter                           #
# Updated 02/07/2017                     #
#----------------------------------------#

# Clear workspace
rm( list = ls() )

# Save current directory
orig_dir = getwd()

# Index
# Lookup - 01:  Load in useful packages
# Lookup - 02:  Define inputs for prior predictive check
# Lookup - 03:  Sample from priors using Stan
# Lookup - 04:  Plot results

###
### Load in useful packages
###
# Lookup - 01

# For geting github packages
# install.packages(devtools)
# library(devtools)

# Miscellanous functions for modeling/plotting
# install_github("rettopnivek/utilityf")
library(utilityf)

# Load in package for Bayesian estimation
library(rstan)
# For parallel processing
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

###
### Define inputs for prior predictive check
###
# Lookup - 02

# Set working directory to location of Stan scripts
setwd('Stan_scripts')

Ns = 48 # Number of subjects
No = 54*Ns # Number of observations
indS = rep( 1:Ns, each = 54 ) # Index for subjects

# Determine means
mu = reverseSoftmax( c( .747, .091, .025, .137 ), 
                     restrict = c( F, F, F, T ) )

# Specify desired priors to test
Priors = cbind( c( 1.700, -.409, -1.700 ),
                c( .3, .3, .3 ),
                c( 2, 2, 2 ),
                c( 8, 8, 8 ) )

###
### Sample from priors using Stan
###
# Lookup - 03

# List of input for Stan
stan_dat = list(
  No = No,
  Ns = Ns,
  K = 4,
  Priors = Priors
)

niter = 1250 # Number of samples to draw
chains = 8 # Number of chains to run in parallel

startTime = Sys.time() # To assess run-time

# Compile model
sm = stan_model(stanc_ret = stanc_builder("SR_prior_check.stan"))

# Draw samples
fit = sampling( sm, data = stan_dat, 
                iter = niter, 
                chains = chains,
                seed = 74893, # For reproducibility
                algorithm = "Fixed_param" )

post = extract(fit)
runTime = Sys.time() - startTime # To assess run-time
print( runTime )

###
### Plot results
###
# Lookup - 04

setwd( orig_dir )
setwd( 'Plots' )

pdf( 'Prior_predictive_check_SR.pdf', width = 12 )

layout( cbind( c(1,3),c(2,4) ) )
hist( post$Targets/No, col = 'grey', border = 'white',
      bty = 'l', xlab = 'P(Target)', freq = T,
      main = 'Prior distribution over targets' )
legend('topright', as.character( round( mean( post$Targets/No ), 2 ) ), 
       bty = 'n' )
hist( post$Competitors/No, col = 'grey', border = 'white',
      bty = 'l', xlab = 'P(Competitors)', freq = T,
      main = 'Prior distribution over competitors' )
legend('topright', as.character( round( mean( post$Competitors/No ), 2 ) ), 
       bty = 'n' )
hist( post$Errors/No, col = 'grey', border = 'white',
      bty = 'l', xlab = 'P(Errors)', freq = T,
      main = 'Prior distribution over errors' )
legend('topright', as.character( round( mean( post$Errors/No ), 2 ) ), 
       bty = 'n' )
hist( post$Unknown/No, col = 'grey', border = 'white',
      bty = 'l', xlab = 'P(Unknown)', freq = T,
      main = 'Prior distribution over unknowns' )
legend('topright', as.character( round( mean( post$Unknown/No ), 2 ) ), 
       bty = 'n' )
dev.off()

setwd( orig_dir )