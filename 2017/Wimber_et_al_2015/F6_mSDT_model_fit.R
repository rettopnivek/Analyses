#--------------------#
# mSDT model fit     #
# Kevin Potter       #
# Updated 01/24/2017 #
#--------------------#

# Clear workspace
rm( list = ls() )

# Save current directory
orig_dir = getwd()

###
### Load in useful packages and data
###

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

# Load in useful functions
source( 'F1_useful_functions.R' )

# Load in data
setwd( 'Data' )
load( 'Recog_mem.RData' )
setwd( orig_dir )

###
### Fit model using Stan
###
setwd( 'Stan_scripts' )

# Define input for Stan
stanDat = list(
  Ns = length( sort( unique( d$S ) ) ), 
  Ni = length( sort( unique( d$IN ) ) ), 
  No = nrow( d ), 
  Y = d$Ch, 
  Co = d$Co, 
  subjIndex = d$S, 
  Cond = cbind( d$B+1, d$IT )
)

warm = 750 # Warm-up
niter = 1250*3 # Number of samples to approximate posterior
chains = 8 # Number of chains to run

startTime = Sys.time() # To assess run-time

# Compile model
sm = stan_model(stanc_ret = stanc_builder("mSDT_model_no_item.stan"))

# Draw samples
fit = sampling( sm, data = stanDat, 
                warmup = warm, iter = warm+niter, 
                chains = chains,
                thin = 3, 
                control = list( adapt_delta = .995,
                                max_treedepth = 14 ) )

post = extract(fit)

setwd( orig_dir )

# Save results to a .RData file
fName = "C:/Users/Kevin/Documents/Posteriors from Stan"
setwd(fName)
setwd( "Wimber_et_al" )
save( stanDat, fit, post, file = 'mSDT_model_no_item_post.RData' )
setwd( orig_dir )
