#-----------------------------#
# nROUSE parameter estimation #
# Kevin Potter                #
# Updated 11/13/2016          #
#-----------------------------#

# Initialize script
source('F5_starting_script.R')

# Load in package for simulating nROUSE model
# install_github("rettopnivek/nROUSE")
library(nROUSE)

# Load in package for adaptive Bayesian estimation
# install.packages('MHadaptive')
library( MHadaptive )

# Load in data
setwd( 'Data' )
load( 'Priming_offset.RData' )
setwd( orig_dir )

# For easy manipulation
d = allDat
colnames( d ) = c('S','PD','PT','O','Co','RT',
                  'Ch','Ac','OT','DT','CDT','Cnd')


###
### Simulation example for one subject
###

# Generating parameters
mu = c( N = .0283, I = 1.2866, Ta = 0.8931, TD = 65 )
sig = c( N = .002, I = .35, Ta = .3, TD = 15 )
# Jitter by some noise
gp = rnorm( 4, mu, sig ); gp[4] = round( gp[4] )

# Create data set to obtain nROUSE model estimates
nDat = cbind( TarDur = gp[4], MaskDur = 500 - gp[4],
              PrimeDur = c(50,2000,50,2000),
              Type = c(2,2,-2,-2), Y = rep(0,4),
              N = rep(0,4) )
nDat = as.data.frame( nDat )
theta = nROUSE_logLik( log(gp[1:3]), nDat, estimate = F )

# Simulate frequency correct
nDat$N = rep( 80, 4 )
nDat$Y = rbinom( 4, nDat$N, theta )

priors = cbind( c(.0283, 1.2866, .8931 ),
                c(.01, .75, .6 ) )

ll_f = function( prm, nDat, priors ) {
  
  ll = nROUSE_logLik( log(prm), nDat )
  ll = ll + sum( dnorm( prm, priors[,1], priors[,2], log = T ) )
  if (is.na(ll)) ll = -Inf
  
  return( ll )
}

st_val = runif( 3, priors[,1] - priors[,2],
                priors[,1] + priors[,2] )


ps = cbind( c(1.1e-6, 3.52e-5, 2.25e-5),
            c(3.52-5,1.6e-2,7.175e-3),
            c(2.25-5,7.12e-3,8.137e-3) )




res = Metro_Hastings( ll_f, st_val, nDat = nDat, priors = priors, 
                      prop_sigma = ps, 
                      par_names = c('N','I','A'),
                      burn_in = 2000, 
                      iterations = 10000 )

# x11()
# hist( p

