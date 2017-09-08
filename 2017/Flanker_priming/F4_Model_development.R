#--------------------#
# Model development  #
# Kevin Potter       #
# Updated 04/30/2017 #
#--------------------#

# Clear workspace
rm( list = ls() )

# Save current directory
orig_dir = getwd()

# Indicate which code segments to run
runCode = c( T, F, F )

# Index
# Lookup - 01: Load in useful packages, functions, and data

###
### Load in useful packages, functions, and data
###
# Lookup - 01

# For geting github packages
# install.packages(devtools)
# library(devtools)

# Integration of C++ and R
# install.packages(Rcpp)
library(Rcpp)

# Miscellanous functions for modeling/plotting
# install_github("rettopnivek/utilityf")
library(utilityf)

# Convenience functions for MLE
# install_github("rettopnivek/mle")
library(mle)

# Distribution functions for assorted sequential sampling models
# install_github("rettopnivek/seqmodels")
library(seqmodels)

# Distribution functions for wiener process
# install.packages('RWiener')
library(RWiener)

# Functions for plotting response time and choice data
# install_github("rettopnivek/rtplots")
library(rtplots)

# Functions for plotting response time and choice data
# install_github("rettopnivek/nROUSE")
library(nROUSE)

# Define useful functions
source('F2_useful_functions.R')

# Load in data
load( 'Data/Flanker_priming.RData' )
  
# For easy manipulation
d = allDat
colnames( d ) = c( 'S', 'Cnd', 'Co', 'Ac', 'Ch', 'RT', 
                   'T', 'P1', 'P3', 'ID', 'FL', 'F', 
                   'PT', 'PTL', 'PD', 'CL', 'RL' )


Nt = data.frame( Co = rep( 0:1, each = 16 ),
                 Cnd = rep( 1:16, 2 ),
                 Nt = 100 )
ms = meta_subject_create( d, Nt = Nt )

# Sort experimental conditions
exp_str = aggregate( d$Cnd, list( d$PD, d$PTL, d$FL, d$Cnd ), unique )
exp_str = exp_str[,-4]
colnames( exp_str ) = c( 'PD', 'PT', 'F','Cnd' )

###
###
###
# Lookup - 02

if ( runCode[1] ) {
  
  # Visualize data to be fitted
  cnd = list( c1 = list( PD = 0, 
                         PT = 'None', 
                         F = 'Identical', 
                         Co = 0,
                         flip = F ),
              c2 = list( PD = 0, 
                         PT = 'None', 
                         F = 'Incompatible', 
                         Co = 0,
                         flip = T ) )
  x11()
  par( mar = c( 4, 1, 3, 5 ) )
  plot_jcdf( ms, cnd, opt = list( yl = c(-1,1) ) )
  text( .2, 1, cnd$c1$F[1], pos = 4 )
  text( .2, -1, cnd$c2$F[1], pos = 4 )
  
  # Extract data to be fitted
  dtbf = quick_cnd_select( ms, 0, "None", 
                           c("Identical","Incompatible"),
                           0:1 )
  
  model = 1
  if ( model == 1 ) {
    # Wald race model
    # Accumulator for consonant (1)
    # Accumulator for vowel (0)
    
    # Model parameters -> matched to k(1), k(0), x(1), x(0), tau
    pNames = c( 'k_c', 'k_v',
                'x_inc_c1','x_id_c1','x_inc_c0','x_id_c0',
                'x_inc_v1','x_id_v1','x_inc_v0','x_id_v0',
                'tau' )
    
    ps = cbind( k1 = 1, k0 = 2, 
                x1 = rep( 3, nrow( dtbf ) ), 
                x0 = rep( 9, nrow( dtbf ) ), 
                tau = rep( 11, nrow( dtbf ) ) )
    ps = as.data.frame( ps )
    ps$x1[ dtbf$FL == 'Identical' & dtbf$CL == 'consonant' ] = 4
    ps$x1[ dtbf$FL == 'Incompatible' & dtbf$CL == 'vowel' ] = 5
    ps$x1[ dtbf$FL == 'Identical' & dtbf$CL == 'vowel' ] = 6
    ps$x0[ dtbf$FL == 'Identical' & dtbf$CL == 'consonant' ] = 10
    ps$x0[ dtbf$FL == 'Incompatible' & dtbf$CL == 'vowel' ] = 7
    ps$x0[ dtbf$FL == 'Identical' & dtbf$CL == 'vowel' ] = 8
    
    prm_type = c( rep( 1, 10 ), 2 )
    
    # Define function for maximum likelihood estimation
    mle_fn1 = function( prm, dat,  
                        sum = T, priors = NULL ) { 
      # Initial setup 
      rt = dat[[1]]$RT
      ch = dat[[1]]$Ch
      ps = dat[[2]]
      prm = tran_par( prm, prm_type )
      prm[11] = prm[11] * min( rt )
      
      # Calculate the log-likelihoods 
      ll = dwaldrace( rt, ch,
                      prm[ ps$k1 ], prm[ ps$x1 ], prm[ ps$tau ],
                      prm[ ps$k0 ], prm[ ps$x0 ], prm[ ps$tau ], ln = 1 )
      if ( !sum ) return( ll )
      
      # Sum the log-likelihoods 
      sll = sum( ll )
      
      # Incorporate priors (Optional) 
      if ( !is.null(priors) ) { 
        sll = sll
      }
      
      # Check for NA values 
      if ( is.na( sll ) ) sll = -Inf 
      
      return( sll ) 
    } 
    
    # Define function to generate dispersed starting values
    start = function() {
      prm = runif( 11, c( rep( .5, 4 ), 0.01, 0.01, .5, .5, .01, .01, -1 ), 
                   c( rep( 2, 4 ), 0.1, 0.1, 2, 2, .1, .1, 1 ) )
      prm = tran_par( prm, c( prm_type[-11], 0 ), reverse = T )
    }
    
    # Create input
    dat = list( dtbf, ps )
    
    # Estimate model
    m1 = MLE( dat, mle_fn1, start, nRep = 5, parNames = pNames )
    
    print( 'Wald race done' )
  }
  
  model = 2
  if ( model == 2 ) {
    
    tst = aggregate( dtbf$Ac, list( dtbf$FL, dtbf$CL ), p_f, val = c("0","1") )
    colnames( tst ) = c( 'Flanker', 'Target', 'Ac' )
    
    tst = rbind(
      cbind( tst[,-3], Ac = 1, P = tst$Ac[,2] ),
      cbind( tst[,-3], Ac = 0, P = tst$Ac[,1] ) )
      
    
    # F - T - R
    # C - C - C
    # V - C - C
    # V - V - V
    # C - V - V
    
    # C - C - V
    # V - C - V
    # V - V - C
    # C - V - C
    
    # Proces model
    # theta(1)     -> Pick target
    # Identical - consonant:    Controlled( kappa(1), xi(3), theta(8), tau(10) )
    # Incompatible - consonant: Controlled( kappa(1), xi(4), theta(8), tau(10) )
    # Identical - vowel:        Controlled( kappa(2), xi(5), theta(8), tau(10) )
    # Incompatible - vowel:     Controlled( kappa(2), xi(6), theta(8), tau(10) )
    # 1 - theta(1) -> Pick flanker
    
    
    # Correct
    # Error
    # Identical - consonant:    Automatic( kappa(1), xi(7), theta(9), tau(10) )
    # Incompatible - consonant: Automatic( kappa(1), xi(7), theta(9), tau(10) )
    # Identical - vowel:        Automatic( kappa(2), xi(7), theta(9), tau(10) )
    # Incompatible - vowel:     Automatic( kappa(2), xi(7), theta(9), tau(10) )
    
    # Shifted wald( kappa, xi, tau )
    
    # Wald (1) -> Correct
    # Wald (0) -> Error
    
    # Mixture of two wald processes
    pNames = c( 'k_inc', 'k_id',
                'x_inc','x_id',
                'tau', 'theta' )
    
    ps = cbind( k = 1, x = rep( 3, nrow( dtbf ) ), 
                tau = 5, theta = 6 )
    ps = as.data.frame( ps )
    ps$k[ dtbf$FL == 'Identical' ] = 4
    ps$x1[ dtbf$FL == 'Incompatible' & dtbf$CL == 'vowel' ] = 5
    ps$x1[ dtbf$FL == 'Identical' & dtbf$CL == 'vowel' ] = 6
    ps$x0[ dtbf$FL == 'Identical' & dtbf$CL == 'consonant' ] = 10
    ps$x0[ dtbf$FL == 'Incompatible' & dtbf$CL == 'vowel' ] = 7
    ps$x0[ dtbf$FL == 'Identical' & dtbf$CL == 'vowel' ] = 8
    
    # Define function for maximum likelihood estimation
    mle_fn1 = function( prm, dat,  
                        sum = T, priors = NULL ) { 
      # Initial setup 
      rt = dat[[1]]$RT
      ch = dat[[1]]$Ch
      ps = dat[[2]]
      prm = tran_par( prm, prm_type )
      prm[10] = prm[10] * min( rt )
      
      # Calculate the log-likelihoods
        
      if ( !sum ) return( ll )
      
      # Sum the log-likelihoods 
      sll = sum( ll )
      
      # Incorporate priors (Optional) 
      if ( !is.null(priors) ) { 
        sll = sll
      }
      
      # Check for NA values 
      if ( is.na( sll ) ) sll = -Inf 
      
      return( sll ) 
    } 
    
  }
  
}

setwd( orig_dir )