#---------------------#
# Data pre-processing #
# Kevin Potter        #
# Updated 05/24/2017  #
#---------------------#

# Clear workspace
rm(list = ls())

# Save current directory
orig_dir = getwd()

# Index
# Lookup - 01:  Initial setup
# Lookup - 02:  Global cut-offs
# Lookup - 03:  Mixture model

###
### Initial setup
###
# Lookup - 01

# Load in .csv file
setwd( 'Data' )
rawDat = read.csv( 'RawData.csv' )
setwd( orig_dir )

# Compute sample size
N = length( unique( rawDat$Subj ) )

# Create a dummy-coded version of accuracy
rawDat$Acc_DC = 0
rawDat$Acc_DC[ rawDat$Acc == 'Correct' ] = 1

# Load in useful packages

# To download packages off of Github
# install.packages( 'devtools' )
# library( devtools )

# Useful functions for plotting & statistics
# install_github("rettopnivek/utilityf")
library( utilityf )

# Define additional useful functions
source( 'F0_Useful_functions.R' )

# Create PDF of trimming results per subject
setwd( 'Plots' )
pdf( 'Data_preprocessing_per_subject.pdf',
     width = 12 )
setwd( orig_dir )

###
### Global cut-offs
###
# Lookup - 02

# Global cut-offs
# Trim response times equal to or less than 0
# Trim response times slower than 10 seconds

# Define logical matrix with values to exclude
exclude = matrix( F, nrow( rawDat ), 3 )
colnames( exclude ) = c( 'Too_fast', 'Too_slow', 'Unlikely' )
# Convert to data frame for easy indexing
exclude = as.data.frame( exclude )

# Identify response times that were too fast
exclude$Too_fast = rawDat$RT1 <= 0
# Identify response times that were too slow
exclude$Too_slow = rawDat$RT1 > 10

###
### Mixture model
###
# Lookup - 03

# Function to plot trimming results
trim_plot = function( rt, ac, ver ) {
  
  # Determine empirical density for observations
  d = densityPoints( rt )
  # Order from fastest to slowest
  o = order( rt )
  
  # Create plot
  par( mar = c(4,4,3,4) )
  plot( d$x, d$y, bty ='u', xlab = ' ',
        ylab = ' ', type = 'l' )
  if ( ver == 1 ) mtext( 'Density', side = 2, line = 2.25 )
  if ( ver == 3 ) mtext( 'Moving accuracy', side = 4, line = 2.25 )
  
  # Identify points that were excluded
  clr = rep( 'black', length( rt ) )
  if ( ver == 1 ) clr[ !keep ] = 'red'
  if ( ver == 2 ) clr[ unlikely ] = 'red'
  
  points( d$x, d$y, pch = 19, col = clr[o] )
  
  # Determine x and y-axis dimensions
  xyd = par("usr")
  
  # Add axes for accuracy
  axis( 4, seq( xyd[3], xyd[4], length = 5 ),
        seq( 0, 1, .25 ) )
  
  # Add in backwards accuracy
  cs = sapply( 1:( length(ac)-1 ),
               function(x) sum( ac[o][ -(1:x) ] ) )
  cs = c( sum(ac[o]), cs )/(length(ac):1)
  h = xyd[4] - xyd[3]
  lines( d$x, h * cs + xyd[3], col = 'grey' )
  # Add in upper boundary for 50% confidence interval on 
  # chance accuracy
  ub = qbinom( .75, length(ac):1, .5 )/( length(ac):1 )
  lines( d$x, ub * h + xyd[3], col = 'grey80', lty = 2 )
}

# Variables to store results
all_bfp = matrix( NA, N, 4 ) # Raw parameter estimates
all_tbfp = matrix( NA, N, 4 ) # Transformed parameter estimates
# Create list to store likelihoods
all_ll = c()
for ( s in 1:N ) all_ll = c( all_ll, list( NULL ) )

# Create progress par to track estimation progress
pb = txtProgressBar( min = 0, max = N, style = 3 )

# Loop over subjects
for ( s in 1:N ) {
  
  # Select subject
  sel = rawDat$Subj == s
  rt = rawDat$RT1[ sel ]
  
  # Apply global cut-offs
  keep = !exclude$Too_fast[sel] & !exclude$Too_slow[sel]
  x = rt[ keep ]
  
  # Starting values for optimization
  starting = c( est_exgauss_param( x ), .95 )
  names( starting )[4] = 'theta'
  if ( is.na( starting[2] ) ) starting[2] = .07
  starting = mixture_tran_par( starting, reverse = T )
  
  # Initial estimation 
  est_1 = optim( starting, mle_mixture_model, x = x, 
                 method = 'Nelder-Mead', 
                 control = list( fnscale = -1 ) )
  
  # Loop through multiple estimation attempts, 
  # jittering initial estimates and rescaling 
  # for better parameter exploration
  nRep = 5
  logLikVal = rep( NA, nRep ) # Track log-likelihood values
  all_est = matrix( NA, nRep, 4 )
  for ( nr in 1:nRep ) {
    starting = est_1$par + runif( 4, -1, 1 )
    
    res = tryCatch(
      optim( starting, mle_mixture_model, x = x, 
             method = 'BFGS', 
             control = list( fnscale = -1,
                             parscale = 1/abs( est_1$par ) ) ),
      error = function(e) return( NULL )
    )
    
    if ( !is.null( res ) ) {
      logLikVal[nr] = res$value
      all_est[nr,] = res$par
    }
    
  }
  
  # Exclude any failed estimation attempts
  succeed = which( !is.na( logLikVal ) )
  mx = which( logLikVal[ succeed ] == max( logLikVal[ succeed ] ) )
  
  # Extract best-fitting parameters (raw)
  bfp = all_est[ succeed[ mx ], ]
  
  # Compute likelihoods for each mixture component 
  # over observations
  ll = mle_mixture_model( x, bfp, sum = F )
  
  # Determine which observations are more likely under the 
  # uniform distribution
  unlikely = which( ll$du/ll$dexg > 1 )
  # Extract unlikely RTs
  x_u = x[ unlikely ]
  # Order from fastest to slowest
  o = order( x_u )
  
  # The mixture model excludes a fair amount of data 
  # in the right tail, so we will combine this technique 
  # with an examination of accuracy
  
  # Extract accuracy for unlikely RTs
  ac = rawDat$Acc_DC[sel][keep]
  ac_u = ac[ unlikely ]
  
  # Exclude response times that are unlikely and have 
  # accuracy within a 50% confidence interval of chance 
  # performance
  ci = qbinom( .75, length(ac_u):1, .5 )
  # Compute reverse cumulative frequency correct for 
  # sorted accuracies
  cs = sapply( 1:( length(ac_u)-1 ),
               function(x) sum( ac_u[o][ -(1:x) ] ) )
  cs = c( sum( ac_u[o] ), cs )
  # Determine which frequencies fall within confidence 
  # interval for chance performance
  trim = which( cs - ci <= 0 )
  # Exclude cases in which response time is less than 1 
  # second
  trim = trim[ x_u[o][ trim ] > 1 ]
  # Take minimum
  trim = min( trim )
  unlikely = unlikely[o][ trim:length( ac_u ) ]
  exclude$Unlikely[sel][keep][unlikely] = T
  
  # Store results
  all_bfp[s,] = bfp
  all_tbfp[s,] = mixture_tran_par( bfp )
  all_ll[[s]] = list( ll )
  
  rt = rawDat$RT1[sel]
  ac = rawDat$Acc_DC[sel]
  rt_k = rt[keep]
  ac_k = ac[keep]
  rt_kl = rt[keep][-unlikely]
  ac_kl = ac[keep][-unlikely]
  
  layout( rbind( 1:3 ) )
  
  trim_plot( rt, ac, 1 )
  trim_plot( rt_k, ac_k, 2 )
  trim_plot( rt_kl, ac_kl, 3 )
  mtext( 'Time (s)', side = 1, outer = T, line = -2 )
  mtext( paste( 'Subject', s ), side = 3, outer = T, line = -2 )
  
  # Update progress bar
  setTxtProgressBar(pb,s)
}
close( pb )

# Close PDf file
dev.off()

# Trim data
trimDat = rawDat[ !exclude$Too_fast & 
                    !exclude$Too_slow & 
                    !exclude$Unlikely, ]

# Save results
setwd( 'Data' )
save( rawDat, trimDat, exclude, N, file = 'FYP_JD.RData' )
setwd( orig_dir )