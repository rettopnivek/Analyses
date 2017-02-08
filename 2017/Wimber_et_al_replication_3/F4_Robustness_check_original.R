#--------------------#
# Robustness check   #
# Kevin Potter       #
# Updated 05/09/2016 #
#--------------------#

# Clear workspace
rm(list = ls())

# Load in data to R workspace
load("Original_Rec_Mem_data.RData")

# Determine the sample size
N = length( unique( OriginalRecMem$Subject ) )

# Load in an R package for fitting mixed effects models
# install.packages(lme4)
library(lme4)

# To examine the robustness of the model, we'll try excluding pairs of 
# subjects and refitting the model.
p_value = matrix( NA, N, N )
effect_size = matrix( NA, N, N )

# Fit a generalized linear model excluding all possible pairs of subjects
pb = txtProgressBar( min = 1, max = N, style = 3 )
for ( n in 1:N ) {
  
  subj = sort( unique( OriginalRecMem$Subject ) )
  
  dS = OriginalRecMem[ OriginalRecMem$Subject != subj[n], ]
  
  for ( n2 in 1:(N-1) ) {
    
    subj2 = sort( unique( dS$Subject ) )
    d = dS[ dS$Subject != subj2[n2], ]
    
    d$Y = d$Accuracy
    d$IT = as.factor( d$ImageType )
    d$SR = as.factor( 1 - d$Baseline )
    d$S = as.factor( d$Subject )
    d$I = as.factor( d$ImageNum )
    TotalTrials = rep(1,nrow(d))
    d$RIF = 0; d$RIF[ d$IT == 2 & d$SR == 1 ] = 1
    
    fit = glmer( Y ~ RIF + (1|S) + (1|I), 
                 data = d, family = binomial(link='logit'),
                 weights = TotalTrials )
    tmp = summary( fit )
    p_value[ subj2[n2], n ] = tmp$coefficients[2,4]
    effect_size[ subj2[n2], n ] = tmp$coefficients[2,1]
    
  }
  
  rm( fit, tmp )
  # Update the progress bar
  setTxtProgressBar(pb,n)
}
close(pb)

print( 'Number of times effect became non-significant' )
print( sum( na.omit( as.vector( p_value ) ) > .05 ) )