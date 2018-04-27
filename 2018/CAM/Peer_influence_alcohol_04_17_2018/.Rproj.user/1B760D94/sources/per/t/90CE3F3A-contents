# Exploratory factor analysis on summed scores
# Written by Kevin Potter
# email: kevin.w.potter@gmail.com
# Please email me directory if you 
# have any questions or comments
# Last updated 2018-04-26

# Table of contents
# 1) Initial setup
#   1.1) lsq_ncp
#   1.2) RMSEA
#   1.3) multiple_EFA
#   1.4) EFA_fit_indices
#   1.5) my_piecewise_reg
#   1.6) rotate_factor_loadings
#   1.7) custom_color_map
#   1.8) horizLines
#   1.9) vertLines
# 2) Exploratory factor analysis
#   2.1) Extract data to fit
#   2.2) Fit common factor models
#   2.3) Extract fit indices for first 14 models
#   2.4) Determine number of factors to consider
#   2.5) Rotations/interpretations

###
### 1) Initial setup
###

# Save current working directory
orig_dir = getwd()

# Navigate to project directory
setwd( '..' ); proj_dir = getwd();

# Load in useful packages

# Collection of useful functions
# devtools::install_github("rettopnivek/utilityf")
library( utilityf )

# Packages for easy manipulating of data frames
# install.packages( 'dplyr' )
library(dplyr)

# Package for rotation of exploratory factor analysis
# install.packages('GPArotation')
library(GPArotation)

# Load in data
setwd( 'Data' )
load( 'Peer_influence_alcohol.RData' )
setwd( proj_dir )

# Specify whether to generate a pdf for figures
savePlot = T
if ( savePlot ) {
  setwd( 'Plots' )
  pdf( 'EFA_results.pdf', width = 12 )
  setwd( proj_dir )
}

# 1.1)
lsq_ncp = function ( ncp, bound, T.ml, fit ) {
  # Purpose:
  # A function to calculate the least squares difference 
  # between the cumulative probability of the test 
  # statistic in a non-central chi-square distribution 
  # given the non-centrality parameter.
  # Arguments:
  # ncp   - The non-centrality parameter for the 
  #         chi-square distribution
  # bound - The desired cumulative probability (e.g. .05, 
  #         .95, etc...)
  # T.ml  - The objective function value
  # fit   - A fit object return by factanal
  # Returns:
  # The least-squares difference between the desired 
  # probability and the generated probabilty 
  # based on the non-centrality parameter.
  
  par = abs(ncp) # Ensure the parameter is positive
  
  # Return the least-squares difference
  out = ( ( bound - 
              pchisq( T.ml, fit$dof, 
                      ncp, lower.tail = F ) )^2 )
  return( out )
}

# 1.2)
RMSEA = function( fit, 
                  N, 
                  conf = .90, 
                  lb = c( 0, 1000 ), 
                  ub = c( 0, 1000 ), 
                  inc = 1, 
                  Bart = 'T', 
                  debug = F ) {
  # Purpose:
  # A function to compute the RMSEA and its confidence interval.
  # Arguments:
  # fit   - A factanal object
  # N     - Sample size
  # conf  - The desired confidence level (e.g. 90%, 95%, 
  #         etc...)
  # lb    - The lower and upper intervals to search for the 
  #         optimal non-centrality parameter for the lower 
  #         boundary of the uncertainty interval
  # ub    - The lower and upper intervals to search for the 
  #         optimal non-centrality parameter for the upper 
  #         boundary of the uncertainty interval
  # inc   - The increments to search over
  # Bart  - Indicates whether Bartlett's correction should be 
  #         used. If it should not, then set 'Bart' equal to the 
  #         observed covariance matrix
  # debug - Indicates whether or not to print out progress of 
  #         the function
  # Returns:
  # Either the estimated RMSEA value, or if a desired confidence 
  # width is specified, a list with the estimated RMSEA value 
  # and the associated uncertainty interval.
  
  if ( is.character(Bart) ) {
    # Extract the negative log of the likelihood
    F.ml = as.numeric( fit$STATISTIC/(N - 1) )
    # Extract the objective function 
    # (i.e. the chi-square statistic)
    T.ml = fit$STATISTIC # (N - 1)*F.ml
  } else {
    # The predicted covariance matrix
    sigma = tcrossprod(fit$loadings) + 
      diag(fit$uniquenesses)
    F.ml = log( det(sigma) ) - 
      log( det(S) ) + 
      sum( diag( (S-sigma) %*% solve(sigma) ) )
    T.ml = F.ml*(N - 1)
    if ( debug ) print(T.ml)
  }
  
  # Correct for bias
  F.o = F.ml - 
    ( fit$dof / ( N - 1 ) )
  # Sets RMSEA to 0 if F.o is negative
  if ( F.o < 0 ) F.o = 0
  # Calculate the RMSEA
  val = sqrt( F.o/fit$dof )
  
  if ( !is.null( conf ) ) {
    
    # Determine the upper and lower boundaries 
    # for the confidence band
    ui = c(
      conf + (1 - conf)/2,
      (1 - conf)/2 )
    
    # Use a brute-force method to determine the 
    # non-centrality parameter that produces the 
    # closest cumulative probabilities to the
    # desired upper and lower boundaries
    
    # Intervals to explore
    ncp_explore = list(
      seq( lb[1], lb[2], inc ),
      seq( ub[1], ub[2], inc ) )
    
    # Lower boundary
    lsq = sapply( ncp_explore[[1]], 
                  lsq_ncp, bound = ui[1], T.ml = T.ml, fit = fit )
    ncp_lb = ncp_explore[[1]][ which.min( lsq ) ]
    if (debug) print(ncp_lb) # For debugging purposes
    
    # Upper boundary
    lsq = sapply( ncp_explore[[2]], 
                  lsq_ncp, bound = ui[2], T.ml = T.ml, fit = fit )
    ncp_ub = ncp_explore[[2]][ which.min( lsq ) ]
    if (debug) print(ncp_ub) # For debugging purposes
    
    # Compute the confidence intervals
    CI.u = sqrt(ncp_ub/((N - 1)*fit$dof))
    CI.l = sqrt(ncp_lb/((N - 1)*fit$dof))
    
    out = list(
      RMSEA = val,
      UI = c( CI.l, CI.u )
    )
    
    names( out$UI ) = 
      paste( '%', round( ui*100, 1 ), sep = '' )
    
  } else {
    out = val
  }
  
  return( out )
}

# 1.3)
multiple_EFA = function( MV, 
                        max_factors = ncol( MV ),
                        ... ) {
  # Purpose:
  # A convenience function that fits a series 
  # of common factor models up to the desired 
  # number of maximum factors.
  # Arguments:
  # MV          - The N x P matrix with the observed 
  #               scores for the P manifest variables
  # max_factors - The maximum number of factors to 
  #               include in the final analysis
  # ...         - Additional parameters to pass to the 
  #               'factanal' function
  # Returns:
  # A list with the 'factanal' output for 
  # the desired number of factors.
  
  # Initialize output
  out = c()
  for ( nf in 1:max_factors ) {
    out = c( out, list( NULL ) )
  }
  names( out ) = paste( 'F', 1:max_factors, sep = '' )
  
  for ( nf in 1:max_factors ) {
    
    est = tryCatch( 
      factanal( x = MV, 
                factors = nf, ... ),
      error = function(e) 'Estimation failed' )
    if ( !is.character( est ) ) {
      out[[nf]] = est
    } else {
      string = paste( 'Estimation failed:',
                      names( out )[nf] )
      warning( string, call. = FALSE )
    }
    
  }
  
  return( out )
}

# 1.4) 
EFA_fit_indices = function( MV, fit, ... ) {
  # Purpose:
  # Computes a set of fit indices for a given 
  # 'factanal' object.
  # Arguments:
  # MV  - A N x P matrix with the observed scores for 
  #       the P manifest variables
  # fit - A 'factanal' object
  # ... - Additional parameters for the RMSEA function
  # Returns:
  # A list with the assorted fit indices for the 
  # exploratory factor analysis.
  
  # Sample size
  N = nrow( MV )
  
  # Observed covariance matrix
  S = cov( MV )
  
  # Number of unique observations in 
  # the covariance matrix
  p = ncol(S)*(ncol(S)+1)/2
  
  # Degrees of freedom
  dof = fit$dof
  
  # Number of free parameters in the model
  k = p - fit$dof
  
  # Residuals
  
  # The predicted covariance matrix
  sigma = tcrossprod(fit$loadings) + 
    diag(fit$uniquenesses)
  
  # The residuals, the difference between the predicted and
  # observed covariance matrix
  res = S - sigma
  
  # Counts the number of residuals greater than .1
  large_res = sum(abs(res) > .1)
  
  # To find the root mean square error 
  # (RMSE; this is not the same as the RMSEA)
  RMSE = sqrt( sum( res[lower.tri( res, diag = F )]^2 ) )
  
  # The formula for the negative log likelihood:
  # F.ml = log(det(sigma)) - 
  #          log(det(S)) + 
  #          sum(diag((S-sigma)%*%solve(sigma)))
  F.ml = as.numeric(fit$STATISTIC/(N - 1))
  
  # Test of perfect fit
  
  # The test statistic
  T.ml = fit$STATISTIC # (N - 1)*F.ml
  
  # The p-value
  p.val = pchisq( T.ml, fit$dof, lower.tail = F ) 
  
  # The test statistic, if the assumptions hold, is
  # distributed as a chi-square distribution. Hence, we
  # can determine the probability of observing a test
  # statistic as high or higher given that the null
  # hypothesis is true (i.e. the model fits in the 
  # population). The function pchisq calculates this
  # probability, where lower.tail = F indicates that 
  # it is a one-tailed test.
  
  # RMSEA
  rmsea = RMSEA( fit, N, ... )
  
  # AIC
  
  # The AIC is defined as: 2*k - 2*ln(L) 
  # where k is the number of parameters 
  # and L is the maximum likelihood. In 
  # factor analysis, the AIC, under the 
  # hypothesis that the model is true 
  # for the population of interest, is 
  # calculated as:
  aic = T.ml - 2*fit$dof
  
  # When comparing models, a AIC with 
  # a smaller value is to be preferred 
  # (If the AIC is negative, a more negative 
  # value is to be preferred).
  
  # Output
  out = list(
    N = N,
    S = S,
    U = p,
    K = k,
    DoF = dof,
    R = res,
    NLR = large_res,
    RMSE = RMSE,
    NLL = F.ml,
    ToPF = c( statistic = T.ml, p = p.val ),
    RMSEA = rmsea,
    AIC = aic
  )
  
  return( out )
}

# 1.5)
my_piecewise_reg = function( x, y, breaks, ... ) {
  # Purpose:
  # Given a predictor and dependent variable, along 
  # with a vector of break points, computes a 
  # piecewise regression.
  # Arguments:
  # x      - A predictor variable
  # y      - A dependent variable
  # breaks - A vector of break points based on the 
  #          predictor variable
  # ...    - Additional options to pass on to the 
  #          'lm' function
  # Returns:
  # A list with the output from the 'lm' function and the 
  # list of data with the dummy codes to construct the 
  # piecewise regression.
  
  # Number of breaks
  nb = length( breaks )
  nbp = nb + 1
  
  # Number of observations
  n = length( y )
  
  # Initialize matrix for x and y values
  M = matrix( NA, n, 2 + nbp*2 )
  M[,1] = x; M[,2] = y;
  
  elbows = data.frame(
    lb = c( min(x) - 1, breaks ),
    ub = c( breaks, max(x) ) )
  
  # Separate intercepts for each segment
  M[ , 1:nbp + 2 ] = apply( elbows, 1, 
                            function(x) 
                              as.numeric( 
                                M[,1] > x[1] & M[,1] <= x[2] ) )
  # Separate slopes for each segment
  M[ , 1:nbp + nbp + 2 ] = 
    M[ , 1:nbp + 2 ] * x
  
  # Define variable labels
  cn = c(
    'x', 'y', 
    paste( 'I', 1:nbp, sep = '' ),
    paste( 'S', 1:nbp, sep = '' ) )
  colnames( M ) = cn
  
  # Convert to data frame
  M = as.data.frame( M )
  
  # Create formula for regression
  piecewise_reg_formula = 
    paste( 'y ~ -1 + ',
           paste( cn[ grep( 'I', cn ) ], collapse = ' + ' ),
           ' + ',
           paste( cn[ grep( 'S', cn ) ], collapse = ' + ' ),
           collapse = '' )
  piecewise_reg_formula = 
    as.formula( piecewise_reg_formula )
  
  # Fit piecewise regression
  fit = lm( piecewise_reg_formula, data = M, ... )
  
  # Output
  out = list(
    lm = lm( piecewise_reg_formula, data = M, ... ),
    data = M )
  
  # To plot unbroken segments, subsequent 
  # segments following the first must 
  # start from the previous value
  
  # Create indices to insert previous values 
  # into data matrix
  break_pos = sapply( breaks, function(x) 
    min( which( M[,1] >= x ) ) )
  break_pos = c( 1, break_pos, n )
  index1 = numeric( n + nb )
  index2 = numeric( n + nb )
  inc = 0
  for ( i in 2:length( break_pos ) ) {
    val2 = break_pos[i-1]:break_pos[i] 
    val = val2
    if ( i > 2 ) val[1] = val2[1] + 1
    index1[1:length(val) + inc] = val
    index2[1:length(val) + inc] = val2
    inc = inc + length(val)
  }
  # Create new matrix of data to use to 
  # generate x and y values for plotting 
  # purposes
  PM = matrix( NA, nrow( M ) + nb,
               ncol( M ) - 2 )
  PM[,1:nbp] = as.matrix( M[ index1, 1:nbp + 2 ] )
  PM[,1:nbp+nbp] = PM[,1:nbp] * x[ index2 ]
  # Generate plotting values
  out$xa = x[ index2 ]
  out$ya = PM %*% cbind( coef( fit ) )
  
  return( out )
}

# 1.6) 
rotate_factor_loadings = function( fit, oblique = T, ... ) {
  # Purpose:
  # ...
  # Arguments:
  # fit     - 
  # oblique - 
  # ...     - 
  # Returns:
  # ...
  
  # Rotate loadings using GPArotation package
  if ( oblique ) {
    rot = GPFoblq( loadings(fit), ... )
  } else {
    rot = GPForth( loadings(fit), ... )
  }
  
  # Extract loadings matrix
  ldf = data.frame(
    Variables = rownames( rot$loadings ),
    Order = 1:nrow( rot$loadings ),
    stringsAsFactors = FALSE
  )
  ldf = cbind( ldf, rot$loadings )
  
  # Round loadings
  tmp = round( rot$loadings * 100 )
  colnames( tmp ) = paste( colnames( tmp ), 'R', sep = '_' )
  ldf = cbind( ldf, tmp )
  rownames( ldf ) = NULL
  
  # Take absolute value of loadings
  tmp = abs( rot$loadings )
  
  # Sort based on factor magnitude
  no = numeric( nrow( ldf ) )
  mx_loading = apply( tmp, 1, function(x) which.max(x) )
  inc = 0
  cutoffs = numeric( nrow( ldf ) )
  for ( fn in 1:ncol( tmp ) ) {
    sel = mx_loading == fn
    cutoffs[ 1:sum(sel) + inc ] = fn
    vrb = ldf$Variables[sel][ order( tmp[ sel, fn ],
                                     decreasing = T ) ]
    no[ 1:sum(sel) + inc ] = vrb
      
    inc = inc + sum(sel)
  }
  o = sapply( no, function(x) which( x == ldf$Variables ) )
  ldf = ldf[o,]
  ldf$Cutoffs = cutoffs
  
  # Cut-offs for highest loadings per factor
  
  
  # Return results
  var_to_keep = c( 'Variables', 'Order', 'Cutoffs', 
                   paste( 'Factor', 
                          1:ncol( tmp ), '_R', sep = '' ),
                   paste( 'Factor', 
                          1:ncol( tmp ), sep = '' ) )
  ldf = ldf[,var_to_keep]
  rownames( ldf ) = as.character( 1:nrow( ldf ) )
  
  out = list(
    Loadings = ldf )
  if ( oblique ) out$Phi = rot$Phi
  
  return( out )
}

# 1.7)
custom_color_map = function( Z, vr = NULL, nc = 11 ) {
  # Purpose:
  # ...
  # Arguments:
  # Z
  # vr
  # nc
  # Returns:
  # ...
  
  # Dimensions
  nx = ncol( Z )
  ny = nrow( Z )
  
  # Range
  if ( is.null( vr ) ) {
    vr[1] = min( Z )
    vr[2] = max( Z )
  }
  
  # Set to odd number
  if ( nc %% 2 != 1 ) {
    nc = nc + 1
  }
  
  # List of plotting variables
  CM = list(
    LM = Z,
    cmap = list(
      z = seq( vr[1], vr[2], length = nc + 1 ),
      clr = c(
        rgb( seq( 0, .8, length = (nc-1)/2 ), 1, 1 ),
        rgb( 1, 1, 1 ),
        rev( rgb( 1, seq( 0, .8, length = (nc-1)/2 ), 1 ) ) ) ),
    x = cbind(
      rep( 0:(nx-1), each = ny ),
      rep( 0:(nx-1), each = ny ),
      rep( 1:nx, each = ny ),
      rep( 1:nx, each = ny ) ),
    y = cbind(
      rep( 0:(ny-1), nx ),
      rep( 1:ny, nx ),
      rep( 1:ny, nx ),
      rep( 0:(ny-1), nx ) )
  )
  # Specify color gradient
  CM$clr = data.frame( clr = rep( rgb( 0,0,0 ), nx * ny ),
                       x = NA, y = NA,
                       z = NA, 
                       stringsAsFactors = F )
  # Loop over matrix cells
  inc = 1
  for ( i in 1:nx ) {
    for ( j in ny:1 ) {
      sel = sum( CM$LM[j,i] >= CM$cmap$z[ -12 ] )
      CM$clr$clr[inc] = 
        CM$cmap$clr[sel]
      CM$clr$x[inc] = i; CM$clr$y[inc] = j
      CM$clr$z[inc] = CM$cmap$z[sel]
      inc = inc + 1
    }
  }
  
  xl = c( 0, nx ); yl = c( 0, ny )
  blankPlot( xl, yl )
  for ( i in 1:(nx*ny) ) {
    polygon( CM$x[i,],
             CM$y[i,],
             border = NA, col = CM$clr$clr[i] )
  }
  CM$xl = xl; CM$yl = yl
  
  return( CM )
}

# 1.8)
horizLines = function( yval, xl, ... ) {
  # Purpose: 
  # Draws a set of horizontal lines onto 
  # an existing plot.
  # Arguments: 
  # yval - A vector of y-axis values at which to 
  #        draw the lines.
  # xl   - A vector giving the lower and upper bounds 
  #        of the x-axis.
  
  l = length( yval )
  
  segments( rep( xl[1], l ), yval,
            rep( xl[2], l ), yval, ... )
  
}

# 1.9)
vertLines = function( xval, yl, ... ) {
  # Purpose: 
  # Draws a set of vertical lines onto 
  # an existing plot.
  # Arguments: 
  # xval - A vector of x-axis values at which to 
  #        draw the lines.
  # yl   - A vector giving the lower and upper bounds 
  #        of the y-axis.
  
  l = length( xval )
  
  segments( xval, rep( yl[1], l ),
            xval, rep( yl[2], l ), ... )
  
}


###
### 2) Exploratory factor analysis
###

# 2.1) Extract data to fit

# Standardize variables
cn = colnames( SC )
sel = c( grep( 'MISS', cn ),
         grep( 'UPPS', cn ),
         grep( 'DMQ', cn ),
         grep( 'BCEOA', cn ), 
         grep( 'Total_days', cn ),
         grep( 'Total_drinks', cn )
)

# Standardize all summary scores
my_standardize = function( x ) {
  out = ( x - mean( x, na.rm = T ) ) / sd( x, na.rm = T )
  return( out )
}

# Standardize the manifest variables
MV = apply( SC[,sel], 2, my_standardize )

# Pairwise deletion for missing data
is_na = apply( MV, 1, function(x) any( is.na(x) ) )
MV = MV[ !is_na, ]

# 2.2) Fit common factor models

# Maximum number of factors to fit
mf = ncol( MV )
all_FA = multiple_EFA( MV, mf, rotation = 'none' )
# Estimation fails after 14 factors

# 2.3) Extract fit indices for first 14 models

fit_details = matrix( NA, 14, 6 )
colnames( fit_details ) = 
  c( 'N_large_residuals',
     'RMSE',
     'RMSEA',
     'RMSEA_UI_L',
     'RMSEA_UI_U',
     'AIC' )
for ( i in 1:14 ) {
  fi = EFA_fit_indices( MV, all_FA[[i]] )
  fit_details[i,1] = fi$NLR
  fit_details[i,2] = fi$RMSE
  fit_details[i,3] = fi$RMSEA$RMSEA
  fit_details[i,4] = fi$RMSEA$UI[1]
  fit_details[i,5] = fi$RMSEA$UI[2]
  fit_details[i,6] = fi$AIC
}

# 2.4) Determine number of factors to consider

if ( !savePlot ) x11( width = 12 )
layout( cbind( 1, 2, 3 ) )

# 2.4.1) Scree plot

# Compute eigenvalues for correlation matrix
ev = eigen( cor( MV ) )

# Create a scree plot

# Plot characteristics
lnSz = 2
ptSz = 1.5
axSz = 1.25
txtSz = 1.5
axPos = -1

xa = 1:length( ev$values )
ya = ev$values

par( mar = c( 4, 5, 3, .5 ) )

# Create a blank plot
xl = c(.5, length(xa ) + .5 )
yl = c( 0, lowerUpper( 1, ev$values )[2] )
blankPlot( xl, yl )
customAxes( xl, yl )

# Generally if eigenvalue is less than 1,
# insufficient additional variance is accounted 
# for by adding more factors
segments( xl[1], 1, xl[2], 1, lwd = lnSz, col = 'grey' )

# Plot eigenvalues
lines( xa, ya, lwd = lnSz )
points( xa, ya, pch = 19, cex = ptSz )

# Add axes and labels
axis( 2, seq( yl[1], yl[2], .5 ), 
      tick = F, cex.axis = axSz, line = axPos )
mtext( 'Eigenvalues', side = 2, line = 2.5, cex = txtSz )
axis( 1, round( seq( xl[1] + .5, xl[2] - .5, 2) ), 
      tick = F, cex.axis = axSz, line = axPos )
mtext( 'Number of factors', side = 1, line = 2.5, cex = txtSz )

# There appear to be elbows at 2 factors and 10 factors
# respectively

# Breaks at 2 and 10 factors
pw1 = my_piecewise_reg( xa, ya, breaks = c( 2, 10 ) )
# Breaks at 3 and 10 factors
pw2 = my_piecewise_reg( xa, ya, breaks = c( 3, 10 ) )
aic = c(
  AIC( pw1$lm ),
  AIC( pw2$lm ) )
aw = round( icWeights( aic ), 2 )
# Indicates that 2 factors is a better cut-off
lines( pw1$xa, pw1$ya, lwd = lnSz, col = 'blue', lty = 2 )

title( 'Scree plot' )

# 2.4.2) RMSEA values

xa = 1:nrow( fit_details )
ya = fit_details[,'RMSEA']

# Error bars
eb = rbind(
  fit_details[,'RMSEA_UI_L'],
  fit_details[,'RMSEA_UI_U'] )

# Create blank plot
xl = c( .5, max( xa ) + .5 )
yl = c( 0, lowerUpper( .1, as.vector( eb ) )[2] )
blankPlot( xl, yl )
customAxes( xl, yl )

# Rule of thumb is an RMSEA less than .1 indicates 
# adequate fit, while an RMSEA less than 0.05 
# indicates excellent fit
segments( rep( xl[1], 2 ), c( .1, .05 ), 
          rep( xl[2], 2 ), c( .1, .05 ), 
          lwd = lnSz, col = 'grey' )

# Add error bars for 90% uncertainty intervals
errorBars( xa, eb, 
           lwd = lnSz, length = .05 )

# Add estimated RMSEA
lines( xa, ya, lwd = lnSz )
points( xa, ya, pch = 19, cex = ptSz )

# Add axes and labels
axis( 2, seq( yl[1], yl[2], .05 ), 
      tick = F, cex.axis = axSz, line = axPos )
mtext( 'RMSEA', side = 2, line = 2.5, cex = txtSz )
axis( 1, round( seq( xl[1] + .5, xl[2] - .5, 2 ) ), 
      tick = F, cex.axis = axSz, line = axPos )
mtext( 'Number of factors', side = 1, line = 2.5, cex = txtSz )

title( 'Measure of discrepancy per degree of freedom' )

# 2.4.3) AIC values

xa = 1:nrow( fit_details )
ya = fit_details[,'AIC']
# Convert to Akaike weights
ya = icWeights( ya )

# Create blank plot
xl = c( .5, max( xa ) + .5 )
yl = c( 0, 1 )
blankPlot( xl, yl )
customAxes( xl, yl )

# Add weights
segments( xa, rep( 0, length(ya) ),
          xa, ya, 
          lwd = lnSz * 1.5 )

# Add axes and labels
axis( 2, seq( 0, 1, .2 ), 
      tick = F, cex.axis = axSz, line = axPos )
mtext( 'Akaike weights', side = 2, line = 2.5, cex = txtSz )
axis( 1, round( seq( xl[1] + .5, xl[2] - .5, 2 ) ), 
      tick = F, cex.axis = axSz, line = axPos )
mtext( 'Number of factors', side = 1, line = 2.5, cex = txtSz )

title( 'Relative probability to fit new data' )

# These figures suggest the number of factors to consider 
# are 2, 4, or 7

# 2.5) Rotations/interpretations

# Variable labels
lbls = rbind(
  c('MISS_CS',      'Consumer suggest. (MISS)' ),  #  1
  c('MISS_P',       'Persuadibility (MISS)' ),     #  2
  c('MISS_PS',      'Phys. suggest. (MISS)' ),     #  3
  c('MISS_PR',      'Phys. react. (MISS)' ),       #  4
  c('MISS_PC',      'Peer conformity (MISS)' ),    #  5
  c('UPPS_NU',      'Negative urgency (UPPS-P)' ), #  6
  c('UPPS_PU',      'Positive urgency (UPPS-P)' ), #  7
  c('UPPS_Pr',      'Premediation (UPPS-P)' ),     #  8
  c('UPPS_Pe',      'Perseverance (UPPS-P)' ),     #  9
  c('UPPS_SS',      'Sensation seek. (UPPS-P)' ),  # 10
  c('DMQ_SM',       'Social motives (DMQ)' ),      # 11
  c('DMQ_Cp',       'Coping (DMQ)' ),              # 12
  c('DMQ_E',       'Enhancement (DMQ)' ),         # 13
  c('DMQ_Cn',      'Conformity (DMQ)' ),          # 14
  c('BCEOA_LC',     'Liquid courage (B-CEOA)' ),   # 15
  c('BCEOA_SP',     'Self-percept. (B-CEOA)' ),    # 16
  c('BCEOA_S',     'Sexuality (B-CEOA)' ),         # 17
  c('BCEOA_TR',     'Tension reduct. (B-CEOA)' ),  # 18
  c('Total_days',   'Total days spent drinking' ), # 19
  c('Total_drinks', 'Total drinks had' )           # 20
)

# Plot characteristics
lnSz = 1
lnSz2 = 3
axSz = 1.25
txtSz = 1.5
axPosX = -1
axPosY = -1

for ( fn in c( 2, 4, 7 ) ) {
  
  # Oblique rotation using quartimin, 
  # Dr. Cudeck's preferred method
  L = rotate_factor_loadings( all_FA[[fn]], 
                              method = 'quartimin', 
                              normalize = T )
  
  # Generate plot of factor loadings
  if ( !savePlot ) x11( width = 12 )
  layout( cbind( 3, 1, 1, 2 ) )
  par( mar = c( 4, 9, 3, 3 ) )
  LM = as.matrix( L$Loadings[,paste('Factor',1:fn,sep='') ] )
  CM = custom_color_map( LM, 
                         c(-1.05,1.05) )
  xl = CM$xl; yl = CM$yl
  
  # Add gridlines
  horizLines( 0:ncol( MV ), c(0,fn), lwd = lnSz )
  vertLines( 0:fn, yl, lwd = lnSz )
  
  # Demark factor structures
  cf = sapply( 2:fn, 
               function(x) 
                 min( which( L$Loadings$Cutoffs == x ) ) )
  horizLines( ncol( MV ) - cf + 1, c( 0, fn ), lwd = lnSz2 )
  
  # Add variable labels
  o = unlist( 
    sapply( lbls[,1], 
            function(x) which( x == L$Loadings$Variables ) ) )
  axis( 2, 1:nrow( L$Loadings ) - .5,
        lbls[o,2],
        tick = F, line = axPosY, cex.axis = axSz, las = 1 )
  
  # Add factor dimensions
  axis( 1, 1:fn - .5, 1:fn, 
        tick = F, line = axPosX, cex.axis = axSz )
  mtext( 'Factor', side = 1, line = 2, cex = txtSz )
  
  title( paste( fn, 'Factor model' ), cex = txtSz )
  
  blankPlot()
  legend( 'left',
          as.character(
            seq( -1, 1, length = length( CM$cmap$clr ) ) ),
          fill = CM$cmap$clr,
          bty = 'n', cex = txtSz )
  blankPlot()
  
}

if ( savePlot ) dev.off()
setwd( orig_dir )