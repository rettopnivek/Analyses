#--------------------#
# Useful functions   #
# Kevin Potter       #
# Updated 05/24/2017 #
#--------------------#

# Index
# Lookup - 01:  dexgauss
# Lookup - 02:  est_exgauss_param
# Lookup - 03:  mixture_tran_par
# Lookup - 04:  mle_mixture_model
# Lookup - 05:  customAxes
# Lookup - 06:  blankPlot
# Lookup - 07:  logistic

# Lookup - 01
dexgauss = function( x, mu, sigma, lambda, log = F ) {
  # Purpose:
  # Computes the density (or log-density) of the 
  # exponentially modified Gaussian distribution
  # (i.e., when a random variable T is the sum of 
  #  a normally distributed random variable D and 
  #  an exponentially distributed random variable R;
  #  T = D + R).
  # Arguments:
  # x      - A vector of response times (can be negative)
  # mu     - The mean of the normally distributed component
  # sigma  - The standard deviation of the normally 
  #          distributed component
  # lambda - The rate parameter for the exponentially 
  #          distributed component
  # log    - Logical; if true, returns the log-likelihood
  # Returns:
  # A vector of likelihoods (or log-likelihoods).
  
  # Vectorize inputs and match vector lengths
  n = sapply( list( x, mu, sigma, lambda ), length )
  mx = max( n )
  x = rep_len( x, mx )
  mu = rep_len( mu, mx )
  sigma = rep_len( sigma, mx )
  lambda = rep_len( lambda, mx )
  
  # Compute density for ex-Gaussian
  sigma2 = sigma^2; # Square standard deviation
  p1 = lambda/2.0;
  p2 = 2.0*mu + lambda*sigma2 - 2.0*x;
  p3 = ( mu + lambda*sigma2 - x )/(sqrt(2)*sigma);
  out = p1*exp(p1*p2)*(1-erf(p3));
  
  # If indicated, compute log-density
  if ( log ) out = log( out );
  
  return( out );
}

# Lookup - 02
est_exgauss_param = function( dat ) {
  # Purpose:
  # Computes the maximum likelihood estimates 
  # for the ex-Gaussian distribution given a 
  # set of response times via the method of 
  # moments.
  # Arguments:
  # rt - A vector of response times
  # Returns:
  # A vector with the estimates for the parameters 
  # mu, sigma, and lambda.
  
  # Compute 1st, 2nd, and 3rd moments ( Terriberry, 2007).
  n = 0; m = 0; m2 = 0; m3 = 0;
  for ( k in 1:length( dat ) ) {
    n1 = n; n = n + 1
    term1 = dat[k] - m; term2 = term1/ n;
    term3 = term1 * term2 * n1
    # Update first moment
    m = m + term2
    # Update third moment
    m3 = m3 + term3 + term2 * (n - 2) - 3 * term2 * m2
    # Update second moment
    m2 = m2 + term3
  }
  # Compute standard deviation of sample
  s = sqrt( m2 / ( n - 1.0 ) )
  # Compute skewness of sample
  y = sqrt( n ) * m3 / ( m2^1.5 );
  # Estimate parameters
  mu_hat = m - s * ( y/2 )^(1/3)
  sigma_hat = sqrt( (s^2) * ( 1 - (y/2)^(2/3) ) )
  lambda_hat = 1/( s * (y/2)^(1/3) )
  return( c( mu = mu_hat, sigma = sigma_hat, 
             lambda = lambda_hat ) )
}

# Lookup - 03
mixture_tran_par = function( prm, reverse = F ) {
  # Purpose: 
  # Transforms the parameter inputs for the 
  # mixture of the uniform and ex-Gauss, 
  # thereby ensuring the values are properly 
  # bounded.
  # Arguments: 
  # prm     - A vector of four parameters, the 
  #           mean, standard deviation, rate, 
  #           and mixture probability of the 
  #           ex-Gaussian
  # reverse - Logical; if true, reverses the 
  #           transformation to produce 
  #           unbounded values
  # Returns: 
  # A vector of four transformed parameters.
  
  if ( !reverse ) {
    prm[2:3] = exp( prm[2:3] )
    prm[4] = logistic( prm[4] )
  } else {
    prm[2:3] = log( prm[2:3] )
    prm[4] = logit( prm[4] )
  }
  
  return( prm )
}

# Lookup - 04
mle_mixture_model = function( x, prm, sum = T ) {
  # Purpose: 
  # A function to compute the sum of the log-likelihoods 
  # for a mixture of the ex-Gaussian and a uniform 
  # distribution. The function may be passed into the 
  # 'optim' function to carry out maximum likelihood 
  # estimation.
  # Arguments: 
  # x   - A vector of response times
  # prm - A vector of four unbounded parameters, 
  #       the mean, raw standard deviation, raw rate, 
  #       and raw mixture probability of the exponential 
  #       Gaussian.
  # sum - Logical; if true, the sum of log-likelihoods is 
  #       returned.
  #       
  # Returns: 
  # Either returns the sum of the log-likelihoods, or 
  # returns a list with the likelihoods for the 
  # ex-Gaussian, uniform, and weighted sum of the 
  # pair of likelihoods.
  
  # Transform parameters so that they are 
  # appropriately bounded
  prm = mixture_tran_par( prm )
  
  # Extract parameters
  mu = prm[1]; sigma = prm[2]; 
  lambda = prm[3]; theta = prm[4];
  
  # Compute likelihoods for ex-Gaussian
  dexg = dexgauss( x, mu, sigma, lambda )
  # Compute likelihoods for the uniform
  du = dunif( x, min(x), max(x) )
  
  # Compute likelihoods for mixture model
  d = theta * dexg + ( 1 - theta ) * du
  
  if ( sum ) {
    sll = sum( log( d ) )
    if ( is.na( sll ) ) sll = -Inf
    return( sll )
  } else {
    return( list( dexg = dexg, du = du, d = d ) )
  }
}

# Lookup - 05
customAxes = function( xl, yl,
                       pos = c( 1, 2 ),
                       lnSz = 2,
                       type = 2,
                       label = NULL,
                       lbPos = c( 2.5, 2.5, 1, 2.5 ),
                       lbSz = 1.5,
                       inc = c(0,0,0,0),
                       axisLabel = NULL,
                       axSz = 1.25,
                       axPos = c( -1, -1, -1, -1 ),
                       prec = 2 ) {
  
  # Check whether labels for axes should be added
  labCheck = F
  if ( !is.null( label ) ) {
    if ( length( label ) != length( pos ) ) {
      warning( "Number of labels does not match number of specified axes" )
    } else labCheck = T
  }
  
  alCheck = F
  if ( !is.null( axisLabel ) ) {
    if ( !is.list( axisLabel ) ) {
      warning( "A list of axis labels is needed" )
    } else if ( length( axisLabel ) != length( pos ) ) {
      warning( "Number of axis labels does not match number of specified axes" )
    } else alCheck = T
  }
  
  for ( i in 1:length( pos ) ) {
    
    # Bottom axis
    if ( pos[i] == 1 | pos[i] == '1' | pos[i] == 'B' |
         pos[i] == 'b' | pos[i] == 'bottom' |
         pos[i] == 'Bottom' ) {
      
      # Draw axis line
      if ( type == 1 | type == 'extend' )
        abline( h = yl[1], lwd = lnSz )
      if ( type == 2 | type == 'truncate' )
        segments( xl[1], yl[1], xl[2], yl[1], lwd = lnSz )
      
      # Add label to axis
      if ( labCheck ) {
        mtext( label[i], side = 1, cex = lbSz, line = lbPos[i] )
      }
      
      # Draw axis tick labels
      
      # Increments
      if ( inc[i] > 0 ) {
        ai = seq( xl[1], xl[2], inc[i] )
        ai = round( ai, prec )
        axis( 1, ai, tick = F, line = axPos[i],
              cex.axis = axSz )
      }
      
      # String labels
      if ( alCheck ) {
        if ( !is.null( axisLabel[[i]] ) ) {
          axis( 1, axisLabel[[i]][[1]], axisLabel[[i]][[2]],
                tick = F, line = axPos[i],
                cex.axis = axSz )
        }
      }
    }
    
    # Left axis
    if ( pos[i] == 2 | pos[i] == '2' | pos[i] == 'L' |
         pos[i] == 'l' | pos[i] == 'left' |
         pos[i] == 'Left' ) {
      
      # Draw axis line
      if ( type == 1 | type == 'extend' )
        abline( v = xl[1], lwd = lnSz )
      if ( type == 2 | type == 'truncate' )
        segments( xl[1], yl[1], xl[1], yl[2], lwd = lnSz )
      
      # Add axis label
      if ( labCheck ) {
        mtext( label[i], side = 2, cex = lbSz, line = lbPos[i] )
      }
      
      # Draw axis tick labels
      
      # Increments
      if ( inc[i] > 0 ) {
        ai = seq( yl[1], yl[2], inc[i] )
        ai = round( ai, prec )
        axis( 2, ai, tick = F, line = axPos[i],
              cex.axis = axSz )
      }
      
      # String labels
      if ( alCheck ) {
        if ( !is.null( axisLabel[[i]] ) ) {
          axis( 2, axisLabel[[i]][[1]], axisLabel[[i]][[2]],
                tick = F, line = axPos[i],
                cex.axis = axSz )
        }
      }
    }
    
    # Top axis
    if ( pos[i] == 3 | pos[i] == '3' | pos[i] == 'T' |
         pos[i] == 't' | pos[i] == 'top' |
         pos[i] == 'Top' ) {
      
      # Draw axis line
      if ( type == 1 | type == 'extend' )
        abline( h = yl[2], lwd = lnSz )
      if ( type == 2 | type == 'truncate' )
        segments( xl[1], yl[2], xl[2], yl[2], lwd = lnSz )
      
      # Add axis label
      if ( labCheck ) {
        mtext( label[i], side = 3, cex = lbSz, line = lbPos[i] )
      }
      
      # Draw axis tick labels
      
      # Increments
      if ( inc[i] > 0 ) {
        ai = seq( xl[1], xl[2], inc[i] )
        ai = round( ai, prec )
        axis( 3, ai, tick = F, line = axPos[i],
              cex.axis = axSz )
      }
      
      # String labels
      if ( alCheck ) {
        if ( !is.null( axisLabel[[i]] ) ) {
          axis( 3, axisLabel[[i]][[1]], axisLabel[[i]][[2]],
                tick = F, line = axPos[i],
                cex.axis = axSz )
        }
      }
      
      
    }
    
    # Right axis
    if ( pos[i] == 4 | pos[i] == '4' | pos[i] == 'R' |
         pos[i] == 'r' | pos[i] == 'right' |
         pos[i] == 'Right' ) {
      
      # Draw axis line
      if ( type == 1 | type == 'extend' )
        abline( v = xl[2], lwd = lnSz )
      if ( type == 2 | type == 'truncate' )
        segments( xl[2], yl[1], xl[2], yl[2], lwd = lnSz )
      
      # Add axis label
      if ( labCheck ) {
        mtext( label[i], side = 4, cex = lbSz, line = lbPos[i] )
      }
      
      # Draw axis tick labels
      
      # Increments
      if ( inc[i] > 0 ) {
        ai = seq( yl[1], yl[2], inc[i] )
        ai = round( ai, prec )
        axis( 4, ai, tick = F, line = axPos[i],
              cex.axis = axSz )
      }
      
      # String labels
      if ( alCheck ) {
        if ( !is.null( axisLabel[[i]] ) ) {
          axis( 1, axisLabel[[i]][[1]], axisLabel[[i]][[2]],
                tick = F, line = axPos[i],
                cex.axis = axSz )
        }
      }
      
    }
    
  }
  
}

# Lookup - 06
blankPlot = function( xDim = c(0,1), yDim = c(0,1) ) {
  
  plot( xDim, yDim, type = 'n', ylab = ' ', xlab = ' ',
        xaxt = 'n', yaxt = 'n', bty = 'n' )
  
}

# Lookup - 07
logistic = function(x) {
  return( 1/(1+exp(-x)) )
}