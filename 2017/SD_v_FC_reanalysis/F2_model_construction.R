#----------------------------------#
# Functions for model construction #
# Kevin Potter                     #
# Updated 07/16/2017               #
#----------------------------------#

# Index
# Lookup - 01:   General purpose functions
# Lookup - 01a:  dot_product
# Lookup - 01b:  p_c
# Lookup - 01c:  extract_cnd
# Lookup - 01d:  qs
# Lookup - 01e:  expand_mat
# Lookup - 01f:  dtnorm
# Lookup - 01g:  bplt
# Lookup - 01h:  compute_rmse
# Lookup - 01i:  find_min_rt
# Lookup - 01j:  resid_bootstrap
# Looku[ - 01k:  estimate

# Lookup - 02:   Diffusion race model functions
# Lookup - 02a:  drm_k
# Lookup - 02b:  drm_x
# Lookup - 02c:  drm_t
# Lookup - 02d:  drm_s
# Lookup - 02e:  drm_pm
# Lookup - 02f:  drm_pm_quick_convert
# Lookup - 02g:  drm_dtbf_create

# Lookup - 03:   Two-boundary wiener process model functions
# Lookup - 03a:  wp_a
# Lookup - 03b:  wp_z
# Lookup - 03c:  wp_v
# Lookup - 03d:  wp_t
# Lookup - 03e:  wp_pm
# Lookup - 03f:  wp_pm_quick_convert
# Lookup - 03g:  wp_dtbf_create

###
### General purpose functions
###
# Lookup - 01

# Lookup - 01a
cppFunction("
NumericVector dot_product( NumericMatrix X, 
                           NumericVector p ) {
  /*
  Purpose: 
  Computes the dot product of a subsets of equal length 
  for vectors x and p.
  Arguments: 
  X  - A design matrix
  p  - A vector of parameters
  Returns: 
  A scalar value.
  */
  
  // Initialize variables
  NumericVector out( X.nrow() );
  
  for ( int r = 0; r < X.nrow(); r++ ) {
    out(r) = sum( X(r,_) * p );
  }
  
  return out;
}
")

# Lookup - 01b
p_c = function( X, prm, ind ) {
  # Purpose: 
  # Given a design matrix and a parameter vector, computes 
  # the parameter subset over observations.
  # Arguments: 
  # X   - A design matrix
  # prm - A vector of parameters
  # ind - A list consisting of...
  #       ix = The subset of columns for the design matrix
  #       ip = The subset of parameters
  #       is = The subset of parameters to rescale
  #       sp = The subset of parameters determining the 
  #            degree of rescaling
  # Returns: 
  # A vector of parameter values. 
  
  if ( ind$sm ) {
    if ( !is.null( ind$is ) ) {
      p = prm[ ind$ip ]
      p[ ind$is ] = p[ ind$is ] * prm[ ind$sp ]
      out = dot_product( X[,ind$ix], p )
    } else {
      out = dot_product( X[,ind$ix], prm[ind$ip] )
    }
  } else out = prm[ ind$ip ][ X[,ind$ix] ]
  
  return( out )
}

# Lookup - 01c
extract_cnd = function( d, task = 'Both' ) {
  # Purpose: 
  # Extracts the combinations of the conditions 
  # along with their labels from the data.
  # Arguments: 
  # d - The data frame with the data
  # Returns: 
  # A list consisting of...
  # Cnd = A data frame listing the unique combinations 
  #       of conditions
  # SH = A data frame with shorthand labels
  
  # Extract conditions
  Cnd = aggregate( d$Ac, list( d$PDL, d$PTL, d$CoL, d$TaL ), 
                   catProp, val = c('0','1') )
  Cnd = Cnd[ , 1:4 ]
  colnames( Cnd ) = c( 'PD', 'PT', 'Co', 'Ta' )
  
  if ( task == 'Forced-choice' ) {
    Cnd = Cnd[ Cnd$Ta == 'Forced-choice', ]
  }
  if ( task == 'Same-different' ) {
    Cnd = Cnd[ Cnd$Ta == 'Same-different', ]
  }
  
  # Create shorthand labels for conditions
  PD_SH = rep( 'S', nrow( Cnd ) )
  PD_SH[ Cnd$PD == .4 ] = 'L'
  PT_SH = rep( 'FP', nrow( Cnd ) )
  PT_SH[ Cnd$PT == 'Target-primed' ] = 'TP'
  Co_SH = rep( 'L', nrow( Cnd ) )
  Co_SH[ Cnd$Co == 'Right' ] = 'R'
  Co_SH[ Cnd$Co == 'Same' ] = 'S'
  Co_SH[ Cnd$Co == 'Different' ] = 'D'
  Ta_SH = rep( 'FC', nrow( Cnd ) )
  Ta_SH[ Cnd$Ta == 'Same-different' ] = 'SD'
  
  mp = cbind( PD = PD_SH, PT = PT_SH, 
              Co = Co_SH, Ta = Ta_SH )
  mp = as.data.frame( mp )
  
  return( list( Cnd = Cnd, SH = mp ) )
}

# Lookup - 01d
qs = function( string, ms ) {
  # Purpose: 
  # A function that allows for quick selection of 
  # unique combinations of conditions via a 
  # shorthand string of labels.
  # Arguments: 
  # string - A string, such as 
  #          'Co_L,Ta_FC'
  #          o Commas separate different variables
  #          o Underscores separate variable names 
  #            from the associated value to select
  # ms     - A list with the data frames containing 
  #          the condition labels.
  # Returns: 
  # A logical vector.
  
  vrb = strsplit( string, split = ',' )[[1]]
  
  K = length( vrb )
  which_true = matrix( NA, nrow( ms$SH ), K )
  for ( k in 1:K ) {
    
    val = strsplit( vrb[k], '_' )[[1]]
    which_true[,k] = ms$SH[,val[1]] == val[2]
  }
  
  out = as.vector( apply( which_true, 1, all ) )
  
  return( out )
}

# Lookup - 01e
expand_mat = function( mat, Cnd, df ) {
  # Purpose: 
  # A function that expands the matrix of parameter 
  # indices per condition into a matrix of parameter 
  # indices for every observation.
  # Arguments: 
  # mat - A matrix of parameter indices per condition
  # Cnd - The data frame of condition labels
  # df  - The data frame with the observations
  # Returns: 
  # An expanded matrix of parameter indices.
  
  # Determine number of observations
  N = nrow( df )
  
  # Create output
  output = matrix( NA, N, ncol( mat ) )
  
  for ( i in 1:nrow( mat ) ) {
    
    sel = 
      df$PDL == Cnd$PD[i] & 
      df$PTL == Cnd$PT[i] & 
      df$CoL == Cnd$Co[i] & 
      df$TaL == Cnd$Ta[i]
    
    for ( k in 1:ncol( mat ) ) {
      output[ sel, k ] = mat[ i, k ]
    }
    
  }
  colnames( output ) = colnames( mat )
  
  return( output )
}

# Lookup - 01f
dtnorm = function( x, m, sd, a, b, log = F ) {
  # Purpose: 
  # Computes the density function for a 
  # truncated normal distribution.
  # Arguments: 
  # x   - A vector of quantiles
  # m   - The mean of the normal distribution
  # sd  - The standard deviation of the normal distribution
  # a   - The lower bound for truncation (can be -Inf)
  # b   - The upper bound for truncation (can be +Inf)
  # log - Logical; if true, the log of the density is returned
  # Returns: 
  # The likelihood or log-likelihood for a truncated normal 
  # distribution.
  
  num = dnorm( x, m = m, sd = sd )
  num[ x < a | x > b ] = 0 # Truncate
  denom = pnorm( b, m = m, sd = sd ) - 
    pnorm( a, m = m, sd = sd ) # Renormalize
  out = num/denom
  
  if ( log ) out = log( out )
  
  return( out )
}

# Lookup - 01g
bplt = function( x, pos = 1, width = .5, 
                 f = NULL, 
                 return = F, 
                 lwd = c(1,1), 
                 lty = c(1,2), 
                 lcol = c('black','black'), 
                 pcol = 'black',
                 bcol = 'white',
                 pch = 19,
                 cex = 1, 
                 ... ) {
  # Purpose: 
  # A function that generates a custom boxplot.
  # Arguments: 
  # x      - A vector of values
  # pos    - The position on the x-axis at which to 
  #          draw the boxplot
  # wdth   - The width of the box
  # f      - An optional function to generate the 5 
  #          y-axis values (the lower end-point, the 
  #          bottom of the box, the midpoint, the top
  #          of the box, and the upper end-point. 
  #          Additional values are treated as extra 
  #          midpoints to draw
  # return - A logical value; if true, the y-axis values 
  #          are returned
  # lwd    - The width of the lines
  # lty    - The type of line
  # lcol   - The color of the line
  # pcol   - The border color of the outlier points
  # bcol   - The background color of the outlier points
  # pch    - The type of outlier point to draw
  # cex    - The size of the outlier points
  # ...    - Additional parameters for the 'polygon' 
  #          function
  # Returns: 
  # The y-axis values.
  
  if ( is.null( f ) ) {
    
    # Summary statistics for x
    sx = list(
      min = min( x ),
      max = max( x ),
      Q1 = quantile( x, .25 ),
      Q2 = quantile( x, .5 ),
      Q3 = quantile( x, .75 ),
      mean = mean( x )
    )
    # Inter-quartile range
    sx$IQR = sx$Q3 - sx$Q1
    # Fences for outliers
    sx$LIF = sx$Q1 - 1.5 * sx$IQR
    sx$UIF = sx$Q3 + 1.5 * sx$IQR
    sx$LOF = sx$Q1 - 3 * sx$IQR
    sx$UOF = sx$Q3 + 3 * sx$IQR
    
    ya = numeric( 6 )
    names( ya ) = c( 'Q0','Q1','Q2','Q3','Q4', 'O1' )
    ya['Q0'] = max( sx$min, sx$LIF )
    ya['Q1'] = sx$Q1
    ya['Q2'] = sx$Q2
    ya['Q3'] = sx$Q3
    ya['Q4'] = min( sx$max, sx$UIF )
    ya['O1'] = sx$mean
    
  } else {
    ya = f( x )
  }
  l = length( ya )
  
  if ( l == 5 )
    names( ya ) = c( 'Q0','Q1','Q2','Q3','Q4' )
  if ( l > 5 ) 
    names( ya ) = c( 'Q0','Q1','Q2','Q3','Q4',
                     paste( 'O', 1:( l - 5 ), 
                                     sep = '' ) )
  
  # Box: (+Q1, -Q1) to (+Q3, -Q3)
  ybox = ya[ c( 'Q1', 'Q3', 'Q3', 'Q1' ) ]
  xbox = pos + c( -width/2, -width/2,
                  width/2, width/2 )
  polygon( xbox, ybox, 
           lwd = lwd[1],
           lty = lty[1],
           border = lcol[1], ... )
  
  # Endpoints: Q0 to Q1, Q3 to Q4
  segments( c( pos, pos ),
            ya[ c('Q0','Q3') ],
            c( pos, pos ),
            ya[ c('Q1','Q4') ], 
            lty = lty[1],
            lwd = lwd[1],
            col = lcol[1] )
  
  # Midpoint: Q2
  segments( pos - width/2, ya['Q2'],
            pos + width/2, ya['Q2'],
            lty = lty[1],
            lwd = lwd[1],
            col = lcol[1] )
  
  # Optional midpoints:
  if ( l > 5 ) {
    for ( j in 6:l ) {
      segments( pos - width/2, ya[j],
                pos + width/2, ya[j],
                lty = lty[2],
                lwd = lwd[2],
                col = lcol[2] )
    }
  }
  
  # Outliers
  sel = x < ya['Q0'] | x > ya['Q4']
  if ( any( sel ) ) {
    points( rep( pos, sum( sel ) ), x[sel],
            pch = pch, col = pcol, bg = bcol,
            cex = cex )
  }
  
  if ( return )
    return( ya )
}

# Lookup - 01h
compute_rmse = function( dtbf, fit, mle_fn ) {
  # Purpose: 
  # Computes the root-mean-square error (RMSE) based on the 
  # predicted and observed joint cumulative distribution 
  # function for a set of data and model fits.
  # Arguments: 
  # dtbf   - The data to be fitted
  # fit    - A fit object from the 'MLE' function
  # mle_fn - The maximum likelihood function (requires a 'predict'
  #          option)
  # Returns: 
  # A list with the computed RMSE and a data frame with 
  # the inputs used to do so.
  
  # Create data frame for output
  rmse = data.frame( rt = dtbf$rt, ch = dtbf$ch,
                     p = 0, phat = 0 )
  
  # Function to compute joint ECDF
  f = function( x ) {
    o1 = order( x$rt[ x$ch == 1] )
    o0 = order( x$rt[ x$ch == 0 ] )
    out = numeric( nrow(x) );
    if ( length( o1 ) > 0 ) {
      p1 = 1:sum( x$ch == 1 )/nrow(x)
      out[ x$ch == 1 ][o1] = p1
    }
    if ( length( o0 ) > 0 ) {
      p0 = 1:sum( x$ch == 0 )/nrow(x)
      out[ x$ch == 0 ][o0] = p0
    }
    return( out )
  }
  
  # Loop over conditions
  for ( i in 1:nrow( dtbf$ms$Cnd ) ) {
    sel = dtbf$PD == dtbf$ms$Cnd$PD[i] & 
      dtbf$PT == dtbf$ms$Cnd$PT[i] & 
      dtbf$Co == dtbf$ms$Cnd$Co[i]
    rmse$p[sel] = f(rmse[sel,])
  }
  
  # Compute predicted joint CDF
  prm = coef( fit )[ which.max( fit$value ), ]
  rmse$phat = mle_fn( prm, dtbf, predict = T )
  
  out = list( rmse = sqrt( mean( (rmse$phat - rmse$p)^2 ) ), 
              input = rmse )
  
  return(out)
}

# Lookup - 01i
find_min_rt = function( mat, ms, dat, choice = F ) {
  # Purpose: 
  # Computes the minimum response time over a desired 
  # set of conditions separately for choices equal to 1
  # and 0.
  # Arguments: 
  # mat    - Matrix to indicate over which conditions to 
  #          compute the minimum response time
  # ms     - A list with two data-frames giving the unique 
  #          combinations of conditions and their shorthand 
  #          labels
  # dat    - The data to be fitted
  # choice - Logical; if true, each choice option is 
  #          assumed to have a separate residual latency
  # Returns: 
  # A list giving the minimum response times either with 
  # matching length to the number of observations or to 
  # the number of unique combinations of conditions.
  
  all_mat = expand_mat( mat, ms$Cnd, dat )
  
  out = c()
  
  for ( j in 1:2 ) {
    
    out = c( out, 
             list( 
               list( 
                  T = matrix( NA, nrow( all_mat ), 1 ),
                  t = matrix( NA, nrow( mat ), 1 ) )
               )
            )
    
    for ( i in unique( mat[,j] ) ) {
      sel1 = all_mat[,j] == i
      sel2 = mat[,j] == i
      out[[j]]$t[sel2,1] = 
        rep( min(dat$RT[ dat$Ch == (j-1) & sel1 ]), sum( sel2 ) )
    }
    out[[j]]$T = expand_mat( out[[j]]$t, ms$Cnd, dat )
  }
  names( out ) = c( 't0','t1' )
  
  if ( !choice ) {
    
    sel = out[[2]]$T < out[[1]]$T
    out[[1]]$T[sel] = out[[2]]$T[sel]
    out[[2]]$T[!sel] = out[[1]]$T[!sel]
    
    sel = out[[2]]$t < out[[1]]$t
    out[[1]]$t[sel] = out[[2]]$t[sel]
    out[[2]]$t[!sel] = out[[1]]$t[!sel]
  }
  
  return( out )
}

# Lookup - 01j
cppFunction("
NumericMatrix resid_bootstrap( NumericMatrix resid, int nRep ) {
  /*
  Purpose:
  A function to that creates a bootstrap sample of the 
  means for a set of residuals via resampling with replacement.
  Arguments:
  resid - A matrix of r residuals over c conditions
  nRep  - The number of bootstrap samples
  Returns:
  A matrix of means over residuals.
  */
  
  // Matrix dimensions
  int R = resid.nrow(); // Number of rows
  int C = resid.ncol(); // Number of columns
  
  // Output
  NumericMatrix out( nRep, C );
  
  // Loop over columns
  for ( int c = 0; c < C; c++ ) {
    
    // Determine if there are any missing values
    int n_na = sum( is_na( resid(_,c) ) );
    // Non-missing sample size
    int N = R - n_na;
    // Vector for computing mean
    NumericVector x( N );
    // Determine non-missing indices
    IntegerVector sample_index( N );
    int inc = 0;
    for ( int r = 0; r < R; r++ ) {
      if ( !NumericVector::is_na(resid(r,c)) ) {
        sample_index( inc ) = r;
        inc += 1;
      }
    }
    
    // Loop over bootstrap iterations
    for ( int nr = 0; nr < nRep; nr++ ) {
      
      IntegerVector ind(N);
      ind = runif(N,0,N-1);
      for ( int n = 0; n < N; n++ ) {
        x(n) = resid( sample_index(ind(n)), c );
      }
      out(nr,c) = mean( x );
    }
  }
  
  return out;
}
")

# Lookup - 01k
estimate = function( ver, subject, task, 
                     data_create_function, 
                     mthd = NULL,
                     n_attempts = 10 ) {
  # Purpose: 
  # A convenience function for maximum likelihood estimation 
  # using the 'optimx' function.
  # Arguments: 
  # ver                  - A vector specifying the structure 
  #                        for each of the parameter types
  # subject              - The subject index for the data to 
  #                        be fitted
  # task                 - The task index for the data to be 
  #                        fitted
  # data_create_function - The function to use to create the 
  #                        list of data to be fitted
  # mthd                 - An optional vector of estimation 
  #                        algorithms to use
  # n_attempts           - The number of iterations to 
  #                        attempt to find parameter estimates
  # Returns: 
  # The estimation output from the 'optimx' function. 
  
  # Data to be fitted ('dtbf' object must already exit)
  dtbf <<- data_create_function( ver, subject, task = task, d = d )
  
  # Model estimation (using multiple estimation methods)
  
  # Default methods to attempt
  if ( is.null( mthd ) ) {
    mthd = c(
      'Nelder-Mead', # SIMPLEX method from Nelder and Mead (1965)
      'BFGS', # Broyden-Fletcher-Goldfarb-Shanno algorithm (1970)
      'nlm', # Non-linear minimization from Schnabel, Koontz,
      # and Weiss (1985)
      'nlminb', # Forthcoming
      'nmkb' # Bounded Nelder-Mead optimization by Varadhan (?)
    )
  }
  
  # Initial estimation
  success = F; attempts = 0
  # Iterate until successful estimation
  while ( !success ) {
    
    # Generate initial starting values
    chk = -Inf; inc = 0
    while ( chk == -Inf | inc < 20 ) {
      start = st_f();
      chk = mle_fn( start, dtbf )
      inc = inc + 1;
    }
    names( start ) = dtbf$pnames
    
    # Attempt to obtain initial estimates
    init_fit = tryCatch(
      optimx( start, mle_fn, dat = dtbf, 
              method = mthd,
              control = list( maximize = T, 
                              maxit = 1000 ) ),
      error = function(e) return( NULL ) )
    
    # Update iteration exit conditions
    if ( !is.null( init_fit ) ) success = T
    attempts = attempts + 1;
    if ( attempts > n_attempts - 1 ) success = T
  }
  
  # Initialize output
  fit = NULL
  
  # If initial estimates were obtained carry out 
  # final estimation attempt
  if ( !is.null( init_fit ) ) {
    
    # Final estimation
    success = F; attempts = 0
    # Iterate until successful estimation
    while ( !success ) {
      
      # Use initial estimates as starting values 
      # and scaling for parameters
      scl = coef( init_fit )[ which.max( init_fit$value ), ]
      
      # Extract initial estimates
      scl = coef( init_fit )[ which.max( init_fit$value ), ]
      
      # Pertube values slightly for new starting values
      chk = -Inf; inc = 0
      while ( chk == -Inf | inc < 20 ) {
        start = start = scl + runif( length( scl ), -.5, .5 )
        chk = mle_fn( start, dtbf )
        inc = inc + 1;
      }
      names( start ) = dtbf$pnames
      
      # Attempt to obtain final estimates
      fit = tryCatch( 
        optimx( start, mle_fn, dat = dtbf, 
                method = mthd,
                control = list( maximize = T, 
                                maxit = 20000,
                                parscale = abs( scl ) ) ),
        error = function(e) return( NULL ) )
      
      # Update iteration exit conditions
      if ( !is.null( fit ) ) success = T
      attempts = attempts + 1;
      if ( attempts >  n_attempts - 1 ) success = T
      
    }
    
  } else {
    # Return error message
    cat( 'Initial estimation failed' )
    stop( '', call. = FALSE )
  }
  
  if ( is.null( fit ) ) {
    # Return error message
    cat( 'Final estimation failed' )
    stop( '', call. = FALSE )
  }
  
  return( fit )
}

###
### Diffusion race model functions
###
# Lookup - 02

# Lookup - 02a
drm_k = function( ver, ms, dat, start ) {
  # Purpose:
  # A function to create different types of parameter 
  # mappings for the threshold parameters of the diffusion 
  # race model.
  # Arguments:
  # ver   - The type of parameter mapping to create
  # ms    - A list with two data frames giving the condition 
  #         labels or their shorthand representation
  # dat   - The data frame with the observations and conditions
  # start - The shift for the parameter and column indices 
  #         respectively.
  # Returns:
  # A list with...
  # X         = A matrix of parameter indices for all observations
  # Xs        = A matrix of parameter indices for the unique 
  #             combination of conditions
  # ind       = A list with...
  #             ip) The subset of parameters for a particular 
  #                 parameter type
  #             ix) The subset of columns for X or Xs associated with 
  #                 a particular parameter type
  #             is) A subset of subsets, the parameters to be 
  #                 rescaled (Optional)
  #             sp) The indices for the rescaling coefficients
  #                 (Optional)
  # nms       = The labels for the parameters
  # end       = The new shift to be applied to subsequent parameter 
  #             and column indices
  # ub        = The upper boundary for the randomized starting points
  # pscl      = The scaling of the parameters for the 'optim' function
  # start_rng = A matrix whose first row gives the lower boundary 
  #             and second row gives the upper boundary for 
  #             use with random starting value generation
  # type      = The type of parameter transformation to apply to 
  #             the parameters
  # priors    = An optional of the prior values to impose on the 
  #             parameters
  
  if ( ver == 'Intercept' ) {
    # Design matrices
    Xs = matrix( 1, nrow( ms$Cnd ), 1 )
    X = expand_mat( Xs, ms$Cnd, dat )
    # Parameter/column indices
    ind = list(
      k1 = list(
        ip = 1 + start[1],
        ix = 1 + start[2],
        sm = F ),
      k0 = list(
        ip = 1 + start[1],
        ix = 1 + start[2],
        sm = F )
    )
    # Transformation type
    type = 1
    # Starting value generation
    start_rng = rbind( log( .5 ), 
                       log( 4 ) )
    # Parameter scaling
    pscl = log( 1.3 )
    # Labels
    nms = c( 'kappa' )
    # Priors
    priors = list()
  }
  
  if ( ver == 'Bias' ) {
    # Design matrices
    Xs = matrix( c(2,1), nrow( ms$Cnd ), 2, byrow = T )
    X = expand_mat( Xs, ms$Cnd, dat )
    # Parameter/column indices
    ind = list(
      k1 = list(
        ip = 1:2 + start[1],
        ix = 1 + start[2],
        sm = F ),
      k0 = list(
        ip = 1:2 + start[1],
        ix = 2 + start[2],
        sm = F )
    )
    # Transformation type
    type = rep( 1, 2 )
    # Starting value generation
    start_rng = rbind( rep( log( .5 ), 2 ), 
                       rep( log( 4 ), 2 ) )
    # Parameter scaling
    pscl = rep( log( 1.3 ), 2 )
    # Labels
    nms = paste( 'kappa', 1:2, sep = '-' )
    # Priors
    priors = list()
  }
  
  if ( ver == 'Prime_match' ) {
    # Design matrices
    Xs = cbind( 
      c( 5, 6, 7, 8, 7, 8, 5, 6 ), # Different
      c( 3, 4, 1, 2, 1, 2, 3, 4 )  # Same
    )
    X = expand_mat( Xs, ms$Cnd, dat )
    # Parameter/column indices
    ind = list(
      k1 = list(
        ip = 1:8 + start[1],
        ix = 1 + start[2],
        sm = F ),
      k0 = list(
        ip = 1:8 + start[1],
        ix = 2 + start[2],
        sm = F )
    )
    # Transformation type
    type = rep( 1, 8 )
    # Starting value generation
    start_rng = rbind( rep( log( .5 ), 8 ), 
                       rep( log( 4 ), 8 ) )
    # Parameter scaling
    pscl = log( 1.3 )
    # Labels
    nms = paste( 'kappa', 1:8, sep = '-' )
    # Priors
    priors = list()
  }
  
  if ( ver == 'Bias_x_duration' ) {
    # Design matrices
    Xs = cbind( rep( 3:4, 4 ), rep( 1:2, 4 ) )
    X = expand_mat( Xs, ms$Cnd, dat )
    # Parameter/column indices
    ind = list(
      k1 = list(
        ip = 1:4 + start[1],
        ix = 1 + start[2],
        sm = F ),
      k0 = list(
        ip = 1:4 + start[1],
        ix = 2 + start[2],
        sm = F )
    )
    # Transformation type
    type = rep( 1, 4 )
    # Starting value generation
    start_rng = rbind( rep( log( .5 ), 4 ), 
                       rep( log( 4 ), 4 ) )
    # Parameter scaling
    pscl = log( 1.3 )
    # Labels
    nms = paste( 'kappa', 1:4, sep = '-' )
    # Priors
    priors = list()
  }
  
  if ( ver == 'Saturated' ) {
    # Design matrices
    Xs = matrix( 1, nrow( ms$Cnd ), 2 )
    Xs[ qs( 'Ta_FC', ms ), 2 ] = 1:8 # Left
    Xs[ qs( 'Ta_FC', ms ), 1 ] = c( 5:8, 1:4 ) + 8 # Right
    Xs[ qs( 'Ta_SD', ms ), 1 ] = 1:8 + 8 # Different
    Xs[ qs( 'Ta_SD', ms ), 2 ] = c( 5:8, 1:4 ) # Same
    X = expand_mat( Xs, ms$Cnd, dat )
    # Parameter/column indices
    ind = list(
      k1 = list(
        ip = 1:16 + start[1],
        ix = 1 + start[2],
        sm = F ),
      k0 = list(
        ip = 1:16 + start[1],
        ix = 2 + start[2],
        sm = F )
    )
    # Transformation type
    type = rep( 1, 16 )
    # Starting value generation
    start_rng = rbind( rep( log( .5 ), 16 ), 
                       rep( log( 4 ), 16 ) )
    # Parameter scaling
    pscl = log( 1.3 )
    # Labels
    nms = paste( 'kappa', 1:16, sep = '-' )
    # Priors
    priors = list()
  }
  
  end = numeric(2)
  end[1] = max( c( ind$k1$ip, ind$k0$ip ) )
  end[2] = max( c( ind$k1$ix, ind$k0$ix ) )  
  
  return( list( X = X, 
                Xs = Xs, 
                ind = ind, 
                nms = nms, 
                end = end,
                start_rng = start_rng,
                type = type, 
                pscl = pscl,
                priors = priors ) )
}

# Lookup - 02b
drm_x = function( ver, ms, dat, start ) {
  # Purpose:
  # A function to create different types of parameter 
  # mappings for the drift rate parameters of the diffusion 
  # race model.
  # Arguments:
  # ver   - The type of parameter mapping to create
  # ms    - A list with two data frames giving the condition 
  #         labels or their shorthand representation
  # dat   - The data frame with the observations and conditions
  # start - The shift for the parameter and column indices 
  #         respectively.
  # Returns:
  # A list with...
  # X         = A matrix of parameter indices for all observations
  # Xs        = A matrix of parameter indices for the unique 
  #             combination of conditions
  # ind       = A list with...
  #             ip) The subset of parameters for a particular 
  #                 parameter type
  #             ix) The subset of columns for X or Xs associated with 
  #                 a particular parameter type
  #             is) A subset of subsets, the parameters to be 
  #                 rescaled (Optional)
  #             sp) The indices for the rescaling coefficients
  #                 (Optional)
  # nms       = The labels for the parameters
  # end       = The new shift to be applied to subsequent parameter 
  #             and column indices
  # ub        = The upper boundary for the randomized starting points
  # pscl      = The scaling of the parameters for the 'optim' function
  # start_rng = A matrix whose first row gives the lower boundary 
  #             and second row gives the upper boundary for 
  #             use with random starting value generation
  # type      = The type of parameter transformation to apply to 
  #             the parameters
  # priors    = An optional of the prior values to impose on the 
  #             parameters
  
  if ( ver == 'Intercept' ) {
    # Design matrices
    Xs = matrix( 1, nrow( ms$Cnd ), 1 )
    X = expand_mat( Xs, ms$Cnd, dat )
    # Parameter/column indices
    ind = list(
      x1 = list(
        ip = 1 + start[1],
        ix = 1 + start[2],
        sm = F ),
      x0 = list(
        ip = 1 + start[1],
        ix = 1 + start[2],
        sm = F )
    )
    # Transformation type
    type = 1
    # Starting value generation
    start_rng = rbind( log( .5 ), 
                       log( 4 ) )
    # Parameter scaling
    pscl = c( .45 )
    # Labels
    nms = paste( 'xi', 1, sep = '-' )
    # Priors
    priors = list()
  }
  
  if ( ver == 'Binary' ) {
    # Design matrices
    Xs = matrix( 1, nrow( ms$Cnd ), 2 )
    Xs[ qs( 'Co_L', ms ), 1 ] = 2
    Xs[ qs( 'Co_S', ms ), 1 ] = 2
    Xs[ qs( 'Co_R', ms ), 2 ] = 2
    Xs[ qs( 'Co_D', ms ), 2 ] = 2
    X = expand_mat( Xs, ms$Cnd, dat )
    # Parameter/column indices
    ind = list(
      x1 = list(
        ip = 1:2 + start[1],
        ix = 1 + start[2],
        sm = F ),
      x0 = list(
        ip = 1:2 + start[1],
        ix = 2 + start[2],
        sm = F )
    )
    # Transformation type
    type = c(1,1)
    # Starting value generation
    start_rng = rbind( rep( log( .5 ), 2 ), 
                       rep( log( 4 ), 2 ) )
    # Parameter scaling
    pscl = c( .7, .2 )
    # Labels
    nms = paste( 'xi', 1:2, sep = '-' )
    # Priors
    priors = list()
  }
  
  if ( ver == 'Saturated' ) {
    # Design matrices
    Xs = matrix( 1, nrow( ms$Cnd ), 2 )
    Xs[ qs( 'Ta_FC', ms ), 2 ] = 1:8 # Left
    Xs[ qs( 'Ta_FC', ms ), 1 ] = c( 5:8, 1:4 ) + 8 # Right
    Xs[ qs( 'Ta_SD', ms ), 1 ] = 1:8 + 8 # Different
    Xs[ qs( 'Ta_SD', ms ), 2 ] = c( 5:8, 1:4 ) # Same
    X = expand_mat( Xs, ms$Cnd, dat )
    # Parameter/column indices
    ind = list(
      x1 = list(
        ip = 1:16 + start[1],
        ix = 1 + start[2],
        sm = F ),
      x0 = list(
        ip = 1:16 + start[1],
        ix = 2 + start[2],
        sm = F )
    )
    # Transformation type
    type = rep(1,16)
    # Starting value generation
    start_rng = rbind( rep( log( .5 ), 16 ), 
                       rep( log( 4 ), 16 ) )
    # Parameter scaling
    pscl = rep( rep( c( .7, .2 ), each = 4 ), 2 )
    # Labels
    nms = paste( 'xi', 1:16, sep = '-' )
    # Priors
    priors = list()
  }
  
  if ( ver == 'Orig (SD)' ) {
    # Design matrices
    Xs = diag( 8 )
    Xs = cbind( Xs, -1 )
    Xs = cbind( -diag( 8 ), 2, Xs )
    X = expand_mat( Xs, ms$Cnd, dat )
    # Parameter/column indices
    ind = list(
      x1 = list(
        ip = c( 5:8, 1:4, 10 ) + start[1],
        ix = 1:9 + start[2],
        sm = T ),
      x0 = list(
        ip = c( 5:8, 1:4, 9 ) + start[1],
        ix = 1:9 + 9 + start[2],
        sm = T )
    )
    # Transformation type
    type = c( rep(1,8), 0, 1 )
    # Starting value generation
    start_rng = rbind( c( rep( log( .5 ), 8 ), 0, log( .5 ) ), 
                       c( rep( log( 4 ), 8 ), 1, log( 4 ) ) )
    # Parameter scaling 
    pscl = c( rep( c( .7, .2 ), each = 4 ), .8, .2 )
    # Labels
    nms = paste( 'xi', 1:10, sep = '-' )
    # Priors
    priors = list()
  }
  
  if ( ver == 'Mirror' ) {
    # Design matrices
    Xs = diag( 8 )
    Xs = cbind( -diag( 8 ), 2, Xs )
    X = expand_mat( Xs, ms$Cnd, dat )
    # Parameter/column indices
    ind = list(
      x1 = list(
        ip = c( 5:8, 1:4, 9 ) + start[1],
        ix = 1:9 + start[2],
        sm = T ),
      x0 = list(
        ip = c( 5:8, 1:4 ) + start[1],
        ix = 1:8 + 9 + start[2],
        sm = T )
    )
    # Transformation type
    type = rep( 1, 9 )
    # Starting value generation
    start_rng = rbind( rep( log( .5 ), 9 ), 
                       rep( log( 4 ), 9 ) )
    # Parameter scaling 
    pscl = c( rep( c( .7, .2 ), each = 4 ), .2 )
    # Labels
    nms = paste( 'xi', 1:9, sep = '-' )
    # Priors
    priors = list()
  }
  
  if ( ver == 'Mirror (Revised)' ) {
    # Design matrices
    Xs = diag( 8 )
    Xs = cbind( Xs, -1 ) # Same
    Xs = cbind( -diag( 8 ), 2, Xs ) # Different + Same
    X = expand_mat( Xs, ms$Cnd, dat )
    # Parameter/column indices
    ind = list(
      x1 = list(
        ip = c( 5:8, 1:4, 10 ) + start[1],
        ix = 1:9 + start[2],
        is = 5:8,
        sp = 11 + start[1],
        sm = T ),
      x0 = list(
        ip = c( 5:8, 1:4, 9 ) + start[1],
        ix = 1:9 + 9 + start[2],
        sm = T )
    )
    # Transformation type
    type = c( rep(1,8), 0, 1, 0 )
    # Starting value generation
    start_rng = rbind( c( rep( log( .5 ), 8 ), 0, log( .5 ), 0 ), 
                       c( rep( log( 4 ), 8 ), 1, log( 4 ), 1 ) )
    # Parameter scaling 
    pscl = c( rep( c( .7, .2 ), each = 4 ), .8, .2, .7 )
    # Labels
    nms = paste( 'xi', 1:11, sep = '-' )
    # Priors
    priors = list()
  }
  
  if ( ver == 'Orig (FC)' ) {
    # Design matrices
    Xs = matrix( 1, nrow( ms$Cnd ), 2 )
    Xs[ qs( 'Ta_FC', ms ), 2 ] = 1:8 # Left
    Xs[ qs( 'Ta_FC', ms ), 1 ] = c( 5:8, 1:4 ) # Right
    X = expand_mat( Xs, ms$Cnd, dat )
    # Parameter/column indices
    ind = list(
      x1 = list(
        ip = 1:8 + start[1],
        ix = 1 + start[2],
        sm = F ),
      x0 = list(
        ip = 1:8 + start[1],
        ix = 2 + start[2],
        sm = F )
    )
    # Transformation type
    type = rep(1,8)
    # Starting value generation
    start_rng = rbind( rep( log( .5 ), 8 ), 
                       rep( log( 4 ), 8 ) )
    # Parameter scaling
    pscl = rep( rep( c( .7, .2 ), each = 4 ), 2 )
    # Labels
    nms = paste( 'xi', 1:8, sep = '-' )
    # Priors
    priors = list()
  }
  
  if ( ver == 'Regression (Orig)' ) {
    # Design matrices
    Xs = matrix( 1, nrow( ms$Cnd ), 4 )
    Xs[,2] = xi_reg[c(5:8,1:4)] # Right
    Xs[,4] = xi_reg[1:8] # Left
    X = expand_mat( Xs, ms$Cnd, dat )
    # Parameter/column indices
    ind = list(
      x1 = list(
        ip = 1:2 + start[1],
        ix = 1:2 + start[2],
        sm = T ),
      x0 = list(
        ip = 3:4 + start[1],
        ix = 3:4 + start[2],
        sm = T )
    )
    # Transformation type
    type = c( 1, 0, 1, 0 )
    # Starting value generation
    start_rng = rbind( c( log( .5 ), 0, log( .5 ), 0 ), 
                       c( log( 4 ), 1, log( 4 ), 1 ) )
    # Parameter scaling
    pscl = rep( 1, 4 )
    # Labels
    nms = paste( 'xi', 1:4, sep = '-' )
    # Priors
    priors = list()
  }
  
  end = numeric(2)
  end[1] = max( c( ind$x1$ip, ind$x0$ip ) )
  end[2] = max( c( ind$x1$ix, ind$x0$ix ) )
  
  return( list( X = X, 
                Xs = Xs, 
                ind = ind, 
                nms = nms, 
                end = end,
                start_rng = start_rng,
                type = type, 
                pscl = pscl,
                priors = priors ) )
}

# Lookup - 02c
drm_t = function( ver, ms, dat, start ) {
  # Purpose:
  # A function to create different types of parameter 
  # mappings for the residual latency parameters of the diffusion 
  # race model.
  # Arguments:
  # ver   - The type of parameter mapping to create
  # ms    - A list with two data frames giving the condition 
  #         labels or their shorthand representation
  # dat   - The data frame with the observations and conditions
  # start - The shift for the parameter and column indices 
  #         respectively.
  # Returns:
  # A list with...
  # X         = A matrix of parameter indices for all observations
  # Xs        = A matrix of parameter indices for the unique 
  #             combination of conditions
  # ind       = A list with...
  #             ip) The subset of parameters for a particular 
  #                 parameter type
  #             ix) The subset of columns for X or Xs associated with 
  #                 a particular parameter type
  #             is) A subset of subsets, the parameters to be 
  #                 rescaled (Optional)
  #             sp) The indices for the rescaling coefficients
  #                 (Optional)
  # nms       = The labels for the parameters
  # end       = The new shift to be applied to subsequent parameter 
  #             and column indices
  # ub        = The upper boundary for the randomized starting points
  # pscl      = The scaling of the parameters for the 'optim' function
  # start_rng = A matrix whose first row gives the lower boundary 
  #             and second row gives the upper boundary for 
  #             use with random starting value generation
  # type      = The type of parameter transformation to apply to 
  #             the parameters
  # priors    = An optional of the prior values to impose on the 
  #             parameters
  
  if ( ver == 'Intercept' ) {
    # Minimum response times
    mat = matrix( 1, nrow( ms$Cnd ), 2 )
    min_rt = find_min_rt( mat, ms, dat, F )
    # Design matrices
    Xs = matrix( 1, nrow( ms$Cnd ), 1 )
    X = expand_mat( Xs, ms$Cnd, dat )
    # Parameter/column indices
    ind = list(
      t1 = list(
        ip = 1 + start[1],
        ix = 1 + start[2],
        sm = F ),
      t0 = list(
        ip = 1 + start[1],
        ix = 1 + start[2],
        sm = F )
    )
    # Transformation type
    type = 2
    # Starting value generation
    start_rng = rbind( logit( .2 ), 
                       logit( .8 ) )
    # Parameter scaling
    pscl = logit( .55 )
    # Labels
    nms = c( 'tau' )
    # Priors
    priors = list()
  }
  
  if ( ver == 'Adaptive' ) {
    # Minimum response times
    mat = matrix( 1:8, nrow( ms$Cnd ), 2 )
    min_rt = find_min_rt( mat, ms, dat, T )
    # Design matrices
    Xs = matrix( 1, nrow( ms$Cnd ), 1 )
    X = expand_mat( Xs, ms$Cnd, dat )
    # Parameter/column indices
    ind = list(
      t1 = list(
        ip = 1 + start[1],
        ix = 1 + start[2],
        sm = F ),
      t0 = list(
        ip = 1 + start[1],
        ix = 1 + start[2],
        sm = F )
    )
    # Transformation type
    type = 2
    # Starting value generation
    start_rng = rbind( logit( .2 ), 
                       logit( .8 ) )
    # Parameter scaling
    pscl = logit( .55 )
    # Labels
    nms = c( 'tau' )
    # Priors
    priors = list()
  }
  
  end = numeric(2)
  end[1] = max( c( ind$t1$ip, ind$t0$ip ) )
  end[2] = max( c( ind$t1$ix, ind$t0$ix ) )
  
  return( list( X = X, 
                Xs = Xs, 
                ind = ind, 
                nms = nms, 
                end = end,
                start_rng = start_rng, 
                type = type, 
                pscl = pscl, 
                min_rt = min_rt,
                priors = priors ) )
}

# Lookup - 02d
drm_s = function( ver, ms, dat, start ) {
  # Purpose:
  # A function to create different types of parameter 
  # mappings for the coefficient of drift parameters of the 
  # diffusion race model.
  # Arguments:
  # ver   - The type of parameter mapping to create
  # ms    - A list with two data frames giving the condition 
  #         labels or their shorthand representation
  # dat   - The data frame with the observations and conditions
  # start - The shift for the parameter and column indices 
  #         respectively.
  # Returns:
  # A list with...
  # X    = A matrix of parameter indices for all observations
  # Xs   = A matrix of parameter indices for the unique 
  #        combination of conditions
  # ind  = A list with...
  #        ip) The subset of parameters for a particular parameter 
  #            type
  #        ix) The subset of columns for X or Xs associated with 
  #            a particular parameter type
  #        is) A subset of subsets, the parameters to be 
  #            rescaled
  #        sp) The indices for the rescaling coefficients
  # lb   = The lower boundary for the randomized starting points
  # ub   = The upper boundary for the randomized starting points
  # pscl = The scaling of the parameters for the 'optim' function
  
  if ( ver == 'Fixed' ) {
    # Design matrices
    Xs = matrix( 1, nrow( ms$Cnd ), 1 )
    X = expand_mat( Xs, ms$Cnd, dat )
    # Parameter/column indices
    ind = list(
      s0 = list(
        ip = 1 + start[1],
        ix = 1 + start[2],
        sm = F )
    )
    # Transformation type
    type = NULL
    # Starting value generation
    start_rng = NULL
    # Parameter scaling
    pscl = NULL
    # Labels
    nms = NULL
    # Priors
    priors = list()
  }
  
  if ( ver == 'Intercept' ) {
    # Design matrices
    Xs = matrix( 1, nrow( ms$Cnd ), 1 )
    X = expand_mat( Xs, ms$Cnd, dat )
    # Parameter/column indices
    ind = list(
      s0 = list(
        ip = 1 + start[1],
        ix = 1 + start[2],
        sm = F )
    )
    # Transformation type
    type = 1
    # Starting value generation
    start_rng = rbind( log( .5 ), 
                       log( 1.5 ) )
    # Parameter scaling
    pscl = log( 1.1 )
    # Labels
    nms = c( 'sigma' )
    # Priors
    priors = list()
  }
  
  end = numeric(2)
  end[1] = max( ind$s0$ip )
  end[2] = max( ind$s0$ix )
  
  return( list( X = X, 
                Xs = Xs, 
                ind = ind, 
                nms = nms, 
                end = end,
                start_rng = start_rng,
                type = type, 
                pscl = pscl,
                priors = priors ) )
}

# Lookup - 02e
drm_pm = function( X, prm, ind ) {
  # Purpose: 
  # A function to create a data frame with the parameter 
  # values associated with each observation.
  # Arguments: 
  # X   - A design matrix with the parameter indices per 
  #       each observation
  # prm - A vector of parameters
  # ind - A list of indices
  # Returns: 
  # A data frame with the paraemter values associated with 
  # each observation per parameter type.
  
  out = matrix( NA, nrow( X ), 7 )
  colnames( out ) = c( 'k1', 'k0',
                       'x1', 'x0',
                       't1', 't0',
                       's0' )
  out = as.data.frame( out )
  
  out$k1 = p_c( X, prm, ind$k1 )
  out$k0 = p_c( X, prm, ind$k0 )
  out$x1 = p_c( X, prm, ind$x1 )
  out$x0 = p_c( X, prm, ind$x0 )
  out$t1 = p_c( X, prm, ind$t1 )
  out$t0 = p_c( X, prm, ind$t0 )
  out$s0 = p_c( X, prm, ind$s0 )
  out$s0[ is.na( out$s0 ) ] = 1
  
  return( out )
}

# Lookup - 02f
drm_pm_quick_convert = function( prm, dtbf ) {
  # Purpose: 
  # A function that creates the parameter matrix associated 
  # with the unique combination of condition levels.
  # Arguments: 
  # prm  - A vector of parameter values
  # dtbf - A list with the data that was fitted
  # Returns: 
  # A parameter matrix for the diffusion race model,
  # giving the thresholds (k1 and k0), drift rates 
  # (x1 and x0), residual latencies (t1 and t0) and 
  # the coefficient of drift (s0) for the unique 
  # combination of the condition levels.
  
  prm = tran_par( prm, dtbf$type )
  out = drm_pm( dtbf$Xs, prm, dtbf$ind )
  out[,'t1'] = out[,'t1'] * dtbf$min_rt$t1$t
  out[,'t0'] = out[,'t0'] * dtbf$min_rt$t0$t
  
  rownames( out ) = apply( dtbf$ms$SH, 1, 
                           function(x) paste( unlist(x), 
                                              collapse = '-' ) )
  as.data.frame( out )
  
  return( out )
}

# Lookup - 02g
drm_dtbf_create = function( ver, subject, task = 'Both', 
                            d = d ) {
  # Purpose: 
  # A function to create the list of data and structures 
  # used to fit the diffusion race model to the observed 
  # data.
  # Arguments: 
  # ver     - A vector of 4 strings, indicating the 
  #           structures for the threshold, drift, 
  #           residual latency, and coefficient of 
  #           drift parameters.
  # subject - The current subject being fitted
  # task    - The task to be fitted
  # d       - The data frame with all of the observations
  # Returns: 
  # A list with...
  # rt = The vector of response times to fit
  # ch = The vector of choices to fit
  # X = The design matrix
  # Xs = The simplified design matrix (for plotting)
  # ind = The list of indices
  # ms = The list of data-frames with the condition combinations
  # start_rng = A matrix with the lower and upper boundaries for 
  #             parameter generation
  # pscl = An optional vector for parameter scaling while fitting
  # min_rt = A list with the minimum response times
  # priors = A list with the parameters for priors
  # pnames = The set of parameter names
  # PD, PT, Co, Ta = The covariates from the experiment

  if ( task == 'Forced-choice' ) {
    sbst = d$S == subject & d$TaL == task
  }
  if ( task == 'Same-different' ) {
    sbst = d$S == subject & d$TaL == task
  }
  if ( task == 'Both' ) {
    sbst = d$S == subject
  }
  
  dat = d[ sbst, ]
  
  ms = extract_cnd( dat, task = task )
  
  # Initialize variables
  X = c(); Xs = c();
  ind = list( k1 = NULL, k0 = NULL,
             x1 = NULL, x0 = NULL,
             t1 = NULL, t0 = NULL,
             s0 = NULL )
  pn = c(); pscl = c(); start_rng = c(); type = c()
  priors = list( mu = NULL, sigma = NULL )
  
  pt = list( c( 'k1','k0' ), c('x1','x0'), c('t1','t0'), 's0' )
  for ( i in 1:4 ) {
    
    if ( i == 1 ) tmp = drm_k( ver[i], ms, dat, c(0,0) )
    if ( i == 2 ) tmp = drm_x( ver[i], ms, dat, tmp$end )
    if ( i == 3 ) tmp = drm_t( ver[i], ms, dat, tmp$end )
    if ( i == 4 ) tmp = drm_s( ver[i], ms, dat, tmp$end )
    
    if ( i == 3 ) min_rt = tmp$min_rt
    
    X = cbind( X, tmp$X ); Xs = cbind( Xs, tmp$Xs )
    for ( j in 1:length( pt[[i]] ) ) {
      ind[[ pt[[i]][j] ]] = tmp$ind[[ pt[[i]][j] ]]
    }
    pn = c( pn, tmp$nms )
    start_rng = cbind( start_rng, tmp$start_rng )
    type = c( type, tmp$type )
    pscl = c( pscl, tmp$pscl )
    priors$mu = c( priors$mu, tmp$priors$mu )
    priors$sigma = c( priors$sigma, tmp$priors$sigma )
  }
  
  priors$mu = as.numeric( unlist( priors$mu ) )
  priors$sigma = as.numeric( unlist( priors$sigma ) )
  
  dtbf = list( rt = dat$RT, ch = dat$Ch,
               X = X,
               Xs = Xs,
               ind = ind,
               ms = ms, 
               start_rng = start_rng, 
               type = type, 
               pscl = pscl,
               min_rt = min_rt, 
               priors = priors, 
               pnames = pn,
               PD = dat$PDL,
               PT = dat$PTL,
               Co = dat$CoL,
               Ta = dat$TaL
  )
  
  return( dtbf )
}

###
### Two-boundary wiener process model functions
###
# Lookup - 03

# Lookup - 03a
wp_a = function( ver, ms, dat, start ) {
  # Purpose:
  # A function to create different types of parameter 
  # mappings for the boundary separation parameters of the 
  # two-boundary wiener process.
  # Arguments:
  # ver   - The type of parameter mapping to create
  # ms    - A list with two data frames giving the condition 
  #         labels or their shorthand representation
  # dat   - The data frame with the observations and conditions
  # start - The shift for the parameter and column indices 
  #         respectively.
  # Returns:
  # A list with...
  # X         = A matrix of parameter indices for all observations
  # Xs        = A matrix of parameter indices for the unique 
  #             combination of conditions
  # ind       = A list with...
  #             ip) The subset of parameters for a particular 
  #                 parameter type
  #             ix) The subset of columns for X or Xs associated with 
  #                 a particular parameter type
  #             is) A subset of subsets, the parameters to be 
  #                 rescaled (Optional)
  #             sp) The indices for the rescaling coefficients
  #                 (Optional)
  # nms       = The labels for the parameters
  # end       = The new shift to be applied to subsequent parameter 
  #             and column indices
  # ub        = The upper boundary for the randomized starting points
  # pscl      = The scaling of the parameters for the 'optim' function
  # start_rng = A matrix whose first row gives the lower boundary 
  #             and second row gives the upper boundary for 
  #             use with random starting value generation
  # type      = The type of parameter transformation to apply to 
  #             the parameters
  # priors    = An optional of the prior values to impose on the 
  #             parameters
  
  if ( ver == 'Intercept' ) {
    # Design matrices
    Xs = matrix( 1, nrow( ms$Cnd ), 1 )
    X = expand_mat( Xs, ms$Cnd, dat )
    # Parameter/column indices
    ind = list(
      a = list(
        ip = 1 + start[1],
        ix = 1 + start[2],
        sm = F )
    )
    # Transformation type
    type = 1
    # Starting value generation
    start_rng = rbind( log( .5 ), 
                       log( 2 ) )
    # Parameter scaling
    pscl = log( 1.3 )
    # Labels
    nms = c( 'alpha' )
    # Priors
    priors = list()
  }
  
  if ( ver == 'Duration' ) {
    # Design matrices
    Xs = cbind( rep( 1:2, 4 ) )
    X = expand_mat( Xs, ms$Cnd, dat )
    # Parameter/column indices
    ind = list(
      a = list(
        ip = 1:2 + start[1],
        ix = 1 + start[2],
        sm = F )
    )
    # Transformation type
    type = rep( 1, 2 )
    # Starting value generation
    start_rng = rbind( rep( log( .5 ), 2 ), 
                       rep( log( 2 ), 2 ) )
    # Parameter scaling
    pscl = rep( log( 1.3 ), 2 )
    # Labels
    nms = paste( 'alpha', 1:2, sep = '-' )
    # Priors
    priors = list()
  }
  
  if ( ver == 'Saturated' ) {
    # Design matrices
    Xs = matrix( 1, nrow( ms$Cnd ), 1 )
    Xs[ qs( 'Ta_FC', ms ), 1 ] = c( 5:8, 1:4 ) # Right
    Xs[ qs( 'Ta_SD', ms ), 1 ] = 1:8 # Different
    X = expand_mat( Xs, ms$Cnd, dat )
    # Parameter/column indices
    ind = list(
      a = list(
        ip = 1:8 + start[1],
        ix = 1 + start[2],
        sm = F )
    )
    # Transformation type
    type = rep( 1, 8 )
    # Starting value generation
    start_rng = rbind( rep( log( .5 ), 8 ), 
                       rep( log( 2 ), 8 ) )
    # Parameter scaling
    pscl = rep( log( 1.3 ), 8 )
    # Labels
    nms = paste( 'alpha', 1:8, sep = '-' )
    # Priors
    priors = list()
  }
  
  end = numeric(2)
  end[1] = max( ind$a$ip )
  end[2] = max( ind$a$ix )
  
  return( list( X = X, 
                Xs = Xs, 
                ind = ind, 
                nms = nms, 
                end = end,
                start_rng = start_rng,
                type = type, 
                pscl = pscl,
                priors = priors ) )
}

# Lookup - 03b
wp_z = function( ver, ms, dat, start ) {
  # Purpose:
  # A function to create different types of parameter 
  # mappings for the starting point parameters of the 
  # two-boundary wiener process.
  # Arguments:
  # ver   - The type of parameter mapping to create
  # ms    - A list with two data frames giving the condition 
  #         labels or their shorthand representation
  # dat   - The data frame with the observations and conditions
  # start - The shift for the parameter and column indices 
  #         respectively.
  # Returns:
  # A list with...
  # X         = A matrix of parameter indices for all observations
  # Xs        = A matrix of parameter indices for the unique 
  #             combination of conditions
  # ind       = A list with...
  #             ip) The subset of parameters for a particular 
  #                 parameter type
  #             ix) The subset of columns for X or Xs associated with 
  #                 a particular parameter type
  #             is) A subset of subsets, the parameters to be 
  #                 rescaled (Optional)
  #             sp) The indices for the rescaling coefficients
  #                 (Optional)
  # nms       = The labels for the parameters
  # end       = The new shift to be applied to subsequent parameter 
  #             and column indices
  # ub        = The upper boundary for the randomized starting points
  # pscl      = The scaling of the parameters for the 'optim' function
  # start_rng = A matrix whose first row gives the lower boundary 
  #             and second row gives the upper boundary for 
  #             use with random starting value generation
  # type      = The type of parameter transformation to apply to 
  #             the parameters
  # priors    = An optional of the prior values to impose on the 
  #             parameters
  
  if ( ver == 'Fixed' ) {
    # Design matrices
    Xs = matrix( 1, nrow( ms$Cnd ), 1 )
    X = expand_mat( Xs, ms$Cnd, dat )
    # Parameter/column indices
    ind = list(
      z = list(
        ip = 1 + start[1],
        ix = 1 + start[2],
        sm = F )
    )
    # Transformation type
    type = NULL
    # Starting value generation
    start_rng = NULL
    # Parameter scaling
    pscl = NULL
    # Labels
    nms = NULL
    # Priors
    priors = list()
  }
  
  if ( ver == 'Intercept' ) {
    # Design matrices
    Xs = matrix( 1, nrow( ms$Cnd ), 1 )
    X = expand_mat( Xs, ms$Cnd, dat )
    # Parameter/column indices
    ind = list(
      z = list(
        ip = 1 + start[1],
        ix = 1 + start[2],
        sm = F )
    )
    # Transformation type
    type = 2
    # Starting value generation
    start_rng = rbind( logit( .2 ), 
                       logit( .8 ) )
    # Parameter scaling
    pscl = logit( .55 )
    # Labels
    nms = c( 'theta' )
    # Priors
    priors = list()
  }
  
  if ( ver == 'Duration' ) {
    # Design matrices
    Xs = cbind( rep( 1:2, 4 ) )
    X = expand_mat( Xs, ms$Cnd, dat )
    # Parameter/column indices
    ind = list(
      z = list(
        ip = 1:2 + start[1],
        ix = 1 + start[2],
        sm = F )
    )
    # Transformation type
    type = rep( 2, 2 )
    # Starting value generation
    start_rng = rbind( rep( logit( .2 ), 2 ), 
                       rep( logit( .8 ), 2 ) )
    # Parameter scaling
    pscl = rep( logit( .55 ), 2 )
    # Labels
    nms = paste( 'theta', 1:2, sep = '-' )
    # Priors
    priors = list()
  }
  
  if ( ver == 'Saturated' ) {
    # Design matrices
    Xs = matrix( 1, nrow( ms$Cnd ), 1 )
    Xs[ qs( 'Ta_FC', ms ), 1 ] = c( 5:8, 1:4 ) # Right
    Xs[ qs( 'Ta_SD', ms ), 1 ] = 1:8 # Different
    X = expand_mat( Xs, ms$Cnd, dat )
    # Parameter/column indices
    ind = list(
      z = list(
        ip = 1:8 + start[1],
        ix = 1 + start[2],
        sm = F )
    )
    # Transformation type
    type = rep( 2, 8 )
    # Starting value generation
    start_rng = rbind( rep( logit( .2 ), 8 ), 
                       rep( logit( .8 ), 8 ) )
    # Parameter scaling
    pscl = rep( logit( .55 ), 8 )
    # Labels
    nms = paste( 'theta', 1:8, sep = '-' )
    # Priors
    priors = list()
  }
  
  end = numeric(2)
  end[1] = max( ind$z$ip )
  end[2] = max( ind$z$ix )
  
  return( list( X = X, 
                Xs = Xs, 
                ind = ind, 
                nms = nms, 
                end = end,
                start_rng = start_rng,
                type = type, 
                pscl = pscl,
                priors = priors ) )
}

# Lookup - 03c
wp_v = function( ver, ms, dat, start ) {
  # Purpose:
  # A function to create different types of parameter 
  # mappings for the drift rate parameters of the 
  # two-boundary wiener process.
  # Arguments:
  # ver   - The type of parameter mapping to create
  # ms    - A list with two data frames giving the condition 
  #         labels or their shorthand representation
  # dat   - The data frame with the observations and conditions
  # start - The shift for the parameter and column indices 
  #         respectively.
  # Returns:
  # A list with...
  # X         = A matrix of parameter indices for all observations
  # Xs        = A matrix of parameter indices for the unique 
  #             combination of conditions
  # ind       = A list with...
  #             ip) The subset of parameters for a particular 
  #                 parameter type
  #             ix) The subset of columns for X or Xs associated with 
  #                 a particular parameter type
  #             is) A subset of subsets, the parameters to be 
  #                 rescaled (Optional)
  #             sp) The indices for the rescaling coefficients
  #                 (Optional)
  # nms       = The labels for the parameters
  # end       = The new shift to be applied to subsequent parameter 
  #             and column indices
  # ub        = The upper boundary for the randomized starting points
  # pscl      = The scaling of the parameters for the 'optim' function
  # start_rng = A matrix whose first row gives the lower boundary 
  #             and second row gives the upper boundary for 
  #             use with random starting value generation
  # type      = The type of parameter transformation to apply to 
  #             the parameters
  # priors    = An optional of the prior values to impose on the 
  #             parameters
  
  if ( ver == 'Intercept' ) {
    # Design matrices
    Xs = matrix( 1, nrow( ms$Cnd ), 1 )
    X = expand_mat( Xs, ms$Cnd, dat )
    # Parameter/column indices
    ind = list(
      v = list(
        ip = 1 + start[1],
        ix = 1 + start[2],
        sm = F )
    )
    # Transformation type
    type = 0
    # Starting value generation
    start_rng = rbind( -2, 
                       2 )
    # Parameter scaling
    pscl = 1
    # Labels
    nms = c( 'xi' )
    # Priors
    priors = list()
  }
  
  if ( ver == 'Saturated' ) {
    # Design matrices
    Xs = matrix( 1, nrow( ms$Cnd ), 1 )
    Xs[ qs( 'Ta_FC', ms ), 1 ] = c( 5:8, 1:4 ) # Right
    Xs[ qs( 'Ta_SD', ms ), 1 ] = 1:8 # Different
    X = expand_mat( Xs, ms$Cnd, dat )
    # Parameter/column indices
    ind = list(
      v = list(
        ip = 1:8 + start[1],
        ix = 1 + start[2],
        sm = F )
    )
    # Transformation type
    type = rep( 0, 8 )
    # Starting value generation
    start_rng = rbind( rep( -2, 8 ), 
                       rep( 2, 8 ) )
    # Parameter scaling
    pscl = rep( 1, 8 )
    # Labels
    nms = paste( 'xi', 1:8, sep = '-' )
    # Priors
    priors = list()
  }
  
  if ( ver == 'Regression' ) {
    # Design matrices
    Xs = matrix( 1, nrow( ms$Cnd ), 2 )
    Xs[,2] = xi_reg[c(5:8,1:4)] - xi_reg[1:8]
    X = expand_mat( Xs, ms$Cnd, dat )
    # Parameter/column indices
    ind = list(
      v = list(
        ip = 1:2 + start[1],
        ix = 1:2 + start[2],
        sm = T )
    )
    # Transformation type
    type = rep( 0, 2 )
    # Starting value generation
    start_rng = rbind( rep( -2, 2 ), 
                       rep( 2, 2 ) )
    # Parameter scaling
    pscl = rep( 1, 2 )
    # Labels
    nms = paste( 'xi', 1:2, sep = '-' )
    # Priors
    priors = list()
  }
  
  end = numeric(2)
  end[1] = max( ind$v$ip )
  end[2] = max( ind$v$ix )
  
  return( list( X = X, 
                Xs = Xs, 
                ind = ind, 
                nms = nms, 
                end = end,
                start_rng = start_rng,
                type = type, 
                pscl = pscl,
                priors = priors ) )
}

# Lookup - 03d
wp_t = function( ver, ms, dat, start ) {
  # Purpose:
  # A function to create different types of parameter 
  # mappings for the residual latency parameters of the 
  # two-boundary wiener process.
  # Arguments:
  # ver   - The type of parameter mapping to create
  # ms    - A list with two data frames giving the condition 
  #         labels or their shorthand representation
  # dat   - The data frame with the observations and conditions
  # start - The shift for the parameter and column indices 
  #         respectively.
  # Returns:
  # A list with...
  # X         = A matrix of parameter indices for all observations
  # Xs        = A matrix of parameter indices for the unique 
  #             combination of conditions
  # ind       = A list with...
  #             ip) The subset of parameters for a particular 
  #                 parameter type
  #             ix) The subset of columns for X or Xs associated with 
  #                 a particular parameter type
  #             is) A subset of subsets, the parameters to be 
  #                 rescaled (Optional)
  #             sp) The indices for the rescaling coefficients
  #                 (Optional)
  # nms       = The labels for the parameters
  # end       = The new shift to be applied to subsequent parameter 
  #             and column indices
  # ub        = The upper boundary for the randomized starting points
  # pscl      = The scaling of the parameters for the 'optim' function
  # start_rng = A matrix whose first row gives the lower boundary 
  #             and second row gives the upper boundary for 
  #             use with random starting value generation
  # type      = The type of parameter transformation to apply to 
  #             the parameters
  # priors    = An optional of the prior values to impose on the 
  #             parameters
  
  if ( ver == 'Intercept' ) {
    # Minimum response times
    mat = matrix( 1, nrow( ms$Cnd ), 2 )
    min_rt = find_min_rt( mat, ms, dat, F )
    # Design matrices
    Xs = matrix( 1, nrow( ms$Cnd ), 1 )
    X = expand_mat( Xs, ms$Cnd, dat )
    # Parameter/column indices
    ind = list(
      t0 = list(
        ip = 1 + start[1],
        ix = 1 + start[2],
        sm = F )
    )
    # Transformation type
    type = 2
    # Starting value generation
    start_rng = rbind( logit( .2 ), 
                       logit( .8 ) )
    # Parameter scaling
    pscl = logit( .55 )
    # Labels
    nms = c( 'tau' )
    # Priors
    priors = list()
  }
  
  if ( ver == 'Adaptive' ) {
    # Minimum response times
    mat = matrix( 1:8, nrow( ms$Cnd ), 2 )
    min_rt = find_min_rt( mat, ms, dat, T )
    # Design matrices
    Xs = matrix( 1, nrow( ms$Cnd ), 1 )
    X = expand_mat( Xs, ms$Cnd, dat )
    # Design matrices
    Xs = matrix( 1, nrow( ms$Cnd ), 1 )
    X = expand_mat( Xs, ms$Cnd, dat )
    # Parameter/column indices
    ind = list(
      t0 = list(
        ip = 1 + start[1],
        ix = 1 + start[2],
        sm = F )
    )
    # Transformation type
    type = 2
    # Starting value generation
    start_rng = rbind( logit( .2 ), 
                       logit( .8 ) )
    # Parameter scaling
    pscl = logit( .55 )
    # Labels
    nms = c( 'tau' )
    # Priors
    priors = list()
  }
  
  end = numeric(2)
  end[1] = max( ind$t0$ip )
  end[2] = max( ind$t0$ix )
  
  return( list( X = X, 
                Xs = Xs, 
                ind = ind, 
                nms = nms, 
                end = end,
                start_rng = start_rng, 
                type = type, 
                pscl = pscl, 
                min_rt = min_rt,
                priors = priors ) )
}

# Lookup - 03e
wp_pm = function( X, prm, ind ) {
  # Purpose: 
  # A function to create a data frame with the 
  # parameter values associated with each observation
  # Arguments: 
  # X   - A design matrix with the parameter indices per 
  #       each observation
  # prm - A vector of parameters
  # ind - A list of indices
  # Returns: 
  # A data frame with the paraemter values associated with 
  # each observation per parameter type.
  
  out = matrix( NA, nrow( X ), 4 )
  colnames( out ) = c( 'a', 'z',
                       'v', 't0' )
  out = as.data.frame( out )
  
  out$a = p_c( X, prm, ind$a )
  out$z = p_c( X, prm, ind$z )
  out$z[ is.na( out$z ) ] = .5
  out$v = p_c( X, prm, ind$v )
  out$t0 = p_c( X, prm, ind$t0 )
  
  return( out )
}

# Lookup - 03f
wp_pm_quick_convert = function( prm, dtbf ) {
  # Purpose: 
  # A function that creates the parameter matrix associated 
  # with the unique combination of condition levels.
  # Arguments: 
  # prm  - A vector of parameter values
  # dtbf - A list with the data that was fitted
  # Returns: 
  # A parameter matrix for the two-boundary wiener process,
  # giving the boundary separation (a), starting point (z),
  # drift rate (v), and residual latency (t0) for the unique 
  # combination of the condition levels.
  
  prm = tran_par( prm, dtbf$type )
  out = wp_pm( dtbf$Xs, prm, dtbf$ind )
  out[,'t0'] = out[,'t0'] * dtbf$min_rt$t0$t
  
  rownames( out ) = apply( dtbf$ms$SH, 1, 
                           function(x) paste( unlist(x), 
                                              collapse = '-' ) )
  as.data.frame( out )
  
  return( out )
}

# Lookup - 03g
wp_dtbf_create = function( ver, subject, task = 'Both', 
                            d = d ) {
  # Purpose: 
  # A function to create the list of data and structures 
  # used to fit the two-boundary wiener process to the 
  # observed data.
  # Arguments: 
  # ver     - A vector of 4 strings, indicating the 
  #           structures for the boundary separation, 
  #           drift rate, residual latency, and the 
  #           starting point parameters.
  # subject - The current subject being fitted
  # task    - The task to be fitted
  # d       - The data frame with all of the observations
  # Returns: 
  # A list with...
  # rt = The vector of response times to fit
  # ch = The vector of choices to fit
  # X = The design matrix
  # Xs = The simplified design matrix (for plotting)
  # ind = The list of indices
  # ms = The list of data-frames with the condition combinations
  # start_rng = A matrix with the lower and upper boundaries for 
  #             parameter generation
  # pscl = An optional vector for parameter scaling while fitting
  # min_rt = A list with the minimum response times
  # priors = A list with the parameters for priors
  # pnames = The set of parameter names
  # PD, PT, Co, Ta = The covariates from the experiment
  
  if ( task == 'Forced-choice' ) {
    sbst = d$S == subject & d$TaL == task
  }
  if ( task == 'Same-different' ) {
    sbst = d$S == subject & d$TaL == task
  }
  if ( task == 'Both' ) {
    sbst = d$S == subject
  }
  
  dat = d[ sbst, ]
  
  ms = extract_cnd( dat, task = task )
  
  # Initialize variables
  X = c(); Xs = c();
  ind = list( a = NULL, 
              z = NULL, 
              v = NULL, 
              t0 = NULL )
  pn = c(); pscl = c(); start_rng = c(); type = c()
  priors = list( mu = NULL, sigma = NULL )
  
  pt = list( c( 'a' ), c('v'), c('t0'), 'z' )
  for ( i in 1:4 ) {
    
    if ( i == 1 ) tmp = wp_a( ver[i], ms, dat, c(0,0) )
    if ( i == 2 ) tmp = wp_v( ver[i], ms, dat, tmp$end )
    if ( i == 3 ) tmp = wp_t( ver[i], ms, dat, tmp$end )
    if ( i == 4 ) tmp = wp_z( ver[i], ms, dat, tmp$end )
    
    if ( i == 3 ) min_rt = tmp$min_rt
    
    X = cbind( X, tmp$X ); Xs = cbind( Xs, tmp$Xs )
    for ( j in 1:length( pt[[i]] ) ) {
      ind[[ pt[[i]][j] ]] = tmp$ind[[ pt[[i]][j] ]]
    }
    pn = c( pn, tmp$nms )
    start_rng = cbind( start_rng, tmp$start_rng )
    type = c( type, tmp$type )
    pscl = c( pscl, tmp$pscl )
    priors$mu = c( priors$mu, tmp$priors$mu )
    priors$sigma = c( priors$sigma, tmp$priors$sigma )
  }
  
  priors$mu = as.numeric( unlist( priors$mu ) )
  priors$sigma = as.numeric( unlist( priors$sigma ) )
  
  dtbf = list( rt = dat$RT, ch = dat$Ch,
               X = X,
               Xs = Xs,
               ind = ind,
               ms = ms, 
               start_rng = start_rng, 
               type = type, 
               pscl = pscl,
               min_rt = min_rt, 
               priors = priors, 
               pnames = pn,
               PD = dat$PDL,
               PT = dat$PTL,
               Co = dat$CoL,
               Ta = dat$TaL
  )
  
  return( dtbf )
}
