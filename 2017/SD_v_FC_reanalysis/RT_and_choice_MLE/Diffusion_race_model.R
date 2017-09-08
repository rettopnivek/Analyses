#---------------------------#
# Diffusion race model fits # 
# Kevin Potter              #
# Updated 07/15/2017        #
#---------------------------#

# Extract arguments
args = commandArgs(trailingOnly=TRUE)

# For direct sourcings
if ( length( args ) == 0 ) {
  args = 7
}

# test if there is at least one argument: if not, return an error
if ( length( args ) == 0 ) {
  stop("At least one argument must be supplied (input file).n", 
       call.=FALSE)
} else if ( length( args ) == 1 ) {
  model_type = as.numeric( args[1] )
}

# Move up a directory
setwd( '..' )

# Initialize script
source('F3_starting_script.R')

# Return to MLE directory
setwd( 'RT_and_choice_MLE' )
cur_dir = getwd()

# Indicate which code to run
runCode = c( F, T )

# Index
# Lookup - 01:   Functions for maximum likelihood estimation
# Lookup - 01a:  mle_fn
# Lookup - 01b:  st_f
# Lookup - 01c:  plot_jcdf_subject
# Lookup - 01d:  model_select
# Lookup - 02:   Example with a single subject
# Lookup - 03:   Estimation over subjects

###
### Functions for maximum likelihood estimation
###
# Lookup - 01

# Lookup - 01a
mle_fn = function( prm, dat, sum = T, priors = NULL, predict = F ) {
  # Arguments:
  # prm     - A vector of parameters
  # dat     - A list that includes...
  #           rt = Response times to be fitted
  #           ch = Choice data to be fitted
  #           X = A design matrix with the indices/values to use  
  #               to create the parameter values by observation
  #           ind = A list of indices for each parameter type
  #           min_rt = A list with the minumum response times 
  #                    to use when constraining the residual 
  #                    latency
  # sum     - Logical; if true, returns the sum of the log-likelihoods
  # priors  - Not used.
  # predict - Logical; if true, returns a vector of the predicted 
  #           joint cumulative probabilities.
  # Returns:
  # The sum or vector of log-likelihoods given the inputted 
  # parameters.
  
  # Initial setup
  prm = tran_par( prm, dat$type )
  ps = drm_pm( dat$X, prm, dat$ind )
  ps[,'t1'] = ps[,'t1'] * dat$min_rt$t1$T
  ps[,'t0'] = ps[,'t0'] * dat$min_rt$t0$T
  
  # If specified, compute the predicted joint 
  # cumulative probabilities
  if ( predict ) {
    prd = pwaldrace( dat$rt, dat$ch, 
                     k1 = ps[,'k1'],
                     x1 = ps[,'x1'], 
                     t1 = ps[,'t1'],
                     k0 = ps[,'k0'], 
                     x0 = ps[,'x0'],
                     t0 = ps[,'t0'],
                     s0 = ps[,'s0'] )
    return( prd )
  }
  
  # Calculate the log-likelihoods 
  ll = dwaldrace( dat$rt, dat$ch, 
                  k1 = ps[,'k1'],
                  x1 = ps[,'x1'], 
                  t1 = ps[,'t1'],
                  k0 = ps[,'k0'], 
                  x0 = ps[,'x0'],
                  t0 = ps[,'t0'],
                  s0 = ps[,'s0'],
                  ln = T )
  if ( !sum ) return( ll )
  
  # Sum the log-likelihoods 
  sll = sum( ll )
  
  # Incorporate priors (Optional) 
  if ( !is.null( priors ) ) {
    # sll = sll + 
  }
  
  # Check for NA values 
  if ( is.na( sll ) ) sll = -Inf 
  
  return( sll ) 
}

# Lookup - 01b
st_f = function() {
  # Purpose: 
  # A function to generate starting values 
  # for the diffusion race model.
  # Details:
  # Requires the variable 'dtbf' to already exist 
  # in the workspace, a list that should include 
  # a matrix with the lower and upper boundaries 
  # to use when generating parameters.
  # Returns: 
  # A vector of raw parameter values.
  
  out = runif( ncol( dtbf$start_rng ),
               dtbf$start_rng[1,],
               dtbf$start_rng[2,] )
  
  return( out )
}

# Lookup - 01c
plot_jcdf_subject = function( res, dtbf, new = T ) {
  # Purpose: 
  # A function to generate a quick plot of the model 
  # predictions against the observed data for a single 
  # subject.
  # Arguments: 
  # res  - A model fit object
  # dtbf - The data frame with the observations that were fitted
  # new  - A logical value; if true, a new plotting window 
  #        is generated
  
  prm = coef(res)[which.max( res$value ),]
  pm = drm_pm_quick_convert( prm, dtbf )
  
  # Loop over tasks
  for ( nw in 1:length( unique( dtbf$ms$Cnd$Ta ) ) ) {
    
    # Create new plotting window
    if ( new ) x11( width = 12 )
    
    layout( matrix( 1:8, 2, 4, byrow = T ) )
    
    # Loop over conditions
    for ( cnd in 1:8 ) {
      
      cnd = cnd + 8 * ( nw - 1 )
      
      # Extract  CDF curve for model predictions
      obj = quickdist( 'wr', type = 'CDF', prm = unlist( pm[cnd,] ) )
      plot( obj, ylim = c(0,1) ) # Generate an empty plot
      # Add curves
      lines( obj ); lines( obj, ch = 0, lty = 2 )
      
      # Add condition label
      title( 
        paste( dtbf$ms$Cnd$PT[cnd], ' for ', dtbf$ms$Cnd$PD[cnd],
               ' s', sep = '' ) )
      
      # Extract model predictions for quantile function
      obj = quickdist( 'wr', type = 'QPE', prm = unlist( pm[cnd,] ) )
      points( obj, pch = 21, bg = 'black' );
      points( obj, ch = 0, pch = 24, bg = 'grey' )
      
      # Extract observed data
      kp = dtbf$PD == dtbf$ms$Cnd$PD[cnd] &
        dtbf$PT == dtbf$ms$Cnd$PT[cnd] & 
        dtbf$Co == dtbf$ms$Cnd$Co[cnd]
      df = data.frame( rt = dtbf$rt, ch = dtbf$ch,
                       subject = subject )
      # Obtain estimated quantile functions
      obs = rtplots( df, label = c( 'rt', 'ch' ), keep = kp, 
                     type = 'QPE', 
                     prob = c( 0, seq(.1,.9,.2), 1 ) )
      points( obs, val = 1, pch = 21, bg = 'white' )
      points( obs, val = 0, pch = 24, bg = 'white' )
      
      # Add additional legends
      if ( cnd == 4 | cnd == 8 ) legend( 'topright',
                                         dtbf$ms$Cnd$Co[cnd],
                                         bty = 'n' )
    }
  }
  
}

# Lookup - 01d
model_select = function( num ) {
  # Purpose: 
  # Generates a model name and a vector of 
  # structure types given a numeric input.
  # Arguments: 
  # A numeric input representing the model to select
  # Returns: 
  # A list with the model name and vector of structure 
  # types.
  
  # Models to fit
  all_mdl = c(
    # Same-different task
    'DRM_SD_M1_Null',   # (1)  3  parameters
    'DRM_SD_M2_Prior',  # (2)  4  parameters
    'DRM_SD_M3_Orig',   # (3)  13 parameters
    'DRM_SD_M4_Thresh', # (4)  19 parameters
    'DRM_SD_M5_Full',   # (5)  21 parameters
    'DRM_SD_M6_Sat',    # (6)  33 parameters
    'DRM_SD_M7_PM',     # (7)  25 parameters
    # Forced-choice task
    'DRM_FC_Null',      # (8)  3  parameters
    'DRM_FC_Prior',     # (9)  4  parameters
    'DRM_FC_Reg_Orig',  # (10) ?  parameters
    'DRM_FC_Orig',      # (11) ?  parameters
    'DRM_FC_Sat'        # (12) ? parameters
  )
  
  # Same-different task
  task = 'Same-different'
  if ( num == 1 ) # Null model
    ver = c( 'Intercept', 'Intercept', 'Intercept', 'Fixed' )
  if ( num == 2 ) # Prior model
    ver = c( 'Intercept', 'Binary', 'Intercept', 'Fixed' )
  if ( num == 3 ) # Original model
    ver = c( 'Bias', 'Orig (SD)', 'Intercept', 'Fixed' )
  if ( num == 4 ) # Threshold model
    ver = c( 'Saturated', 'Binary', 'Adaptive', 'Fixed' )
  if ( num == 5 ) # Full model
    ver = c( 'Bias_x_duration', 'Saturated', 'Adaptive', 'Fixed' )
  if ( num == 6 ) # Saturated model
    ver = c( 'Saturated', 'Saturated', 'Adaptive', 'Fixed' )
  if ( num == 7 ) # Prime matching model
    ver = c( 'Prime_match', 'Saturated', 'Adaptive', 'Fixed' )
  
  # Forced-choice task
  if ( num > 7 ) task = 'Forced-choice'
  if ( num == 8 )
    ver = c( 'Intercept', 'Intercept', 'Intercept', 'Fixed' )
  if ( num == 9 )
    ver = c( 'Intercept', 'Binary', 'Intercept', 'Fixed' )
  if ( num == 10 )
    ver = c( 'Intercept', 'Orig (FC)', 'Intercept', 'Fixed' )
  if ( num == 11 )
    ver = c( 'Intercept', 'Regression (Orig)', 'Adaptive', 'Fixed' )
  if ( num == 12 )
    ver = c( 'Bias_x_duration', 'Saturated', 'Adaptive', 'Fixed' )
  
  return( list( ver = ver, mn = all_mdl[num], task = task ) )
}

###
### Example with a single subject
###
# Lookup - 02

if ( runCode[1] ) {
  
  subject = 2
  tmp = model_select( model_type )
  ver = tmp$ver; mn = tmp$mn; task = tmp$task; rm( tmp )
  if ( model_type == 9 ) {
    setwd( 'Estimation_results' )
    load( 'DRM_SD_Orig_xi.RData' )
    setwd( '..' )
    xi_reg = scale( xi_DRM_SD_Orig[ subject, ] )
  }
  
  dtbf = NULL
  fit = estimate( ver, subject, task,
                  drm_dtbf_create )
  plot_jcdf_subject( fit, dtbf )
  rmse = compute_rmse( dtbf, fit, mle_fn )
  
}

###
### Estimation over subjects
###
# Lookup - 03

if ( runCode[2] ) {
  
  setwd( 'Estimation_results' )
  tmp = model_select( model_type )
  ver = tmp$ver; mn = tmp$mn; task = tmp$task; rm( tmp )
  fname = paste( mn, 'results.RData', sep = '_' )
  
  # Create a log file
  sink( paste( mn, '_log.txt', sep = '' ), split = TRUE )
  
  # Loop over subjects
  startTime = Sys.time()
  print( Sys.time() )
  print( 'Estimation started' )
  for (s in 1:N) {
    
    # Subject index
    subject = s
    
    if ( model_type == 9 ) {
      load( 'DRM_SD_Orig_xi.RData' )
      xi_reg = scale( xi_DRM_SD_Orig[ subject, ] )
    }
    
    # Estimate parameters
    dtbf = NULL
    fit = estimate( ver, subject, task,
                    drm_dtbf_create )
    
    # Extract relevant parts of estimation process
    cur_prm = coef(fit)[ which.max( fit$value ), ]
    pm = as.matrix( drm_pm_quick_convert( cur_prm, dtbf ) )
    rmse = compute_rmse( dtbf, fit, mle_fn )
    
    # Create storage for results
    if ( s == 1 ) {
      raw_prm = matrix( NA, N, length( cur_prm ) )
      colnames( raw_prm ) = dtbf$pnames
      rownames( raw_prm ) = 1:N
      all_prm = raw_prm
      all_pm = array( NA, dim = c( N, nrow( pm ), ncol( pm ) ) )
      pred = array( NA, dim = c( N, nrow( pm ), 11 ) )
      obs = array( NA, dim = c( N, nrow( pm ), 11 ) )
      mdl_fit = matrix( NA, N, 7 )
      colnames( mdl_fit ) = c( 'K', 'DF', 'LogLik', 'AIC', 'BIC',
                               'RMSE', 'Converged' )
      par_type = dtbf$type
    }
    
    # Number of observations
    n = length( dtbf$rt )
    # Number of parameters
    k = length( cur_prm )
    
    # Store parameter estimates
    raw_prm[s,] = cur_prm
    all_prm[s,] = tran_par( cur_prm, par_type )
    all_pm[s,,] = pm
    # Extract model diagnostics
    mdl_fit[s,1] = length( cur_prm )
    mdl_fit[s,2] = n - k
    mdl_fit[s,3] = fit$value[ which.max( fit$value ) ]
    mdl_fit[s,4] = -2 * mdl_fit[s,3] + 2*k # AIC
    mdl_fit[s,5] = -2 * mdl_fit[s,3] + log(n)*k # BIC
    mdl_fit[s,6] = rmse$rmse # Root mean square error
    # Convergence
    if ( fit$convcode[ which.max( fit$value ) ] == 0 ) 
      mdl_fit[s,7] = 1 else mdl_fit[s,7] = NA
    
    # Extract observed and predicted values
    for ( cnd in 1:nrow( dtbf$ms$Cnd ) ) {
      obj = quickdist( 'wr', type = 'QPE', 
                       prm = unlist( pm[cnd,] ) )
      pred[s,cnd,1:10] = c( obj$q_pv$x1, obj$q_pv$x0 )
      pred[s,cnd,11] = pwaldrace( Inf, 1, 
                                  k1 = pm[cnd,'k1'], 
                                  x1 = pm[cnd,'x1'], 
                                  t1 = pm[cnd,'t1'], 
                                  k0 = pm[cnd,'k0'],
                                  x0 = pm[cnd,'x0'], 
                                  t0 = pm[cnd,'t0'],
                                  s0 = pm[cnd,'s0'] )
      
      kp = dtbf$PD == dtbf$ms$Cnd$PD[cnd] &
        dtbf$PT == dtbf$ms$Cnd$PT[cnd] & 
        dtbf$Co == dtbf$ms$Cnd$Co[cnd]
      df = data.frame( rt = dtbf$rt, ch = dtbf$ch )
      cur_obs = by( df$rt[kp], list( df$ch[kp] ), 
                    quantile, prob = seq(.1,.9,.2) )
      if ( !is.null( cur_obs$`0` ) ) 
        obs[s,cnd,6:10] = cur_obs$`0`
      if ( !is.null( cur_obs$`1` ) ) 
        obs[s,cnd,1:5] = cur_obs$`1`
      obs[s,cnd,11] = mean( df$ch[kp] == 1 )
    }
    
    # Track progress
    print( paste( s, '/', N, ' done', sep = '' ) )
    
  }
  mdl_fit = as.data.frame( mdl_fit )
  
  ###
  ### Estimate sampling error via a bootstrap method
  ###
  
  # Test statistic
  T_x = median
  
  # Aggregate parameters
  gp = apply( all_pm, c(2,3), T_x )
  colnames( gp ) = c( 'k1', 'k0',
                           'x1','x0',
                           't1','t0',
                           's0' )
  gp = as.data.frame( gp )
  
  pred_bs = matrix( NA, 8, 11 )
  for ( i in 1:8 ) {
    for ( j in 1:0 ) {
      sel = 1:5 + 5*abs(j-1)
    pred_bs[i,sel] = qwaldrace( seq( .1, .9, .2 ), j,
                                k1 = gp$k1[i],
                                k0 = gp$k0[i],
                                x1 = gp$x1[i],
                                x0 = gp$x0[i],
                                t1 = gp$t1[i],
                                t0 = gp$t0[i],
                                s0 = gp$s0[i] )
    }
    pred_bs[i,11] = pwaldrace( Inf, 1,
                                k1 = gp$k1[i],
                                k0 = gp$k0[i],
                                x1 = gp$x1[i],
                                x0 = gp$x0[i],
                                t1 = gp$t1[i],
                                t0 = gp$t0[i],
                                s0 = gp$s0[i] )
  }
  rownames( pred_bs ) = 
    apply( dtbf$ms$SH, 1, paste, collapse = '-' )
  
  f = function( x ) {
    
    out = numeric( 11 )
    
    out[ 1:5 ] = tryCatch( 
      quantile( x$rt[ x$ch == 1 ], seq( .1, .9, .2 ) ),
      error = function(e) rep( NA, 5 ) )
    out[ 6:10 ] = tryCatch( 
      quantile( x$rt[ x$ch == 0 ], seq( .1, .9, .2 ) ),
      error = function(e) rep( NA, 5 ) )
    out[11] = mean( x$ch == 1 )
    
    return( out )
  }
  
  nRep = 10000 # Number of iterations to use for bootstrap
  prog = seq( 0, 1, .1 ) * nRep
  
  bootstrap_resid = array( NA, dim = c( nRep, 8, 11 ) )
  
  inc = 1
  for ( nr in 1:nRep ) {
    
    sim = rwaldrace( 640,
                     k1 = rep( gp$k1, each = 80 ),
                     k0 = rep( gp$k0, each = 80 ),
                     x1 = rep( gp$x1, each = 80 ),
                     x0 = rep( gp$x0, each = 80 ),
                     t1 = rep( gp$t1, each = 80 ),
                     t0 = rep( gp$t0, each = 80 ),
                     s0 = rep( gp$s0, each = 80 )
    )
    
    sim_bs = by( sim, list( rep(1:8,each=80) ), f )
    sim_bs = matrix( unlist(sim_bs), 8, 11, byrow = T )
    bootstrap_resid[nr,,] = sim_bs - pred_bs
    
    if ( nr > prog[inc] ) {
      string = paste( round( 100 * prog[inc]/nRep ),
                      '% done', sep = '' )
      print( string )
      inc = inc + 1
    }
  }
  
  sampling_error = list(
    pred_bs = pred_bs,
    bootstrap_resid = bootstrap_resid,
    ui = apply( bootstrap_resid, c(2,3), 
                quantile, prob = c(.025,.975),
                na.rm = T )
  )
  
  ###
  ### Save results
  ###
  
  save( raw_prm, mdl_fit, pred, obs, task, ver, mn,
        all_pm, all_prm, par_type, sampling_error, file = fname )
  print( 'Estimation finished' )
  runTime = Sys.time() - startTime
  print( runTime )
  
  sink()
}

setwd( cur_dir )