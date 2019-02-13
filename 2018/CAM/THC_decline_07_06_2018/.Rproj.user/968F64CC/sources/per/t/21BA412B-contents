# Exponential decay functions
# Written by Kevin Potter
# email: kevin.w.potter@gmail.com
# Please email me directly if you 
# have any questions or comments
# Last updated 2018-11-20

# Table of contents
# 1) Initial setup
# 2) Define functions
#   2.1) ed_model
#   2.2) mle_ed_tobit
#   2.3) estimate_ed_lm
#   2.4) plot_decay
#   2.5) estimate_ed_brms
#   2.6) extract_brm_res
#   2.7) compute_coverage_prob
#   2.8) compute_R2

###
### 1) Initial setup
###

# Package for Bayesian estimation
Sys.setenv(USE_CXX14 = 1) # For Rstan to work
my_package_load( 'brms' )

# Make sure 'loo' package is installed
# install.packages( 'loo' )

###
### 2) Define functions
###

# 2.1)
ed_model = function( param, x, log = F ) {
  # Purpose:
  # Computes the predicted THCCOOH values 
  # using an exponential decay model.
  # Arguments:
  # param - A vector with the starting level of 
  #         THCCOOH and the elimation rate
  # x     - A vector with the days since THC 
  #         ingestion
  # log   - Logical; if TRUE, returns the 
  #         log THCCOOH values instead
  # Returns:
  # The predicted log THCCOOH levels per each day.
  
  yhat = log( param[1] ) - param[2] * x
  if ( !log ) yhat = exp( yhat )
  
  return( yhat )
}

# 2.2)
mle_ed_tobit = function( param, dat ) {
  # Purpose:
  # Computes the negative of the summed log-likelihoods 
  # for a tobit regression approach to fit an 
  # exponential decay model with log-normal error.
  # Arguments:
  # param - The starting value, elimination rate, and the 
  #         residual standard deviation parameters.
  # dat   - A list with the days since THC ingestion, 
  #         the THCCOOH levels, and the cut-off below 
  #         which values are rounded down to 0.
  # Returns:
  # The negative of the summed log-likelihoods.
  
  # Extract data
  x = dat$x # Independent variable
  y = dat$y # Dependent variable
  a = dat$a # Cut-off for censored data
  
  # Indicator function for censored data
  I_a = y < a
  y_star = y; y_star[ I_a ] = a
  y_star = log( y_star )
  
  # Restrict slope to be positive
  param[2] = exp( param[2] )
  
  # Compute mean
  mu = ed_model( param[1:2], x, log = T )
  # Extract standard devation for residuals
  sigma = param[3]
  
  # Carry out tobit regression
  p1 = I_a * pnorm( ( y_star - mu )/sigma, log.p = T )
  p2 = (1-I_a) * (dnorm( (y_star - mu)/sigma, log = T )- log( sigma ) )
  sll = sum( p1 + p2 )
  if ( is.na( sll ) ) sll = -Inf
  
  return( -sll )
}

# 2.3)
estimate_ed_lm = function( x, y ) {
  # Purpose:
  # Fits an exponential decay model via 
  # tobit log-linear regression.
  # Arguments:
  # x    - The independent variable
  # y    - The dependent variable
  # Returns:
  # A list with the parameter estimates,
  # the predicted THCCOOH levels, 
  # and the residals.
  
  # Straightforward estimation of the 
  # power law and exponential decay models 
  # requires log transformations, so we 
  # first check for zeroes in the dependent 
  # variable
  is_zero = y == 0
  
  # If there are no zero values, we can simply 
  # use a log-linear model to estimate the 
  # 3 parameters in the model
  if ( !any( is_zero ) ) {
    
    lmf = lm( log( y ) ~ 1 + x )
    est = c( coef( lmf ), sigma( lmf ) )
    est[1] = exp( est[1] )
    est[2] = -est[2]
    
  } else {
    
    # If there are zeros, we can treat them 
    # as censored data (i.e., the actual 
    # values are actually non-zero, but 
    # any value less than a cut-off was 
    # rounded down to zero). Therefore, 
    # we can use tobit regression to 
    # fit the log-linear model
    
    # First, we use the naive estimator 
    # to determine the starting values 
    # for the optimization routine
    sel = y > 0
    ne = lm( log( y[sel] ) ~ 1 + x[sel] )
    st_val = c( coef( ne ), sigma( ne ) )
    st_val[1] = exp( st_val[1] )
    st_val[2] = -st_val[2]
    
    # Set up data to analyze
    dat = list(
      x = x,
      y = y,
      a = 1
    )
    # We can then use the tobit regression likelihood 
    # to estimate the parameters of interest despite 
    # the censored data
    tbtf = tryCatch( {
      st_val[2] = log( abs( st_val[2] ) ); st_val[3] = 1;
      optim( st_val, mle_ed_tobit, dat = dat,
             control = list( maxit = 10000 ) )
      },
      error = function(e) NULL
    )
    
    # If the MLE worked
    if ( !is.null( tbtf ) ) {
      
      est = tbtf$par
      est[2] = exp( est[2] )
      
    } else {
      est = st_val
      warning( 'Tobit regression failed' )
    }
    
  }
  names( est ) = c( 'Baseline', 'Elimination', 'Error' )
  
  # Generate predicted THCCOOH levels
  pred = ed_model( est, x )
  # Compute residuals
  resid_val = y - pred
  
  return( list( est = est, 
                x = x, 
                pred = pred, 
                resid = resid_val ) )
  
}

# 2.4) 
plot_decay = function( dat, new = T, estimate = T ) {
  # Purpose:
  # Plots the decay in THCCOOH and log THCCOOH levels 
  # over days since ingestion.
  # Arguments:
  # dat      - A data frame with the Time, THCCOOH,
  #            and ID variables. If there are 
  #            multiple subjects, the average 
  #            will also be plotted
  # new      - Logical; if true, a new plotting 
  #            pane is generated
  # estimate - Logical; if true, the 
  #            exponential decay model is fit 
  #            to the data using a tobit regression 
  #            approach.
  
  # Extract data
  x = dat$Time
  y = dat$THCCOOH
  no_na = !is.na( x ) & !is.na( y )
  x = x[ no_na ]
  y = y[ no_na ]
  y_log = log( y )
  id = dat$ID
  Ns = length( unique( id ) )
  
  clr = rep( 'black', length( y ) )
  clr[ y == 0 ] = 'grey'
  
  if ( Ns > 1 ) {
    
    plt = dat %>% 
      group_by( Time ) %>% 
      summarise(
        R = mean( THCCOOH )
      )
    plt2 = dat %>%
      filter( THCCOOH > 0 ) %>% 
      group_by( Time ) %>% 
      summarise(
        R = mean( log( THCCOOH ) )
      )
    
    clr = rep( 'grey', length( y ) )
  }
  
  if ( estimate ) {
    est = estimate_ed_lm( x, y )
    o = order( est$x )
    x_est = est$x[o]
    y_hat = est$pred[o]
  }
  
  # Create two panes
  if ( new ) x11( width = 12 )
  layout( cbind( 1, 2 ) )
  
  ### Raw data
  
  # Create a blank plot
  xl = c( -1, 35 )
  yl = lowerUpper( 100, y )
  yl[1] = -7.5
  blankPlot( xl, yl )
  
  # Add observations
  points( x, y, pch = 19, col = clr )
  
  # Add mean
  if ( Ns > 1 ) lines( plt$Time, plt$R, lwd = 2 )
  # Add estimates
  if ( estimate ) lines( x_est, y_hat, lwd = 2, col = 'blue' )
  
  # Add axes and labels
  customAxes( xl, yl )
  
  # x-axis
  axis( 1, seq( 0, 32, 4 ),
        tick = F, line = -1, cex.axis = 1.25 )
  mtext( 'Days since THC ingested',
         side = 1, line = 2, cex = 1.25 )
  
  # y-axis
  axis( 2, round( seq( 0, yl[2], length = 5 ) ),
        tick = F, line = -1, cex.axis = 1.25 )
  mtext( 'THCCOOH - ng/mL',
         side = 2, line = 2, cex = 1.25 )
  
  ### Log-transformed data
  
  # Create a blank plot
  sel = y_log > -Inf & y_log < Inf
  
  if ( sum(sel) > 1 ) {
    yl = lowerUpper( 1, y_log[sel] )
  } else {
    yl = lowerUpper( 1, c( 0, y_log[sel] ) )
  }
  if ( estimate ) {
    yl = lowerUpper( 1, log( y_hat[ y_hat > 0 ] ) )
  }
  blankPlot( xl, yl )
  
  horizLines( log(1), xl, lwd = 2, lty = 2, col = 'grey' )
  
  customAxes( xl, yl,
              label = 
                c( 'Days since THC ingestion', 'THCCOOH - log(ng/mL)' ),
              axSz = 1.25 )
  axis( 1, seq( 0, 32, 4 ),
        tick = F, line = -1, cex.axis = 1.25 )
  axis( 2, round( seq( yl[1], yl[2], length = 5 ), 1 ),
        tick = F, line = -1.5, cex.axis = 1.25 )
  
  points( x, y_log, pch = 19, col = clr )
  
  if ( estimate & Ns == 1 ) {
    if ( any( y_hat < 1 ) ) {
      points( x[ y_hat < 1 ], log( y_hat[ y_hat < 1 ] ),
              pch = 19, col = 'grey' )
    }
  }
  
  # Add mean
  if ( Ns > 1 ) lines( plt2$Time, plt2$R, lwd = 2 )
  # Add estimates
  if ( estimate ) {
    sel = log( y_hat ) > yl[1]
    lines( x_est[sel], log( y_hat[sel] ), lwd = 2, col = 'blue' )
  }
  
}

# 2.5) 
estimate_ed_brms = function( dtbf, 
                             covariates,
                             algorithm = 
                               list(
                                 wup = 1000,
                                 d_iter = 10000/4,
                                 thin = 1,
                                 seed = NA ),
                             horseshoe = F ) {
  # Purpose:
  # Fits an exponential decay model to 
  # THCCOOH levels over days with both 
  # population and subject-level effects.
  # Arguments:
  # dtbf       - A data frame with the 
  #              log THCCOOH levels, an 
  #              intercept vector 'Start_point',
  #              days coded as 'Elimination_Rate',
  #              and the standardized predictors
  #              with the tags '.SP' or '.ER' 
  #              to denote with non-linear 
  #              parameter they are associated with
  # covariates - An optional list with two vectors,
  #              'SP', and 'ER', indicating the 
  #              predictors to include for 
  #              the start point and elimination
  #              rate parameters
  # algorithm  - A list with the warmup iterations,
  #              number of desired samples for each 
  #              chain, the thinning value, and 
  #              the seed for the random number 
  #              generator
  # horseshoe  - Logical; if true, instead fits 
  #              a model with all predictors 
  #              included and a horseshoe prior 
  #              to shrink small effects without 
  #              overly impacted larger effects
  # Returns:
  # A list with the brms output,
  # the estimated LOO-CV, the 
  # estimation run time, and 
  # the formula used to fit 
  # the data.
  
  # Estimate run time
  tick()
  
  # Extract settings for estimation 
  # algorithm
  
  # Warm up period
  wup = algorithm$wup
  # Desired number of iterations
  d_iter = algorithm$d_iter
  # Thinning value
  thin = algorithm$thin
  n_iter = d_iter * thin + wup
  seed = algorithm$seed
  
  # If not using a horseshoe prior
  if ( !horseshoe ) {
    
    ### Create estimation formula
    
    # Create brms formula
    lhs = 'log_THCCOOH ~ 0 +'
    
    # Check for additional predictors of 
    # the start point
    if ( !is.null( covariates$SP ) ) {
      
      SP = paste( 'Start_point +',
                  paste( 
                    covariates$SP, 
                    collapse = ' + ' ),
                  '+' )
    } else {
      SP = 'Start_point +'
    }
    
    # Check for additional predictors of 
    # slope
    if ( !is.null( covariates$ER ) ) {
      
      # Create formula
      ER = paste( 'Elimination_rate +',
                  paste( 
                    covariates$ER, 
                    collapse = ' + ' ),
                  '+ ' )
    } else {
      ER = 'Elimination_rate +'
    }
    
    # Subject-level effects
    SLE = '(0 + Start_point + Elimination_rate | ID)'
    
    frm = paste( lhs, SP, ER, SLE )
    
    ### Specify priors
    
    # Weakly informative prior on residual error
    prior = set_prior( 'student_t( 19, .33, .22 )', 
                       class = 'sigma' )
    
    # Prior on starting point
    prior = c( prior,
               set_prior( 'normal( 4.6, 1.1 )', 
                          class = 'b', coef = 'Start_point' ) )
    
    # Prior on elimination rate
    prior = c( prior,
               set_prior( 'student_t( 15, .23, .16 )', 
                          class = 'b', coef = 'Elimination_rate' ) )
    
    # Prior on correlation between random effects
    prior = c( prior,
               set_prior( 'lkj_corr_cholesky(1)', 
                          class = 'L' ) )
    
    if ( !is.null( covariates$SP ) ) {
      # Set mild regularizing priors on predictors
      for ( i in 1:length( covariates$SP ) ) {
        prior = c( prior,
                   set_prior( 'normal( 0, 1 )', 
                              class = 'b', 
                              coef = covariates$SP[i] ) )
      }
    }
    
    if ( !is.null( covariates$ER ) ) {
      # Set mild regularizing priors on predictors
      for ( i in 1:length( covariates$ER ) ) {
        prior = c( prior,
                   set_prior( 'normal( 0, 1 )', 
                              class = 'b', 
                              coef = covariates$ER[i] ) )
      }
    }
    
    # Fit model using brms
    mf = brm( as.formula( frm ),
              # Data and priors
              data = dtbf, prior = prior, 
              # Estimation settings
              warmup = wup,
              iter = n_iter,
              chains = 4,
              cores = 4,
              thin = thin, 
              control = list( adapt_delta = 0.995,
                              max_treedepth = 15 ) )
    
    
    
  } else {
    
    ### Horshoe prior example
    
    # Weakly informative prior on residual error
    prior = set_prior( 'student_t( 10, .33, .5 )', 
                       class = 'sigma' )
    # Prior on starting point
    prior = c( prior,
               set_prior( 'normal( 4.6, 2.5 )', 
                          class = 'b', coef = 'Start_point',
                          nlpar = 'eta1' ) )
    
    # Prior on elimination rate
    prior = c( prior,
               set_prior( 'student_t(15, .23, .16 )', 
                          class = 'b', coef = 'Elimination_rate',
                          nlpar = 'eta1' ) )
    
    # Prior on correlation between random effects
    prior = c( prior,
               set_prior( 'lkj_corr_cholesky(1)', 
                          class = 'L' ) )
    
    # Set a horseshoe prior on predictors to shrink 
    # small effects and increase large effects
    prior = c( prior,
               prior( horseshoe(1), nlpar = 'eta2' ) )
    
    frm = bf( log_THCCOOH ~ eta1 + eta2,
              eta1 ~ 0 + Start_point + Elimination_rate + 
                (0 + Start_point + Elimination_rate|ID),
              eta2 ~ 0 + 
                # Predictors for start point
                Sex.SP + 
                zBMI.SP + 
                zYears_of_MJ_use.SP + 
                zLevel_of_MJ_use.SP + 
                Race.SP + 
                # Predictors for elimination rate
                Sex.ER + 
                zBMI.ER + 
                zYears_of_MJ_use.ER + 
                zLevel_of_MJ_use.ER + 
                Race.ER,
              nl = TRUE )
    
    mf = brm( frm, 
              # Data and priors
              data = dtbf, prior = prior, 
              # Estimation settings
              warmup = wup,
              iter = n_iter,
              chains = 4,
              cores = 4,
              thin = thin, 
              seed = seed, 
              control = list( adapt_delta = 0.9999,
                              max_treedepth = 15 ) )
    
  }
  
  # Compute approximation to leave-one-out
  # cross-validation (LOOCV)
  loocv = loo::loo( mf )
  
  tock()
  
  out = list(
    model_fit = mf,
    loocv = loocv,
    run_time = run_time,
    formula = frm
  )
  
  return( out )
}

# 2.6)
extract_brm_res = function( mf ) {
  # Purpose:
  # A function to extract the estimated 
  # effect size, posterior standard deviation,
  # 95% credible interval, and posterior 
  # p-value for the regression coefficients,
  # subject-level variances, and measurement 
  # error variability.
  # Arguments:
  # mf - The brms output object
  # Returns:
  # A matrix with the effect size, 
  # posterior standard deviation, 
  # 95% credible intervals, and 
  # one-sided posterior p-values for 
  # a subset of parameters.
  
  # Function to compute posterior p-values
  pv = function( x ) {
    if ( mean(x) > 0 ) {
      p = sum( x < 0 )/length( x )
    } else {
      p = sum( x > 0 )/length( x )
    }
    return( p )
  }
  
  # Extract summary of results
  sm = summary( mf )
  
  # Combine results of interest
  res = rbind(
    sm$fixed[,1:4],
    sm$random$ID[,1:4],
    Measurement_error = sm$spec_pars[,1:4]
  )
  colnames( res ) = c(
    'Effect_size',
    'SD',
    'Lower_CI_95',
    'Upper_CI_95'
  )
  
  # Compute posterior p-values
  out = cbind( res, p_value = NA )
  
  # Extract posterior samples
  pst = as.matrix( mf )
  # Compute posterior p-values
  out[,'p_value'] = apply( pst[,1:nrow(out)], 2, pv )
  
  out = data.frame(
    Variable = c(
      rownames( sm$fixed ),
      'sd(Start_point)',
      'sd(Elimination_rate)',
      'cor(Start_point,Elimination_rate)',
      'Measurement_error'
    ),
    X = out,
    stringsAsFactors = F
  )
  
  return( out )
}

# 2.7)
compute_coverage_prob = function( mf, test ) {
  # Purpose:
  # Computes the proportion of observations 
  # that fall within the 95% credible 
  # prediction intervals.
  # Arguments:
  # mf   - A brms output object
  # test - A data frame with the holdout observations
  # Returns:
  # The proportion of observations that fell 
  # within the 95% credible prediction intervals.
  
  # Generate predictions based on 
  # holdout sample
  pred = predict( mf$model_fit, 
                  newdata = test, allow_new_levels = T )
  # Compute coverage probabilities
  coverage_prob = 
    sum( pred[,3] < test$log_THCCOOH & 
           pred[,4] > test$log_THCCOOH )/nrow( test )
  return( coverage_prob )
}

# 2.8)
compute_R2 = function( mf, test ) {
  # Purpose:
  # Computes the R-squared coefficient 
  # between the model predictions and the 
  # observed data for a holdout test 
  # sample.
  # Arguments:
  # mf   - A brms output object
  # test - A data frame with the holdout observations
  # Returns:
  # The R-squared coefficient, its standard deviation, 
  # and its 95% credible interval.
  
  pred = predict( mf$model_fit, 
                  newdata = test, 
                  allow_new_levels = T,
                  summary = F )
  all_r2 = apply( pred, 1, function(x) cor( x, test$log_THCCOOH )^2 )
  
  out = c(
    R2 = mean( all_r2 ),
    SD = sd( all_r2 ),
    CI_2.5 = quantile( all_r2, .025 ),
    CI_97.5 = quantile( all_r2, .975 )
  )
  
  return( out )
}


setwd( R_dir )
