# Estimation functions
# Written by Kevin Potter
# email: kevin.w.potter@gmail.com
# Please email me directory if you 
# have any questions or comments
# Last updated 2018-04-23

# Table of contents
# 1) Initial setup
# 2) Functions (General purpose)
#   2.1) convergenceExtract
#   2.2) findDec
#   2.3) plotConvergence
#   2.4) quick_CI_plot
# 3) Functions (Mixed effects)
#   3.1) estimate_sdt
#   3.1) Compile C++ code for posterior predictive checks
# 4) Functions (fNIRS)
#   4.1) estimate_fNIRS

###
### 1) Initial setup
###

setwd( proj_dir )
setwd( 'Stan_scripts' )
stan_dir = getwd()

# Compile stan script for estimating SDT 
# model with mixed effects
if ( !exists('sdt_mix') ) {
  message( 'Estimation script' )
  comp_time = Sys.time()
  sdt_mix = stan_model(stanc_ret = 
                         stanc_builder("SDT_mixed_effects_fit.stan"))
  comp_time = Sys.time() - comp_time
  message(
    paste( 'Compilation time:', round(comp_time,2), 
           attributes(comp_time)$units ) )
  rm( comp_time )
}

# Compile stan script for estimating SDT 
# model with mixed effects
if ( !exists('fNIRS_fit') ) {
  message( 'Estimation script' )
  comp_time = Sys.time()
  fNIRS_fit = stan_model(stanc_ret = 
                           stanc_builder("fNIRS_fit.stan"))
  comp_time = Sys.time() - comp_time
  message(
    paste( 'Compilation time:', round(comp_time,2), 
           attributes(comp_time)$units ) )
  rm( comp_time )
}

setwd( proj_dir )

###
### 2) Functions (General purpose)
###

# 2.1)
convergenceExtract = function( fit, parName = NULL ) {
  # Purpose:
  # Extract convergence diagnostics from a Stan fit object.
  # Arguments:
  # fit     - A Stan fit object
  # parName - An optional string giving the final parameter label 
  #           of the subset of the output to include
  # Notes:
  # Extracting the convergence statistics can be slow, especially 
  # when a large number of parameters were stored.
  # Returns:
  # A list with the Gelman-Rubin convergence statistic, the 
  # effective number of samples, and the total number of samples.
  
  Rhat = summary(fit)$summary[,"Rhat"]
  n_eff = summary(fit)$summary[,"n_eff"]
  totSampSize = length(extract(fit, pars = "lp__")[[1]])
  # We're only interested in a subset of parameters
  if ( length( parName ) == 0 ) 
    parName = names( Rhat )[ 
      which( names( Rhat ) == "logLik[1]" ) - 1 ]
  sel = 1:max( grep( parName, names(Rhat) ) )
  Rhat = Rhat[sel]; n_eff = n_eff[sel];
  
  return( list( Rhat = Rhat, n_eff = n_eff, 
                totSampSize = totSampSize ) )
}

# 2.2)
findDec = function( x, spacing = 10 ) {
  # Purpose:
  # Determines the rounded leading digit and the 
  # number of trailing zeros for a number or the 
  # number of decimal places.
  # Arguments:
  # x       - A vector of values
  # spacing - The value whose exponent should be increased
  # Returns:
  # A vector giving the leading digit, the number of 
  # trailing/leading zeros, the same but in scientific 
  # notation, and 1 if it's trailing zeros, -1 if it's 
  # decimal places.
  
  mx = max( x )
  rnd = mx
  
  if ( round( mx ) > 1 ) {
    
    inc = 0;
    
    while ( rnd > 1 ) {
      inc = inc + 1;
      rnd = round( mx/( spacing^inc ) )
    }
    
    v = round( mx/spacing^(inc-1) )
    f = spacing^(inc-1)
    out = c( v,f,inc-1,1)
  }
  
  if ( round( mx ) == 1 ) {
    
    out = c( 1, 1, 1, 0 )
    
  }
  
  if ( round( mx ) == 0 ) {
    
    inc = 0;
    
    while ( rnd < 1 ) {
      inc = inc + 1;
      rnd = round( mx*( spacing^inc ) )
    }
    
    v = round( mx*spacing^(inc) )
    f = spacing^(inc)
    out = c( v,f,inc,-1)
    
  }
  
  return( out )
}

# 2.3)
plotConvergence = function( conv, savePlot ) {
  # Purpose:
  # Generates a plot of the Gelman-Rubin statistics 
  # and the effective number of samples for the 
  # marginal posterior samples of the parameters.
  # Arguments:
  # conv     - Output from the 'convergenceExtract' 
  #            function
  # savePlot - A logical value, which when false generates a 
  #            new plotting window

  # Remove NaN for correlation matrix
  sel = which( is.na( conv$Rhat ) )
  if ( length( sel ) == 0 ) sel = -(1:length(conv$Rhat))
  conv$Rhat = conv$Rhat[ -sel ]
  conv$n_eff = conv$n_eff[ -sel ]
  
  if (!savePlot) x11(width=12)
  layout( cbind(1,2) )
  
  # Plot a histogram of the Gelman-Rubin statistics for the 
  # marginal posterior samples of the parameters
  
  tmp = hist( conv$Rhat, plot = F )
  scl = findDec( tmp$density )
  if ( scl[4] == 1 ) scl = scl[1]*(scl[2]/10) else 
    scl = scl[1]/(scl[2])
  
  yl = lowerUpper( scl, tmp$density )
  yl[1] = 0
  
  xl = lowerUpper( .1, conv$Rhat )
  xl[2] = max( xl[2], 1.12 )
  xl[1] = min( xl[1], .98 )
  
  plot( xl, yl, type = 'n', cex.axis = 1.5, cex.lab = 1.5,
        xlab = expression( hat(R) ), ylab = 'Density',
        bty = 'l', main = 'Gelman-Rubin Statistic' )
  
  segments( tmp$mids, rep(0,length(tmp$mids)),
            tmp$mids, tmp$density, lwd = 3,
            col = 'grey' )
  abline( v = 1.1, lty = 2, lwd = 2 )
  
  # Plot a histogram of the effective sample size for the 
  # set of parameters
  
  tmp = hist( conv$n_eff, plot = F )
  scl = findDec( tmp$density )
  if ( scl[4] == 1 ) scl = scl[1]*(scl[2]/10) else 
    scl = scl[1]/(scl[2])
  
  yl = lowerUpper( scl, tmp$density )
  yl[1] = 0
  
  xl=c(0,conv$totSampSize)
  plot( xl, yl, type = 'n', cex.axis = 1.5, cex.lab = 1.5,
        xlab = expression(N[sample]), ylab = 'Density',
        bty = 'l', main = 'Effective sample size' )
  
  segments( tmp$mids, rep(0,length(tmp$mids)),
            tmp$mids, tmp$density, lwd = 3,
            col = 'grey' )
  
}

# 2.4)
quick_CI_plot = function( pos, i, w, CI_clr, pred, center ) {
  # Purpose:
  # A convenience function to add 95% and 68% credible 
  # intervals and a desired center to an existing plot.
  # Arguments:
  # pos    - The x-axis position at which to draw the CI
  # i      - The index in the 
  # w      - The half-width of the interval
  # CI_clr - The colors to use for the two CIs
  # pred   - The data frame with the descriptive stats 
  #          for the posterior predictive check
  # center - The type of central measure to use 
  #          (either 'Mean', 'Mode', or 'Median')
  
  xb = pos + c( -w, -w, w, w )
  yb = pred[ i, c( 'CI_2.5', 'CI_97.5',
                   'CI_97.5', 'CI_2.5' ) ]
  polygon( xb, yb, col = CI_clr[1], border = NA )
  
  xb = pos + c( -w, -w, w, w )
  yb = pred[ i, c( 'CI_16', 'CI_84',
                   'CI_84', 'CI_16' ) ]
  polygon( xb, yb, col = CI_clr[2], border = NA )
  
  horizLines( pred[ i, center ], pos + c(-w,w),
              lwd = lnSz, col = 'white' )
  
}

###
### 3) Functions (Mixed effects)
###

# 3.1)
estimate_sdt = function( df, 
                         X, 
                         random, 
                         Priors, 
                         var_names = c( 'Counts', 
                                        'Trials',
                                        'Target',
                                        'Subject' ),
                         niter = 1667,
                         chains = 6,
                         warm = 1000,
                         adapt_delta = .995,
                         max_treedepth = 14,
                         debug = F, 
                         ... ) {
  # Purpose:
  # Given a data frame of observations and 
  # a list of two design matrices, estimates 
  # the posterior for a basic SDT model 
  # with a set of random effects over subjects.
  # Arguments:
  # df            - The data frame with all 
  #                 observations to be modeled
  # X             - A list with the design matrices 
  #                 for the d' and criterion parameters,
  #                 respectively
  # random        - A vector of the positions in the 
  #                 'group_param' variable (the combined 
  #                 set of coefficients for the d' and 
  #                 criterion parameters) for the 
  #                 random variables
  # var_names     - A vector specifying the name  
  #                 of the observed counts, total 
  #                 number of trials, the type of trial,
  #                 and the subject index in the data 
  #                 frame
  # n_iter        - The number of post-warm-up iterations
  #                 to sample for each chain
  # chains        - The number of chains to run
  # warmup        - The number of iterations for the 
  #                warm
  # adapt_delta   - A control parameter for Stan's 
  #                 estimation algorithm controlling 
  #                 the resolution of the sampler
  # max_treedepth - A control parameter for Stan's 
  #                 estimation algorithm controlling 
  #                 the maximum tree depth for the 
  #                 N-U-Turn-Sampler
  # Returns:
  # A list where...
  #   post     = A list of the posterior estimates 
  #              returned by Stan;
  #   conv     = A list with the convergence diagnostics 
  #              obtained from Stan;
  #   random   = A vector indicating the positions of 
  #              the random variables in the overall 
  #              set of coefficients;
  #   run_time = How long the code took to run.
  
  # Initialize variable for run time duration
  run_time = Sys.time()
  
  # Input for stan script
  stan_dat = list(
    # Number of observations
    No = nrow( df ),
    # Number of subjects
    Ns = max( df[ ,var_names[4] ] ),
    # Number of group-level d' parameters
    Nd = ncol( X[[1]] ),
    # Number of group-level criterion parameters
    Nc = ncol( X[[2]] ),
    # Number of subject-level parameters
    Nsp = length( random ), 
    # Observed counts
    Y = df[ ,var_names[1] ], 
    # Total number of trials per observation
    Nt = df[ ,var_names[2] ],
    # Presence of target (1 = yes)
    Co = df[ ,var_names[3] ],
    # Subject index
    subjIndex = df[ ,var_names[4] ],
    # Design matrix for d'
    Xd = X[[1]],
    # Design matrix for criterion
    Xc = X[[2]],
    # Index for coefficients with random effects
    eta_pos = random, 
    # Matrix with parameters governing priors
    Priors = Priors
  )
  
  # Run stan script
  fit = sampling( sdt_mix, 
                  data = stan_dat, 
                  warmup = warm, 
                  iter = warm + niter, 
                  chains = chains, 
                  control = list( adapt_delta = adapt_delta,
                                  max_treedepth = max_treedepth ),
                  ... )
  
  # Extract convergence diagnostics
  conv = convergenceExtract( fit, 'Tau' )
  
  # Extract posterior samples
  post = rstan::extract( fit )
  
  # Compute run time duration
  run_time = Sys.time() - run_time
  
  # Return output
  out = list(
    post = post,
    convergence = conv,
    random = random, 
    runtime = run_time )
  
  if ( debug )
    out$fit = fit
  
  return( out )
}

# 3.2)
# Compile C++ code for posterior predictive checks
setwd( proj_dir )
setwd( 'src' )
Rcpp::sourceCpp( file = 
  'Functions_for_posterior_predictive_simulations.cpp' )
setwd( proj_dir )

###
### 4) Functions (fNIRS)
###

# 4.1)
estimate_fNIRS = function( df, X, Priors,
                           niter = 1667,
                           chains = 6,
                           warm = 1000,
                           adapt_delta = .995,
                           max_treedepth = 14,
                           debug = F, 
                           ... ) {
  # Purpose:
  # ...
  # Arguments:
  # ...
  # Returns:
  # ...
  
  # Initialize variable for run time duration
  run_time = Sys.time()
  
  # Input for stan script
  stan_dat = list(
    No = nrow( df ),
    Ns = max( df$Subject ),
    Nroi = max( df$indexROI ), 
    Nc = ncol( X ),
    indexSubject = df$Subject,
    indexROI = df$indexROI,
    X = X,
    zHbO = df$zHbO,
    nu = df$nu,
    priors_beta = Priors$priors_beta,
    priors_Omega = Priors$priors_Omega,
    priors_Tau = Priors$priors_Tau
  )
  
  # Run stan script
  fit = sampling( fNIRS_fit, 
                  data = stan_dat, 
                  warmup = warm, 
                  iter = warm + niter, 
                  chains = chains, 
                  control = list( adapt_delta = adapt_delta,
                                  max_treedepth = max_treedepth ),
                  ... )
  
  # Extract convergence diagnostics
  conv = convergenceExtract( fit, 'Tau' )
  
  # Extract posterior samples
  post = rstan::extract( fit )
  
  # Compute run time duration
  run_time = Sys.time() - run_time
  
  # Return output
  out = list(
    post = post,
    convergence = conv,
    runtime = run_time )
  
  if ( debug )
    out$fit = fit
  
  return( out )
}

