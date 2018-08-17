# Estimating SDT models using probit regression
# Written by Kevin Potter
# email: kevin.w.potter@gmail.com
# Please email me directly if you 
# have any questions or comments
# Last updated 2018-05-01

# Table of contents
# 1) Initial setup
#   1.1) link_function
# 2) Simple mixed effects probit regression example
# 3) Mixed effects probit regression example
# 4) SDT model with random intercepts for bias/d'
# 5) Realistic example based on experimental design

###
### 1) Initial setup
###

# Save current working directory
orig_dir = getwd()

# Package for manipulating data frames
# install.packages( 'dplyr' )
library( dplyr )

# Package for GLM mixed effects estimation
# install.packages( 'lme4' )
library( lme4 )

# Package for GLM mixed effects estimation (Bayesian)
# install.packages( 'rstanarm' )

# Indicators for which code segments to run
run_code = c(
  F, # Logistic/Probit regression (Random intercept)
  F, # Logistic/Probit regression (Random slope)
  F, # SDT model with random intercepts for bias/d'
  T  # Realistic example based on experimental design
)

# Indicators for estimation type
est_type = c(
  F, # lme4
  F, # Prior predictive check
  T, # Bayesian estimation
  F  # Posterior predictive check
)

# 1.1)
link_function = function( x, type = 'probit', reverse = F ) {
  # Purpose:
  # Either 1) converts an unbounded value to a probability
  # bounded between 0 and 1 using one of two types of 
  # link functions, or 2) converts a probability to 
  # an unbounded value with the corresponding 
  # inverse link functions.
  # Arguments:
  # x       - A vector of numeric values
  # type    - The type of link function to use, either
  #           'probit' = The standard normal cumulative distribution
  #           'logistic' = The logistic function
  # reverse - If true, applies the inverse of the specified 
  #           link function
  # Returns:
  # The converted values after applying the link function.
  
  out = NULL
  
  if ( reverse ) {
    if ( type == 'logistic' ) {
      out = log( x/(1-x) )
    }
    
    if ( type == 'probit' ) {
      out = qnorm( x )
    }
  } else {
    if ( type == 'logistic' ) {
      out = 1/( 1 + exp( -x ) )
    }
    
    if ( type == 'probit' ) {
      out = pnorm( x )
    }
  }
  
  return( out )
}

###
### 2) Logistic/Probit regression (Random intercept)
###

if ( run_code[1] ) {
  
  lf = 'logistic'
  
  ### Design
  
  Ns = 100
  # Design
  design = data.frame(
    I = 1, 
    A = rep( 0:1, each = 2 ),
    B = rep( 0:1, 2 ),
    N = rep( 50, 4 )
  )
  design$AxB = design$A * design$B
  # Data frame to store simulated data
  sim = data.frame(
    Subject = rep( 1:Ns, each = nrow( design ) ), 
    I = 1,
    A = rep( design$A, Ns ),
    B = rep( design$B, Ns ),
    AxB = rep( design$AxB, Ns ),
    Y = NA,
    N = rep( design$N, Ns )
  )
  
  ### Generating parameters
  
  # Regression coefficients
  beta = runif( 4, -.5, .5 )
  # Subject-level variability
  sigma = runif( 1, .25, 1 )
  # Subject-level deviations
  eta = rnorm( Ns, 0, sigma )
  
  # Group-level means
  X = as.matrix( sim[,c('I','A','B','AxB')] )
  sim$X_beta =  as.numeric( X %*% cbind( beta ) ) + 
    eta[ sim$Subject ]
  sim$Theta = link_function( sim$X_beta, type = lf )
  sim$Y = rbinom( nrow( sim ), sim$N, sim$Theta )
  sim$Z = sim$N - sim$Y
  
  ### Parameter recovery
  if ( lf == 'logistic' ) {
    est = lme4::glmer( cbind( Y, Z ) ~ -1 + I + A + B + AxB + 
                         ( -1 + I | Subject ),
                       family = binomial(link = "logit"),
                       data = sim )
  } else {
    est = lme4::glmer( cbind( Y, Z ) ~ -1 + I + A + B + AxB + 
                         ( -1 + I | Subject ),
                       family = binomial(link = "probit"),
                       data = sim )
  }
  
  # Plot relationship between generating 
  # and estimated parameters
  fe = fixef( est )
  re = as.numeric( ranef( est )[[1]][,1] )
  
  x11( width = 12 );
  layout( cbind( 1, 2 ) )
  plot( beta, fe, pch = 19,
        bty = 'n', xlab = 'Generating',
        ylab = 'Estimated' )
  abline( a = 0, b = 1 )
  plot( eta, re, pch = 19,
        bty = 'n', xlab = 'Generating',
        ylab = 'Estimated' )
  abline( a = 0, b = 1 )
  
}

###
### 3) Logistic/Probit regression (Random slope)
###

if ( run_code[2] ) {
  
  lf = 'probit'
  
  ### Design
  
  Ns = 100
  # Design
  design = data.frame(
    I = 1, 
    A = rep( 0:1, each = 2 ),
    B = rep( 0:1, 2 ),
    N = rep( 50, 4 )
  )
  design$AxB = design$A * design$B
  # Data frame to store simulated data
  sim = data.frame(
    Subject = rep( 1:Ns, each = nrow( design ) ), 
    I = 1,
    A = rep( design$A, Ns ),
    B = rep( design$B, Ns ),
    AxB = rep( design$AxB, Ns ),
    Y = NA,
    N = rep( design$N, Ns )
  )
  
  ### Generating parameters
  
  # Regression coefficients
  beta = runif( 4, -.5, .5 )
  # Correlation matrix
  Omega = diag( 2 )
  Omega[2,1] = Omega[1,2] = runif( 1, -.5, .5 )
  # Standard deviations
  Tau = runif( 2, .25, 1 )
  # Covariance matrix
  Sigma = diag( Tau ) %*% Omega %*% diag( Tau )
  # Subject-level deviations
  eta = MASS::mvrnorm( Ns, cbind( rep( 0, 2 ) ),
                       Sigma )
  
  # Group-level means
  X = as.matrix( sim[,c('I','A','B','AxB')] )
  sim$X_beta = X[,1] * ( beta[1] + eta[ sim$Subject, 1 ] ) + 
    X[,2] * ( beta[2] + eta[ sim$Subject, 2 ] ) + 
    X[,3] * beta[3] + X[,4] * beta[4]
  sim$Theta = link_function( sim$X_beta, type = lf )
  sim$Y = rbinom( nrow( sim ), sim$N, sim$Theta )
  sim$Z = sim$N - sim$Y
  
  ### Parameter recovery
  if ( lf == 'logistic' ) {
    est = lme4::glmer( cbind( Y, Z ) ~ -1 + I + A + B + AxB + 
                         ( -1 + I + A | Subject ),
                       family = binomial(link = "logit"),
                       data = sim )
  } else {
    est = lme4::glmer( cbind( Y, Z ) ~ -1 + I + A + B + AxB + 
                         ( -1 + I + A | Subject ),
                       family = binomial(link = "probit"),
                       data = sim )
  }
  
  # Plot relationship between generating 
  # and estimated parameters
  fe = fixef( est )
  re = ranef( est )$Subject
  
  x11( width = 12 );
  layout( cbind( 1, 2, 3 ) )
  plot( beta, fe, pch = 19,
        bty = 'n', xlab = 'Generating',
        ylab = 'Estimated' )
  abline( a = 0, b = 1 )
  plot( eta[,1], re[,1], pch = 19,
        bty = 'n', xlab = 'Generating',
        ylab = 'Estimated' )
  abline( a = 0, b = 1 )
  plot( eta[,2], re[,2], pch = 19,
        bty = 'n', xlab = 'Generating',
        ylab = 'Estimated' )
  abline( a = 0, b = 1 )
  
}

###
### 4) SDT model with random intercepts for bias/d'
###

if ( run_code[3] ) {
  
  ### Design
  
  Ns = 100
  # Design
  design = data.frame(
    Bias = 1, 
    Type = 0:1,
    Type_label = c( 'Noise', 'Signal' ), 
    N = rep( 50, 2 ),
    stringsAsFactors = FALSE
  )
  # Data frame to store simulated data
  sim = data.frame(
    Subject = rep( 1:Ns, each = nrow( design ) ), 
    Bias = 1,
    Type = rep( design$Type, Ns ),
    Type_label = rep( design$Type_label, Ns ),
    Y = NA,
    N = rep( design$N, Ns )
  )
  
  ### Generating parameters
  
  # Regression coefficients
  beta = runif( 2, c( -.25, .5 ), c( .25, 2.5 ) )
  # Correlation matrix
  Omega = diag( 2 )
  Omega[2,1] = Omega[1,2] = runif( 1, -.5, .5 )
  # Standard deviations
  Tau = runif( 2, c( .1, .25 ), c( .2, .75 ) )
  # Covariance matrix
  Sigma = diag( Tau ) %*% Omega %*% diag( Tau )
  # Subject-level deviations
  eta = MASS::mvrnorm( Ns, cbind( rep( 0, 2 ) ),
                       Sigma )
  
  # Group-level means
  X = as.matrix( sim[,c('Bias','Type')] )
  sim$X_beta = X[,1] * ( beta[1] + eta[ sim$Subject, 1 ] ) + 
    X[,2] * ( beta[2] + eta[ sim$Subject, 2 ] )
  sim$Theta = link_function( sim$X_beta, type = 'probit' )
  sim$Y = rbinom( nrow( sim ), sim$N, sim$Theta )
  sim$Z = sim$N - sim$Y
  
  # Estimation
  
  # Formula
  glmer_formula = paste(
    'cbind( Y, Z ) ~ ',
    '-1 + Bias + Type + ',
    '( -1 + Bias + Type',
    ' | Subject)', sep = '' )
  glmer_formula = as.formula( glmer_formula )
  
  if ( est_type[1] ) {
    
    est = lme4::glmer( glmer_formula, 
                       family = binomial(link = "probit"), 
                       data = sim,
                       control = 
                         glmerControl( optimizer="bobyqa", 
                                       optCtrl = list( maxfun=2e5 ) )
    )
    
  }
  
  if ( est_type[2] ) {
    
    prior_pc = rstanarm::stan_glmer( glmer_formula,
                                     family = 
                                       binomial(link = "probit"),
                                     data = sim,
                                     prior_PD = T )
  }
  
  if ( est_type[3] ) {
    
    est = rstanarm::stan_glmer( glmer_formula,
                                family = 
                                  binomial(link = "probit"),
                                data = sim )
  }
  
  if ( all( est_type[3:4] ) ) {
    
    # Posterior predictive check
    post_pc = posterior_predict(est)
    fn = function(x) by( x, list( sim$Type ), mean )[1:2]
    ppc = apply( apply( post_pc, 1, fn ), 1, 
                 quantile, prob = c( .05, .95 ) )
    
  }
  
  if ( exists( 'est' ) ) {
    
    # Plot relationship between generating 
    # and estimated parameters
    if ( est_type[1] ) {
      fe = fixef( est )
      re = ranef( est )$Subject
    }
    if ( est_type[3] ) {
      fe = rstanarm::fixef( est )
      re = rstanarm::ranef( est )[[1]]
    }
    
    x11( width = 12 );
    layout( cbind( 1, 2, 3 ) )
    plot( beta, fe, pch = 19,
          bty = 'n', xlab = 'Generating',
          ylab = 'Estimated' )
    abline( a = 0, b = 1 )
    plot( eta[,1], re[,1], pch = 19,
          bty = 'n', xlab = 'Generating',
          ylab = 'Estimated' )
    abline( a = 0, b = 1 )
    plot( eta[,2], re[,2], pch = 19,
          bty = 'n', xlab = 'Generating',
          ylab = 'Estimated' )
    abline( a = 0, b = 1 )
    
  }
  
}

###
### 5) Realistic example based on experimental design
###

if ( run_code[4] ) {
  
  ### Design
  
  Ns = 66 # Number of subjects
  
  # Design
  design = data.frame(
    Condition = rep(
      rep( c( 'Placebo', 'Drug' ), each = 4 ), 2 ), 
    Timepoints = rep( 
      rep( c( 'T1_Pre_drug', 'T2_Post_drug' ), each = 2 ), 4 ), 
    Task = rep( c( 'NBack_0', 'NBack_2' ), each = 8 ), 
    Response_type = rep( c('False_alarms','Hits'), 8 ), 
    Counts = NA,
    Trials = c( rep( c(26,22), 4 ), rep( c(53,19), 4 ) ), 
    stringsAsFactors = FALSE
  )
  
  # Numeric coding for categorical variables
  X = matrix( 0, nrow( design ), 8 )
  X[ design$Task == 'NBack_0', 1 ] = 1
  X[ design$Task == 'NBack_2', 2 ] = 1
  all_cnd = rbind(
    c( 'Drug', 'T2_Post_drug' ),
    c( 'Placebo', 'T1_Pre_drug' ),
    c( 'Placebo', 'T2_Post_drug' ) )
  all_cnd = cbind( rep( unique( design$Task ), each = 3 ),
                   rbind( all_cnd, all_cnd ) )
  for ( i in 1:nrow( all_cnd ) ) {
    sel = design$Task == all_cnd[i,1] & 
      design$Condition == all_cnd[i,2] & 
      design$Timepoints == all_cnd[i,3]
    X[ sel, i + 2 ] = 1
  }
  colnames( X ) = c(
    'Baseline_0',
    'Baseline_2',
    'Drug_0_2',
    'Placebo_0_1',
    'Placebo_0_2',
    'Drug_2_2',
    'Placebo_2_1',
    'Placebo_2_2' )
  
  Xc = X
  colnames( Xc ) = paste( 'crt', colnames( X ), sep = '_' )
  Xd = X * .5
  Xd[ design$Response_type == 'Hits' ] = 
    Xd[ design$Response_type == 'Hits' ] * -1
  colnames( Xd ) = paste( 'dp', colnames( X ), sep = '_' )
  X = cbind( Xc, Xd )
  
  # Data frame for simulated data
  sim = data.frame(
    Subject = rep( 1:Ns, each = nrow( design ) ),
    Condition = rep( design$Condition, Ns ),
    Timepoints = rep( design$Timepoints, Ns ),
    Task = rep( design$Task, Ns ),
    Response_type = rep( design$Response_type, Ns ),
    Counts = NA,
    Trials = rep( design$Trials, Ns ),
    stringsAsFactors = FALSE )
  
  # Expand design matrix
  DM = matrix( NA, nrow( sim ), ncol( Xc ) + ncol( Xd ) )
  for ( i in 1:ncol( X ) ) {
    DM[,i] = rep( X[,i], Ns )
  }
  colnames( DM ) = colnames( X )
  
  ### Generating parameters
  
  # Regression coefficients
  beta = runif( 16, 
                rep( c( -.25, .5 ), each = 8 ),
                rep( c( .25, 2.5 ), each = 8 ) )
  
  # Define number of random effects
  Nef = 4
  # Specify position of random effects
  ref_pos = c( 1:2, 1:2 + ncol(Xc) )
  
  # Generate a covariance matrix
  set.seed( 938 ) # For reproducibility
  A = matrix( runif(4^2)*2-1, ncol = Nef )
  Sigma = t(A) %*% A
  
  # Extract the correlation matrix
  Omega = cov2cor( Sigma )
  # Extract the variances
  Tau = diag( Sigma )
  
  # Scale down variances
  set.seed( 397 ) # For reproducibility
  Tau = runif( Nef, 
               c( .1, .1, .25, .25 ),
               c( .2, .2, .75, .75 ) )
  # New covariance matrix
  Sigma = diag( Tau ) %*% Omega %*% diag( Tau )
  # Subject-level deviations
  eta = MASS::mvrnorm( Ns, cbind( rep( 0, 4 ) ),
                       Sigma )
  
  # Group-level means
  sim$X_beta = numeric( nrow( sim ) )
  tmp = matrix( 0, nrow( sim ), ncol( DM ) )
  tmp[,ref_pos] = eta[ sim$Subject, ]
  for ( i in 1:ncol( DM ) ) {
    sim$X_beta = sim$X_beta + ( DM[,i] * 
                                  ( beta[i] + tmp[,i] ) )
  }
  sim$Theta = link_function( sim$X_beta, type = 'probit' )
  sim$Z = rbinom( nrow( sim ), sim$Trials, sim$Theta )
  sim$Y = sim$Trials - sim$Z
  sim$Count = sim$Y
  
  ### Parameter recovery
  
  dtbf = sim %>% 
    select( Task, Condition, Timepoints, Response_type,
            Y, Z, Subject )
  dtbf = cbind( dtbf, DM )
  
  # Estimation
  
  # Formula
  glmer_formula = paste(
    'cbind( Z, Y ) ~ -1 + ',
    paste( colnames( DM ), collapse = ' + ' ),
    ' + ( -1 + ',
    paste( colnames( DM )[ ref_pos ], collapse = ' + ' ),
    ' | Subject)', sep = '' )
  glmer_formula = as.formula( glmer_formula )
  
  if ( est_type[1] ) {
    
    est = lme4::glmer( glmer_formula, 
                       family = binomial(link = "probit"), 
                       data = dtbf,
                       control = 
                         glmerControl( optimizer="bobyqa", 
                                       optCtrl = list( maxfun=2e5 ) )
    )
    
  }
  
  if ( est_type[2] ) {
    
    prior_pc = rstanarm::stan_glmer( glmer_formula,
                                     family = 
                                       binomial(link = "probit"),
                                     data = dtbf,
                                     prior_PD = T )
    # Prior predictive check
    prior_pc = posterior_predict(prior_pc)
    fn = function(x) by( x, list( dtbf$Condition,
                                  dtbf$Timepoints,
                                  dtbf$Response_type ), mean )[1:8]
    ppc = apply( apply( prior_pc, 1, fn ), 1, 
                 quantile, prob = c( .05, .95 ) )
  }
  
  if ( est_type[3] ) {
    
    est = rstanarm::stan_glmer( glmer_formula,
                                family = 
                                  binomial(link = "probit"),
                                data = dtbf )
  }
  
  if ( all( est_type[3:4] ) ) {
    
    # Posterior predictive check
    post_pc = posterior_predict(est)
    fn = function(x) by( x, list( sim$Type ), mean )[1:2]
    ppc = apply( apply( post_pc, 1, fn ), 1, 
                 quantile, prob = c( .05, .95 ) )
    
  }
  
  if ( exists( 'est' ) ) {
    
    # Plot relationship between generating 
    # and estimated parameters
    if ( est_type[1] ) {
      fe = fixef( est )
      re = ranef( est )$Subject
    }
    if ( est_type[3] ) {
      fe = rstanarm::fixef( est )
      re = rstanarm::ranef( est )[[1]]
    }
    
    x11( width = 12 );
    layout( matrix( 1:8, 2, 4, byrow = T ) )
    
    plot( beta, fe, pch = 19,
          bty = 'n', xlab = 'Generating',
          ylab = 'Estimated' )
    abline( a = 0, b = 1 )
    title( 'Group-level parameters' )
    for ( i in 1:4 ) {
      plot( eta[,i], re[,i], pch = 19,
            bty = 'n', xlab = 'Generating',
            ylab = 'Estimated' )
      abline( a = 0, b = 1 )
    }
    
  }
  
}

setwd( orig_dir )