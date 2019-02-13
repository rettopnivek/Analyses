# Modeling results
# Written by Kevin Potter
# email: kevin.w.potter@gmail.com
# Please email me directly if you 
# have any questions or comments
# Last updated 2019-01-07

# Table of contents
# 1) Initial setup
# 2) Define functions
#   2.1)  qp
#   2.2)  qpp
#   2.3)  ln_m
#   2.4)  ln_sd
#   2.5)  pv
#   2.6)  num2str.1
#   2.7)  num2str.2
#   2.8)  mvrnormArma
#   2.9)  pred_int_variance
#   2.10) beta_eta
#   2.11) blm
#   2.12) find_y_val
# 3) Estimate overall model results
# 4) Primary summary of results
# 5) Abstract
# 6) Exclusion criteria
# 7) Analytic approach
# 8) Sample description
# 9) CN-THCCOOH concentrations
# 10) Rate of CN-THCCOOH decay
# 11) PK model validation
# 12) Verifying abstinence
#   12.1) Determine priors for subsequent bayesian regressions
#   12.2) Example of formula
#   12.3) Minimum decrease for necessary confidence
#   12.4) Figure 3
#   12.5) Summary
# 13) Limitations

###
### 1) Initial setup
###

source( 'S01_Folder_paths.R' )

# Indicate code segment to run
run_code = c(
  F, # Estimate overall model results
  F, # Abstract
  F, # Exclusion criteria
  F, # Analytic approach
  F, # Sample description
  F, # CN-THCCOOH concentrations
  F, # Rate of CN-THCCOOH decay
  F, # PK model validation
  F, # Limitations
  F  # Marginal posteriors for subsequent priors
)

# Load in useful packages

# Collection of useful functions
my_package_load( 'utilityf', github = T )

# Package for working with data frames
my_package_load( 'dplyr' )

# Package for Bayesian estimation
Sys.setenv(USE_CXX14 = 1) # For Rstan to work
my_package_load( 'brms' )

# Package for creating Office documents
my_package_load( 'officer' )

# Package for creating nice tables
my_package_load( 'flextable' )

# Package for inverse gamma distribution
my_package_load( 'invgamma' )

# Load in packages for exponential decay model
source( 'S03_Exponential_decay_functions.R' )

# Load in data
setwd( dat_dir )
load( 'THC_decay.RData' )
# Load in data with imputted values
load( 'Imputted_data.RData' )

# If present, load in 
# best-fitting model 
# results
if ( 'Best_fitting_model.RData' %in% dir() ) {
  load( 'Best_fitting_model.RData' )
}

# Initialize a word document with a 
# summary of the results
SoR =  read_docx() 

###
### 2) Define functions
###

# 2.1)
qp = function( lst ) {
  # Purpose:
  # A convenience function for pasting together 
  # disparate elements.
  # Arguments:
  # lst - A list of elements that 
  #       can be converted to a character 
  #       string
  # Returns:
  # A single character string.
  
  out = paste( unlist( lst ), collapse = '' )
  
  return( out )
}

# 2.2) 
qpp = function( vec ) {
  # Purpose:
  # A convenience function for pasting together 
  # multiple character strings
  # Arguments:
  # vec - A vector of character strings
  # Returns:
  # A single character string.
  
  out = paste( vec, collapse = '' )
  
  return( out )
}

# 2.3) 
ln_m = function( mu, sigma ) {
  # Purpose:
  # Computes the mean for a log-normal distribution.
  # Arguments:
  # mu    - Mean parameter for log-normal distribution
  # sigma - Standard deviation parameter for log-normal 
  #         distribution
  # Returns:
  # The mean for a log-normal distribution.
  
  s2 = pow( sigma, 2 )
  out = exp( mu + s2/2 )
  return( out )
}

# 2.4) 
ln_sd = function( mu, sigma ) {
  # Purpose:
  # Computes the standard deviation for a 
  # log-normal distribution.
  # Arguments:
  # mu    - Mean parameter for log-normal distribution
  # sigma - Standard deviation parameter for log-normal 
  #         distribution
  # Returns:
  # The standard deviation for a log-normal distribution.
  
  s2 = pow( sigma, 2 )
  out = ( exp( s2 ) - 1 ) * exp( 2*mu + s2 )
  return( sqrt( out ) )
}

# 2.5) 
pv = function( x ) {
  # Purpose:
  # Function to compute posterior p-values from 
  # a MCMC sample.
  # Arguments:
  # x - A vector of values
  # Returns:
  # A posterior p-value.
  
  if ( mean(x) > 0 ) out = sum( x < 0 )/length(x)
  if ( mean(x) <= 0 ) out = sum( x > 0 )/length(x)
  
  return( out )
}

# 2.6) 
num2str.1 = function( val ) {
  # Purpose:
  # Converts a single digit into its associated 
  # character string.
  # Arguments:
  # val - A single digit number
  # Returns:
  # A character string.
  
  nmb = as.character( val )
  
  if ( nmb == '0' ) out = 'zero'
  if ( nmb == '1' ) out = 'one'
  if ( nmb == '2' ) out = 'two'
  if ( nmb == '3' ) out = 'three'
  if ( nmb == '4' ) out = 'four'
  if ( nmb == '5' ) out = 'five'
  if ( nmb == '6' ) out = 'six'
  if ( nmb == '7' ) out = 'seven'
  if ( nmb == '8' ) out = 'eight'
  if ( nmb == '9' ) out = 'nine'
  
  return( out )
}

# 2.7) 
num2str.2 = function( val ) {
  # Purpose:
  # Converts a single or double digit number into 
  # its associated character string.
  # Arguments:
  # val - A single or double digit number
  # Returns:
  # A character string.
  
  val = as.character( val )
  if ( nchar( val ) == 1 ) {
    out = num2str.1( val )
  }
  
  if ( nchar( val ) == 2 ) {
    
    nmb = strsplit( val, split = "" )[[1]]
    
    digit_2 = num2str.1( nmb[2] )
    
    if ( nmb[1] == '1' ) {
      if ( nmb[2] == '0' ) out = 'ten'
      if ( nmb[2] == '1' ) out = 'eleven'
      if ( nmb[2] == '2' ) out = 'twelve'
      if ( nmb[2] == '3' ) out = 'thirteen'
      if ( nmb[2] == '4' ) out = 'fourteen'
      if ( nmb[2] == '5' ) out = 'fifteen'
      if ( nmb[2] == '6' ) out = 'sixteen'
      if ( nmb[2] == '7' ) out = 'seventeen'
      if ( nmb[2] == '8' ) out = 'eighteen'
      if ( nmb[2] == '9' ) out = 'nineteen'
    } else {
      
      vals = c(
        'twenty',
        'thirty',
        'forty',
        'fifty',
        'sixty',
        'seventy',
        'eighty',
        'ninety'
      )
      for ( i in 2:9 ) {
        if ( nmb[1] == as.character(i) ) {
          if ( nmb[2] == '0' ) {
            out = vals[i-1]
          } else {
            out = paste( vals[i-1], digit_2, sep = '-' )
          }
        }
      }
      
    }
    
  }
  
  # Capitalize first letter
  tmp = strsplit( out, split = "" )[[1]]
  tmp[1] = LETTERS[ letters %in% tmp[1] ]
  out = paste( tmp, collapse = "" )
  
  return( out )
}

# 2.8)
cpp_code = 
  "arma::mat mvrnormArma(int n, arma::vec mu, arma::mat sigma) {
     int ncols = sigma.n_cols;
     arma::mat Y = arma::randn(n, ncols);
     return arma::repmat(mu, 1, n).t() + Y * arma::chol(sigma);
   }"
Rcpp::cppFunction(
  # Purpose:
  # Generates values from a multivariate normal distribution.
  # Arguments:
  # n     - Number of samples to generate
  # mu    - Mean vector
  # sigma - Covariance matrix
  # Returns:
  # The specified n draws from a multivariate normal 
  # distribution with mean vector mu and covariance 
  # matrix sigma.
  code = cpp_code,
  depends = "RcppArmadillo"
)

# 2.9)
cpp_code2 = 
"arma::mat pred_int_variance( arma::mat pfsfm, 
                              arma::vec mds, 
                              double time ) {
  
  // Variable declaration
  int S = pfsfm.n_rows; // Number of posterior samples
  
  arma::vec mu(2); // Mean vector for random effects
  mu(0) = 0; mu(1) = 0;
  // Covariance matrix with fixed population parameters
  arma::mat Sigma_fixed( 2, 2, arma::fill::zeros );
  // Covariance matrix with uncertain population parameters
  arma::mat Sigma( 2, 2, arma::fill::zeros );
  
  // Starting point estimate
  double SP = 0.0;
  // Elimination rate estimate
  double ER = 0.0;
  
  // Initialize output
  arma::mat out( 6, S, arma::fill::zeros );
  
  // Fill out covariance matrix with fixed population parameters
  Sigma_fixed(0,0) = mds(2)*mds(2);
  Sigma_fixed(1,1) = mds(3)*mds(3);
  Sigma_fixed(0,1) = mds(4)*mds(3)*mds(2);
  Sigma_fixed(1,0) = mds(4)*mds(3)*mds(2);
  
  // Loop over posterior samples
  for ( int s = 0; s < S; s++ ) {
    
    // Simulate random effects (Fixed)
    arma::mat fix_eta = mvrnormArma( 1, mu, Sigma_fixed );
    
    // Fill out covariance matrix with uncertain population parameters
    Sigma(0,0) = pfsfm(s,2)*pfsfm(s,2);
    Sigma(1,1) = pfsfm(s,3)*pfsfm(s,3);
    Sigma(0,1) = pfsfm(s,4)*pfsfm(s,3)*pfsfm(s,2);
    Sigma(1,0) = pfsfm(s,4)*pfsfm(s,3)*pfsfm(s,2);
    
    // Simulate random effects (Uncertain)
    arma::mat cur_eta = mvrnormArma( 1, mu, Sigma );
    
    // Compute new log(CN-THCCOOH) using fixed population parameters
    SP = mds(0) + fix_eta(0);
    ER = (-time)*( mds(1) + fix_eta(1) );
    out(0,s) = SP + ER + R::rnorm( 0, mds(5) );
    // Variance based on subject-level effects
    SP = mds(0) + fix_eta(0);
    ER = (-time)*( mds(1) + fix_eta(1) );
    out(1,s) = SP + ER;
    // Variance based on residual
    SP = mds(0);
    ER = (-time)*( mds(1) );
    out(2,s) = SP + ER + R::rnorm( 0, mds(5) );
    
    // Compute new log(CN-THCCOOH) using uncertain population parameters
    SP = pfsfm(s,0) + cur_eta(0);
    ER = (-time)*( pfsfm(s,1) + cur_eta(1) );
    out(3,s) = SP + ER + R::rnorm( 0, pfsfm(s,5) );
    // Variance based on subject-level effects
    SP = pfsfm(s,0) + cur_eta(0);
    ER = (-time)*( pfsfm(s,1) + cur_eta(1) );
    out(4,s) = SP + ER;
    // Variance based on residual
    SP = pfsfm(s,0);
    ER = (-time)*( pfsfm(s,1) );
    out(5,s) = SP + ER + R::rnorm( 0, pfsfm(s,5) );
  }
  
  return out;
}"
cppFunction( cpp_code2,
             # Purpose:
             # ...
             # Arguments:
             # pfsfm
             # mds
             # time
             # Returns:
             # ...
             depends = "RcppArmadillo",
             includes = cpp_code
)

# 2.10)
cpp_code2 = 
"arma::mat beta_eta( arma::mat beta_post, 
                    arma::mat Sigma_post ) {
  
  // Number of posterior samples
  int S = beta_post.n_rows;
  // Number of coefficients
  int Nc = beta_post.n_cols;
  
  arma::vec mu(2); // Mean vector for random effects
  mu(0) = 0; mu(1) = 0;
  // Covariance matrix with uncertain population parameters
  arma::mat Sigma( 2, 2, arma::fill::zeros );
  
  // Initialize output
  arma::mat out( S, Nc, arma::fill::zeros );
  
  // Loop over posterior samples
  for ( int s = 0; s < S; s++ ) {
    
    // Fill out covariance matrix with uncertain population parameters
    Sigma(0,0) = Sigma_post(s,0)*Sigma_post(s,0);
    Sigma(1,1) = Sigma_post(s,1)*Sigma_post(s,1);
    Sigma(0,1) = Sigma_post(s,2)*Sigma_post(s,1)*Sigma_post(s,0);
    Sigma(1,0) = Sigma_post(s,2)*Sigma_post(s,1)*Sigma_post(s,0);
    
    // Simulate random effects
    arma::mat cur_eta = mvrnormArma( 1, mu, Sigma );
    
    // Loop over columns
    for ( int nc = 0; nc < Nc; nc++ ) {
      if ( nc < 2 ) {
        out(s,nc) = beta_post(s,nc) + cur_eta(nc);
      } else {
        out(s,nc) = beta_post(s,nc);
      }
    }
    
  }
  
  return out;
}"
cppFunction( cpp_code2,
             # Purpose:
             # ...
             # Arguments:
             # pfsfm
             # mds
             # time
             # Returns:
             # ...
             depends = "RcppArmadillo",
             includes = cpp_code
)

# Clean up workspace
rm( cpp_code, cpp_code2 )

# 2.11)
blm = function( y, X, mu_0, Lambda_0, a_0, b_0 ) {
  # Purpose:
  # Fits a Bayesian linear regression model to a 
  # set of observations given a matrix of predictors 
  # and priors for the regression coefficients and 
  # residual standard deviation.
  # Arguments:
  # y        - A column vector of observed values to fit
  # X        - A design matrix with the predictors
  # mu_0     - A column vector with the means for the 
  #            conditional normal prior on the coefficients
  # Lambda_0 - A precision matrix (inverse of the covariance 
  #            matrix for the coefficients)
  # a_0      - The shape parameter for the inverse gamma 
  #            prior on the residual standard deviation
  #            (equivalent to the shape parameter for 
  #            the gamma prior on the precision)
  # b_0      - The scale parameter for the inverse gamma 
  #            prior on the residual standard deviation 
  #            (equivalent to  1 over the rate or inverse scale 
  #            parameter for the gamma prior on precision)
  # Returns:
  # 
  
  # Ensure dependent variable is a column vector
  y = as.matrix( y )
  # Ensure predictors are in matrix form
  if ( is.null( dim( X ) ) ) {
    X = rbind( X )
  }
  
  # Sample size
  n = nrow( y )
  
  # Check that inputs are correctly specified
  issues = c(
    paste( 'Number of prior means for coefficients', 
           'must match number of predictors' ),
    paste( 'Prior precision matrix must be',
           'positive-definite square matrix' ),
    paste( 'Number of predictors must',
           'match dimensions of prior precision matrix' ),
    paste( 'Shape and scale parameters must be positive' )
  )
  
  checks = c(
    nrow( mu_0 ) == ncol( X ),
    nrow( Lambda_0 ) == ncol( Lambda_0 ) & 
      det( Lambda_0 ) > 0,
    nrow( Lambda_0 ) == ncol( X ),
    a_0 > 0 & b_0 > 0
  )
  
  # If not, return error message
  if ( !all( checks ) ) {
    
    val = sum( !checks )
    if ( val > 1 ) {
      string = paste(
        issues[ !check ][ -val ],
        collapse = ', ' )
      string = paste( string, ', and ', 
                      issues[ !check ][val], sep = '' )
    } else {
      string = issues[ !check ]
    }
    
    stop( string, call. = F )
  }
  
  # Update prior parameters based on current data
  
  # Posterior precision matrix for coefficients
  tXX = t( X ) %*% X
  Lambda_n = tXX + Lambda_0
  # New covariance matrix
  Sigma_n = solve( Lambda_n )
  
  # Posterior means for coefficients
  mu_n = solve( Lambda_n ) %*% ( Lambda_0 %*% mu_0 + t(X) %*% y )
  
  # Posterior shape parameter
  a_n = a_0 + n/2
  # Posterior scale parameter
  b_n = ( t(y) %*% y ) + ( t(mu_0) %*% Lambda_0 %*% mu_0 )
  b_n = b_n - ( t( mu_n ) %*% Lambda_n %*% mu_n )
  b_n = b_0 + .5*b_n
  
  # Create output
  out = list(
    mu_n = mu_n,
    Lambda_n = Lambda_n,
    Sigma_n = Sigma_n,
    a_n = a_n,
    b_n = b_n
  )
  
  return( out )
}

# 2.12) 
find_y_val = function( y_0, 
                       day_0, 
                       day_1, 
                       priors, 
                       intervals = c( .8, .9, .95, .99 ), 
                       level_of_use = 0, 
                       min_val = 0, 
                       M_level_of_use = 20.2,
                       SD_level_of_use = 8.21381, 
                       inc = .0014 ) {
  # Purpose:
  # For a given set of days and the value of CN-THCCOOH 
  # for the earlier day, determines the associated log 
  # CN-THCCOOH value at the later day corresponding with 
  # an elimination rate greater than 0 for a specified 
  # set of credible intervals.
  # Arguments:
  # y_0             - log( CN-THCCOOH ) for earlier time point
  # day_0           - Day of earlier time point
  # day_1           - Day of later time point
  # priors          - Named list of priors, where...
  #                   mu_0 = The prior means for 
  #                          the coefficients describing 
  #                          log( CN-THCCOOH ) at baseline, 
  #                          change in baseline( CN-THCCOOH ) 
  #                          due to standardized level of use,
  #                          and elimination rate
  #                   Lambda_0 = The precision matrix 
  #                              (inverse of the covariance matrix) 
  #                              for the three aforementioned 
  #                              regression coefficients
  #                   a_0 = The shape parameter for the inverse 
  #                         gamma prior on the residual variance
  #                   b_0 = The scale parameter for the inverse 
  #                         gamma prior on the residual variance
  # intervals       - Width of credible intervals
  # level_of_use    - The number of days spent using cannabis 
  #                   in the last month
  # M_level_of_use  - The mean of level of use, used to 
  #                   define z-scores
  # SD_level_of_use - The standard deviation of level of use, 
  #                   used to define z-scores
  # min_val         - The minimum log( CN-THCCOOH ) value 
  #                   to consider when carrying out 
  #                   the search routine
  # inc             - The increments of log( CN-THCCOOH ) 
  #                   to search over
  # Returns:
  # The minimum value of log( CN-THCCOOH ) needed at the later 
  # time point (day_1) such that the specified credible 
  # interval does not overlap with zero.
  
  # Create column vector with the log of the CN-THCCOOH values
  log_CN_THCCOOH = matrix( y_0, 2, 1 )
  
  # z-score for level of use
  zLoU = ( level_of_use - M_level_of_use )/( SD_level_of_use )
  
  # Create design matrix for Bayesian linear regression
  X = rbind( c( 1, zLoU, -day_0 ), c( 1, zLoU, -day_1 ) )
  
  # Define function to fit a bayesian linear regression 
  # given a log( CN-THCCOOH ) value for day_1 and 
  # compute the posterior probability that the elimination 
  # rate is below zero (i.e., no elimination)
  f = function( y ) {
    log_CN_THCCOOH[2] = y
    est = blm( log_CN_THCCOOH,
               X,
               priors$mu_0,
               priors$Lambda_0,
               priors$a_0,
               priors$b_0
    )
    out = pnorm( 0,
                 est$mu_n[3],
                 sqrt( diag( est$Sigma_n ) )[3] )
    return( out )
  }
  
  # Define sequence of log( CN-THCCOOH ) values to 
  # explore
  input = seq( y_0, min_val, -inc )
  # Compute posterior p-values over inputs
  ppv = sapply( input, f )
  
  out = sapply( (1 - intervals)/2,
                function(x) input[ min( which( ppv < x ) ) ] )
  names( out ) = paste(
    'CI',
    intervals*100,
    sep = '_' )
  
  return( out )
}

###
### 3) Estimate overall model results
###

if ( run_code[1] ) {
  
  # Best-fitting model
  # On average, level of MJ use significantly 
  # predicted start point
  cvrts = list(
    SP = 'zLevel_of_MJ_use.SP',
    ER = NULL
  )
  
  # Estimation settings
  algorithm = list(
    wup = 1000,
    d_iter = 10000/4,
    thin = 1,
    seed = 2991
  )
  
  # Estimate model
  bfm = estimate_ed_brms( dtbf, cvrts,
                          algorithm = algorithm )
  
  # Model with all predictors
  
  # Estimation settings
  algorithm_fm = list(
    wup = 1000,
    d_iter = 10000/4,
    thin = 1,
    seed = 9348
  )
  
  # Use horseshoe prior
  cvrts_fm = list(
    SP = NULL,
    ER = NULL
  )
  fm = estimate_ed_brms( dtbf, cvrts_fm,
                         algorithm = algorithm_fm,
                         horseshoe = T )
  
  # Save posterior estimates
  setwd( dat_dir )
  save( bfm, fm,
        file = 'Best_fitting_model.RData' )
  setwd( R_dir )
  
}

###
### 4) Primary summary of results
###

# Number of subjects
Ns = length( unique( dtbf$ID ) )

# Extract posterior estimates for population-level 
# parameters and subject-level variability
post = as.matrix( bfm$model_fit )

# Extract columns for subject-level posteriors
pn = colnames( post )
sp_subj = grep( 'Start_point', pn )[ -(1:3) ]
er_subj = grep( 'Elimination_rate', pn )[ -(1:3) ]

# Subject-level elimination rates
er = colMeans( post[,'b_Elimination_rate'] + post[,er_subj] )
# Subject-level start points
sp = colMeans( post[,'b_Start_point'] + post[,sp_subj] )

# Initialize list for result summaries
all_summary_results = list()

# 4.1) Half-life

# Half-life = log(2) / Elimination_rate
all_summary_results$Half_life = list(
  # Average half-life
  Mean = round( mean( log(2)/post[,'b_Elimination_rate'] ), 1 ),
  # Standard deviations based on subject-level half-lives
  SD = round( sd( log(2)/er ), 1 ),
  # 95% confidence interval
  CI_95 = round( quantile( log(2)/log(2)/post[,'b_Elimination_rate'],
                           c( .025, .975 ) ), 1 ),
  # Range
  Range = round( range( log(2)/er ) )
)

# 4.2) Window of detection

# Limit of quantitation = 5 ng/mL
# log(5) = log( a ) - b * x
# Window of detection = ( log(5) - log( a ) )/-b
all_summary_results$Window_of_detection = list(
  Mean = round( mean( ( log(5) - post[,'b_Start_point'] ) / 
                        -post[,'b_Elimination_rate'] ) ),
  CI_95 = round( quantile( ( log(5) - post[,'b_Start_point'] ) / 
                             -post[,'b_Elimination_rate'], 
                           c( .025, .975 ) ) ),
  Range = round( range( ( log(5) - sp )/-er ) )
)

# Prediction interval for window of detection
mx_days = 160
nd = dtbf[ 1:(mx_days+1), ]
nd$ID = '999_NEW'
# Set values to averages
nd$zLevel_of_MJ_use.SP = 0
nd$Elimination_rate = -(0:mx_days)
f = function( x ) {
  x < log(5)
}
ex = predict( bfm$model_fit,
              newdata = nd,
              transform = f,
              allow_new_levels = T,
              summary = F )
pi = apply( ex, 1, function(x) min( which(x) )+1 )

all_summary_results$Window_of_detection$PI_95 = 
  round( quantile( pi, c( .025, .975 ) ) )

# Clean up workspace
rm( nd, f, ex, pi )

# 4.3) Starting level of CN-THCCOOH 

all_summary_results$Starting_level = list(
  Mean = round( mean( ln_m( post[,'b_Start_point'],
                            post[,'sd_ID__Start_point'] ) ) ),
  SD = round( mean( ln_sd( post[,'b_Start_point'],
                           post[,'sd_ID__Start_point'] ) ) ),
  Range = 
    round( 
      range( apply(  post[,'b_Start_point'] + 
                       post[,sp_subj], 2, function(x) 
                         mean( ln_m( x, 
                                     post[,'sd_ID__Start_point'] ) ) ) ) )
)

# Clean up workspace
rm( post, pn, sp_subj, er_subj )

###
### 5) Abstract
###

if ( run_code[2] ) {
  
  print( 'Abstract' )
  
  # Average number of days since baseline measurement 
  # for each visit
  avg_days_measured = dtbf %>% 
    group_by( Visit_number ) %>% 
    summarise( M = round( mean( Days_from_baseline ) ) )
  avg_days_measured = 
    paste( paste( avg_days_measured$M[ -c(1,nrow(avg_days_measured)) ], 
                  collapse = ', ' ),
           'and',
           avg_days_measured$M[ nrow( avg_days_measured ) ],
           collapse = '' )
  
  # Average number of days spent abstinent
  tmp = dtbf %>% 
    group_by( ID ) %>% 
    summarize(
      Days_of_abst = max( Days_from_baseline )
    )
  days_of_abst = qp( list( 
    round( mean( tmp$Days_of_abst ), 1 ),
    ' days of abstinence (SD = ',
    round( sd( tmp$Days_of_abst, 1 ) ),
    '). '
  ) )
  
  # Detectable THCCOOH by last visit
  tmp = dtbf %>%
    group_by( ID ) %>% 
    summarize(
      Duration = max( Time ),
      Y = THCCOOH_orig[ Time == max( Time ) ],
      D = Dipstick[ Time == max( Time ) ],
      O = THCCOOH_no_CN[ Time == max( Time ) ]
    )
  # Subjects with detectable amounts of THCCOOH 
  # by final visit
  a_25 = tmp$Duration >= 25
  N_abst_25 = sum( a_25 ) # Number of subjects abstinent for 25+ days
  N_detect = sum( tmp$Y > 0 & a_25 )
  Per_detect = round( 100 * N_detect/sum( a_25 ), 1 )
  N_pos = sum( tmp$O > 15 & !is.na( tmp$O ) & tmp$D == 1 )
  Per_pos = round( 100 * N_pos/sum( a_25 ), 1 )
  
  # Model estimation results
  avg_start_lev = all_summary_results$Starting_level$Mean
  avg_half_life = all_summary_results$Half_life$Mean
  avg_win_det = all_summary_results$Window_of_detection$Mean
  win_det_rng = all_summary_results$Window_of_detection$Range
  
  # Template for abstract:
  abstract = qpp( c(
    'Background: Despite adolescents being the most ',
    'frequent users of cannabis, all information on ',
    'cannabis drug testing interpretation is based on data ',
    'from adults. The objective of this study was ',
    'to define the time course of urinary ',
    '11-nor-9-carboxy-THC (THCCOOH) excretion using a ',
    'population pharmacokinetics (PK) model among ',
    Ns,
    ' frequent, adolescent, non-treatment seeking ',
    'cannabis users during one month of ',
    'biochemically-verified cannabis abstinence. ',
    'Methods: Urine specimens were collected at ',
    'non-abstinent baseline and after, on average, ',
    avg_days_measured,
    ' days of abstinence. Per federal drug testing ',
    'guidelines, specimens were considered "positive" ',
    'for THCCOOH when analyte levels exceeded both 50 ',
    'ng/mL on a "rapid" immunoassay drug test and ',
    '15 ng/mL on a confirmatory assay via liquid ',
    'chromatography-tandem mass spectrometry (LC/MS/MS), ',
    'with a 5 ng/mL limit of quantitation (LOQ). ',
    'Results: Participants had an average of ',
    days_of_abst,
    'The average initial creatinine-adjusted ',
    'THCCOOH concentration (CN-THCCOOH) was ',
    qp( list( avg_start_lev, ' ng/mg, ' ) ),
    'and the average half-life was ',
    qp( list( avg_half_life, ' days. ' ) ),
    'The average window of urinary ',
    'CN-THCCOOH detection was ',
    qp( list( avg_win_det,  'days, ' ) ),
    'with an estimated range of ',
    qp( list( win_det_rng[1], ' to ', 
              win_det_rng[2], ' days. ' ) ),
    'At the final timepoint and among those with > 25 ',
    'days of abstinence ',
    qp( list( '(n = ', N_abst_25, '), ') ),
    qp( list( Per_detect, '% ' ) ),
    qp( list( '(n = ', N_detect, ') ' ) ), 
    'had THCCOOH concentrations > 5 ng/mL (i.e., ',
    'detectable on the confirmatory assay) and ',
    qp( list( Per_pos, '% ' ) ),
    qp( list( '(n = ', N_pos, ') ' ) ),
    'were "positive" per federal drug testing ',
    'guidelines. More frequent past month cannabis ',
    'use predicted greater starting levels of ',
    'CN-THCCOOH but not rate of elimination. ',
    'Sex, BMI, race, and years of cannabis use ',
    'did not predict starting concentrations ',
    'or rate of elimination. ',
    'Nested 5-fold cross-validation suggested ',
    'high reliability and predictive validity ',
    'for the final statistical model. ',
    'Conclusions: A population PK model can be ',
    'used to robustly estimate the rate of ',
    'elimination of cannabis metabolites among ',
    'adolescents, and findings provide guidelines ',
    'for differentiating new use from residual ',
    'drug excretion in this population. Findings ',
    'underscore the fact that, as with adults, the ',
    'presence of detectable cannabinoid metabolites ',
    'does not necessarily indicate recent use in ',
    'adolescents.'
  ) )
  
  # Add to word document
  SoR = SoR %>% 
    body_add_par(value = 'Abstract', 
                 style = "centered") %>% 
    body_add_par( ' ', style = "Normal") %>% 
    body_add_par( abstract, style = "Normal") %>% 
    body_add_par( ' ', style = "Normal")
  
  # Clean up workspace
  rm( a_25, N_detect, avg_days_measured,
      N_abst_25, Per_pos, Per_detect,
      days_of_abst, tmp, abstract,
      avg_win_det,
      avg_half_life,
      avg_start_lev )
  
}

###
### 6) Exclusion criteria
###

if ( run_code[3] ) {
  
  print( 'Study sample' )
  
  tmp = dtbf %>% 
    group_by( ID ) %>% 
    summarize(
      V = unique( Years_of_MJ_use )
    )
  
  tmp = paste(
    '(mean years used = ',
    round( mean( tmp$V ) ),
    ', SD = ',
    round( sd( tmp$V ) ),
    '), ', sep = ''
  )
  
  string = qpp( c(
    'Urine cannabinoid concentrations were ',
    'determined from ',
    Ns, 
    ' non-treatment seeking adolescents between ',
    'the ages of ',
    floor( min( dtbf$Age ) ), ' and ',
    floor( max( dtbf$Age ) ), 
    ' years, drawn from a larger Boston-based ',
    'study on cognition during cannabis ',
    'abstinence (Schuster et al., 2016). ',
    'Participants were recruited via peer ',
    'referral, community advertisements, and via ',
    'advertising at a public high school in a ',
    'northwest Boston suburb. Participants ',
    'included in this sub-study used cannabis ',
    'at least weekly ',
    tmp,
    'used within two days of the baseline visit, ',
    'and had at least two detectable values of THCCOOH. ',
    'All participants in this sub-study were ',
    'incentivized to complete 30 days of complete cannabis ',
    'abstinence as part of the parent project ',
    '(Schuster et al., 2016).'
  ) )
  
  # Clean up workspace
  rm( tmp )
  
  # Add to word document
  SoR = SoR %>% 
    body_add_par(value = 'Study sample', 
                 style = "centered") %>% 
    body_add_par( ' ', style = "Normal") %>% 
    body_add_par( string, style = "Normal") %>% 
    body_add_par( ' ', style = "Normal")
  
  print( 'Specimen collection and analysis' )
  
  # Average number of days since baseline measurement 
  # for each visit
  avg_days_measured = dtbf %>% 
    group_by( Visit_number ) %>% 
    summarise( M = round( mean( Days_from_baseline ) ) )
  avg_days_measured = 
    paste( paste( avg_days_measured$M[ -c(1,nrow(avg_days_measured)) ], 
                  collapse = ', ' ),
           'and',
           avg_days_measured$M[ nrow( avg_days_measured ) ],
           collapse = '' )
  
  string = qpp( c(
    'Urine specimens were obtained at baseline ',
    'and again after, on average, ',
    avg_days_measured,
    ' days of abstinence.'
  ) )
  
  
  # Add to word document
  SoR = SoR %>% 
    body_add_par(value = 'Specimen collection and analysis', 
                 style = "centered") %>% 
    body_add_par( ' ', style = "Normal") %>% 
    body_add_par( string, style = "Normal") %>% 
    body_add_par( ' ', style = "Normal")
  
  # Clean up workspace
  rm( avg_days_measured )
  
  print( 'Exclusion criteria' )
  
  # Load in details for missing data
  setwd( dat_dir )
  load( 'Summary_of_missing_data.RData' )
  
  tst = dtbf %>% group_by( ID ) %>% summarize( Ns = length( ID ) )
  sodi = list(
    No = nrow( dtbf ),
    Ms = round( mean( tst$Ns ), 1 ),
    Rs = range( tst$Ns )
  )
  
  string = qpp( c(
    'Among the ', Ns, ' study participants, ',
    missing_data$DC[1], 
    ' data points (across ', 
    missing_data$DC[2], 
    ' participants) were unavailable due to participants who became ',
    'non-compliant or lost to follow-up, ',
    missing_data$LV[1],
    ' data points (across ',
    missing_data$LV[2], 
    ' participants) were unavailable because samples leaked in ',
    ' transit and the quantitative assay could not be performed, ',
    missing_data$ND[1],
    ' data points (across ',
    missing_data$ND[2], 
    ' participants) were unavailable because they ',
    'exceeded the linear range of the assay and ',
    'sample dilutions were not performed, and ',
    'finally data for ',
    missing_data$NR[1],
    ' visits were not recorded (across ',
    missing_data$NR[2],
    ' participants). ',
    'This resulted in a total of ',
    sodi$No,
    ' usable urine specimens (M = ',
    sodi$Ms, 
    ' samples per participant, range: ',
    sodi$Rs[1],
    ' - ',
    sodi$Rs[2],
    '). '
  ) )
  
  # Add to word document
  SoR = SoR %>% 
    body_add_par(value = 'Data exclusions', 
                 style = "centered") %>% 
    body_add_par( ' ', style = "Normal") %>% 
    body_add_par( string, style = "Normal") %>% 
    body_add_par( ' ', style = "Normal")
  
  # Clean up workspace
  rm( sodi, string, tst )
  
}

###
### 7) Analytic approach
###

if ( run_code[4] ) {
  
  print( 'Analytic approach' )
  
  string = qpp( c(
    'A nested 5-fold cross-validation approach was ',
    'conducted to avoid over-fitting and to minimize ',
    'the risk of false positives (Cawley & Talbot, 2010). ',
    'The cross-validation procedure first consisted of ',
    'fitting three models to a subset of the data (',
    4*Ns/5, ' randomly selected subjects, the training ',
    'sample): 1) a null model with no predictors, 2) a ',
    'full model with coefficients for all predictors ',
    'weighted to reduce weak effects to zero while ',
    'preserving strong effects (e.g., Carvalho, Polson, ',
    '& Scott, 2010), and 3) a final model that only ',
    'included significant predictors. Leave-one-out ',
    'cross-validation was applied to select the most ',
    'predictive model, and the selected model was then ',
    'validated against the remaining ',
    Ns/5, ' subjects (the test sample). This approach was ',
    'repeated 5 times such that every participant was part ',
    'of the test sample one time. ',
    'Finally, because cross-validation results can ',
    'vary based on how the data are partitioned, we ',
    'repeated this two-step validation procedure 30 times, ',
    'splitting the data into 5 new partitions each time. ',
    'This three-step approach allowed determination of the ',
    'set of predictors that robustly influenced starting ',
    'levels and elimination rate for CN-THCCOOH. Results ',
    'are based on the most reliable and predictive model ',
    'as applied to the full data set following the ',
    'cross-validation procedure. '
  ) )
  
  # Add to word document
  SoR = SoR %>% 
    body_add_par(value = 'Analytic approach', 
                 style = "centered") %>% 
    body_add_par( ' ', style = "Normal") %>% 
    body_add_par( string, style = "Normal") %>% 
    body_add_par( ' ', style = "Normal")
  
  # Extract versions for packages
  tst = sessionInfo()
  
  string = qpp( c(
    'Analyses were conducted in R (version ',
    paste( tst$R.version$major, tst$R$minor, sep = '.' ), 
    '; R Core Team, 2017). ',
    'Data were prepared using the "dplyr" package (version ',
    '0.7.1',
    '; Wickham, Francois, Henry, & Muller, 2017), ',
    'and the mixed effects PK model was fit using the R ',
    'package "brms" (version 2.3.1; Bürkner, 2017). ',
    'The R code used for the full set of analyses can be ',
    'found here: ',
    paste( 'https://github.com/rettopnivek/Analyses/',
           'tree/master/2018/CAM/THC_decline_07_06_2018. ',
           sep = '' )
  ) )
  
  # Add to word document
  SoR = SoR %>% 
    body_add_par( string, style = "Normal") %>% 
    body_add_par( ' ', style = "Normal")
  
  # Clean up workspace
  rm( string, tst )
  
}

###
### 8) Sample description
###

if ( run_code[5] ) {
  
  print( 'Sample description' )
  
  # Initialize list with demographics information
  dmg = list()
  
  # Sex
  tmp = dtbf %>% 
    group_by( ID ) %>% 
    summarize(
      S = unique( Sex )
    )
  dmg$N_males = sum( tmp$S == -1 )
  dmg$N_females = sum( tmp$S != -1 )
  
  # Age
  tmp = dtbf %>% 
    summarize(
      M = round( mean( Age ), 1 ),
      SD = round( sd( Age ), 1 ),
      Min = round( min( Age ) ),
      Max = round( max( Age ) )
    )
  dmg$Age_m = tmp$M
  dmg$Age_sd = tmp$SD
  dmg$Age_r = c( tmp$Min, tmp$Max )
  
  # Number of days in last month spent using MJ 
  tmp = raw_dat %>% 
    filter( !is.na( tlfb_mj_14b ) & 
              id %in% dtbf$ID ) %>% 
    select( id, tlfb_mj_14b )
  dmg$Days_mj_m = round( mean( tmp$tlfb_mj_14b ) )
  dmg$Days_mj_sd = round( sd( tmp$tlfb_mj_14b ) )
  dmg$Days_mj_r = range( tmp$tlfb_mj_14b )
  
  # Number of times spent smoking per day
  tmp = raw_dat %>% 
    filter( !is.na( tlfb_mj_6 ) & 
              id %in% dtbf$ID ) %>% 
    select( id, tlfb_mj_6 )
  dmg$Smoke_m = round( mean( tmp$tlfb_mj_6 ), 1 )
  dmg$Smoke_sd = round( sd( tmp$tlfb_mj_6 ), 1 )
  dmg$Smoke_r = range( tmp$tlfb_mj_6 )
  
  # Age first used MJ
  tmp = dtbf %>% 
    group_by( ID ) %>% 
    summarize(
      Y = unique( Age_first_used_MJ )
    )
  dmg$AFU_m = round( mean( tmp$Y ) )
  dmg$AFU_sd = round( sd( tmp$Y ) )
  dmg$AFU_r = round( range( tmp$Y ) )
  
  # Years spent using MJ
  tmp = dtbf %>% 
    group_by( ID ) %>% 
    summarize(
      Y = unique( Years_of_MJ_use )
    )
  dmg$YU_m = round( mean( tmp$Y ) )
  dmg$YU_sd = round( sd( tmp$Y ) )
  dmg$YU_r = round( range( tmp$Y ), 1 )
  
  # Recency of MJ use
  tmp = dtbf %>% 
    group_by( Recency_of_MJ_use ) %>% 
    summarize( 
      Ns = length( unique( ID ) )
    )
  dmg$U_0d = tmp$Ns[1]
  dmg$U_1d = tmp$Ns[2]
  dmg$U_2d = tmp$Ns[3]
  
  # CUDIT summed scores
  tmp = dtbf %>% 
    group_by( ID ) %>% 
    summarize( 
      CDT = unique( CUDIT.SS )
    )
  dmg$Hz_n = sum( tmp$CDT >= 8 )
  dmg$Hz = round( 100*sum( tmp$CDT >= 8 )/nrow( tmp ) )
  dmg$CU = round( 100*sum( tmp$CDT >= 12 )/nrow(tmp) )
  
  # Duration of cannabis abstinence
  tmp = dtbf %>% 
    group_by( ID ) %>% 
    summarize(
      Y = max( Time )
    )
  dmg$CA_m = round( mean( tmp$Y ), 1 )
  dmg$CA_sd = round( sd( tmp$Y ), 1 )
  dmg$CA_r = round( range( tmp$Y ) )
  dmg$CA_25 = sum( tmp$Y >= 25 )
  
  # 
  tmp = dtbf %>% 
    filter( Time >= 25 ) %>% 
    group_by( ID ) %>% 
    summarize(
      Y = THCCOOH_orig[ which.max( Time ) ] > 0,
      D = Positive_test_fed[ which.max( Time ) ] == 1
    )
  dmg$CA_ndet = sum( tmp$Y )
  dmg$CA_det = round( 100*sum( tmp$Y )/nrow( tmp ), 1 )
  dmg$CA_nfed = sum( tmp$D )
  dmg$CA_fed = round( 100*sum( tmp$D )/nrow( tmp ), 1 )
  
  # Paragraph 1
  
  string = qpp( c(
    'Study participants (N = ', Ns, '; ',
    dmg$N_males, ' males and ', dmg$N_females, ' females) ',
    'had a mean age of ', dmg$Age_m, ' years (SD = ', dmg$Age_sd, 
    '; range: ', dmg$Age_r[1], ' - ', dmg$Age_r[2], ' years) ',
    'and used cannabis on average ', dmg$Days_mj_m, 
    ' days in the last month (SD = ', dmg$Days_mj_sd, 
    '; range: ', dmg$Days_mj_r[1], ' - ', dmg$Days_mj_r[2], 
    ' days) and ', dmg$Smoke_m, ' times per smoking day ',
    '(SD = ', dmg$Smoke_sd, 
    '; range: ', dmg$Smoke_r[1], ' - ', dmg$Smoke_r[2], 
    ' times). The average age of first cannabis ',
    'exposure was ', dmg$AFU_m, 
    ' years (SD = ', dmg$AFU_sd,
    '; range: ', dmg$AFU_r[1], ' - ', dmg$AFU_r[2], 
    ' years), and the mean number of years of use was ',
    dmg$YU_m, 
    ' years (SD = ', dmg$YU_sd, 
    '; range: ', dmg$YU_r[1], ' - ', dmg$YU_r[2], ' years). ',
    num2str.2(dmg$U_0d), ' ', 
    'participants used cannabis on the first day of the study, ',
    dmg$U_1d, ' participants last used one day prior to ',
    'baseline, and ', dmg$U_2d, 
    ' participants last used two days prior to baseline. ',
    'Approximately ', dmg$Hz, '% (n = ', dmg$Hz_n, ') ',
    'used cannabis at potentially hazardous levels and ',
    dmg$CU, '% met criteria for possible cannabis ',
    'use disorder based on responses on the CUDIT-R. ',
    'See Table 1 for additional summary ',
    'participant characteristics. '
  ))
  
  # Add to word document
  SoR = SoR %>% 
    body_add_par(value = 'Sample description', 
                 style = "centered") %>% 
    body_add_par( ' ', style = "Normal") %>% 
    body_add_par( string, style = "Normal") %>% 
    body_add_par( ' ', style = "Normal")
  
  # Paragraph 2
  
  string = qpp( c( 
    'Table 2 summarizes CN-THCCOOH concentrations ',
    'for all study participants (see supplementary ',
    'appendix 2 for CN-THCCOOH concentrations by ',
    'study visit). Large inter-subject variability, ',
    'both in terms of starting CN-THCCOOH concentrations ',
    'and decay rates was observed and is illustrated ',
    'in Figure 1. Supplementary appendix 1 presents a ',
    'visualization of the data imputation results for the ',
    sum( dtbf$THCCOOH_orig == 0 ),
    ' samples (',
    round(100*sum( dtbf$THCCOOH_orig == 0 )/nrow(dtbf)),
    '% of the total collected) that fell below the LOQ for THCCOOH.'
  ))
  
  # Add to word document
  SoR = SoR %>% 
    body_add_par( string, style = "Normal") %>% 
    body_add_par( ' ', style = "Normal")
  
  # Paragraph 3
  
  string = qpp( c(
    'Rates of continuous abstinence, with 95% certainty ',
    '(Schwilke et al., 2011), were high in this sample. ',
    'Participants had an average of ',
    dmg$CA_m, ' days of cannabis abstinence ',
    '(SD = ', dmg$CA_sd, '; ',
    'range: ', dmg$CA_r[1], ' - ', dmg$CA_r[2], ' days), ',
    'tabulated from the estimated day of last cannabis ',
    'use (Table 2). At the final measurement timepoint ',
    'and among those who had at least 25 days of cannabis ',
    'abstinence (n = ', dmg$CA_25, '), ',
    dmg$CA_det, '% (n = ', dmg$CA_ndet, ') ',
    'still had detectable THCCOOH concentrations ',
    '(i.e., > 5 ng/mL), and ',
    dmg$CA_fed, '% (n = ', dmg$CA_nfed, ') ',
    'were "positive" by federal drug testing. '
  ))
  
  # Add to word document
  SoR = SoR %>% 
    body_add_par( string, style = "Normal") %>% 
    body_add_par( ' ', style = "Normal")
  
  # Clean up workspace
  rm( tmp, dmg, string )
}

###
### 9) CN-THCCOOH concentrations
###

if ( run_code[6] ) {
  
  print( 'CN-THCCOOH concentrations' )
  
  # Extract posterior samples from 
  # fit of full model example
  post = as.matrix( fm$model_fit )
  sel = grepl( 'b_eta', colnames( post ) )
  
  # Initialize details about output
  rslt = list()
  
  # Starting point
  rslt$SP_m = all_summary_results$Starting_level$Mean
  rslt$SP_sd = all_summary_results$Starting_level$SD
  rslt$SP_r = all_summary_results$Starting_level$Range
  
  # Predictors
  rslt$beta = data.frame(
    Variable = colnames( post[,sel] ),
    Beta = round( colMeans( post[,sel] ), 2 ),
    p_val = round( apply( post[,sel], 2, pv ), 3 ),
    stringsAsFactors = F
  )
  rownames( rslt$beta ) = 1:nrow( rslt$beta )
  # Remove intercepts
  rslt$beta = rslt$beta[ -(1:2), ]
  
  f = function( p_val ) {
    
    if ( p_val > 0 ) {
      out = paste( 'p =', p_val )
    } else {
      out = 'p < 0.001'
    }
    
    return( out )
  }
  rslt$beta$p_string = sapply( rslt$beta$p_val, f )
  
  string = qpp( c(
    'On the day of self-reported last cannabis ',
    'exposure, the average estimated initial ',
    'CN-THCCOOH concentration based on the ',
    'PK model was ',
    rslt$SP_m, ' ng/mg ',
    '(SD = ', rslt$SP_sd, ' ng/mg; ',
    'range: ', rslt$SP_r[1], ' - ', rslt$SP_r[2], ' ng/mg). ',
    'Frequency of past month cannabis use predicted ',
    'starting levels of CN-THCCOOH, with more ',
    'frequent users having higher starting urine ',
    'CN-THCCOOH concentrations ',
    '(beta = ', rslt$beta$Beta[4], ', ', 
    rslt$beta$p_string[4], '). ',
    'In contrast, sex, BMI, race, and years of ',
    'cannabis use did not account for variability ',
    'in starting levels (i.e., ',
    'beta_sex = ', rslt$beta$Beta[1], ', ', 
    rslt$beta$p_string[1], ', ',
    'beta_BMI = ', rslt$beta$Beta[2], ', ', 
    rslt$beta$p_string[2], ', ',
    'beta_race = ', rslt$beta$Beta[5], ', ', 
    rslt$beta$p_string[5], ', ',
    'and beta_years = ', rslt$beta$Beta[3], ', ', 
    rslt$beta$p_string[3], ') .'
  ) )
  
  # Add to word document
  SoR = SoR %>% 
    body_add_par(value = 'CN-THCCOOH concentrations', 
                 style = "centered") %>% 
    body_add_par( ' ', style = "Normal") %>% 
    body_add_par( string, style = "Normal") %>% 
    body_add_par( ' ', style = "Normal")
  
  # Clean up workspace
  rm( post, sel, rslt, f, string )
}

###
### 10) Rate of CN-THCCOOH decay
###

if ( run_code[7] ) {
  
  print( 'Rate of CN-THCCOOH decay' )
  
  # Extract posterior samples from 
  # fit of full model example
  post = as.matrix( bfm$model_fit )
  
  # Initialize details about output
  rslt = list()
  
  # Function to convert p-values 
  # into a string
  f = function( p_val ) {
    
    if ( p_val > 0 ) {
      out = paste( 'p =', p_val )
    } else {
      out = 'p < 0.001'
    }
  }
  
  # Correlation between start point and 
  # elimination rate
  sel = grepl( 'cor_', colnames( post ) )
  rslt$R_m = round( mean( post[,sel] ), 2 )
  rslt$R_p = f( round( pv( post[,sel] ), 3 ) )
  
  string = qpp( c(
    'Based on the PK model, the average half-life ',
    'of CN-THCCOOH was ',
    all_summary_results$Half_life$Mean,
    ' days (SD = ', all_summary_results$Half_life$SD, 
    ' days). The average window of urinary CN-THCCOOH ',
    'detection (i.e., the duration until THCCOOH dropped ',
    'below the LOQ of 5 ng/mL) was ',
    all_summary_results$Window_of_detection$Mean,
    ' days (95% CI: ',
    all_summary_results$Window_of_detection$CI_95[1],
    ' - ',
    all_summary_results$Window_of_detection$CI_95[2],
    ' days; 95% prediction interval: ',
    all_summary_results$Window_of_detection$PI_95[1],
    ' - ',
    all_summary_results$Window_of_detection$PI_95[2],
    ' days). ',
    'Variabilility among individuals was much higher, ',
    'with an estimated range of ',
    all_summary_results$Window_of_detection$Range[1],
    ' to ',
    all_summary_results$Window_of_detection$Range[2],
    ' days based on the PK model (Figure 2). ',
    'Estimated initial CN-THCCOOH concentration was ',
    'correlated with decay rate (R = ',
    rslt$R_m, ', ', rslt$R_p, '). '
  ) )
  
  # Add to word document
  SoR = SoR %>% 
    body_add_par(value = 'Rate of CN-THCCOOH decay', 
                 style = "centered") %>% 
    body_add_par( ' ', style = "Normal") %>% 
    body_add_par( string, style = "Normal") %>% 
    body_add_par( ' ', style = "Normal")
  
  # Extract posterior samples from 
  # fit of full model example
  post = as.matrix( fm$model_fit )
  sel = grepl( 'b_eta', colnames( post ) )
  
  # Initialize details about output
  rslt = list()
  
  # Predictors
  rslt$beta = data.frame(
    Variable = colnames( post[,sel] ),
    Beta = round( colMeans( post[,sel] ), 2 ),
    p_val = round( apply( post[,sel], 2, pv ), 3 ),
    stringsAsFactors = F
  )
  rownames( rslt$beta ) = 1:nrow( rslt$beta )
  # Remove intercepts
  rslt$beta = rslt$beta[ -(1:2), ]
  # Adjust p-value for cases where p < 0.001
  rslt$beta$p_string = sapply( rslt$beta$p_val, f )
  
  string = qpp( c(
    'Sex, BMI, years of cannabis use, race (white versus ',
    'non-white), and frequency of past month cannabis use ',
    'were not associated with decay rate (i.e., ',
    'beta_sex = ', rslt$beta$Beta[6], ', ', 
    rslt$beta$p_string[6], ', ', 
    'beta_BMI = ', rslt$beta$Beta[7], ', ', 
    rslt$beta$p_string[7], ', ', 
    'beta_years = ', rslt$beta$Beta[8], ', ',
    rslt$beta$p_string[8], ', ', 
    'beta_frequency = ', rslt$beta$Beta[9], ', ',
    rslt$beta$p_string[9], ', and ',
    'beta_race = ', rslt$beta$Beta[10], ', ',
    rslt$beta$p_string[10], '. '
  ) )
  
  # Add to word document
  SoR = SoR %>% 
    body_add_par( string, style = "Normal") %>% 
    body_add_par( ' ', style = "Normal")
  
  # Clean up workspace
  rm( post, sel, rslt, f, string )
}

###
### 11) PK model validation
###

if ( run_code[8] ) {
  
  print( 'PK model validation' )
  
  setwd( dat_dir )
  setwd( 'Posterior_estimates' )
  
  all_res = data.frame(
    Sample = rep( 1:30, each = 5 ),
    Fold = rep( 1:5, 30 ),
    Model = NA,
    PPC = NA,
    R2 = NA,
    Coef = 'None'
  )
  
  for ( it in 1:30 ) {
    
    # it = 1
    
    load( dir()[ it ] )
    
    # Which model fit the best on average
    K = length( output )
    param = matrix( NA, K, 3 )
    which_coef = c()
    for ( k in 1:K ) which_coef = c( which_coef, list( NULL ) )
    names( which_coef ) = paste( 'Fold', 1:K, sep = '_' )
    
    f = function( x ) {
      paste( x, collapse = ' --- ' )
    }
    
    for ( k in 1:K ) {
      sel = which.max( as.numeric(  output[[k]]$model_comparisons ) )
      param[k,1] = sel - 1
      param[k,2] = unlist( output[[k]]$CV_coverage_prob[sel] )
      param[k,3] = unlist( output[[k]]$R_squared[sel][1] )[1]
      
      # Determine coefficients included in model 2
      sel = output[[k]]$results$Model == 'Model_2'
      vrb = output[[k]]$results$Variable[sel]
      which_coef[[k]] = vrb[ c( grep( '.SP', vrb ), grep( '.ER', vrb ) ) ]
    }
    
    all_res$Model[ all_res$Sample == it ] = param[,1]
    all_res$PPC[ all_res$Sample == it ] = param[,2]
    all_res$R2[ all_res$Sample == it ] = param[,3]
    all_res$Coef = sapply( 1:K, function(x) f( which_coef[[x]] ) )
    
  }
  
  # Percentage out of predictors where 
  # level of MJ use predicted starting levels 
  f = function(x) {
    if ( !is.null( unlist(x) ) ) {
      val =  unlist(x) == 'zLevel_of_MJ_use.SP'
      out = 100*sum(val)/length( unlist(x) )
    } else {
      out = 0
    }
  }
  all_res$Level_of_use = sapply( all_res$Coef, f )
  
  # Summarize results
  
  # Which model was preferred
  mdl = all_res %>% 
    group_by( Sample ) %>% 
    summarize( 
      M0 = sum( Model == 0 )/5, 
      M1 = sum( Model == 1 )/5, 
      M2 = sum( Model == 2 )/5 )
  mdl_res = mdl %>% 
    summarize(
      Mean_M0 = round( mean( M0 )*100 ),
      Min_M0 = round( min( M0 )*100 ),
      Max_M0 = round( max( M0 )*100 ),
      Mean_M1 = round( mean( M1 )*100 ),
      Min_M1 = round( min( M1 )*100 ),
      Max_M1 = round( max( M1 )*100 ),
      Mean_M2 = round( mean( M2 )*100 ),
      Min_M2 = round( min( M2 )*100 ),
      Max_M2 = round( max( M2 )*100 )
    )
  
  # Determine significant predictors
  sig_pred = 
    100*table( 
      unlist( strsplit( all_res$Coef, split = '---' ) ) )/150
  
  # Posterior predictive checks
  all_res$Cur_vrb = all_res$Model
  sm = all_res %>% 
    group_by( Sample ) %>% 
    summarize( 
      Mean = mean( Cur_vrb ), 
      SD = sd( Cur_vrb ), 
      Min = min( Cur_vrb ), 
      Max = max( Cur_vrb ) )
  
  # Predictor selection
  string = qpp( c(
    'As stated, the most predictive variables were ',
    'determined via nested 5-fold cross-validation, ',
    'repeated 30 times to control for variability ', 
    'inherent to the cross-validation approach. ',
    'Recall that we compared performance for the null model ',
    '(no predictors), the full model (all predictors), ',
    'and finally the reduced model (significant ',
    'predictors only). Over the 5 folds and 30 repetitions, ',
    'the null model ',
    qp( list( 'was preferred ', mdl_res$Mean_M0, '% ' ) ), 
    'of the time, while the full model ',
    qp( list( 'was preferred ', mdl_res$Mean_M1, '% ' ) ), 
    'of the time, and the reduced model ',
    qp( list( 'was preferred ', mdl_res$Mean_M2, '% ' ) ), 
    'of the time. Hence, predictive models were preferred ',
    qp( list( mdl_res$Mean_M1 + mdl_res$Mean_M2, '% ' ) ),
    'of the time. Across the models with predictors, ',
    'level of marijuana use predicted starting levels of ', 
    'CN-THCCOOH ',
    qp( list( round( mean( all_res$Level_of_use > 0 )*100 ), '% ' )),
    'of the time, and in fact was the only significant predictor ',
    qp( list( round( mean( all_res$Level_of_use == 100 )*100 ), '% ' )),
    'of the time.'
  ) )
  
  # Add to word document
  SoR = SoR %>% 
    body_add_par(value = 'PK model validation', 
                 style = "centered") %>% 
    body_add_par( ' ', style = "Normal") %>% 
    body_add_par( string, style = "Normal") %>% 
    body_add_par( ' ', style = "Normal")
  
  # Model performance
  string = qpp( c(
    'Model performance was good. Over the 5 folds and ',
    '30 repetitions, 95% predictive coverage intervals ',
    'for the best-fitting model captured the withheld ',
    qp( list( 'sample of ', length( unique( dtbf$ID ) )/5 ) ),
    qp( list( ' subjects ', round( 100*mean( all_res$PPC ) ), '% ' )),
    'of the time ',
    qp( list( '(range: ', round( 100*min( all_res$PPC ) ),
        ' - ', round( 100*max( all_res$PPC ) ), '%). ' )),
    'Furthermore, the best-fitting models on average captured ',
    qp( list( round( mean( all_res$R2 )*100 ), '% ' ) ),
    'of the variance in CN-THCCOOH ',
    qp( list( '(range: ', round( 100*min( all_res$R2 ) ),
        ' - ', round( 100*max( all_res$R2 ) ), '%). ' ))
  ) )
  
  # Add to word document
  SoR = SoR %>% 
    body_add_par( string, style = "Normal" ) %>% 
    body_add_par( ' ', style = "Normal" )
  
  # Clean up workspace
  rm( string, all_res, sig_pred, sm, mdl, mdl_res, 
      f, param, which_coef, K, k )
  
}

###
### 12) Verifying abstinence
###

if ( run_code[9] ) {
  
  print( 'Verifying abstinence' )
  
  # 12.1) Determine priors for subsequent bayesian regressions
  
  # Isolate posterior samples for the population level parameters
  post = as.matrix( bfm$model_fit )
  pfsfm = data.frame(
    b_SP = post[,"b_Start_point"],
    b_LMJU_SP = post[,"b_zLevel_of_MJ_use.SP"],
    b_ER = post[,"b_Elimination_rate"],
    sigma_SP = post[,"sd_ID__Start_point"],
    sigma_ER = post[,"sd_ID__Elimination_rate"],
    R_SP_ER = post[,"cor_ID__Start_point__Elimination_rate"],
    sigma = post[,"sigma"]
  )
  
  # Extract priors based on marginal posteriors
  
  # Prior means for beta
  mu_0 = apply( pfsfm[,grep('b_',colnames(pfsfm))], 2, mean )
  
  # Prior covariance for beta
  tst = beta_eta( as.matrix( pfsfm[,c(1,3,2)] ),
                  as.matrix( pfsfm[,4:6] ) )
  tst = tst[,c(1,3,2)]
  Sigma_0 = var( tst )
  # Prior precision matrix for beta
  Lambda_0 = solve( Sigma_0 )
  
  qh = function( x, xlb = ' ' ) {
    hist( x, 
          freq = F, 
          col = 'grey', border = 'white',
          breaks = 'FD', bty = 'l',
          xlab = xlb,
          ylab = 'Density' )
  }
  
  # Determine best-fitting gamma parameters
  lsd = function(prm,dn) {
    
    ey = dn$y
    py = dinvgamma( dn$x, shape = prm[1], scale = prm[2] )
    
    out = sum( pow( py-ey, 2 ) )
    return(out)
  }
  
  dn = density( pow( pfsfm$sigma, 2 ) )
  est = optim( c(1,1), lsd, dn = dn )
  
  b("
    x11( width = 12 );
    layout( cbind( 1, 2 ) )
    qh( pfsfm$sigma, xlb = 'Residual standard deviation' )
    qh( pow( pfsfm$sigma, 2 ), xlb = 'Residual variance' )
    lines( dn$x, dinvgamma( dn$x, shape = est$par[1],
    scale = est$par[2] ), lty = 2, col = 'red' )
    ")
  
  a_0 = est$par[1]
  b_0 = 1/est$par[2]
  
  # Save list with priors
  priors = list(
    # Prior means for coefficients
    mu_0 = mu_0,
    # Prior precision matrix for coefficients
    # (Inverse of covariance matrix)
    Lambda_0 = Lambda_0,
    # Prior covariance matrix for coefficients
    Sigma_0 = Sigma_0,
    # Prior shape parameter for residual variance
    # (Inverse gamma distribution)
    a_0 = a_0,
    # Prior scale parameter for residual variance
    # (Inverse gamma distribution)
    b_0 = b_0
  )
  
  # Clean up workspace
  rm( est, lsd, qh, tst, post, mu_0, Lambda_0, Sigma_0,
      a_0, b_0 )
  
  # 12.2) Example of formula
  # (Not run)
  if ( FALSE ) {
    
    # Extract details on level of use
    tmp = dtbf %>% 
      group_by( ID ) %>% 
      summarize(
        LoU = unique( Level_of_MJ_use )
      )
    # Compute mean and standard deviation
    # Mean = 20.2
    # SD = 8.21381
    mLoU = 20.2
    sLoU = 8.21381
    
    # Extract data for first subject as example
    ex_y = dtbf$THCCOOH_orig[ dtbf$ID == dtbf$ID[1] ]
    ex_days = dtbf$Time[ dtbf$ID == dtbf$ID[1] ]
    ex_days = ex_days[ ex_y > 0 ]
    ex_y = ex_y[ ex_y > 0 ]
    
    # Standardize level of use
    zLoU = (number_of_days_spent_using - mLoU)/sLoU
    # Design matrix
    X = cbind( 1, zLoU, -ex_days )
    colnames( X ) = c( 'Log_SL', 'zLoU', 'ER' )
    # Input
    y = cbind( log( ex_y ) )
    
    # Bayesian regression
    est = blm(
      y,
      X,
      priors$mu_0,
      priors$Lambda_0,
      priors$a_0,
      priors$b_0
    )
    
    # Confidence intervals
    cbind(
      qnorm( .025, est$mu_n, sqrt( diag( est$Sigma_n ) ) ),
      qnorm( .975, est$mu_n, sqrt( diag( est$Sigma_n ) ) )
    )
    
  }
  
  # 12.3) Minimum decrease for necessary confidence
  
  # Specify range of CN-THCCOOH values at baseline
  starting_level = seq( 25, 700, length = 100 )
  # Specify day of second time point
  day_2nd_time_point = 16
  # Compute highest value of log( CN-THCCOOH ) 
  # where probability of non-zero elimination rate 
  # is above specified cut-off (80%, 90%, 95%, and 99% 
  # credible intervals)
  tick()
  max_y = sapply( log( starting_level ), 
                function(x) find_y_val( x, 
                                        0, 
                                        day_2nd_time_point,
                                        priors ) )
  tock()
  print( run_time )
  
  # For average baseline level and level of use 
  # compute highest value of log( CN-THCCOOH ) 
  # where probability of non-zero elimination rate 
  # is above specified cut-off (80%, 90%, 95%, and 99% 
  # credible intervals)
  tick()
  days = 2:31
  max_y_days = sapply( days, 
                       function(x) find_y_val( priors$mu_0[1], 
                                               0, 
                                               x,
                                               priors ) )
  tock()
  print( run_time )
  
  # 12.4) Figure 3
  
  # x11()
  setwd( fig_dir )
  png( 'Figure_3.png', 
       width = 7, height = 7, units = 'in',
       res = 200 )
  
  clr = colortools::tetradic( 'orange' )
  
  xa = log( c( 5, 10, 25, 50, 100, 200 ) )
  ya = log( c( 25, 50, 100, 200, 400, 700 ) )
  
  xl = c( 1.5, 5.5 )
  yl = c( 3, 7 )
  par( mar = c( 4, 4, 4, 4 ) )
  blankPlot( xl, yl )
  
  # Grid lines
  b("
  horizLines( seq( yl[1], yl[2], 100 ), xl,
              col = 'grey' )
  vertLines( c( xa, xl[2] ), yl, col = 'grey' )
  ")
  
  # Average at baseline
  horizLines( priors$mu_0[1],
              xl, lty = 2, col = 'grey' )
  text( log(125), priors$mu_0[1],
        'Average over sample',
        pos = 1 )
  
  for ( i in 1:nrow( max_y ) ) {
    lines( log( starting_level ) - max_y[i,],
           log( starting_level ),
           lwd = 2, col = clr[i] )
  }
  
  lbl = paste( 'Minimum necessary decrease by day',
               day_2nd_time_point,
               '- ng/mg' )
  customAxes( xl, yl,
              pos = 1:4, 
              label = c( lbl,
                         'Starting CN-THCCOOH - ng/mg',
                         'Decrease - log(ng/mg)',
                         'Start - log(ng/mg)'
                         ) )
  axis( 2, ya, 
        round( exp( ya ) ), 
        tick = F, line = -1.5, cex.axis = 1.25 )
  axis( 1, xa,
        round( exp( xa ) ),
        tick = F, line = -1.5, cex.axis = 1.25 )
  axis( 4, round( ya, 1 ), 
        tick = F, line = -1.5, cex.axis = 1.25 )
  axis( 3, round( xa, 1 ),
        tick = F, line = -1.5, cex.axis = 1.25 )
  
  legend( 2, yl[2] + diff(yl)*.2,
          paste( 
            sapply( 
              strsplit( rownames( max_y ), split = '_' ), 
              function(x) unlist(x)[2] ), 
            '% CI', sep = '' ),
          fill = clr,
          horiz = T, 
          bty = 'n',
          xpd = T )
  
  dev.off()
  setwd( R_dir )
  
  # 12.5) Summary
  
  string = qpp( c( 
  'Like Schwilke et al. (2010), our results can be extended ',
  'to examine data from new individuals and estimate ',
  'whether their elimination rate is non-zero. ',
  'Specifically, the marginal posterior distributions from the PK model ',
  'can be converted into conjugate prior distributions for a ',
  "standard Bayesian regression model predicting an individual's ",
  'log(CN-THCCOOH) from their level of use and the days at ',
  'which specimens were collected. Close-formed solutions ',
  'exist for the posterior distribution (e.g., Walter & Augustin, ',
  '2009) of the regression model given an inverse gamma prior ',
  'on the residual variance and a conditional multivariate ',
  'normal prior on the regression coefficients. ',
  'Based on the PK model, the distribution for the residual ',
  'variance is nicely approximated by an inverse gamma with ',
  'a shape parameter of ', round( priors$a_0, 1 ), ' and ',
  'a scale parameter of ', round( priors$b_0, 1 ), '. ', 
  'In turn, the joint distribution for the three regression ',
  'coefficients (log of the baseline level of CN-THCCOOH, ',
  'standardized level of use, and elimination rate) ',
  'is nicely approximated by a multivariate normal with ',
  'mean vector ', 
  paste( '[', paste( round( priors$mu_0, 2 ),
                     collapse = ', ' ), ']', sep = '' ), ' and ',
  'precision matrix (inverse of the covariance matrix) ',
  paste( '[', paste( round( priors$Lambda_0[1,], 2 ),
                     collapse = ', ' ), '; ', sep = '' ),
  paste( paste( round( priors$Lambda_0[2,], 2 ),
                     collapse = ', ' ), '; ', sep = '' ),
  paste( paste( round( priors$Lambda_0[3,], 2 ),
                collapse = ', ' ), ']. ', sep = '' ),
  'Note the precision matrix incorporates variability ',
  'due both to individual differences and uncertainty ',
  'in population-level estimates. '
  ) )
  
  # Add to word document
  SoR = SoR %>% 
    body_add_par(value = 'Abstinence verification', 
                 style = "centered") %>% 
    body_add_par( ' ', style = "Normal") %>% 
    body_add_par( string, style = "Normal") %>% 
    body_add_par( ' ', style = "Normal")
  
  string = qpp( c(
    'Figure 3 shows an example application of the Bayesian ',
    'regression approach with the priors based on the PK ',
    'model results. The figure reports the minimum necessary ',
    'decrease needed to be observed at day 16, in log(mg/ng), ',
    'to have a 80%, 90%, ', 
    '95%, and 99% probability (conditioned on the prior and ',
    'specified data) of a non-zero elimination rate, ',
    'over a range of log(CN-THCCOOH) levels at ingestion. ',
    'For convenience, we fixed level of use to the mean. ',
    'We chose day 16 for the example because this is the day ',
    'by which the estimated population mean drops below the LOQ. ',
    'For convenience, on the standard x and y-axes, we report ',
    'the raw CN-THCCOOH concentrations corresponding to the ',
    'log units used in the figure. ',
    'At low starting concentrations of CN-THCCOOH, ',
    'the minimum necessary decrease can be quite small. ',
    'However, due to the nature of exponential decay, ',
    'higher starting concentrations of CN-THCCOOH ', 
    'require much larger decreases in order to have ',
    'high confidence of non-zero elimination. '
  ) )
  
  # Add to word document
  SoR = SoR %>% 
    body_add_par( string, style = "Normal" ) %>% 
    body_add_par( ' ', style = "Normal" )
  
  string = qpp( c(
    'Figure 3 is based on two specimens, one at day 0 and one ',
    'at day 16. The algorithm is easily extended to other days ',
    'and additional measurements. If necessary, even a single ',
    'specimen can be used. However, we emphasize that with ',
    'fewer specimens, results will be more biased towards the priors, ',
    'the population-level values, leading to more conservative ',
    'estimates. For example, consider the minimum necessary decrease ',
    'to have a 95% probability of non-zero elimination on ',
    'the third day of abstinence. A subject with an average level of ',
    'CN-THCCOOH at ingestion (101 ng/mg) and average level of use ',
    'still needs to see a decrease of ',
    round( exp( priors$mu_0[1] ) - 
             exp( find_y_val( y_0 = 4.6, day_0 = 0, 
                              day_1 = 3, priors = priors ) ) )[3],
    ' ng/mg. '
  ) )
  
  # Add to word document
  SoR = SoR %>% 
    body_add_par( string, style = "Normal" ) %>% 
    body_add_par( ' ', style = "Normal" )
  
}

###
### 13) Limitations
###

if ( run_code[10] ) {
  
  # Determine the proportion of the variance in the 
  # prediction interval that is dependent on
  # 1) Uncertainty around population parameters
  # 2) Subject-level variance
  # 3) Residual variance
  
  # Isolate posterior samples for the population level parameters
  post = as.matrix( bfm$model_fit )
  pfsfm = data.frame(
    b_SP = post[,"b_Start_point"],
    # b_LMJU_SP = post[,"b_zLevel_of_MJ_use.SP"],
    b_ER = post[,"b_Elimination_rate"],
    sigma_SP = post[,"sd_ID__Start_point"],
    sigma_ER = post[,"sd_ID__Elimination_rate"],
    R_SP_ER = post[,"cor_ID__Start_point__Elimination_rate"],
    sigma = post[,"sigma"]
  )
  
  pfsfm = as.matrix( pfsfm )
  mds = apply( pfsfm, 2, findMode )
  
  nRep = 100
  days = 0:31
  V = array( NA, dim = c(nRep,length(days),6) )
  for ( j in 1:length( days ) ) {
    print( paste( 'Day', days[j] ) )
    for ( i in 1:nRep ) {
      tst = pred_int_variance( pfsfm, mds, days[j] )
      V[i,j,] = apply( tst, 1, var )
    }
  }
  
  # Extract proportions 
  f = function(v) {
    tmp = apply( v, 2, mean )
    out = numeric(3)
    out[1] = 1 - (tmp[1]/tmp[4])
    out[2] = tmp[2]/tmp[4]
    out[3] = tmp[3]/tmp[4]
    names(out) = c( 'Sample_size', 'Subject', 'Residual' )
    return(out)
  }
  tmp = apply( V, 2, f )
  
  string = qpp( c(
    'One limitation to note is that estimates for individuals ',
    'were only based off of, at most, seven specimens. This ',
    'could result in a higher degree of uncertainty around ',
    'subject-level estimates, and is reflected, for instance, ',
    'in the wide range for the window of detection across the ',
    'prediction interval and individual subjects. With future ',
    'studies, collection of more specimens will likely reduce ',
    'this range. '
  ) )
  
  # Add to word document
  SoR = SoR %>% 
    body_add_par(value = 'Limitations', 
                 style = "centered") %>% 
    body_add_par( ' ', style = "Normal") %>% 
    body_add_par( string, style = "Normal") %>% 
    body_add_par( ' ', style = "Normal")
  
  string = qpp( c(
  'Fortunately, a benefit of a multilevel model is that ',
  'estimates for individual subjects are biased slightly ',
  'towards the population-level estimates (shrinkage), ',
  'which gives a beneficial reduction in variance. ',
  'One way to quantify this is to look at the influence on ',
  'the variance of the prediction interval from the ',
  'the uncertainty in the population-level estimates, ',
  'subject-level variance, and unexplained (residual) variance. ',
  'At the day of ingestion, about ',
  round( 100*tmp[1,1] ), '% of the ',
  'variance reflected uncertainty in the population parameters, ',
  round( 100*tmp[2,1] ), '% ',
  'reflected subject-level variance, and ',
  100 - round( 100*(tmp[1,1]+tmp[2,1]) ), '% ',
  'reflected unexplained variance. At later days, such as 16 ',
  'days following ingestion (the point at which, on average, ',
  'CN-THCCOOH would drop below detectable levels) ',
  'residual variance had less influence, most likely due to ',
  'imputted data, with ',
  round( 100*tmp[1,17] ), '% ',
  'reflecting uncertainty in population-level estimates',
  round( 100*tmp[2,17] ), '% ',
  'reflecting subject-level variance, and only ',
  100 - round( 100*(tmp[1,17]+tmp[2,17]) ), '% ',
  'reflecting residual variance. '
  ) )
  
  # Add to word document
  SoR = SoR %>% 
    body_add_par( string, style = "Normal" ) %>% 
    body_add_par( ' ', style = "Normal" )
  
  string = qpp( c(
  'Critically, little of the variance in predictions reflected ',
  'uncertainty in population-level estimates, a product of the ',
  'larger sample size of ', Ns, ' subjects. This also further ',
  'emphasizes the high degree of individual differences in ',
  'the use and metabolism of cannabis.'
  ) )
  
  # Add to word document
  SoR = SoR %>% 
    body_add_par( string, style = "Normal" ) %>% 
    body_add_par( ' ', style = "Normal" )
  
  string = qpp( c(
    'Note another limitation is that the data imputation ',
    'could lead to a slight and artificial reduction in ',
    'residual variance for later days. Fortunately, the ',
    'impact of this is minor, as ',
    'evidenced by the lack of change in the proportion ',
    'of variance due to uncertainty in population-level ',
    'estimates (', 
    round( 100*tmp[1,1] ), '% versus ',
    round( 100*tmp[1,17] ), '%), ',
    'though future studies can consider using ',
    'multiple imputation techniques instead. '
  ) )
  
  # Add to word document
  SoR = SoR %>% 
    body_add_par( string, style = "Normal") %>% 
    body_add_par( ' ', style = "Normal")
  
}

# Output Word document
setwd( proj_dir )
setwd( 'Documents' )
print( SoR, target = 'Summary_of_results.docx' )

setwd( R_dir )
