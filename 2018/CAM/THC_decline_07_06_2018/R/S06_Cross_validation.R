# Nested cross-validation
# Written by Kevin Potter
# email: kevin.w.potter@gmail.com
# Please email me directly if you 
# have any questions or comments
# Last updated 2018-11-20

# Table of contents
# 1) Initial setup
# 2) RNG seeds for reproducibility
# 3) Nested cross-validation
#   3.1) Fit models to training data
#     3.1.1) Model 0 - Null model
#     3.1.2) Model 1 - All predictors with horseshoe prior
#     3.1.3) Model 2 - Significant predictors only
#   3.2) Cross-validation for inner fold
#   3.3) Save results for subset

###
### 1) Initial setup
###

source( 'S01_Folder_paths.R' )

# Indicate code segment to run
run_code = c(
  T
)

# Load in useful packages

# Collection of useful functions
my_package_load( 'utilityf', github = T )

# Package for working with data frames
my_package_load( 'dplyr' )

# Package for Bayesian estimation
Sys.setenv(USE_CXX14 = 1) # For Rstan to work
my_package_load( 'brms' )

# Load in packages for exponential decay model
source( 'S03_Exponential_decay_functions.R' )

# Load in data
setwd( dat_dir )
load( 'THC_decay.RData' )
# Load in data with imputted values
load( 'Imputted_data.RData' )

###
### 2) RNG seeds for reproducibility
###

create_seeds = F
if ( create_seeds ) {
  
  # Specify number of folds
  K = 5
  
  # Number of times to run nested K-fold CV
  n_cv = 30
  
  # Number of models
  n_m = 3
  
  # Total number of seeds to generate
  n_seeds = n_cv * (K*n_m + 1)
  
  chk = rep( T, n_seeds )
  
  while( any( chk ) ) {
    
    # Specify seeds for random number generator 
    # to ensure reproducibility
    rng_seeds = data.frame(
      Iteration = 1:n_cv
    )
    tmp = matrix( round( runif(n_seeds)*100000 ),
                  n_cv, (K*n_m) + 1 )
    
    # Check if any values are the same
    for ( i in 1:n_seeds ) chk[i] = tmp[i] %in% tmp[-i]
    
  }
  
  colnames( tmp ) = c( 
    'Splits', 
    paste( 'Fold_', rep( 1:K, each = n_m ), 
           '_m', rep( 0:(n_m-1), K ), sep = '' )
  )
  
  # Save RNG seeds
  rng_seeds = cbind( rng_seeds, tmp )
  setwd( dat_dir )
  save( K, n_cv, rng_seeds, file = 'RNG_seeds.RData' )
  setwd( R_dir )
  
  # Clean up workspace
  rm( n_seeds, chk, tmp, i )
  
} else {
  # Load in previous RNG seeds
  setwd( dat_dir )
  load( 'RNG_seeds.RData' )
  setwd( R_dir )
}

###
### 3) Nested cross-validation
###

# Estimation settings
algorithm = list(
  # Iterations for warm-up
  wup = 1000,
  # Desired number of iterations
  d_iter = 10000/4,
  # Thinning
  thin = 1,
  # RNG seed set within loop
  seed = NA
)

if ( run_code[1] ) {
  
  # Loop over multiple cases of 
  # K-fold cross-validation
  
  for ( cur_cv_iter in 26:nrow( rng_seeds ) ) {
    
    # Track progress
    print( paste( ' --- Start iteration', 
                  cur_cv_iter,
                  '---' ) )
    
    # Specify iteration
    # cur_cv_iter = 1
    
    # Estimate total run time
    tick()
    total_run_time = run_time
    
    # Use function from R package 'loo' to 
    # randomly split subjects into 
    # subsets
    subj = unique( dtbf$ID )
    
    # For reproducibility
    set.seed( rng_seeds$Splits[ cur_cv_iter ] )
    K_fold = loo::kfold_split_random( 
      K, 
      N = length( unique( dtbf$ID ) ) )
    
    # Initialize output for each fold
    output = c()
    for ( k in 1:K ) output = c( output, list(NULL) )
    names( output ) = paste( 'Fold', 1:K, sep = '_' )
    
    # Loop over folds
    for ( k in 1:K ) {
      
      # Current fold
      cur_fold = paste( 'Fold', k, sep = '_' )
      print( cur_fold )
      
      # Initialize lists to store details 
      # for model fit
      coverage_prob = list(
        m0 = NULL,
        m1 = NULL,
        m2 = NULL
      )
      
      R2_values = list(
        m0 = NULL,
        m1 = NULL,
        m2 = NULL
      )
      
      # Specify training and test data
      test_sel = subj[ K_fold == k ]
      test = dtbf[ dtbf$ID %in% test_sel, ]
      train_sel = subj[ K_fold != k ]
      train = dtbf[ dtbf$ID %in% train_sel, ]
      
      # 3.1) Fit models to training data
      
      # 3.1.1) Model 0 - Null model
      
      # Set seed
      seed_type = paste( cur_fold, 'm0', sep = '_' )
      algorithm$seed = rng_seeds[ cur_cv_iter, seed_type ]
      
      m0 = estimate_ed_brms( train, list(),
                             algorithm = algorithm )
      
      # Extract and save results
      res0 = extract_brm_res( m0$model_fit )
      coverage_prob$m0 = compute_coverage_prob( m0, test )
      R2_values$m0 = compute_R2( m0, test )
      
      # 3.1.2) Model 1 - All predictors with horseshoe prior
      
      # Set seed
      seed_type = paste( cur_fold, 'm1', sep = '_' )
      algorithm$seed = rng_seeds[ cur_cv_iter, seed_type ]
      
      cvrts = list(
        SP = NULL,
        ER = NULL
      )
      m1 = estimate_ed_brms( train, cvrts,
                             algorithm = algorithm,
                             horseshoe = T )
      
      # Extract and save results
      res1 = extract_brm_res( m1$model_fit )
      coverage_prob$m1 = compute_coverage_prob( m1, test )
      R2_values$m1 = compute_R2( m1, test )
      
      # Determine significant coefficients
      val = which( 
        res1$Variable %in% c( 'eta1_Start_point', 
                              'eta1_Elimination_rate',
                              'sd(Start_point)',
                              'sd(Elimination_rate)', 
                              'cor(Start_point,Elimination_rate)',
                              'Measurement_error' ) )
      
      # Extract significant predictors
      # Sampling error around .05
      # qbinom( c( .025, .975 ), 10000, .05 )/10000
      which_sig = which( res1$X.p_value[-val] < .0458 )
      
      # 3.1.3) Model 2 - Significant predictors only
      
      if ( length( which_sig ) > 0 ) {
        
        # Determine which predictors should be included
        sig_pred = res1$Variable[-val][ which_sig ]
        sig_pred = unlist( strsplit( sig_pred, split = 'eta2_' ) )
        sig_pred = sig_pred[ sig_pred != '' ]
        
        # Fit model with significant predictors only
        cvrts = list(
          SP = NULL,
          ER = NULL
        )
        cvrts$SP = sig_pred[ grep( 'SP', sig_pred ) ]
        cvrts$ER = sig_pred[ grep( 'ER', sig_pred ) ]
        if ( length( cvrts$SP ) == 0 ) cvrts$SP = NULL
        if ( length( cvrts$ER ) == 0 ) cvrts$ER = NULL
        
        # Set seed
        seed_type = paste( cur_fold, 'm2', sep = '_' )
        algorithm$seed = rng_seeds[ cur_cv_iter, seed_type ]
        
        m2 = estimate_ed_brms( train, cvrts,
                               algorithm = algorithm )
        
        # Extract and save results
        res2 = extract_brm_res( m2$model_fit )
        coverage_prob$m2 = compute_coverage_prob( m2, test )
        R2_values$m2 = compute_R2( m2, test )
        
      }
      
      # 3.2) Cross-validation for inner fold
      
      if ( length( which_sig ) > 0 ) {
        
        # Determine the model with the best LOO-CV
        model_comp = loo::loo_model_weights( list( 
          m0 = m0$loocv, 
          m1 = m1$loocv, 
          m2 = m2$loocv ) )
        
        # Combined results
        all_res = data.frame(
          Model = c( rep( 'Model_0', nrow( res0 ) ), 
                     rep( 'Model_1', nrow( res1 ) ), 
                     rep( 'Model_2', nrow( res2 ) ) ),
          stringsAsFactors = F
        )
        all_res = cbind( all_res, rbind( res0, res1, res2 ) )
        
      } else {
        
        # Remove last model
        coverage_prob$m2 = NULL
        
        # Determine the model with the best LOO-CV
        model_comp = loo::loo_model_weights( list( 
          m0 = m0$loocv, 
          m1 = m1$loocv ) )
        
        # Combine results
        all_res = data.frame(
          Model = c( rep( 'Model_0', nrow( res0 ) ), 
                     rep( 'Model_1', nrow( res1 ) ) ),
          stringsAsFactors = F
        )
        all_res = cbind( all_res, rbind( res0, res1 ) )
        
      }
      
      # 3.3) Save results for subset
      
      output[[k]] = list(
        results = all_res,
        CV_coverage_prob = coverage_prob,
        R_squared = R2_values,
        model_comparisons = model_comp
      )
      
    }
    
    # Total run time
    tick()
    total_run_time = run_time - total_run_time
    print( total_run_time )
    
    # Create file name
    fname = paste(
      cur_cv_iter, '_',
      K, '_',
      'fold_cross_validation_results.RData',
      sep = ''
    )
    
    # Save results
    setwd( dat_dir )
    setwd( 'Posterior_estimates' )
    save( output, 
          file = fname )
    setwd( R_dir )
    
    # Clean up workspace
    rm( coverage_prob,
        cvrts,
        m0, m1, m2,
        output,
        R2_values,
        res0, res1, res2,
        test, train,
        cur_fold,
        fname, k, K_fold,
        test_sel,
        train_sel,
        val,
        which_sig
    )
    
    print( paste( ' --- End iteration', 
                  cur_cv_iter,
                  '---' ) )
    
  }
  
}

setwd( R_dir )

