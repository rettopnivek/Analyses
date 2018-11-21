# Bayesian population-level model (Example)
# Written by Kevin Potter
# email: kevin.w.potter@gmail.com
# Please email me directly if you 
# have any questions or comments
# Last updated 2018-11-20

# Table of contents
# 1) Initial setup
# 2) Basic regression example
# 3) Horseshoe prior example

###
### 1) Initial setup
###

source( 'S01_Folder_paths.R' )

# Indicate code segments to run
run_code = c(
  F,
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
### 2) Basic regression example
###

if ( run_code[1] ) {
  
  ### Specify empirical bayes priors for start point/elimination rate
  
  ebp = all_est %>% 
    summarize(
      M_LSP = mean( Log_start_point ),
      SD_LSP = sd( Log_start_point ),
      SEM_LSP = sem( Log_start_point ),
      M_ER = mean( Elimination_rate ),
      SD_ER = sd( Elimination_rate ),
      SEM_ER = sem( Elimination_rate ),
      M_RSD = mean( Residual_SD ),
      SD_RSD = sd( Residual_SD ),
      SEM_RSD = sem( Residual_SD )
    )
  # round( ebp, 2 )
  
  # Weakly informative prior on residual error
  prior = set_prior( 'student_t( 19, .34, .23 )', 
                     class = 'sigma' )
  
  # Prior on starting point
  prior = c( prior,
             set_prior( 'normal( 4.7, 1.1 )', 
                        class = 'b', coef = 'Start_point' ) )
  
  # Prior on elimination rate
  prior = c( prior,
             set_prior( 'student_t( 15, .23, .16 )', 
                        class = 'b', coef = 'Elimination_rate' ) )
  
  # Prior on correlation between random effects
  prior = c( prior,
             set_prior( 'lkj_corr_cholesky(1)', 
                        class = 'L' ) )
  
  # Warm up period
  wup = 1000
  # Desired number of iterations
  d_iter = 10000/4
  # Thinning value
  thin = 1
  # Total samples to draw
  n_iter = d_iter * thin + wup
  
  # Model estimation
  tick()
  m0 = brm( log_THCCOOH ~ 0 + # Remove standard intercept
              # Population effects (Start point)
              Start_point + 
              # Population effects (Elimination rate)
              Elimination_rate + 
              # Subject-level effects
              ( 0 + Start_point + Elimination_rate | ID ),
            # Data and priors
            data = dtbf, 
            prior = prior, 
            # Estimation settings
            warmup = wup,
            iter = n_iter,
            chains = 4,
            cores = 4,
            thin = thin, 
            control = list( adapt_delta = 0.995,
                            max_treedepth = 15 ) )
  tock()
  print( run_time )
  
  # Modeling results
  tick()
  
  # Observed versus predicted distribution
  x11(); pp_check( m0 )
  
  x11(); plot( marginal_effects( m0, effects = 'Elimination_rate' ), 
               points = T )
  
  Ns = length( unique( dtbf$ID ) )
  for ( i in 1:ceiling( Ns/9 ) ) {
    
    x11()
    ind = 1:9 + 9 * (i-1)
    if ( max( ind ) > Ns ) ind = min(ind):Ns
    conditions = data.frame( 
      ID = unique(dtbf$ID)[ind], stringsAsFactors = F )
    rownames(conditions) = unique(dtbf$ID)[ind]
    Fit_by_subject = marginal_effects( m0, 
                                       conditions = conditions,
                                       effects = 'Elimination_rate', 
                                       re_formula = NULL, 
                                       method = "predict")
    plot( Fit_by_subject, nrow = 3, ncol = 3, points = TRUE)
    
  }
  tock()
  print( run_time )
  
}

###
### 3) Horseshoe prior example
###

if ( run_code[2] ) {
  
  # Weakly informative prior on residual error
  prior = set_prior( 'student_t( 19, .34, .23 )', 
                     class = 'sigma' )
  
  # Prior on starting point
  prior = c( prior,
             set_prior( 'normal( 4.7, 1.1 )', 
                        class = 'b', coef = 'Start_point',
                        nlpar = 'eta1' ) )
  
  # Prior on elimination rate
  prior = c( prior,
             set_prior( 'student_t( 15, .23, .16 )', 
                        class = 'b', coef = 'Elimination_rate',
                        nlpar = 'eta1' ) )
  
  # Prior on correlation between random effects
  prior = c( prior,
             set_prior( 'lkj_corr_cholesky(1)', 
                        class = 'L' ) )
  
  # Horseshoe prior on remaining coefficients
  prior = c( prior,
             prior( horseshoe(1), nlpar = 'eta2' ) )
  
  bformula = bf( log_THCCOOH ~ eta1 + eta2,
                 eta1 ~ 0 + 
                   Start_point + 
                   Elimination_rate + 
                   (0 + Start_point + Elimination_rate|ID),
                 eta2 ~ 0 + 
                   # Predictors for start-point
                   Sex.SP + 
                   Race.SP + 
                   zBMI.SP + 
                   zYears_of_MJ_use.SP + 
                   zLevel_of_MJ_use.SP + 
                   # Predictors for elimination rate
                   Sex.ER + 
                   Sex.ER + 
                   zBMI.ER + 
                   zYears_of_MJ_use.ER + 
                   zLevel_of_MJ_use.ER,
                 nl = TRUE)
  
  # Warm up period
  wup = 1000
  # Desired number of iterations
  d_iter = 10000/4
  # Thinning value
  thin = 1
  n_iter = d_iter * thin + wup
  
  m1 = brm( bformula, 
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
  
}

setwd( R_dir )

