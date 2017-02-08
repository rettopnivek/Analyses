functions {
  #include "Wald_race_stan_functions.stan"
}
data {
  int No; // Total number of observations
  int Ns; // Number of subjects
  int Ni; // Number of items
  int Nx; // Number of conditions for drift rate
  int Nk; // Number of conditions for threshold
  vector[2] Y[No]; // Matrix of response times and choice
  int Co[No]; // Position of correct response (-1 = left, 1 = right)
  int subjIndex[No];
  int itemIndex[No];
  row_vector[Nx] X_x_R[No]; // Design matrices
  row_vector[Nx] X_x_L[No]; 
  row_vector[Nk] X_k_R[No]; 
  row_vector[Nk] X_k_L[No]; 
  real<lower=0.0> min_RT[Ns]; // Smallest response time for each subject
  matrix[ Nx + Nk + 5, 2 ]Priors; // Matrix of parameters for priors
}
parameters {
  vector[Nx] beta_x; // Group-level effects for drift rate
  vector[Nk] beta_k; // Group-level effects for threshold
  vector[3] eta[Ns]; // Subject-specific effects on drift 
                     // and threshold
  cholesky_factor_corr[3] L_Omega; // Correlations for subject 
                                   // effects
  vector<lower=0.0>[3] tau; // Scale parameters for subject effects
  real zeta[Ni]; // Item-specific effects on drift
  real<lower=0.0> sigma_zeta; // Standard deviation for item effects
  real<lower=0.0,upper=1.0> theta[Ns]; // Proportion for residual
                                       // latency
  real<lower=0.0,upper=1.0> phi; // Mean for residual latency
                                 // proportion
  real<lower=0.1> lambda; // Total counts for residual latency
                          // proportion
}
transformed parameters {
  // Variable declaration
  vector[8] param[No];
  
  // Loop over observations
  for ( no in 1:No ) {
    
    // Calculate parameters for wald race model
    
    // Threshold for picking right
    param[ no, 1 ] = X_k_R[ no ] * beta_k + 
      eta[ subjIndex[ no ], 3 ];
    
    // Drift rate for picking right
    param[ no, 2 ] = 
      X_x_R[ no ] * beta_x + // Group-level effects
      eta[ subjIndex[ no ], 1 ] + // Base
      Co[ no ] * eta[ subjIndex[ no ], 2 ] + // Subject
      Co[ no ] * zeta[ itemIndex[ no ] ] // Item
    ;
    
    // Coefficient of drift
    param[ no, 3 ] = 1.0;
    
    // Residual latency
    param[ no, 4 ] = theta[ subjIndex[ no ] ] * 
      min_RT[ subjIndex[ no ] ];
    
    // Threshold for picking left
    param[ no, 5 ] = X_k_L[ no ] * beta_k + 
      eta[ subjIndex[ no ], 3 ];
    
    // Drift rate for picking left
    param[ no, 6 ] = 
      X_x_L[ no ] * beta_x + // Group-level effects
      eta[ subjIndex[ no ], 1 ] + // Base
      -Co[ no ] * eta[ subjIndex[ no ], 2 ] + // Subject
      -Co[ no ] * zeta[ itemIndex[ no ] ] // Item
    ;
    
    // Coefficient of drift
    param[ no, 7 ] = 1.0;
    
    // Residual latency
    param[ no, 8 ] = theta[ subjIndex[ no ] ] * 
      min_RT[ subjIndex[ no ] ];
  }
  
}
model {
  // Variable declaration
  vector[ No ] summands; // For custom likelihood
  vector[3] Mu; // Declare vector for location of subject effects
  matrix[3, 3] Sigma; // Covariance matrix
  real alpha_rl; // Standard parameters for beta distribution
  real beta_rl;
  int inc;
  
  // Fix location of subject effects to 0
  Mu[1] = 0.0; Mu[2] = 0.0; Mu[3] = 0.0;
  // Calculate cholesky decomposition of covariance matrix
  Sigma = diag_pre_multiply(tau, L_Omega);
  // Transform values for residual latency hierarchy
  alpha_rl = lambda * phi;
  beta_rl = lambda * (1.0 - phi);

  // Priors
  inc = 1;
  for (i in 1:Nx) {
    beta_x[i] ~ normal( Priors[inc,1], Priors[inc,2] );
    inc = inc + 1;
  }
  for (i in 1:Nk) {
    beta_k[i] ~ normal( Priors[inc,1], Priors[inc,2] );
    inc = inc + 1;
  }
  phi ~ beta( Priors[inc,1], Priors[inc,2] );
  inc = inc + 1;
  lambda ~ pareto( Priors[inc,1], Priors[inc,2] );
  inc = inc + 1;
  sigma_zeta ~ gamma( Priors[inc,1], Priors[inc,2] );
  inc = inc + 1;
  tau ~ gamma( Priors[inc,1], Priors[inc,2] );
  inc = inc + 1;
  L_Omega ~ lkj_corr_cholesky( Priors[inc,1] );
  
  // Hierarchy
  eta ~ multi_normal_cholesky(Mu,Sigma);
  zeta ~ normal( 0.0, sigma_zeta );
  theta ~ beta( alpha_rl, beta_rl );
  
  // Likelihood
  inc = 1;
  for ( no in 1:No ) {
    summands[ no ] = waldrace_lpdf( Y[ no ] | param[ no ] );
    if ( summands[ no ] == log( 0 ) && inc == 1 ) {
      print( no );
      print( param[ no ] );
      print( beta_x )
      print( eta[ subjIndex[ no ] ] );
      print( zeta[ itemIndex[ no ] ] );
      inc = inc + 1;
    };
  }
  
  // Call to the sampler
  target += sum( summands );
}
generated quantities {
  vector[No] logLik;
  corr_matrix[3] Omega;
  
  Omega = multiply_lower_tri_self_transpose(L_Omega);
  
  for (no in 1:No) logLik[no] = 
    waldrace_lpdf( Y[no] | param[ no ] );
}
