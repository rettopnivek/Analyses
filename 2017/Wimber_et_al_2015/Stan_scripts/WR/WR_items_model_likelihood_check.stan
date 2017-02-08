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
model {
  // Included so that script compiles
}
generated quantities {
  // Variable declaration
  vector[No] logLik;
  corr_matrix[3] Omega;
  vector[8] param[No];
  
  Omega = multiply_lower_tri_self_transpose(L_Omega);
  
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
      
    logLik[no] = waldrace_lpdf( Y[no] | param[ no ] );
  }
  
}
