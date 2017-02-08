functions {
  #include "Wald_race_stan_functions.stan"
}
data {
  int No; // Total number of observations
  int Ns; // Number of subjects
  int Ni; // Number of items
  int Nx; // Number of conditions for drift rate
  int Nk; // Number of conditions for threshold
  int Co[No]; // Position of correct response (-1 = left, 1 = right)
  int subjIndex[No];
  int itemIndex[No];
  row_vector[Nx] X_x_R[No]; // Design matrices
  row_vector[Nx] X_x_L[No]; 
  row_vector[Nk] X_k_R[No]; 
  row_vector[Nk] X_k_L[No]; 
  real<lower=0.0> min_RT[Ns]; // Fastest RTs
  // Parameters
  vector[Nx] beta_x; // Group-level effects for drift rate
  vector[Nk] beta_k; // Group-level effects for threshold
  corr_matrix[3] Omega; // Correlations for subject effects
  vector<lower=0.0>[3] tau; // Scale parameters for subject effects
  real<lower=0.0> sigma_zeta; // Standard deviation for item effects
  real<lower=0.0,upper=1.0> phi; // Mean for residual latency proportion
  real<lower=0.1> lambda; // Total counts for residual latency proportion
}
model {
  // Included so script will compile
}
generated quantities {
  vector[2] Y[No]; // Observed responses
  real zeta [Ni]; // Item effects
  vector[3] eta[Ns]; // Subject effects
  real theta[Ns]; // Proportion for minimum RT
  
  // Create local block
  {
    vector[3] Mu;
    matrix[3,3] Sigma;
    vector[8] param;
    real alpha_rl;
    real beta_rl;
    
    // Generate subject and item effects
    for (ni in 1:Ni) 
      zeta[ni] = normal_rng( 0.0, sigma_zeta );
    for (r in 1:3)
      Sigma[r, r] = tau[r] * tau[r] * Omega[r, r];
    for ( r in 2:3 ) {
      for ( c in 1:(4 - r) ) {
        Sigma[ r + (c - 1), r - 1 ] = 
          tau[r] * tau[r - 1] * Omega[ r + (c - 1), r - 1];
        Sigma[ r - 1, r + (c - 1) ] = Sigma[ r + (c - 1), r - 1 ];
      }
    }
    Mu[1] = 0.0; Mu[2] = 0.0; Mu[3] = 0.0;
    alpha_rl = lambda * phi;
    beta_rl = lambda * (1.0 - phi);
    for (ns in 1:Ns) {
      eta[ns] = multi_normal_rng( Mu, Sigma );
      theta[ns] = beta_rng( alpha_rl, beta_rl );
    }
    
    // Loop over observations
    for ( no in 1:No ) {
      
      // Calculate parameters for wald race model
    
    // Threshold for picking right
    param[ 1 ] = exp( X_k_R[ no ] * beta_k + 
      eta[ subjIndex[ no ], 3 ] );
    
    // Drift rate for picking right
    param[ 2 ] = exp( 
      X_x_R[ no ] * beta_x + // Group-level effects
      eta[ subjIndex[ no ], 1 ] + // Base
      Co[ no ] * eta[ subjIndex[ no ], 2 ] + // Subject
      Co[ no ] * zeta[ itemIndex[ no ] ] // Item
    );
    
    // Coefficient of drift
    param[ 3 ] = 1.0;
    
    // Residual latency
    param[ 4 ] = theta[ subjIndex[ no ] ] * 
      min_RT[ subjIndex[ no ] ];
    
    // Threshold for picking left
    param[ 5 ] = exp( X_k_L[ no ] * beta_k + 
      eta[ subjIndex[ no ], 3 ] );
    
    // Drift rate for picking left
    param[ 6 ] = exp( 
      X_x_L[ no ] * beta_x + // Group-level effects
      eta[ subjIndex[ no ], 1 ] + // Base
      -Co[ no ] * eta[ subjIndex[ no ], 2 ] + // Subject
      -Co[ no ] * zeta[ itemIndex[ no ] ] // Item
    );
    
    // Coefficient of drift
    param[ 7 ] = 1.0;
    
    // Residual latency
    param[ 8 ] = theta[ subjIndex[ no ] ] * 
      min_RT[ subjIndex[ no ] ];
          
      Y[no] = waldrace_rng( param );
    }
  }
  
}
