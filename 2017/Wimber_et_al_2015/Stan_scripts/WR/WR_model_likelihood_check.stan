functions {
  #include "Wald_race_stan_functions.stan"
}
data {
  int No; // Total number of observations
  int Ns; // Number of subjects
  int Nx; // Number of conditions for drift rate
  int Nk; // Number of conditions for threshold
  vector[2] Y[No]; // Matrix of response times and choice
  int subjIndex[No];
  int condIndex[No,4];
  real<lower=0.0> min_RT[Ns]; // Smallest response time for 
                              // each subject
  // Subject-level parameters
  matrix[ Nx + Nk, Ns ] raw_prm; // Standardized subject parameters
  real<lower=0.0,upper=1.0> theta[Ns]; // Proportions for residual
                                       // latency
  // Hyper-parameters
  vector[ Nx + Nk ] mu;
  vector<lower=0.0>[ Nx + Nk ] sigma; 
  cholesky_factor_corr[ Nx + Nk ] L_Omega;
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
  vector[Nx] xi[Ns]; // Drift rate
  vector[Nk] kappa[Ns]; // Threshold
  vector[No] logLik;
  vector[8] param[No];
  
  Omega = multiply_lower_tri_self_transpose(L_Omega);
  
  // Declare a local block
  {
    matrix[ Nx + Nk, Ns ] prm_scl;
    vector[ Nx + Nk ] prm_hat;
    
    // Re-scale raw parameters
    prm_scl = (diag_pre_multiply(sigma, L_Omega) * raw_prm);
    
    // Loop over subjects
    for ( ns in 1:Ns) {
      
      // Re-center raw parameters
      prm_hat = exp( mu + col(prm_scl,ns) );
      
      xi[ns] = prm_hat[ 1:Nx ];
      kappa[ns] = prm_hat[ (Nx+1):(Nx+Nk) ];
      
    }
    
  }
  
  // Loop over observations
  for ( no in 1:No ) {
    
    // Calculate parameters for wald race model
    
    // Threshold for picking right
    param[ no, 1 ] = kappa[ subjIndex[ no ], condIndex[ no, 1 ] ];
    
    // Drift rate for picking right
    param[ no, 2 ] = xi[ subjIndex[ no ], condIndex[ no, 3 ] ];
    
    // Coefficient of drift
    param[ no, 3 ] = 1.0;
    
    // Residual latency
    param[ no, 4 ] = theta[ subjIndex[ no ] ] * 
      min_RT[ subjIndex[ no ] ];
    
    // Threshold for picking right
    param[ no, 5 ] = kappa[ subjIndex[ no ], condIndex[ no, 2 ] ];
    
    // Drift rate for picking right
    param[ no, 6 ] = xi[ subjIndex[ no ], condIndex[ no, 4 ] ];
    
    // Coefficient of drift
    param[ no, 7 ] = 1.0;
    
    // Residual latency
    param[ no, 8 ] = theta[ subjIndex[ no ] ] * 
      min_RT[ subjIndex[ no ] ];
      
    // Calculate log-density
    logLik[ no ] = waldrace_lpdf( Y[ no ] | param[ no ] );
  }
  
}
