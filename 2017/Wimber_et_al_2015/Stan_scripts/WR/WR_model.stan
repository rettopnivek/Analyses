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
  matrix[ Nx + Nk + 4, 2 ]Priors; // Parameters for priors
}
parameters {
  // Subject-level parameters
  vector<lower=0.0>[Nx] xi[Ns]; // Drift rates
  vector<lower=0.0>[Nk] kappa[Ns]; // Thresholds
  real<lower=0.0,upper=1.0> theta[Ns]; // Proportions for residual
                                       // latency
  // Hyper-parameters
  vector<lower=0.0>[ Nx + Nk ] Mu; // Location parameters
  vector<lower=0.0>[ Nx + Nk ] Tau; // Scale parameters
  cholesky_factor_corr[ Nx + Nk ] L_Omega; // Correlations
  real<lower=0.0,upper=1.0> phi; // Mean for residual latency
                                 // proportion
  real<lower=0.1> lambda; // Total counts for residual latency
                          // proportion
}
transformed parameters {
  vector[ Nx + Nk ] all_prm[Ns];
  
  for (ns in 1:Ns) {
    all_prm[ ns, 1:Nx ] = xi[ns];
    all_prm[ ns, (Nx+1):(Nx+Nk) ] = kappa[ns];
  }
  
}
model {
  // Variable declaration
  vector[ No ] summands; // For custom likelihood
  vector[8] param;
  matrix[ Nx + Nk, Nx + Nk ] Sigma;
  real alpha_rl; // Standard parameters for beta distribution
  real beta_rl;
  int inc;
  
  // Transform values for residual latency hierarchy
  alpha_rl = lambda * phi;
  beta_rl = lambda * (1.0 - phi);
  // Calculate cholesky decomposition of covariance matrix
  Sigma = diag_pre_multiply(Tau, L_Omega);
  
  // Priors
  inc = 1;
  for (i in 1:(Nx+Nk) ) {
    Mu[i] ~ normal( Priors[inc,1], Priors[inc,2] );
    inc = inc + 1;
  }
  phi ~ beta( Priors[inc,1], Priors[inc,2] );
  inc = inc + 1;
  # lambda ~ normal( Priors[inc,1], Priors[inc,2] );
  inc = inc + 1;
  # sigma ~ gamma( Priors[inc,1], Priors[inc,2] );
  inc = inc + 1;
  L_Omega ~ lkj_corr_cholesky( Priors[inc,1] );
  
  // Hierarchy
  all_prm ~ multi_normal_cholesky(Mu,Sigma);
  theta ~ beta( alpha_rl, beta_rl );
  
  // Likelihood
  
  // Loop over observations
  inc = 1;
  for ( no in 1:No ) {
    
    // Calculate parameters for wald race model
    
    // Threshold for picking right
    param[ 1 ] = kappa[ subjIndex[ no ], condIndex[ no, 1 ] ];
    
    // Drift rate for picking right
    param[ 2 ] = xi[ subjIndex[ no ], condIndex[ no, 3 ] ];
    
    // Coefficient of drift
    param[ 3 ] = 1.0;
    
    // Residual latency
    param[ 4 ] = theta[ subjIndex[ no ] ] * 
      min_RT[ subjIndex[ no ] ];
    
    // Threshold for picking right
    param[ 5 ] = kappa[ subjIndex[ no ], condIndex[ no, 2 ] ];
    
    // Drift rate for picking right
    param[ 6 ] = xi[ subjIndex[ no ], condIndex[ no, 4 ] ];
    
    // Coefficient of drift
    param[ 7 ] = 1.0;
    
    // Residual latency
    param[ 8 ] = theta[ subjIndex[ no ] ] * 
      min_RT[ subjIndex[ no ] ];
    
    // Calculate log-density
    summands[ no ] = waldrace_lpdf( Y[ no ] | param );
  }
  
  // Call to the sampler
  target += sum( summands );
}
generated quantities {
  vector[No] logLik;
  vector[8] param[No];
  corr_matrix[Nx+Nk] Omega;
  
  Omega = multiply_lower_tri_self_transpose(L_Omega);
  
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
