functions {
  #include "RvSDT_prob.stan"
}
data {
  int No; // Number of total observations
  int Ns; // Number of subjects
  int Ni; // Number of items
  int Nd; // Number of conditions for d'
  int Nc; // Number of conditions for criterion
  int Co[No]; // Position of correct response (0 = left, 1 = right)
  int subjIndex[No];
  int itemIndex[No];
  row_vector[Nd] Xd[No]; // Design matrix for d'
  row_vector[Nc] Xc[No]; // Design matrix for criterion
  vector[Nd] beta_dp; // Group-level effects on d'
  vector[Nc] beta_c; // Group-level effects on criterion;
  corr_matrix[2] Omega; // Correlation matrix for subject effects
  vector<lower=0.0>[2] tau; // Scale parameters for subject effects
  real mu_lambda; // Mean probability of recall
  real<lower=0.0> sigma_lambda; // Variability for recall probabilities
}
model {
  // Included so script will compile
}
generated quantities {
  real theta[No]; // Probability of picking right
  vector[2] eta[Ns]; // Subject-specific effects on d' and criterion
  vector[Ni] lambda; // Logit values for recall probabilities
  int Y[No]; // Observed responses
  
  // Define probabilities for bernoulli distribution
  {
    vector[2] Mu;
    matrix[2,2] Sigma;
    real dp;
    real crt;
    real lmb;
    
    // Generate subject and item effects
    for (ni in 1:Ni) 
      lambda[ni] = normal_rng( 0.0, sigma_lambda );
    for (r in 1:2)
      Sigma[r, r] = tau[r] * tau[r] * Omega[r, r];
    Sigma[1, 2] = tau[1] * tau[2] * Omega[1, 2];
    Sigma[2, 1] = Sigma[1, 2];
    Mu[1] = 0.0; Mu[2] = 0.0;
    for (ns in 1:Ns)
      eta[ns] = multi_normal_rng( Mu, Sigma );
    
    // Loop over observations and calculate probability of picking right
    for ( no in 1:No ) {
      
      // Determine parameters for comparison process
      dp = Xd[no] * beta_dp + eta[ subjIndex[ no ], 1 ];
      crt = Xc[no] * beta_c + eta[ subjIndex[ no ], 2 ];
      // Determine probability of recall
      lmb = inv_logit( mu_lambda + lambda[ itemIndex[ no ] ] );
      // Calculate probability of picking right
      theta[ no ] = RvSDT_prob( lmb, dp, crt, Co[no] );
      // Simulate observations
      Y[ no ] = bernoulli_rng( theta[ no ] );

    }
  }
  
}

