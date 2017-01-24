functions {
  #include "SDT_prob.stan"
}
data {
  int No; // Number of total observations
  int Ns; // Number of subjects
  int Ni; // Number of items
  int Nd; // Number of conditions for d'
  int Nc; // Number of conditions for criterion
  int<lower=0,upper=1> Co[No]; // Position of correct response (0 = left)
  int subjIndex[No];
  int itemIndex[No];
  matrix[No,Nd] Xd; // Design matrix for d'
  matrix[No,Nc] Xc; // Design matrix for criterion
  vector[Nd] beta_dp; // Group-level effects on d'
  vector[Nc] beta_c; // Group-level effects on criterion
  corr_matrix[2] Omega; // Correlation matrix for subject effects
  vector<lower=0.0>[2] tau; // Scale parameters for subject effects
  real<lower=0.0> sigma_zeta;
}
model {
  // Included so script will compile
}
generated quantities {
  real theta[No]; // Probability of picking right
  int Y[No]; // Observed responses
  real zeta [Ni]; // Item effects
  vector[2] eta[Ns]; // Subject effects
  
  // Create local block
  {
    vector[2] Mu;
    matrix[2,2] Sigma;
    vector[No] mu_d;
    vector[No] mu_c;
    real dp;
    real crt;
    
    // Generate subject and item effects
    for (ni in 1:Ni) 
      zeta[ni] = normal_rng( 0.0, sigma_zeta );
    for (r in 1:2)
      Sigma[r, r] = tau[r] * tau[r] * Omega[r, r];
    Sigma[1, 2] = tau[1] * tau[2] * Omega[1, 2];
    Sigma[2, 1] = Sigma[1, 2];
    Mu[1] = 0.0; Mu[2] = 0.0;
    for (ns in 1:Ns)
      eta[ns] = multi_normal_rng( Mu, Sigma );
    
    // Calculate group-level effects
    mu_d = Xd * beta_dp;
    mu_c = Xc * beta_c;
    
    // Loop over observatinos
    for ( no in 1:No ) {
      // Calculate d' as sum of group, subject, and item effects
      dp =  mu_d[no] + eta[ subjIndex[ no ], 1 ] + 
        zeta[ itemIndex[ no ] ];
      // Calculate criterion based on subject effects
      crt = mu_c[no] + eta[ subjIndex[ no ], 2 ];
      // Calculate probability of picking right
      theta[ no ] = SDT_prob( dp, crt, Co[no] );
      // Simulate observations
      Y[ no ] = bernoulli_rng( theta[ no ] );
    }
  }
  
}

