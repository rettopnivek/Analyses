functions {
  #include "mSDT_prob.stan"
}
data {
  int No; // Number of total observations
  int Ns; // Number of subjects
  int Ni; // Number of items
  int Co[No]; // Position of correct response (0 = left, 1 = right)
  int subjIndex[No]; 
  int itemIndex[No];
  int Cond[No,2]; // Indices for conditions
  real beta[2]; // Group-level effects on d'
  corr_matrix[2] Omega; // Correlation matrix for subject effects
  vector<lower=0.0>[2] tau; // Scale parameters for subject effects
  real<lower=0.0> sigma_zeta; // Variability of item effects
  real<lower=0.0,upper=1.0> phi[2]; // Mean for mixture probabilities
  real<lower=0.1> nu[2]; // Total counts for mixture probabilities
}
model {
  // Included so script will compile
}
generated quantities {
  real theta[No]; // Probability of picking right
  int Y[No]; // Observed responses
  real zeta [Ni]; // Item effects
  vector[2] eta[Ns]; // Subject effects
  matrix[Ns,2] lambda; // Mixture probabilities
  
  {
    matrix[2,2] Sigma;
    vector[2] Mu;
    real dp;
    real crt;
    real lmb;
    
    // Generate subject and item effects
    for (ni in 1:Ni) 
      zeta[ni] = normal_rng( 0.0, sigma_zeta );
    for (r in 1:2)
      Sigma[r, r] = tau[r] * tau[r] * Omega[r, r];
    Sigma[1, 2] = tau[1] * tau[2] * Omega[1, 2];
    Sigma[2, 1] = Sigma[1, 2];
    Mu[1] = 0.0; Mu[2] = 0.0;
    for (ns in 1:Ns) {
      eta[ns] = multi_normal_rng( Mu, Sigma );
      lambda[ns,1] = beta_rng( nu[1] * phi[1], nu[1] * (1.0 - phi[1]) );
      lambda[ns,2] = beta_rng( nu[2] * phi[2], nu[2] * (1.0 - phi[2]) );
    }
    // Define probabilities for bernoulli distribution
    for ( no in 1:No ) {
      dp = exp( beta[ Cond[ no, 1 ] ] + eta[ subjIndex[ no ], 1 ] + 
        zeta[ itemIndex[ no ] ] );
      crt = eta[ subjIndex[ no ], 2 ];
      lmb = lambda[ subjIndex[ no ], Cond[ no, 2 ] ];
      // Calculate probability of picking right
      theta[ no ] = mSDT_prob( dp, crt, lmb, Co[no] );
      Y[ no ] = bernoulli_rng( theta[ no ] );
    }
  }
}

