functions {
  #include "SDT_prob.stan"
}
data {
  int No; // Number of total observations
  int Ns; // Number of subjects
  int Ni; // Number of items
  int Nd; // Number of conditions for d'
  int Nc; // Number of conditions for criterion
  int Y[No]; // Observed responses
  int Co[No]; // Position of correct response (0 = left, 1 = right)
  int subjIndex[No]; 
  int itemIndex[No];
  row_vector[Nd] Xd[No]; // Design matrix for d'
  row_vector[Nc] Xc[No]; // Design matrix for criterion
  matrix[ Nd + Nc + 3, 2 ] Priors; // Matrix of prior values
}
parameters {
  vector[Nd] eta_dp[Ns]; // Subject-specific effects on d'
  vector[Nc] eta_c[Ns]; // Subject-specific effects on criterion
  vector[Nd+Nc] Mu; // Group-level means
  cholesky_factor_corr[Nd+Nc] L_Omega; // Correlations for subject effects
  vector<lower=0.0>[Nd+Nc] tau; // Scale parameters for subject effects
  real zeta[Ni]; // Item-specific effects on d'
  real<lower=0.0> sigma_zeta;
}
transformed parameters {
  real theta[No]; // Probability of picking right
  
  // Define probabilities for bernoulli distribution
  {
    real dp;
    real crt;
    
    for ( no in 1:No ) {
      // Calculate d' as sum of subject and item effects
      dp = Xd[ no ] * eta_dp[ subjIndex[ no ] ] + 
        zeta[ itemIndex[ no ] ];
      // Calculate criterion based on subject effects
      crt = Xc[ no ] * eta_c[ subjIndex[ no ] ];
      // Calculate probability of picking right
      theta[ no ] = SDT_prob( dp, crt, Co[no] );
    }
  }
  
}
model {
  matrix[ Nd+Nc, Nd+Nc ] Sigma;
  vector[ Nd + Nc ] eta[ Ns ];
  int inc;
  
  // Fill in vector
  for (ns in 1:Ns) {
    eta[ ns, 1:Nd ] = eta_dp[ ns ];
    eta[ ns, (Nd+1):(Nd+Nc) ] = eta_c[ ns ];
  }
  
  // Calculate cholesky decomposition of covariance matrix
  Sigma = diag_pre_multiply(tau, L_Omega);
  
  // Priors
  inc = 1;
  for (i in 1:Nd) {
    Mu[i] ~ normal( Priors[inc,1], Priors[inc,2] );
    inc = inc + 1;
  }
  for (i in 1:Nc) {
    Mu[i] ~ normal( Priors[inc,1], Priors[inc,2] );
    inc = inc + 1;
  }
  sigma_zeta ~ gamma( Priors[inc,1], Priors[inc,2] );
  inc = inc + 1;
  tau ~ gamma( Priors[inc,1], Priors[inc,2] );
  inc = inc + 1;
  L_Omega ~ lkj_corr_cholesky( Priors[inc,1] );
  
  // Hierarchy
  eta ~ multi_normal_cholesky(Mu,Sigma);
  zeta ~ normal( 0.0, sigma_zeta );
  
  // Likelihood
  Y ~ bernoulli(theta);
}
generated quantities {
  vector[No] logLik;
  corr_matrix[Nd+Nc] Omega;
  
  Omega = multiply_lower_tri_self_transpose(L_Omega);
  
  for (no in 1:No) logLik[no] = bernoulli_lpmf( Y[no] | theta[no] );
}
