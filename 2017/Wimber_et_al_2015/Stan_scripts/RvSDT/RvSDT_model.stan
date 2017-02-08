functions {
  #include "RvSDT_prob.stan"
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
  matrix[ Nd + Nc + 4, 2 ] Priors; // Matrix of prior values
}
parameters {
  vector[Nd] beta_dp; // Group-level effects on d'
  vector[Nc] beta_c; // Group-level effects on criterion;
  vector[2] eta[Ns]; // Subject-specific effects on d' and criterion
  cholesky_factor_corr[2] L_Omega; // Correlations for subject effects
  vector<lower=0.0>[2] tau; // Scale parameters for subject effects
  vector[Ni] lambda; // Logit values for recall probabilities
  real mu_lambda; // Mean probability of recall
  real<lower=0.0> sigma_lambda; // Variability for recall probabilities
}
transformed parameters {
  real theta[No]; // Probability of picking right
  
  // Define probabilities for bernoulli distribution
  {
    real dp;
    real crt;
    real lmb;

    // Loop over observations and calculate probability of picking right
    for ( no in 1:No ) {
      
      dp = Xd[no] * beta_dp + eta[ subjIndex[ no ], 1 ];
      crt = Xc[no] * beta_c + eta[ subjIndex[ no ], 2 ];
      
      lmb = inv_logit( lambda[ itemIndex[ no ] ] );
      theta[ no ] = RvSDT_prob( lmb, dp, crt, Co[no] );
      
    }
  }
  
}
model {
  vector[2] Mu; // Declare vector for location of subject effects
  matrix[2, 2] Sigma;
  int inc;
  
  // Fix location of subject effects to 0
  Mu[1] = 0.0; Mu[2] = 0.0;
  // Calculate cholesky decomposition of covariance matrix
  Sigma = diag_pre_multiply(tau, L_Omega);
  
  // Priors
  inc = 1;
  for (i in 1:Nd) {
    beta_dp[i] ~ normal( Priors[inc,1], Priors[inc,2] );
    inc = inc + 1;
  }
  for (i in 1:Nc) {
    beta_c[i] ~ normal( Priors[inc,1], Priors[inc,2] );
    inc = inc + 1;
  }
  mu_lambda ~ normal( Priors[inc,1], Priors[inc,2] );
  inc = inc + 1;
  sigma_lambda ~ gamma( Priors[inc,1], Priors[inc,2] );
  inc = inc + 1;
  tau ~ gamma( Priors[inc,1], Priors[inc,2] );
  inc = inc + 1;
  L_Omega ~ lkj_corr_cholesky( Priors[inc,1] );
  
  // Hierarchy
  eta ~ multi_normal_cholesky(Mu,Sigma);
  lambda ~ normal( mu_lambda, sigma_lambda );
  
  // Likelihood
  Y ~ bernoulli(theta);
}
generated quantities {
  vector[No] logLik;
  corr_matrix[2] Omega;
  
  Omega = multiply_lower_tri_self_transpose(L_Omega);
  
  for (no in 1:No) logLik[no] = bernoulli_lpmf( Y[no] | theta[no] );
}

