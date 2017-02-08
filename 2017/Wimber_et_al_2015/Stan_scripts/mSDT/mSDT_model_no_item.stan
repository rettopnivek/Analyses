functions {
  #include "mSDT_prob.stan"
}
data {
  int No; // Number of total observations
  int Ns; // Number of subjects
  int Nd; // Total number of conditions for d' values
  int Nc; // Total number of conditions for criterion values
  int Nm; // Total number of conditions for mixture probabilities
  int Y[No]; // Observed responses
  int Co[No]; // Position of correct response (0 = left, 1 = right)
  int subjIndex[No];
  matrix[No,Nd] Xd; // Design matrix for d'
  matrix[No,Nc] Xc; // Design matrix for criterion
  matrix[No,Nm] Xm; // Design matrix for mixture probability
  matrix[ Nd + Nc + Nm + 3, 2 ] Priors; // Matrix of prior values
}
parameters {
  vector[Nd] beta_dp; // Group-level effects for d'
  vector[Nc] beta_c; // Group-level effects for criterion
  vector[2] eta[Ns]; // Subject-specific effects on d' and criterion
  cholesky_factor_corr[2] L_Omega; // Correlations for subject effects
  vector<lower=0.0>[2] tau; // Scale parameters for subject effects
  vector[Nm] beta_m; // Group-level effects for mixture probabilities
  vector[Ns] lambda; // Mixture probabilities
  real<lower=0.0> sigma_lambda; // Variability for mixture probabilities
}
transformed parameters {
  real theta[No]; // Probability of picking right
  
  // Define probabilities for bernoulli distribution
  {
    real dp;
    real crt;
    real lmb;
    vector[No] mu_d;
    vector[No] mu_c;
    vector[No] mu_m;
    
    // Calculate group-level effects
    mu_d = Xd * beta_dp;
    mu_c = Xc * beta_c;
    mu_m = Xm * beta_m;
    
    // Loop over observations and calculate probability of picking right
    for ( no in 1:No ) {
      
      dp = mu_d[ no ] + eta[ subjIndex[ no ], 1 ];
      crt = mu_c[ no ] + eta[ subjIndex[ no ], 2 ];
      lmb = inv_logit( mu_m[ no ] + lambda[ subjIndex[ no ] ] );
      theta[ no ] = mSDT_prob( dp, crt, lmb, Co[no] );
      
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
  for (i in 1:Nm) {
    beta_m[i] ~ normal( Priors[inc,1], Priors[inc,2] );
    inc = inc + 1;
  }
  sigma_lambda ~ gamma( Priors[inc,1], Priors[inc,2] );
  inc = inc + 1;
  tau ~ gamma( Priors[inc,1], Priors[inc,2] );
  inc = inc + 1;
  L_Omega ~ lkj_corr_cholesky( Priors[inc,1] );
  
  // Hierarchy
  eta ~ multi_normal_cholesky(Mu,Sigma);
  lambda ~ normal( 0.0, sigma_lambda );
  
  // Likelihood
  Y ~ bernoulli(theta);
}
generated quantities {
  vector[No] logLik;
  corr_matrix[2] Omega;
  
  Omega = multiply_lower_tri_self_transpose(L_Omega);
  
  for (no in 1:No) logLik[no] = bernoulli_lpmf( Y[no] | theta[no] );
}


