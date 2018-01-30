data {
  int No; // Number of observations
  int Ns; // Number of subjects
  int K; // Number of IVs
  row_vector[K] X[No]; // Design matrix
  int subj_index[No]; // Index for subjects
  real Y[No]; // Observed data
}
parameters {
  // Standardized subject-level coefficients
  vector[K] beta_raw[Ns];
  // Parameters for multivariate normal hierarchy
  vector[K] Mu;
  cholesky_factor_corr[K] L_Omega;
  vector<lower=0.0>[K] Tau;
  // Error term for regression
  real<lower=0.0> sigma;
}
transformed parameters {
  vector[K] beta[Ns]; // Unstandardized coefficients
  vector[No] mu;
  
  // Loop over subjects
  for ( ns in 1:Ns ) {
    beta[ns] = Mu + Tau .* (L_Omega * beta_raw[ns]);
  }
  
  // Loop over observations
  for ( no in 1:No ) {
    mu[no] = dot_product( X[no], beta[ subj_index[no] ] );
  }
  
}
model {
  
  // Priors
  Mu ~ normal( 0.0, 1.0 );
  L_Omega ~ lkj_corr_cholesky(2.0);
  Tau ~ cauchy( 0.0, 5.0 );
  sigma ~ cauchy( 0.0, 5.0 );
  
  // Hierarchy
  for ( ns in 1:Ns ) 
  beta_raw[ns] ~ normal( 0.0, 1.0 );
  
  // Log-likelihood
  Y ~ normal( mu, sigma );
}
generated quantities {
  // Reconstruct correlation matrix
  corr_matrix[K] Omega;
  
  Omega = multiply_lower_tri_self_transpose(L_Omega);
}

