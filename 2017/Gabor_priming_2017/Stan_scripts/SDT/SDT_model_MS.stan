functions {
  #include "SDT_prob.stan"
}
data {
  int Ns; // Number of subjects
  int No; // Total number of observations
  int Nd; // Number of d' values to estimate
  int Nc; // Number of criterion values to estimate
  row_vector[Nd] X_d[No]; // Design matrix for d'
  row_vector[Nc] X_crt[No]; // Design matrix for criterion
  int subj_index[No]; // Index for subjects
  real Co[No]; // Position of correct answer (0 = left, 1 = right)
  int Y[No]; // Observed responses (0 = left, 1 = right)
  matrix[Nd+Nc+1,4] Priors;
}
transformed data {
  int K; // Total number of variables

  K = Nd + Nc;
}
parameters {
  vector[K] beta_raw[Ns]; // Standardized coefficients
  // Parameters for multivariate normal hierarchy
  vector[K] Mu;
  cholesky_factor_corr[K] L_Omega;
  vector<lower=0.0>[K] Tau;
}
transformed parameters {
  real theta[No]; // Probability of picking right
  vector[K] beta[Ns]; // Unstandardized coefficients
  
  // Local block
  {
    // Variable declarations
    real cur_dp;
    real cur_crt;
    int s;
    
    // Compute unstandardized coefficients per subject
    for ( ns in 1:Ns ) {
      beta[ns] = Mu + Tau .* (L_Omega * beta_raw[ns]);
    }
    
    // Calculate probability of picking right
    for (no in 1:No) {
      
      s = subj_index[ no ];
      cur_dp = dot_product( X_d[no], beta[s,1:Nd] );
      cur_crt = dot_product( X_crt[no], beta[s,(Nd+1):K] );
      
      theta[no] = SDT_prob( cur_dp, cur_crt, Co[no] );
    }
  }
  
}
model {
  
  // Priors
  for ( k in 1:K ) {
    Mu[k] ~ normal( Priors[k,1], Priors[k,2] );
    Tau[k] ~ gamma( Priors[k,3], Priors[k,4] );
  }
  L_Omega ~ lkj_corr_cholesky( Priors[K+1,1] );
  
  // Hierarchy
  for ( ns in 1:Ns ) 
    beta_raw[ns] ~ normal( 0.0, 1.0 );
  
  // Likelihood
  Y ~ bernoulli(theta);
}
generated quantities {
  // Reconstruct correlation matrix
  corr_matrix[K] Omega;
  
  Omega = multiply_lower_tri_self_transpose(L_Omega);
}

