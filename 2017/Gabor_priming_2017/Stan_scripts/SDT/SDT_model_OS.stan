functions {
  #include "SDT_prob.stan"
}
data {
  int No; // Number of observations
  int Nd; // Number of d' values to estimate
  int Nc; // Number of criterion values to estimate
  row_vector[Nd] X_d[No]; // Design matrix for d'
  row_vector[Nc] X_crt[No]; // Design matrix for criterion
  real Co[No]; // Position of correct answer (0 = left, 1 = right)
  int Y[No]; // Observed responses (0 = left, 1 = right)
  matrix[Nd+Nc,2] Priors;
}
parameters {
  vector[Nd+Nc] beta_raw; // Standardized coefficients
}
transformed parameters {
  real theta[No]; // Probability of picking right
  vector[Nd+Nc] beta; // Unstandardized coefficients
  
  // Local block
  {
    real cur_dp;
    real cur_crt;
    
    // Compute unstandardized coefficients
    for (i in 1:(Nd+Nc)) {
      beta[i] = beta_raw[i] * Priors[i,2] + Priors[i,1];
    }
    
    // Calculate probability of picking right
    for (no in 1:No) {
      
      cur_dp = dot_product( X_d[no], beta[ 1:Nd ] );
      cur_crt = dot_product( X_crt[no], beta[ (Nd+1):(Nd+Nc) ] );
      
      theta[no] = SDT_prob( cur_dp, cur_crt, Co[no] );
    }
  }
  
}
model {
  
  // Priors
  beta_raw ~ normal( 0.0, 1.0 );
  
  // Likelihood
  Y ~ bernoulli(theta);
}

