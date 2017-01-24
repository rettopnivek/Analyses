functions {
  #include "SDT_prob.stan"
}
data {
  int No; // Number of observations
  int Cond[No]; // Condition index
  int Y[No]; // Observed responses (0 = left, 1 = right)
  real Co[No]; // Position of correct answer (0 = left, 1 = right)
}
parameters {
  real<lower=0.0> dp[4]; // d' values
  real crt; // Criterion values
}
transformed parameters {
  real theta[No]; // Probability of picking right
  
  // Calculate probability of picking right
  for (no in 1:No) {
    theta[no] = SDT_prob( dp[ Cond[no] ], crt, Co[no] );
  }
}
model {
  
  // Priors
  dp ~ normal( 0.0, 1.0 );
  crt ~ normal( 0.0, 1.0 );
  
  // Likelihood
  Y ~ bernoulli(theta);
}

