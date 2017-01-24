functions {
  #include "mSDT_prob.stan"
}
data {
  int No; // Number of observations
  int Cond[No,2]; // Condition index
  int Y[No]; // Observed responses (0 = left, 1 = right)
  real Co[No]; // Position of correct answer (0 = left, 1 = right)
}
parameters {
  real<lower=0.0> dp[2];
  real crt;
  real<lower=0.0,upper=1.0> lmb[2];
}
transformed parameters {
  real theta[No];
  
  // Calculate probability of picking right
  for (no in 1:No) {
    theta[no] = mSDT_prob( dp[ Cond[no,1] ], crt, 
      lmb[ Cond[no,2] ], Co[no] );
  }
}
model {
  
  // Priors
  dp ~ normal( 0.0, 1.0 );
  crt ~ normal( 0.0, 1.0 );
  lmb ~ beta( 4.0, 1.0 );
  
  // Likelihood
  Y ~ bernoulli(theta);
}

