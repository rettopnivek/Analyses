functions {
  #include "mSDT_prob.stan"
}
data {
  int No; // Number of total observations
  int Ns; // Number of subjects
  int Ni; // Number of items
  int Y[No]; // Observed responses
  int Co[No]; // Position of correct response (0 = left, 1 = right)
  int subjIndex[No]; 
  int itemIndex[No];
  int Cond[No,2]; // 1st column: ?, 2nd column, ?
}
parameters {
  real beta[2]; // Group-level effects on d'
  vector[2] eta[Ns]; // Subject-specific effects on d' and criterion
  cholesky_factor_corr[2] L_Omega; // Correlations for subject effects
  vector<lower=0.0>[2] tau; // Scale parameters for subject effects
  real zeta[Ni]; // Item-specific effects on d'
  real<lower=0.0> sigma_zeta; // Variability of item effects
  matrix<lower=0.0,upper=1.0>[Ns,2] lambda; // Mixture probabilities
  real<lower=0.0,upper=1.0> phi[2]; // Mean for mixture probabilities
  real<lower=0.1> nu[2]; // Total counts for mixture probabilities
}
transformed parameters {
  real theta[No]; // Probability of picking right
  real<lower=0.0> alpha[2]; // Positive count parameter
  real<lower=0.0> kappa[2]; // Positive count parameter
  
  // Define probabilities for bernoulli distribution
  {
    real dp;
    real crt;
    real lmb;
    for ( no in 1:No ) {
      dp = exp( beta[ Cond[ no, 1 ] ] + eta[ subjIndex[ no ], 1 ] + 
        zeta[ itemIndex[ no ] ] );
      crt = eta[ subjIndex[ no ], 2 ];
      lmb = lambda[ subjIndex[ no ], Cond[ no, 2 ] ];
      theta[ no ] = mSDT_prob( dp, crt, lmb, Co[no] );
    }
  }
  
  // Reparameterize beta distribution parameters
  for (i in 1:2) {
    alpha[i] = nu[i] * phi[i];
    kappa[i] = nu[i] * (1.0 - phi[i]);
  }
  
}
model {
  vector[2] Mu; // Declare vector for location of subject effects
  matrix[2, 2] Sigma;
  
  // Fix location of subject effects to 0
  Mu[1] = 0.0; Mu[2] = 0.0;
  // Calculate cholesky decomposition of covariance matrix
  Sigma = diag_pre_multiply(tau, L_Omega);
  
  // Priors
  beta ~ normal( 0.0, 1.0 );
  L_Omega ~ lkj_corr_cholesky(2.0);
  tau ~ cauchy(0, 2.5);
  
  // Hierarchy
  eta ~ multi_normal_cholesky(Mu,Sigma);
  zeta ~ normal( 0.0, sigma_zeta );
  col(lambda,1) ~ beta( alpha[1], kappa[1] );
  col(lambda,2) ~ beta( alpha[2], kappa[2] );
  
  // Likelihood
  Y ~ bernoulli(theta);
}
generated quantities {
  vector[No] logLik;
  corr_matrix[2] Omega;
  
  Omega = multiply_lower_tri_self_transpose(L_Omega);
  
  for (no in 1:No) logLik[no] = bernoulli_lpmf( Y[no] | theta[no] );
}

