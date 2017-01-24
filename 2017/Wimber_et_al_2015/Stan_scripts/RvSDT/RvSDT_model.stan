functions {
  #include "RvSDT_lpmf.stan"
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
  real<lower=0.0,upper=1.0> theta[Ni]; // Probability of recall
  real<lower=0.0,upper=1.0> phi; // Mean for recall probability
  real<lower=0.1> nu; // Total counts for recall probability
}
transformed parameters {
  real alpha;
  real kappa;
  
  alpha = nu * phi;
  kappa = nu * (1.0 - phi);
}
model {
  vector[2] Mu; // Declare vector for location of subject effects
  matrix[2, 2] Sigma;
  real dp;
  real crt;
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
  phi ~ beta( Priors[inc,1], Priors[inc,2] );
  inc = inc + 1;
  nu ~ normal( Priors[inc,1], Priors[inc,2] );
  inc = inc + 1;
  tau ~ gamma( Priors[inc,1], Priors[inc,2] );
  inc = inc + 1;
  L_Omega ~ lkj_corr_cholesky( Priors[inc,1] );
  
  // Hierarchy
  eta ~ multi_normal_cholesky(Mu,Sigma);
  theta ~ beta( alpha, kappa );
  
  // Likelihood
  for (no in 1:No) {
    
    dp = Xd[no] * beta_dp + eta[ subjIndex[ no ], 1 ];
    crt = Xc[no] * beta_c + eta[ subjIndex[ no ], 2 ];
    
    target += RvSDT_lpmf( Y[no] | 
                          theta[ itemIndex[ no ] ],
                          dp,
                          crt,
                          Co[ no ] );
  }
}
generated quantities {
  vector[No] logLik;
  corr_matrix[2] Omega;
  
  Omega = multiply_lower_tri_self_transpose(L_Omega);
  
  {
    real dp;
    real crt;
    
    // Likelihood
    for (no in 1:No) {
      
      dp = Xd[no] * beta_dp + eta[ subjIndex[ no ], 1 ];
      crt = Xc[no] * beta_c + eta[ subjIndex[ no ], 2 ];
      
      logLik[no] = RvSDT_lpmf( Y[no] | 
                               theta[ itemIndex[ no ] ],
                               dp,
                               crt,
                               Co[ no ] );
    }
  }
}
