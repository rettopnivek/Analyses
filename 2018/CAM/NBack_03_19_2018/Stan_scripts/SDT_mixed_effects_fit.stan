functions {
  #include "SDT_prob.stan"
}
data {
  int No; // Number of total observations
  int Ns; // Number of subjects
  int Nd; // Number of coefficients for d'
  int Nc; // Number of coefficients for criterion
  // Number of subject-level dimensions
  int<lower=2,upper=Nd+Nc> Nsp;
  int Y[No]; // Observed counts
  int Nt[No]; // Total number of trials per observations
  int Co[No]; // Presence of target (1 = yes, 0 = no)
  int<lower=1,upper=Ns> subjIndex[No]; // Index for subjects
  row_vector[Nd] Xd[No]; // Design matrix for d'
  row_vector[Nc] Xc[No]; // Design matrix for criterion
  // Design matrix position for subject-level parameters
  int<lower=1,upper=Nd+Nc> eta_pos[Nsp];
  matrix[ Nd+Nc+2, 2 ] Priors; // Parameters governing priors
}
transformed data {
  int dp_sel[Nd];
  int crt_sel[Nc];
  vector[Nsp] Mu;
  
  // Create indices for separating 
  // d' and criterion parameters
  for ( i in 1:Nd ) dp_sel[i] = i;
  for (i in 1:Nc) crt_sel[i] = i + Nd;
  
  // Create a zero-vector
  Mu = rep_vector( 0.0, Nsp );
}
parameters {
  // Standardized subject-level coefficients
  vector[Nsp] subj_param[Ns];
  // Group-level coefficients
  vector[Nd+Nc] group_param;
  // Parameters for multivariate normal distribution 
  // for subject-level coefficients
  corr_matrix[Nsp] Omega;
  vector<lower=0.0>[Nsp] Tau;
}
transformed parameters {
  real theta[No]; // Probability of saying the target is present
  
  // Define probabilities for binomial distribution
  {
    int is;
    real dp;
    real crt;
    // Matrix of all possible subject effects
    vector[Nd+Nc] eta;
    
    // Initialize eta with all zeroes
    eta = rep_vector( 0.0, Nd+Nc);
    
    for ( no in 1:No ) {
      
      eta[ eta_pos ] = subj_param[ subjIndex[no] ];
      
      // Compute d' value from weighted combinations
      dp = dot_product( Xd[no], 
        group_param[ dp_sel ] + eta[ dp_sel ] );
      // Compute criterion value from weighted combinations
      crt = dot_product( Xc[no], 
        group_param[ crt_sel ] + eta[ crt_sel ] );
      
      // Compute probability
      theta[no] = SDT_prob( dp, crt, Co[no] );
    }
  }
  
}
model {
  // Variable declarations
  int inc;
  
  // Priors
  inc = 1;
  for ( i in 1:(Nd+Nc) ) {
    group_param[i] ~ normal( Priors[i,1], Priors[i,2] );
    inc = inc + 1;
  }
  Omega ~ lkj_corr(Priors[inc,1]);
  inc = inc + 1;
  Tau ~ cauchy( Priors[inc,1], Priors[inc,2] );
  
  // Hierarchy
  for ( ns in 1:Ns ) 
    subj_param[ns] ~ multi_normal( Mu, quad_form_diag(Omega,Tau) );
    
  // Likelihood
  Y ~ binomial( Nt, theta );
}
generated quantities {
  vector[No] logLik;
  
  for (no in 1:No) 
    logLik[no] = binomial_lpmf( Y[no] | Nt[no], theta[no] );
}
