functions {
  #include "SDT_prob.stan"
}
data {
  int No; // Number of total observations
  int Nd; // Number of coefficients for d'
  int Nc; // Number of coefficients for criterion
  // Number of subject-level dimensions
  int<lower=2,upper=Nd+Nc> Nsp;
  int Nt[No]; // Total number of trials per observations
  int Co[No]; // Presence of target (1 = yes, 0 = no)
  row_vector[Nd] Xd[No]; // Design matrix for d'
  row_vector[Nc] Xc[No]; // Design matrix for criterion
  // Design matrix position for subject-level parameters
  int<lower=1,upper=Nd+Nc> eta_pos[Nsp];
  // Posterior samples for group-level
  int S; // Number of posterior samples
  vector[Nd+Nc] group_param[S];
  corr_matrix[Nsp] Omega[S];
  vector<lower=0.0>[Nsp] Tau[S];
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
model {
  // Necessary for script to compile
}
generated quantities {
  matrix[No,S] P; // Simulated proportions
  
  // Define probabilities for binomial distribution
  {
    real theta;
    int is;
    real dp;
    real crt;
    // Matrix of all possible subject effects
    vector[Nd+Nc] eta;
    // Variables for generating subject-level parameters
    matrix[Nsp,Nsp] L_Omega;
    matrix[Nsp,Nsp] L_Sigma;
    vector[Nsp] subj_param;
    
    // Initialize eta with all zeroes
    eta = rep_vector( 0.0, Nd+Nc);
    
    // Loop posterior samples
    for ( s in 1:S ) {
      
      // Generate subject-level parameters
      L_Omega = cholesky_decompose( Omega[s] );
      L_Sigma = diag_pre_multiply(Tau[s], L_Omega);
      subj_param = multi_normal_cholesky_rng(Mu,L_Sigma);
        
      // Fill in non-zero values for eta
      eta[ eta_pos ] = subj_param;
      
      // Loop over observations
      for ( no in 1:No ) {
        
        // Compute d' value from weighted combinations
        dp = dot_product( Xd[no], 
          group_param[ s, dp_sel ] + eta[ dp_sel ] );
        // Compute criterion value from weighted combinations
        crt = dot_product( Xc[no], 
          group_param[ s, crt_sel ] + eta[ crt_sel ] );
        
        // Compute probability
        theta = SDT_prob( dp, crt, Co[no] );
        
        // Simulate an observed count
        P[no,s] = binomial_rng( Nt[no], theta );
        // Convert to proportion
        P[no,s] = P[no,s] / Nt[no];
      }
    }
  }
}

