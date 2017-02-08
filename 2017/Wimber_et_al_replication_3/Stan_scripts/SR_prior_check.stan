data {
  int<lower=1> No; // Total number of observations
  int<lower=1> Ns; // Total number of subjects
  int<lower=2> K; // Number of response categories
  int indS[No]; // Subject index per trial
  matrix[K-1,4] Priors; // Matrix of parameter values for priors
}
model {
  // Included so script will compile
}
generated quantities {
  // Variable declarations
  matrix[ Ns, K-1 ] beta_raw; // Subject-level coefficients
  real mu_beta[K-1]; // Group-level means
  real<lower=0> sigma_beta[K-1]; // Group-level standard deviations
  matrix[ Ns, K ] beta; // Coefficients for all categories
  int<lower=1,upper=K> Y[No]; // Simulated responses
  int Targets;
  int Competitors;
  int Errors;
  int Unknown;
  vector[Ns] zeros;
  
  // Generate parameters from priors
  
  for (i in 1:(K-1)) {
    mu_beta[i] = normal_rng( Priors[i,1], Priors[i,2] );
    sigma_beta[i] = gamma_rng( Priors[i,3], Priors[i,4] );
  }
  
  // Hierarchy
  for (i in 1:3) {
    for (ns in 1:Ns) {
      beta_raw[ns,i] = normal_rng( mu_beta[i], sigma_beta[i] );
    }
  }
  zeros = rep_vector(0, Ns);
  beta = append_col( beta_raw, zeros );
  
  // Simulate observations
  Targets = 0;
  Competitors = 0;
  Errors = 0;
  Unknown = 0;
  for ( no in 1:No ) {
    Y[no] = categorical_rng( softmax( to_vector( beta[ indS[no] ] ) ) );
    # Category frequencies
    if ( Y[no] == 1 ) Targets = Targets + 1;
    if ( Y[no] == 2 ) Competitors = Competitors + 1;
    if ( Y[no] == 3 ) Errors = Errors + 1;
    if ( Y[no] == 4 ) Unknown = Unknown + 1;
  }
  
}

