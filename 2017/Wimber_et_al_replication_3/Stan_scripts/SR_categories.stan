data {
  int<lower=1> No; // Total number of observations
  int<lower=1> Ns; // Total number of subjects
  int<lower=2> K; // Number of response categories
  int<lower=1,upper=K> Y[No]; // Observed responses
  int indS[No]; // Subject index per trial
  matrix[K-1,4] Priors; // Matrix of parameter values for priors
}
parameters{
  matrix[ Ns, K-1 ] beta_raw; // Subject-level coefficients
  real mu_beta[K-1]; // Group-level means
  real<lower=0> sigma_beta[K-1]; // Group-level standard deviations
}
transformed parameters {
  vector[Ns] zeros;
  matrix[ Ns, K ] beta; // Coefficients for all categories
  simplex[ K ] theta[ Ns ]; // Category probabilities
  
  zeros = rep_vector(0, Ns);
  beta = append_col( beta_raw, zeros ); // Fix 'unknown' category to 0
  for (ns in 1:Ns ) theta[ ns ] = softmax( to_vector( beta[ ns ] ) );
}
model {
  
  // Priors
  mu_beta ~ normal( col( Priors, 1 ), col( Priors, 2 ) );
  sigma_beta ~ gamma( col( Priors, 3 ), col( Priors, 4 ) );
  
  // Hierarchy
  for (i in 1:3) {
    col( beta_raw, i ) ~ normal( mu_beta[i], sigma_beta[i] );
  }
  
  // Likelihood
  for ( no in 1:No ) {
    Y[no] ~ categorical( theta[ indS[no] ] );
  }
  
}
  
