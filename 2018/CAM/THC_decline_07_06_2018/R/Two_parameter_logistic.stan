data {
  int<lower=1> Ni; // Number of items
  int<lower=1> Ns; // Number of subjects
  int<lower=1> No; // Number of observations
  int<lower=1, upper=Ni> item[No]; // Index for items
  int<lower=1, upper=Ns> subject[No]; // Index for subjects
  int<lower=0, upper=1> y[No]; // Observed binary response
  matrix[3,2] priors; // Prior values for alpha, beta, and theta
}
parameters {
  vector<lower=0>[Ni] alpha; // Forthcoming
  vector[Ni] beta; // Forthcoming
  vector[Ns] theta; // Forthcoming
}
model {
  // Variable declarations
  vector[No] eta;
  
  // Informative priors
  alpha ~ lognormal( priors[1,1], priors[1,2] );
  beta ~ normal( priors[2,1], priors[2,2] );
  theta ~ normal( priors[3,1], priors[3,2] );
  
  // Loop over observations
  for (n in 1:No)
    eta[n] = alpha[item[n]] * ( theta[subject[n]] - beta[item[n]] );
  // Log-likelihood
  y ~ bernoulli_logit(eta);
}
