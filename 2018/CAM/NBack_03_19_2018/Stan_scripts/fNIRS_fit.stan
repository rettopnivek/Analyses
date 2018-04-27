data {
  int No; // Number of observations
  int Ns; // Number of subjects
  int Nroi; // Number of ROI
  int Nc; // Number of regression coefficients
  int indexSubject[No]; // Subject index over observations
  int indexROI[No]; // ROI index for channels over observations
  matrix[Nroi,Nc] X; // Design matrix
  real zHbO[No]; // Standardized BOLD estimates
  real nu[No]; // Degrees of freedom for ROI values
  vector[Nc] priors_beta[2];
  real priors_Omega;
  vector[Nroi] priors_Tau[2];
}
parameters {
  // ROI estimates per subject
  row_vector[Nroi] ROI_subjects[Ns];
  real<lower=0.0> sigma_channels[Ns];
  // Parameters for multivariate normal distribution 
  // for ROI over subjects
  vector[Nc] beta; // Regression coefficients for mean
  corr_matrix[Nroi] Omega;
  vector<lower=0.0>[Nroi] Tau;
}
model {
  // Variable declarations
  real mu_zHbO[No];
  real sigma_zHbO[No];
  
  // Loop over observations
  for (no in 1:No ) {
    // Mean for student-t distribution
    mu_zHbO[no] = 
      ROI_subjects[ indexSubject[ no ], indexROI[ no ] ];
    // Standard deviation for student-t distribution
    sigma_zHbO[no] = sigma_channels[ indexSubject[ no ] ];
  }
  
  // Priors
  beta ~ normal( priors_beta[1], priors_beta[2] );
  Omega ~ lkj_corr( priors_Omega );
  Tau ~ cauchy( priors_Tau[1], priors_Tau[2] );
  
  // Hierarchy
  for ( ns in 1:Ns ) 
    ROI_subjects[ns] ~ multi_normal( X * beta, 
      quad_form_diag(Omega,Tau) );
  
  // Likelihood
  
  // Individual ROI follow a student t-distribution
  // for robustness
  zHbO ~ student_t( nu, mu_zHbO, sigma_zHbO );
  
}

