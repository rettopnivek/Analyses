real mSDT_prob( real dp, real crt, real lmb, real Co ) {
  // Variable declarations
  real theta;
  real P_RR;
  real P_RL;
  
  // Calculate the probability of picking right given the 
  // correct answer is on the right.
  P_RR = lmb*( 1.0 - normal_cdf( crt, dp/2.0, 1.0 ) ) + 
    ( 1.0 - lmb )*( 1.0 - normal_cdf( crt, 0.0, 1.0 ) );
  
  // Calculate the probability of picking right given the 
  // correct answer is on the left.
  P_RL = lmb*( 1.0 - normal_cdf( crt, -dp/2.0, 1.0 ) ) + 
    ( 1.0 - lmb )*( 1.0 - normal_cdf( crt, 0.0, 1.0 ) );
  
  # Select the appropriate probability for the given trial
  theta = Co*P_RR + (1-Co)*P_RL;
  
  return( theta );
}


