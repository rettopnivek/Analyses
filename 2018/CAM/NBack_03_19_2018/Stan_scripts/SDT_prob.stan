real SDT_prob( real dp, real crt, real Co ) {
  // Purpose:
  // Calculates the probability of pressing a key 
  // to indicate the presence of a target in an 
  // N-back task (i.e., the probability of 
  // making a hit or false alarm).
  // Arguments:
  // dp  - A d' value
  // crt - A criterion value (positive values indicate 
  //       a bias to the left)
  // Co  - Whether the target is present (1) or not (0).
  // Returns:
  // The probability of picking right conditioned on the 
  // position of the correct choice.
  
  // Variable declarations
  real theta;
  real P_H;
  real P_FA;
  
  // Calculate the probability of pressing a key when 
  // the target is present
  P_H = 1.0 - normal_cdf( crt, dp/2.0, 1.0 );
  
  // Calculate the probability of pressing a key when 
  // the target is absent
  P_FA = 1.0 - normal_cdf( crt, -dp/2.0, 1.0 );
  
  // Select the appropriate probability for the given trial
  theta = Co*P_H + (1-Co)*P_FA;
  
  return( theta );
}

