real SDT_prob( real dp, real crt, real Co ) {
  // Purpose:
  // Calculates the probability of picking right given a 
  // signal detection theory comparison process.
  // Arguments:
  // dp  - A d' value
  // crt - A criterion value (positive values indicate 
  //       a bias to the left)
  // Co  - The position of the correct choice ( 0 = left)
  // Returns:
  // The probability of picking right conditioned on the 
  // position of the correct choice.
  
  // Variable declarations
  real theta;
  real P_RR;
  real P_RL;
  
  // Calculate the probability of picking right given the 
  // correct answer is on the right.
  P_RR = 1.0 - normal_cdf( crt, dp/2.0, 1.0 );
  
  // Calculate the probability of picking right given the 
  // correct answer is on the left.
  P_RL = 1.0 - normal_cdf( crt, -dp/2.0, 1.0 );
  
  // Select the appropriate probability for the given trial
  theta = Co*P_RR + (1-Co)*P_RL;
  
  return( theta );
}


