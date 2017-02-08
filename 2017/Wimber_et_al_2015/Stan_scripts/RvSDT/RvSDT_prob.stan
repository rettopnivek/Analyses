real RvSDT_prob( real theta, real crt, real dp, real Co ) {
  // Purpose:
  // Calculates the probability of picking right given a mixture 
  // between a recall process and a SDT comparison process, where 
  // the recall process operates ony for the left choice.
  // Arguments:
  // theta - The probability of recall
  // crt   - The criterion for the SDT process
  // dp    - The d' value
  // Co    - The position of the correct response
  // Returns:
  // The probability of picking right.
  
  // Variable declarations
  real sdt;
  real out;
  real pr;
  
  // Define probability for comparison and recall processes
  sdt = 0.0;
  pr = 0.0;
  
  # If correct answer is on the right
  if ( Co == 1.0 ) {
    pr = 0.0;
    sdt = 1.0 - normal_cdf( crt, dp/2.0, 1.0 );
  }
  if ( Co == 0.0 ) {
    pr = theta;
    sdt = 1.0 - normal_cdf( crt, -dp/2.0, 1.0 );
  }
  
  // Return log-density for mixture
  out = pr + ( 1.0 - pr )*sdt;
  
  return( out );
}
