real RvSDT_lpmf( int y, real theta, real crt, real dp, real Co ) {
  // Purpose:
  // Calculates a mixture between a recall process and a SDT 
  // comparison process, where the recall process operates ony 
  // for the left choice.
  // Arguments:
  // y     - An observed response (0 = left, 1 = right)
  // theta - The probability of recall
  // crt   - The criterion for the SDT process
  // dp    - The d' value
  // Co    - The position of the correct response
  // Returns:
  // The log-density.
  
  // Variable declarations
  real PR;
  vector[2] lcd;
  real sdt;
  real out;
  
  // Define probability of recall
  PR = (1.0 - Co) * theta;
  
  // Define log-density for recall process
  lcd[1] = log( PR ) + bernoulli_lpmf( y | 1.0 );
  
  // Define log-density for comparison process
  sdt = Co * log( 1.0 - normal_cdf( crt | dp/2.0, 1.0 ) ) + 
    (1.0 - Co) * log( 1.0 - normal_cdf( crt | -dp/2.0, 1.0 ) );
  lcd[2] = log( 1.0 - PR ) + sdt;
  
  // Return log-density for mixture
  out = log_sum_exp( lcd );
  
  return( out );
}
