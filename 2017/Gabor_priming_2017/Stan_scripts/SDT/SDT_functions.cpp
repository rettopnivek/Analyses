#include <Rcpp.h>

void my_fill_vector( Rcpp::NumericVector v_to, 
                     Rcpp::NumericVector v_from, 
                     int start, int end ) {
  /*
  Purpose:
  Fills a vector with values from another vector 
  within a given range
  Arguments:
  v_to   - The vector to fill
  v_from - The vector to fill from
  start  - The index to start from for v_from
  end    - The index to end at for v_from
  Returns:
  A side-effect, where v_to is updated with new 
  values.
  */
  
  int inc = 0;
  for ( int i = start - 1; i < end; i++ ) {
    v_to[inc] = v_from[i];
    inc++;
  }
  
}

double my_dot_product( Rcpp::NumericVector x, Rcpp::NumericVector y ) {
  /*
  Purpose:
  A function to compute the dot product for two vectors of matching
  length.
  Arguments:
  x - The row vector
  y - The column vector
  Returns:
  The sum of the element-wise products of the two vectors.
  */
  
  double out;
  
  out = sum( x * y );
  
  return out;
}

// [[Rcpp::export]]
double SDT_prob( double dp, double crt, double Co ) {
  /*
  Purpose:
  Calculates the probability of picking right given a 
  signal detection theory comparison process.
  Arguments:
  dp  - A d' value
  crt - A criterion value (positive values indicate 
        a bias to the left)
  Co  - The position of the correct choice ( 0 = left)
  Returns:
  The probability of picking right conditioned on the 
  position of the correct choice.
  */
  
  // Variable declarations
  double theta;
  double P_RR;
  double P_RL;
  
  // Calculate the probability of picking right given the 
  // correct answer is on the right.
  P_RR = 1.0 - R::pnorm( crt, dp/2.0, 1.0, 1, 0 );
  
  // Calculate the probability of picking right given the 
  // correct answer is on the left.
  P_RL = 1.0 - R::pnorm( crt, -dp/2.0, 1.0, 1, 0 );
  
  // Select the appropriate probability for the given trial
  theta = Co*P_RR + (1.0-Co)*P_RL;
  
  return( theta );
}

// [[Rcpp::export]]
Rcpp::NumericMatrix SDT_rng( Rcpp::List L, Rcpp::NumericMatrix B ) {
  /*
  Purpose:
  
  Arguments:
  
  Returns:
  
  */
  
  // Dimensions
  int No = L[0]; // Number of observations
  int Nd = L[1]; // Number of d' values to estimate
  int Nc = L[2]; // Number of criterion values to estimate
  
  // Design matrix for d'
  Rcpp::NumericMatrix X_d = L[3];
  // Design matrix for criterion
  Rcpp::NumericMatrix X_crt = L[4];
  // Position of correct answer (0 = left, 1 = right)
  Rcpp::NumericVector Co = L[5];
  
  // Number of iterates
  int Ni = B.nrow();
  
  // Additional variable declarations
  double theta;
  double cur_dp;
  double cur_crt;
  Rcpp::NumericVector beta_dp(Nd);
  Rcpp::NumericVector beta_crt(Nc);
  Rcpp::NumericVector X_d_row(Nd);
  Rcpp::NumericVector X_crt_row(Nc);
  
  // Initialize output
  Rcpp::NumericMatrix Y(No,Ni);
  
  // Loop over iterates
  for ( int ni = 0; ni < Ni; ni++ ) {
    
    // Loop over trials
    for ( int no = 0; no < No; no++ ) {
      
      my_fill_vector( X_d_row, X_d(no, Rcpp::_ ), 1, Nd );
      my_fill_vector( X_crt_row, X_crt(no,Rcpp::_), 1, Nc );
      my_fill_vector( beta_dp, B(ni,Rcpp::_), 1, Nd );
      my_fill_vector( beta_crt, B(ni,Rcpp::_), Nd+1, Nd+Nc );
      
      cur_dp = my_dot_product( X_d_row, beta_dp );
      cur_crt = my_dot_product( X_crt_row, beta_crt );
      
      theta = SDT_prob( cur_dp, cur_crt, Co(no) );
      Y(no,ni) = R::rbinom( 1, theta );
      
    }
    
  }
  
  return Y;
}
