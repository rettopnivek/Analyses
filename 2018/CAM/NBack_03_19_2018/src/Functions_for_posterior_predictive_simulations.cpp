#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

/*
 Title
 Written by Kevin Potter
 email: kevin.w.potter@gmail.com
 Please email me directory if you have any questions or comments
 Last updated 2018-03-25
 
 Table of contents
 1) sdt_prob_cpp
 2) mvrnormArma
 3) sdt_mix_rng
 
*/

// 1)

double sdt_prob_cpp( double dp, double crt, double co ) {
  /*
   Purpose:
   Computes the probability of indicating the 
   presence of a target based on whether 
   a target is present or absent.
   Arguments:
   dp  - Discriminability parameter (higher values 
         indicate a greater likelihood of hits 
         and correct rejections)
   crt - Criterion parameter (positive values 
         indicate a bias against saying the 
         target is present)
   co  - Whether the target is present or not
   Returns:
   A vector of probabilities.
   */
  
  // Probability of a hit
  double prob_H = 1.0 - R::pnorm( crt, dp/2.0, 1.0, 1, 0 );
  // Probability of a false alarm
  double prob_FA = 1.0 - R::pnorm( crt, -dp/2.0, 1.0, 1, 0 );
  // Probability of indicating the presence of a target
  double theta = co * prob_H + (1.0-co) * prob_FA;
  
  return theta;
}

// 2)

arma::mat mvrnormArma(int N, arma::mat Omega, arma::vec Tau ) {
  /*
   Purpose:
   Given a correlation matrix and vector of variances, 
   generate random draws from a multivariate normal with a 
   mean vector of 0 for N subjects.
   Arguments:
   N     - The number of subjects
   Omega - The correlation matrix
   Tau   - The vector of variances
   References:
   Adapted from
   http://gallery.rcpp.org/articles/simulate-multivariate-normal/
   Returns:
   A matrix with draws from the multivariate normal distribution.
  */
  
  // Number of variables
  int Nc = Omega.n_cols;
  
  // Reconstruct covariance matrix
  arma::mat Tau_diag(Nc,Nc,arma::fill::zeros);
  Tau_diag.diag() = Tau;
  arma::mat Sigma = Tau_diag * Omega * Tau_diag;
  
  // Simulate from standard normal
  arma::mat Y = arma::randn(N, Nc);
  
  // Convert to multivariate normal draws
  return Y * arma::chol(Sigma);
}


// 3)
// [[Rcpp::export]]
Rcpp::NumericMatrix sdt_mix_rng( Rcpp::List post,
                                 int No,
                                 int Ns,
                                 Rcpp::IntegerVector Nt,
                                 Rcpp::NumericVector Co,
                                 Rcpp::IntegerVector subjIndex,
                                 Rcpp::List X,
                                 Rcpp::IntegerVector eta_pos ) {
  /*
   Purpose:
   ...
   Arguments:
   ...
   Returns:
   ...
  */
  
  // Extract group-level posterior samples
  arma::mat group_param = post[1];
  arma::cube Omega = post[2];
  arma::mat Tau = post[3];
  
  // Extract design matrices
  arma::mat Xd = X[0];
  arma::mat Xc = X[1];
  
  // Dimensions
  int Nd  = Xd.n_cols;
  int Nc = Xc.n_cols;
  
  // Number of random effects
  int Nef = eta_pos.size();
  
  // Determine number of posterior samples
  int S = group_param.n_rows;
  
  // Initialize output
  Rcpp::NumericMatrix out(No,S);
  
  // Loop over posterior samples
  for (int s = 0; s < S; s++) {
    
    // Extract correlation matrix and variances
    arma::mat Omega_cur( Nef, Nef );
    arma::vec Tau_cur(Nef);
    for ( int o_r = 0; o_r < Nef; o_r++ ) {
      Tau_cur(o_r) = Tau(s,o_r);
      for ( int o_c = 0; o_c < Nef; o_c++ ) {
        Omega_cur(o_r,o_c) = Omega(s,o_r,o_c);
      }
    }
    
    // Simulate subject level effects
    arma::mat eta = mvrnormArma( Ns, 
                                 Omega_cur, 
                                 Tau_cur );
    
    // Initialize column vectors
    arma::vec beta_dp(Nd);
    arma::vec eta_dp(Nd); eta_dp.fill(0.0);
    arma::vec beta_crt(Nc);
    arma::vec eta_crt(Nc); eta_crt.fill(0.0);
    
    // Loop over observations
    for (int no = 0; no < No; no++) {
      
      // Fill in column vectors
      int inc = 0;
      for ( int i = 0; i < Nd; i++ ) {
        beta_dp(i) = group_param(s,i);
        if ( i == eta_pos(inc) - 1 ) {
          eta_dp(i) = eta( subjIndex(no) - 1, inc );
          inc += 1;
        }
      }
      
      for ( int i = 0; i < Nc; i++ ) {
        beta_crt(i) = group_param(s,i+Nd);
        if ( i == eta_pos(inc) - 1 ) {
          eta_crt(i) = eta( subjIndex(no) - 1, inc );
          inc += 1;
        }
      }
      
      arma::vec dp = Xd.row(no) * ( beta_dp + eta_dp );
      
      arma::vec crt = Xc.row(no) * ( beta_crt + eta_crt );
      
      double theta = sdt_prob_cpp( dp(0), crt(0), Co(no) );
      
      // Simulate proportion
      out(no,s) = R::rbinom( Nt(no), theta ) / Nt(no );
      
    }
    
  }
  
  return out;
}
