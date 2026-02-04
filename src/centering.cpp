#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// Double-centering function for matrices
// [[Rcpp::export]]
arma::mat centering(const arma::mat& A) {
  int n = A.n_rows;
  int m = A.n_cols;

  arma::rowvec col_means = mean(A, 0);
  arma::colvec row_means = mean(A, 1);
  double grand_mean = accu(A) / (n * m);

  arma::mat result = A;
  result.each_col() -= row_means;
  result.each_row() -= col_means;
  result += grand_mean;

  return result;
}
