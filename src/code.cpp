#define ARMA_64BIT_WORD 1
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;
//' Multiply two sparse matrices
//'
//' @param A Matrix A.
//' @param B Matrix B
//'
//' @export
// [[Rcpp::export]]
arma::mat multiply(arma::sp_mat A, arma::mat B) {
  return A * B;
}

//' transpose a sparse matrix
//'
//' @param X Matrix A.
//'
//' @export
// [[Rcpp::export]]
arma::sp_mat trans_(arma::sp_mat X) {
  return trans(X);
}

// [[Rcpp::export]]
int trace_(arma::sp_mat X) {
  return trace(X);
}
