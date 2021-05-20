#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;

//' Diffuse
//'
//' @param A Matrix A.
//' @param B Matrix B
//' @param niter
//'
// [[Rcpp::export]]
arma::sp_mat diffuse_arma(arma::sp_mat A, arma::sp_mat B, int niter) {

  for (int i = 0; i < niter; i++) {
    B = A * B;
  }

  return B;
}

//' transpose a sparse matrix
//'
//' @param X Matrix A.
//'
// [[Rcpp::export]]
arma::sp_mat trans_(arma::sp_mat X) {
  return trans(X);
}
