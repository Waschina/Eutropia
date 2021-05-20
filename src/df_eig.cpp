#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;
using Eigen::MappedSparseMatrix;
using Eigen::SparseMatrix;

//' Diffuse Eigen
//'
//' @param A Matrix A.
//' @param B Matrix B
//' @param niter
//'
// [[Rcpp::export]]
Eigen::SparseMatrix<double> diffuse_eigen(Eigen::SparseMatrix<double>& A,
                                          Eigen::SparseMatrix<double>& B,
                                          int niter) {
  //SparseMatrix<double> AA(A);
  //SparseMatrix<double> BB(B);

  for (int i = 0; i < niter; i++) {
    B = A * B;
  }

  return B;
}
