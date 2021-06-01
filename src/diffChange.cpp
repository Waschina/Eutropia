#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::mat diffChange(arma::mat adjmat, arma::mat nneighbors, arma::mat conc, int niter) {

  arma::mat ychg(conc.n_rows, conc.n_cols, arma::fill::zeros);

  for(int k = 0; k < niter; k++) {
    ychg = arma::mat(conc.n_rows, conc.n_cols, arma::fill::zeros);

    // It might be possible to parallelize the following with columns on individual workers

    // Inflow
    for(unsigned j = 0; j < adjmat.n_rows; j++) {
      int fn = (int) adjmat(j,0);
      int tn = (int) adjmat(j,1);
      ychg.row(tn-1) += conc.row(fn-1)/13;
    }

    // Outflow
    for(unsigned j = 0; j < nneighbors.n_rows; j++) {
      int fn = (int) nneighbors(j,0); // index of field
      double nn = (double) nneighbors(j,1); // n neighbors (mostly 13, with self-count)
      ychg.row(fn-1) -= conc.row(fn-1)*((nn-1)/13);
    }

    conc += ychg;
  }

  return conc;
}

