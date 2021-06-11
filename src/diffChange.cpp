#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::mat diffChange(arma::mat adjmat, arma::mat nneighbors, arma::mat conc, arma::vec niter) {

  const double b_div = 1.0 / 13;

  arma::vec ychg(conc.n_rows, arma::fill::zeros);

  for(unsigned i = 0; i < niter.size(); i++) {

    for(int k = 0; k < niter(i); k++) {

      ychg = ychg * 0;

      // Inflow
      for(unsigned j = 0; j < adjmat.n_rows; j++) {
        ychg(((int) adjmat(j,1))-1) += conc(((int) adjmat(j,0))-1, i);
      }
      ychg = ychg * b_div;

      // Outflow
      for(unsigned j = 0; j < nneighbors.n_rows; j++) {
        ychg(((int) nneighbors(j,0))-1) -= conc(((int) nneighbors(j,0))-1, i)*((((double) nneighbors(j,1))-1)*b_div);
      }

      conc.col(i) += ychg;
    }
  }


  return conc;
}

// [[Rcpp::export]]
arma::mat diffChangeVec(arma::mat adjmat, arma::mat nneighbors, arma::Row<double> conc, int niter) {



  for(int k = 0; k < niter; k++) {
    arma::Row<double> ychg(conc.size(), arma::fill::zeros);

    // It might be possible to parallelize the following with columns on individual workers

    // Inflow
    for(unsigned j = 0; j < adjmat.n_rows; j++) {
      ychg(((int) adjmat(j,1))-1) += conc(((int) adjmat(j,0))-1)/13;
    }

    // Outflow
    for(unsigned j = 0; j < nneighbors.n_rows; j++) {
      ychg(((int) nneighbors(j,0))-1) -= conc(((int) nneighbors(j,0))-1)*((((double) nneighbors(j,1))-1)/13);
    }

    conc += ychg;
  }

  return conc;
}


