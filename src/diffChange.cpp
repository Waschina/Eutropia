#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <cmath>

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
arma::mat diffChangePar(arma::Mat<int> adjmat, arma::Mat<int> nneighbors,
                        arma::Mat<double> conc,
                        arma::Col<int> niter, arma::Col<int> indvar, int ncores) {

  const double b_div = 1.0 / 13;

  #pragma omp parallel for num_threads(ncores)
  for(unsigned i = 0; i < indvar.size(); i++) {

    arma::vec ychg(conc.n_rows, arma::fill::zeros);
    for(int k = 0; k < niter(i); k++) {

      ychg = ychg * 0;

      // Inflow
      for(unsigned j = 0; j < adjmat.n_rows; j++) {
        ychg((adjmat(j,1))-1) += conc((adjmat(j,0))-1, indvar(i)-1);
      }
      ychg = ychg * b_div;

      // Outflow
      for(unsigned j = 0; j < nneighbors.n_rows; j++) {
        ychg((nneighbors(j,0))-1) -= conc((nneighbors(j,0))-1, indvar(i)-1)*((((double) nneighbors(j,1))-1)*b_div);
      }

      conc.col(indvar(i)-1) += ychg;
    }
  }


  return conc;
}
