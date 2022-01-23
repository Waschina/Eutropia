## usethis namespace: start
#' @importFrom Rcpp sourceCpp
#' @useDynLib Eutropia, .registration = TRUE
## usethis namespace: end
NULL

#' @importFrom utils installed.packages
.onAttach <- function(libname, pkgname) {
  packageStartupMessage("Welcome to Eutropia.")

  # choose 'best' solver
  supported_solvers <- c("cplexAPI","glpkAPI","clpAPI")
  supported_solvers <- supported_solvers[supported_solvers %in% rownames(installed.packages())]
  SYBIL_SETTINGS("SOLVER", supported_solvers[1])
}
