## usethis namespace: start
#' @importFrom Rcpp sourceCpp
#' @useDynLib Eutropia, .registration = TRUE
## usethis namespace: end
NULL


.onAttach <- function(libname, pkgname) {
  packageStartupMessage("Welcome to Eutropia.")
}

.onLoad <- function(libname, pkgname) {
  if("cplexAPI" %in% installed.packages()) {
    sybil::SYBIL_SETTINGS("SOLVER","cplexAPI")
    cat("Solver: cplexAPI\n")
  } else {
    sybil::SYBIL_SETTINGS("SOLVER","glpkAPI")
    cat("Solver: glpkAPI\n")
  }
}
