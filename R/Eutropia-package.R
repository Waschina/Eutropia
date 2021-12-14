## usethis namespace: start
#' @importFrom Rcpp sourceCpp
#' @useDynLib Eutropia, .registration = TRUE
## usethis namespace: end
NULL


.onAttach <- function(libname, pkgname) {
  packageStartupMessage("Welcome to Eutropia.")
}

.onLoad <- function(libname, pkgname) {
  if("cplexAPI" %in% utils::installed.packages()) {
    sybil::SYBIL_SETTINGS("SOLVER","cplexAPI")
    packageStartupMessage("Solver: cplexAPI")
  } else {
    sybil::SYBIL_SETTINGS("SOLVER","glpkAPI")
    packageStartupMessage("Solver: glpkAPI")
  }
}
