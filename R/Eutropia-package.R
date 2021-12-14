## usethis namespace: start
#' @importFrom Rcpp sourceCpp
#' @useDynLib Eutropia, .registration = TRUE
## usethis namespace: end
NULL


.onAttach <- function(libname, pkgname) {
  packageStartupMessage("Welcome to Eutropia.")
}
