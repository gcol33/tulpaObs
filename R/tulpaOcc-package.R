#' @keywords internal
"_PACKAGE"

#' @useDynLib tulpaOcc, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL

.onLoad <- function(libname, pkgname) {
  # Ensure Matrix's CHOLMOD stubs are available before tulpa's DLL init
  requireNamespace("Matrix", quietly = TRUE)
}
