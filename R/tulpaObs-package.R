#' @keywords internal
"_PACKAGE"

#' @useDynLib tulpaObs, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL

.onLoad <- function(libname, pkgname) {
  # Ensure Matrix's CHOLMOD stubs are available before tulpa's DLL init
  requireNamespace("Matrix", quietly = TRUE)
  # Ensure tulpa's DLL (with the CellCouplingSpec registry + the
  # `tulpa_register_cell_coupling` R_GetCCallable entry) is loaded before
  # registering OccuCoverLognormalCoupling against it.
  requireNamespace("tulpa", quietly = TRUE)
  cpp_register_occu_cover_lognormal_coupling()
}
