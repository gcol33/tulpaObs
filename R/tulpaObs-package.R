#' @keywords internal
"_PACKAGE"

#' @useDynLib tulpaObs, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom stats nobs simulate update coef fitted predict residuals
#' @importFrom stats dnorm plogis pnorm rbinom rnorm binomial cor glm ks.test
#' @importFrom stats model.matrix optim quantile runif sd setNames var
#' @importFrom methods as
#' @importFrom utils modifyList
NULL

# `.data` is the rlang/dplyr pronoun used inside conditional dplyr/ggplot2
# helpers; register it so R CMD check's global-variable analysis stays quiet.
utils::globalVariables(".data")

.onLoad <- function(libname, pkgname) {
  # Ensure Matrix's CHOLMOD stubs are available before tulpa's DLL init
  requireNamespace("Matrix", quietly = TRUE)
  # Ensure tulpa's DLL (with the CellCouplingSpec registry + the
  # `tulpa_register_cell_coupling` R_GetCCallable entry) is loaded before
  # registering the OccuCover specs against it.
  requireNamespace("tulpa", quietly = TRUE)
  cpp_register_occu_cover_lognormal_coupling()
  cpp_register_occu_cover_beta_coupling()
  cpp_register_occu_cover_lognormal_agg_coupling()
  cpp_register_occu_cover_beta_agg_coupling()
  cpp_register_occu_only_coupling()
}
