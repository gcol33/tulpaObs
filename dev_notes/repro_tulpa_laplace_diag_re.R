# =============================================================================
# repro_tulpa_laplace_diag_re.R
#
# Upstream gap (gcol33/tulpa): the Laplace random-effect path parameterises a
# multi-coefficient RE block by a *per-coefficient vector of marginal sigmas*
# (a diagonal covariance across coefficients). There is no Cholesky factor /
# off-diagonal covariance, so a correlated random intercept + slope
# `(1 + x | g)` (a 2x2 group covariance with a correlation term) is not
# representable under engine = "laplace".
#
# Evidence in tulpa source (R/fit_laplace.R, non-spatial path):
#   re_sigma_list <- lapply(re_list, function(r) {
#     s <- r$sigma; nc <- r$n_coefs %||% 1L
#     if (length(s) == 1L && nc > 1L) rep(s, nc) else s   # <- per-coef sigma
#   })
#   result <- cpp_laplace_fit_multi_re(..., re_sigma_list = re_sigma_list,
#                                      re_ncoefs = re_ncoefs, ...)
#   # ^ no correlation / Cholesky argument
#
# Downstream consequence (tulpaObs): correlated slopes error out of
# .validate_re_laplace() with a pointer to NUTS. `(x || g)` (uncorrelated) and
# `(0 + x | g)` (slope-only) work because they ARE diagonal; `(1 + x | g)` and
# `(1 + x + z | g)` are NUTS-only.
# =============================================================================

library(tulpa)

# ---- API gap: no correlation / covariance input for RE blocks ---------------
# cpp_laplace_fit_multi_re is an internal Rcpp wrapper (not exported).
cpp_fit  <- getFromNamespace("cpp_laplace_fit_multi_re", "tulpa")
lap_args <- names(formals(tulpa_laplace))
cpp_args <- names(formals(cpp_fit))

cat("tulpa_laplace() formals:\n  ", paste(lap_args, collapse = ", "), "\n")
cat("cpp_laplace_fit_multi_re() formals:\n  ",
    paste(cpp_args, collapse = ", "), "\n\n")

pat <- "chol|corr|cov|covar|sigma_mat|L_re"
stopifnot(!any(grepl(pat, cpp_args, ignore.case = TRUE)))
cat("=> The multi-RE Laplace fitter takes re_sigma_list (per-coefficient\n",
    "   marginal sigmas) + re_ncoefs, but no cross-coefficient covariance.\n",
    "   A correlated (1 + x | g) block is therefore not representable.\n",
    sep = "")

# ---- Requested API -----------------------------------------------------------
# Allow re_list[[k]] to carry a lower-triangular Cholesky factor L (n_coefs x
# n_coefs) instead of / alongside the marginal-sigma vector, and add the
# corresponding term to the joint Hessian. With that, engine = "laplace" could
# fit correlated random slopes (PQL-style variance bias aside), matching the
# bar forms the NUTS sampler already supports.
