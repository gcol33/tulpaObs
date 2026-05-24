# =============================================================================
# repro_tulpa_laplace_no_prior.R
#
# Upstream gap (gcol33/tulpa): the Laplace block fitter `tulpa_laplace()` and
# the EM driver `tulpa_em_laplace()` expose no fixed-effect prior / penalty on
# beta. The MI/Gibbs correction phases (.mi_correction / .gibbs_correction)
# refit each block by calling `tulpa_laplace()` via the `m_step_encode`
# callback, so they inherit the *unpenalised* MAP.
#
# Downstream consequence (tulpaObs): a penalized MI/Gibbs occupancy fit is
# impossible. tulpaObs's weakly-informative prior (occu_priors(), which breaks
# the psi-p identifiability ridge at small J) is applied through a *forked*
# penalized IRLS loop (.tobs_em_laplace_penalized) precisely because there is
# no place to inject it into tulpa's block fit. Hence tobs(method =
# "laplace_gibbs" / "laplace_mi") must disable the prior.
# =============================================================================

library(tulpa)

# ---- 1. API gap: no prior argument on either entry point --------------------
lap_args <- names(formals(tulpa_laplace))
em_args  <- names(formals(tulpa_em_laplace))

cat("tulpa_laplace() formals:\n  ",  paste(lap_args, collapse = ", "), "\n")
cat("tulpa_em_laplace() formals:\n  ", paste(em_args,  collapse = ", "), "\n\n")

pat <- "prior|penal|ridge|lambda|beta_mean|beta_sd"
stopifnot(!any(grepl(pat, lap_args, ignore.case = TRUE)))
stopifnot(!any(grepl(pat, em_args,  ignore.case = TRUE)))
cat("=> Neither the Laplace block fitter nor the EM driver accepts a\n",
    "   fixed-effect prior / penalty. The correction phases call\n",
    "   tulpa_laplace() to refit blocks, so MI/Gibbs cannot be penalized.\n",
    sep = "")

# ---- 2. Requested API --------------------------------------------------------
# A minimal fix would be an optional Gaussian prior on the fixed effects,
# threaded from tulpa_em_laplace(...) down to tulpa_laplace(... , beta_prior =)
# and added to the IRLS objective as  + sum((beta - mu)^2 / (2 sd^2)).
# With that hook, the MI/Gibbs correction phases would inherit the same
# regularization for free (they reuse m_step_encode / tulpa_laplace).
