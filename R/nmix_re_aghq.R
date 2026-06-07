# =============================================================================
# Single-species N-mixture random effects via the shared AGHQ engine
# (gcol33/tulpaObs#13).
#
# The community-N-mixture path (R/nmix_laplace_re.R) routes a native compiled
# oracle (NMixCommunityOracle) through tulpa::tulpa_re_aghq(): the grouping
# factor is species, the per-group RE is the FULL coefficient vector with a
# Gaussian community prior, and theta is the community mean. That is the WRONG
# shape for plain GLMM-style RE on a single species: the grouping factor is
# something else (station, observer-per-site, site cluster), the per-group RE
# is a SUBSET of coefficients on ONE arm (intercept-only or one slope), and
# there is no "community mean" -- theta is the plain fixed-effect vector.
#
# This file builds the SINGLE-SPECIES grouped-RE oracle (NMixGroupedOracle,
# src/nmix_re_oracle.{h,cpp}) and routes tulpa::tulpa_re_aghq() through it as a
# prebuilt native oracle. The oracle owns:
#   * Site-level RE on either the abundance arm (lambda) or the detection arm
#     (p). On the p arm Z is still site-level: the engine-supplied per-site
#     scalar offset is applied UNIFORMLY to all of that site's visits, leaving
#     the per-site marginal a function of one scalar offset per site -- the
#     per-row separability the engine contract needs. Visit-level p RE where
#     observers vary WITHIN a site would couple multiple b's through the
#     shared latent N and breaks that factorization; it would need a
#     different engine, not just a different oracle. The data path is
#     unreachable from the public API today (det_visit_formula doesn't parse
#     structured terms), so the gate is the structural restriction.
#   * Poisson or negative binomial. Under NB the global dispersion log_r is
#     the (p_lambda + p_p + 1)-th theta entry; the per-site kernel already
#     accepts r, so rebind(theta) picks it up and the engine threads its
#     theta-gradient through the oracle's theta_score.
#
# Scope. One shared grouping factor across all RE terms (the per-group integral
# factorizes only then; the engine returns NULL otherwise and the caller keeps
# the no-RE fit). RE total dimension per group <= 3 (the engine grid is
# n_quad^dim). RE arm split across BOTH arms is rejected upstream (split-arm
# fits would need two simultaneous oracles).
#
# Production speed. The native oracle owns the per-site marginal and the
# per-group site list, so each (g, b) call touches only that group's sites
# and stays in C++. n_quad = 1 is the joint Laplace (glmer nAGQ = 1); n_quad
# > 1 is the AGHQ debias -- the same shape as the occupancy / detection RE
# path. At the recovery suite scale (N = 100, J = 4, 10 groups) a Poisson
# fit takes ~2 s, NB ~6-10 s (probe: dev_notes/probe_nmix_re_oracle_speedup.R).
# =============================================================================


# AGHQ refinement (and n_quad = 1 fit) of a single-species COUNT model with
# grouped random effects. Family-agnostic: `make_oracle` is the native
# REGroupOracle factory (cpp_nmix_grouped_oracle / cpp_removal_grouped_oracle,
# identical signatures); everything else -- the applicability gate, the
# covariance-only re_terms spec, the site-level Z assembly, the AGHQ call, and
# the BLUP re-pack -- is shared across families. `model` is the autoscaled
# `tobs_model` whose `X_processes` are (X_lambda, X_p) and whose `y_long` /
# `site_idx` are the long-form counts (N-mixture visits or removal passes);
# `design` the per-term RE design (.tobs_re_design output, each tagged with its
# `arm`: "lambda" if the random effect enters the abundance predictor, "p" if it
# enters detection); `beta_lambda` / `beta_p` / `Sigma_list` / `b` the warm
# starts; `mixture` the abundance mixing distribution code; `K_max` the
# marginal-sum truncation (resolved by the family wrapper). Returns a list of
# refined estimates or NULL when the engine declines (caller falls back to the
# no-RE fit).
.tobs_count_re_aghq <- function(make_oracle, model, design, beta_lambda, beta_p,
                                Sigma_list, b, mixture = "P", r_init = 10,
                                K_max, n_quad = 1L, lkj_eta = 1.5,
                                theta_prior_sd = 100, max_iter = 200L,
                                verbose = FALSE) {
  # ---- applicability: single arm, one shared grouping factor, RE dim <= 3 --
  arm <- unique(vapply(design, function(d) d$arm %||% "lambda", character(1)))
  if (length(arm) != 1L) return(NULL)
  idx1 <- as.integer(design[[1]]$idx)
  ng   <- as.integer(design[[1]]$n_groups)
  one_group <- all(vapply(design, function(d)
    identical(as.integer(d$idx), idx1) &&
      identical(as.integer(d$n_groups), ng), logical(1)))
  if (!one_group) return(NULL)
  dtot <- sum(vapply(design, function(d) as.integer(d$n_coefs), integer(1)))
  if (dtot > 3L) return(NULL)

  X_lambda <- model$X_processes[[1]]
  X_p      <- model$X_processes[[2]]
  y_long   <- as.integer(model$y_long)
  site_idx <- as.integer(model$site_idx)
  N <- nrow(X_lambda)
  p_lam <- ncol(X_lambda); p_p <- ncol(X_p)

  # Engine wants n_obs = number of "rows". Rows are sites here: the per-site
  # marginal is per-row in the engine's sense, with each site belonging to one
  # group via idx (length N).
  if (any(vapply(design, function(d) length(d$idx) != N, logical(1))) ||
      any(vapply(design, function(d) nrow(d$Z) != N, logical(1)))) {
    return(NULL)  # non-site-level design: visit-level p RE etc. not factorizable
  }

  K_max <- as.integer(K_max)
  is_nb <- identical(mixture, "NB")

  # ---- per-term RE spec for the engine -----------------------------------
  # On the oracle path re_terms is COVARIANCE-ONLY (idx = Z = NULL): the per-
  # row design lives in the oracle, and the engine reads only the Sigma block
  # structure (n_coefs / correlated / n_groups) plus the shared-grouping-factor
  # check (identical(NULL, NULL) -> TRUE, so only n_groups has to match across
  # terms). nested_laplace_re_cov.R:90 documents this form.
  re_terms <- lapply(design, function(d) list(
    n_groups = as.integer(d$n_groups),
    n_coefs = as.integer(d$n_coefs),
    correlated = isTRUE(d$correlated)))

  # Site-level RE design (N x dtot): cbind of each term's Z, with intercept-only
  # blocks expanded to a column of 1s (matching .re_cov_block_layout's default).
  Z_site <- do.call(cbind, lapply(design, function(d) {
    if (is.null(d$Z)) matrix(1, N, as.integer(d$n_coefs))
    else              as.matrix(d$Z)
  }))
  if (ncol(Z_site) != dtot) {
    stop(sprintf(".tobs_nmix_re_aghq: assembled Z_site has %d columns, ",
                 "expected dtot = %d.", ncol(Z_site), dtot), call. = FALSE)
  }

  orc <- make_oracle(
    arm = if (arm == "lambda") 0L else 1L,
    y = y_long, site_idx = site_idx,
    X_lambda = X_lambda, X_p = X_p,
    Z_site = Z_site, site_group = idx1,
    n_sites = N, n_groups = ng, K_max = K_max,
    nb = is_nb)

  theta0 <- c(as.numeric(beta_lambda), as.numeric(beta_p))
  if (is_nb) theta0 <- c(theta0, log(r_init))

  ref <- tulpa::tulpa_re_aghq(
    theta0 = theta0, re_terms = re_terms, Sigma0 = Sigma_list,
    oracle = orc,
    n_quad = as.integer(n_quad), lkj_eta = lkj_eta,
    theta_prior_sd = theta_prior_sd, maxit = as.integer(max_iter))
  if (is.null(ref)) return(NULL)

  beta_lambda_ref <- ref$theta[seq_len(p_lam)]
  beta_p_ref      <- ref$theta[p_lam + seq_len(p_p)]
  r_ref           <- if (is_nb) exp(ref$theta[p_lam + p_p + 1L]) else NA_real_
  log_r_ref       <- if (is_nb) ref$theta[p_lam + p_p + 1L]     else NA_real_

  # Re-pack the per-term BLUP matrices (n_groups x n_coefs) into the term-major,
  # group-major (byrow) layout .tobs_re_param_block() / .tobs_re_offset() use.
  b_out    <- unlist(lapply(ref$blup,     function(M) as.numeric(t(M))),
                     use.names = FALSE)
  bvar_out <- unlist(lapply(ref$blup_var, function(M) as.numeric(t(M))),
                     use.names = FALSE)

  list(
    ok            = TRUE,
    arm           = arm,
    mixture       = mixture,
    beta_lambda   = beta_lambda_ref,
    beta_p        = beta_p_ref,
    log_r         = log_r_ref,
    r             = r_ref,
    vcov          = ref$theta_cov,
    theta_se      = ref$theta_se,
    Sigma_list    = ref$Sigma_list,
    b             = b_out,
    b_var         = bvar_out,
    log_marginal  = ref$log_marginal,
    n_quad        = ref$n_quad,
    lkj_eta       = ref$lkj_eta,
    converged     = ref$converged
  )
}


# N-mixture wrapper: the native nmix grouped oracle + the N-mixture K_max default
# (largest single count + 100). Signature kept stable for .tobs_fit_nmix_re.
.tobs_nmix_re_aghq <- function(model, design, beta_lambda, beta_p,
                               Sigma_list, b, mixture = "P", r_init = 10,
                               K_max = NULL, n_quad = 1L, lkj_eta = 1.5,
                               theta_prior_sd = 100, max_iter = 200L,
                               verbose = FALSE) {
  if (is.null(K_max)) K_max <- as.integer(max(as.integer(model$y_long)) + 100L)
  .tobs_count_re_aghq(cpp_nmix_grouped_oracle, model, design,
                      beta_lambda, beta_p, Sigma_list, b, mixture, r_init,
                      K_max, n_quad, lkj_eta, theta_prior_sd, max_iter, verbose)
}


# Removal wrapper: the native removal grouped oracle + the removal K_max default
# (largest per-site removal TOTAL + 100; depletion sums over passes, so the
# truncation must clear the per-site total, not the per-pass max). The per-site
# kernel depletes the available count per pass; the RE enters one arm exactly as
# in the N-mixture, so the shared AGHQ helper applies verbatim.
.tobs_removal_re_aghq <- function(model, design, beta_lambda, beta_p,
                                  Sigma_list, b, mixture = "P", r_init = 10,
                                  K_max = NULL, n_quad = 1L, lkj_eta = 1.5,
                                  theta_prior_sd = 100, max_iter = 200L,
                                  verbose = FALSE) {
  if (is.null(K_max)) {
    site_tot <- tapply(as.integer(model$y_long),
                       factor(as.integer(model$site_idx),
                              levels = seq_len(model$n_sites)), sum)
    site_tot[is.na(site_tot)] <- 0L
    K_max <- as.integer(max(as.integer(site_tot)) + 100L)
  }
  .tobs_count_re_aghq(cpp_removal_grouped_oracle, model, design,
                      beta_lambda, beta_p, Sigma_list, b, mixture, r_init,
                      K_max, n_quad, lkj_eta, theta_prior_sd, max_iter, verbose)
}


# Shared grouped-RE fit for a single-species count model (N-mixture / removal /
# ...): warm-start the betas with the no-RE Laplace fit (`warm_fun`), seed Sigma
# at a diagonal 0.25 per term, refine through the family AGHQ helper (`aghq_fun`,
# one of .tobs_nmix_re_aghq / .tobs_removal_re_aghq), and assemble the `raw`
# object build_nmix_fit consumes. `family_label` only flavours the error text;
# the structure is identical because the per-site marginals share the
# (lambda, p[, log_r]) coefficient layout and the REGroupOracle interface.
.tobs_fit_count_re <- function(model, re, warm_fun, aghq_fun, family_label,
                               mixture = "P", K_max = NULL,
                               max_iter = 100L, tol = 1e-6, verbose = TRUE,
                               n_quad = 1L, lkj_eta = 1.5, theta_prior_sd = 100) {
  if (inherits(re, "tobs_re")) re <- list(re)
  arms <- .tobs_re_split_two_arms(
    re, model, "lambda", "p",
    sprintf(paste0("A random effect shared across the abundance (lambda) and ",
                   "detection (p) arms is not supported on the %s AGHQ path."),
            family_label))
  # v1: single arm. Both populated -> reject (cross-arm AGHQ needs a joint
  # two-arm oracle, a separate engine path).
  if (length(arms$lambda) && length(arms$p)) {
    stop(sprintf(paste0("Random effects on BOTH the abundance and detection ",
                        "arms in one %s fit are not yet supported; the AGHQ ",
                        "path integrates one arm at a time. Put the RE on ",
                        "lambda OR p, not both."), family_label), call. = FALSE)
  }
  design <- if (length(arms$lambda)) arms$lambda else arms$p

  X_lambda <- model$X_processes[[1]]
  X_p      <- model$X_processes[[2]]
  y_long   <- model$y_long

  warm <- tryCatch(
    warm_fun(model, mixture = mixture, K_max = K_max,
             max_iter = as.integer(max_iter), tol = as.numeric(tol)),
    error = function(e) NULL)
  beta_lambda_init <- if (!is.null(warm)) warm$beta_lambda
                      else c(log(max(mean(y_long), 0.1)), rep(0, ncol(X_lambda) - 1L))
  beta_p_init      <- if (!is.null(warm)) warm$beta_p else rep(0, ncol(X_p))
  r_init <- if (!is.null(warm) && identical(mixture, "NB") &&
                is.finite(warm$r %||% NA_real_)) as.numeric(warm$r) else 10

  Sigma_init <- lapply(design, function(d) diag(0.25, d$n_coefs))
  b_init <- numeric(sum(vapply(design,
                               function(d) as.integer(d$n_groups * d$n_coefs),
                               integer(1))))

  ref <- aghq_fun(model, design,
                  beta_lambda = beta_lambda_init, beta_p = beta_p_init,
                  Sigma_list = Sigma_init, b = b_init,
                  mixture = mixture, r_init = r_init, K_max = K_max,
                  n_quad = as.integer(n_quad), lkj_eta = lkj_eta,
                  theta_prior_sd = theta_prior_sd, max_iter = as.integer(max_iter),
                  verbose = isTRUE(verbose))
  if (is.null(ref) || !isTRUE(ref$ok)) {
    stop(sprintf(paste0("%s AGHQ random-effect refinement did not produce a ",
                        "usable fit (singular marginal Hessian or non-finite ",
                        "optimum). Try a different K_max or simplify the RE ",
                        "structure."), family_label), call. = FALSE)
  }

  raw <- list(
    mixture     = mixture,
    beta_lambda = ref$beta_lambda,
    beta_p      = ref$beta_p,
    log_r       = ref$log_r,
    r           = ref$r,
    vcov        = ref$vcov,
    log_lik     = ref$log_marginal,
    converged   = ref$converged,
    K_max       = K_max
  )
  re_post <- list(arm = ref$arm, design = design,
                  Sigma_list = ref$Sigma_list,
                  b = ref$b, b_var = ref$b_var,
                  n_quad = ref$n_quad, lkj_eta = ref$lkj_eta)
  build_nmix_fit(raw, model, spatial = NULL, re_post = re_post)
}
