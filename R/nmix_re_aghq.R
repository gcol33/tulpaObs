# =============================================================================
# Single-species N-mixture random effects via the shared AGHQ engine.
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
# grouped random effects. Family-agnostic: `make_oracle` is a CLOSURE
# `function(arm_code, Z_site, K_max)` that returns a native REGroupOracle XPtr,
# capturing whatever per-site data its family kernel needs (the N-mixture /
# removal long-form counts, the distance per-site bin counts + quadrature, ...).
# Everything else -- the applicability gate, the covariance-only re_terms spec,
# the site-level Z assembly, the AGHQ call, and the BLUP re-pack -- is shared
# across families and reads ONLY the RE `design` (per-term .tobs_re_design
# output, each tagged with its `arm`: "lambda" if the random effect enters the
# abundance predictor, "p"/the detection arm otherwise), the warm-start `beta_lambda`
# / `beta_p` / `Sigma_list` / `b`, the mixing-distribution code `mixture`, and the
# marginal-sum truncation `K_max` (resolved by the family wrapper). The fixed-effect
# coefficient counts come from the warm-start lengths, the number of sites from the
# design, so no model layout is assumed here. Returns a list of refined estimates
# or NULL when the engine declines (caller falls back to the no-RE fit).
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

  p_lam <- length(beta_lambda); p_p <- length(beta_p)
  N <- length(idx1)

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
    stop(sprintf(paste0(".tobs_nmix_re_aghq: assembled Z_site has %d columns, ",
                        "expected dtot = %d."), ncol(Z_site), dtot), call. = FALSE)
  }

  orc <- make_oracle(if (arm == "lambda") 0L else 1L, Z_site, K_max)

  theta0 <- c(as.numeric(beta_lambda), as.numeric(beta_p))
  if (is_nb) theta0 <- c(theta0, log(r_init))

  ref <- tulpa::tulpa_re_aghq(
    theta0 = theta0, re_terms = re_terms, Sigma0 = Sigma_list,
    oracle = orc,
    n_quad = as.integer(n_quad), lkj_eta = lkj_eta,
    theta_prior_sd = theta_prior_sd, max_iter = as.integer(max_iter))
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
  # A group the engine could not solve carries NA BLUPs into `b` / `b_var`, so
  # the fit is not reported as converged and says which groups failed
  # (gcol33/tulpaObs#281).
  gstat <- .tobs_aghq_group_status(ref)

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
    converged     = .tobs_aghq_converged(ref, gstat),
    group_ok      = gstat$group_ok,
    groups_failed = gstat$failed,
    n_iter        = .tobs_aghq_n_iter(ref),
    K_max         = K_max
  )
}


# N-mixture wrapper: the native nmix grouped oracle (long-form visit counts) +
# the N-mixture K_max default (largest single count + 100). Signature kept stable
# for .tobs_fit_nmix_re.
.tobs_nmix_re_aghq <- function(model, design, beta_lambda, beta_p,
                               Sigma_list, b, mixture = "P", r_init = 10,
                               K_max = NULL, n_quad = 1L, lkj_eta = 1.5,
                               theta_prior_sd = 100, max_iter = 200L,
                               verbose = FALSE) {
  if (is.null(K_max)) K_max <- as.integer(max(as.integer(model$y_long)) + 100L)
  is_nb <- identical(mixture, "NB")
  idx1  <- as.integer(design[[1]]$idx)
  make_oracle <- function(arm_code, Z_site, K_max)
    cpp_nmix_grouped_oracle(
      arm = arm_code, y = as.integer(model$y_long),
      site_idx = as.integer(model$site_idx),
      X_lambda = model$X_processes[[1]], X_p = model$X_processes[[2]],
      Z_site = Z_site, site_group = idx1,
      n_sites = length(idx1), n_groups = as.integer(design[[1]]$n_groups),
      K_max = as.integer(K_max), nb = is_nb)
  .tobs_count_re_aghq(make_oracle, model, design,
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
  is_nb <- identical(mixture, "NB")
  idx1  <- as.integer(design[[1]]$idx)
  make_oracle <- function(arm_code, Z_site, K_max)
    cpp_removal_grouped_oracle(
      arm = arm_code, y = as.integer(model$y_long),
      site_idx = as.integer(model$site_idx),
      X_lambda = model$X_processes[[1]], X_p = model$X_processes[[2]],
      Z_site = Z_site, site_group = idx1,
      n_sites = length(idx1), n_groups = as.integer(design[[1]]$n_groups),
      K_max = as.integer(K_max), nb = is_nb)
  .tobs_count_re_aghq(make_oracle, model, design,
                      beta_lambda, beta_p, Sigma_list, b, mixture, r_init,
                      K_max, n_quad, lkj_eta, theta_prior_sd, max_iter, verbose)
}


# Distance wrapper: the native distance grouped oracle (per-site bin counts +
# the per-fit detection quadrature) + the distance K_max default (3 * largest
# per-site detected total + 100; the latent total spans the undetected
# individuals, so the truncation needs the same multiplicative margin the
# distance Laplace fit uses). Abundance-arm RE only, half-normal key only -- the
# hazard-rate key carries a global scalar shape coordinate that is not a per-site
# design column, so it is not expressible in the count-family theta layout; the
# detection arm couples a site's bins through the shared latent N, so it does not
# factorize into the per-site scalar offset the base assumes. Both are gated in
# .tobs_fit_distance_re; this wrapper assumes them.
# Insert a log_shape row/col into a [lambda | sigma | log_r?] AGHQ vcov at the
# position build_distance_fit expects (after the sigma block, before log_r), with
# the profile variance on its diagonal and zero off-diagonals (the profile-
# likelihood approximation: the shape uncertainty is treated as independent of the
# beta block, which the outer profile does not couple). Keeps the vcov dimension
# in step with the coefficient vector so build_distance_fit does not blank it.
.tobs_distance_insert_shape_vcov <- function(V, p_lam, p_sig, is_nb, var_shape) {
  V <- as.matrix(V); n <- nrow(V); pos <- p_lam + p_sig; ins <- pos + 1L
  old_to_new <- c(seq_len(pos), if (n > pos) (pos + 2L):(n + 1L))
  Vn <- matrix(0, n + 1L, n + 1L)
  Vn[old_to_new, old_to_new] <- V
  Vn[ins, ins] <- if (is.finite(var_shape)) var_shape else NA_real_
  Vn
}

.tobs_distance_re_aghq <- function(model, design, beta_lambda, beta_p,
                                   Sigma_list, b, mixture = "P", r_init = 10,
                                   K_max = NULL, n_quad = 1L, lkj_eta = 1.5,
                                   theta_prior_sd = 100, max_iter = 200L,
                                   verbose = FALSE) {
  if (is.null(K_max)) {
    R_max <- if (length(model$y)) max(rowSums(model$y)) else 0L
    K_max <- as.integer(3L * R_max + 100L)
  }
  is_nb    <- identical(mixture, "NB")
  idx1     <- as.integer(design[[1]]$idx)
  y_bins   <- matrix(as.integer(model$y), nrow(model$y), ncol(model$y))
  hazard   <- identical(model$key, "hazard")
  key_code <- .dist_key_code(model$key)
  p_lam    <- ncol(model$X_processes[[1]]); p_sig <- ncol(model$X_processes[[2]])
  # Built ONCE for the whole profile -- every candidate eta_b's oracle reuses it
  # instead of rebuilding the Gauss-Legendre quadrature per candidate
  # (cutpoints/transect/quad_order don't vary with eta_b).
  quad_xptr <- cpp_distance_build_quad(as.numeric(model$cutpoints),
                                       .dist_transect_code(model$transect),
                                       as.integer(model$quad_order))

  # A make_oracle closure at a FIXED hazard shape eta_b (0 for the half-normal
  # key). Profiling rebuilds the closure per candidate eta_b.
  make_oracle_eb <- function(eta_b) function(arm_code, Z_site, K_max)
    cpp_distance_grouped_oracle(
      arm = arm_code, y_bins = y_bins,
      X_lambda = model$X_processes[[1]], X_sigma = model$X_processes[[2]],
      Z_site = Z_site, site_group = idx1,
      n_sites = length(idx1), n_groups = as.integer(design[[1]]$n_groups),
      quad_xptr = quad_xptr,
      K_max = as.integer(K_max), nb = is_nb, key = key_code, eta_b = eta_b)
  fit_at <- function(eta_b)
    .tobs_count_re_aghq(make_oracle_eb(eta_b), model, design,
                        beta_lambda, beta_p, Sigma_list, b, mixture, r_init,
                        K_max, n_quad, lkj_eta, theta_prior_sd, max_iter, verbose)

  if (!hazard) return(fit_at(0))

  # Hazard-rate key: the scalar log-shape eta_b is not a per-site design column, so
  # rather than a second global theta slot it is PROFILED over the AGHQ log-
  # marginal (each candidate is a full AGHQ fit at a fixed shape). Warm the shape
  # from a no-RE hazard Laplace fit, then optimise it in a window around that.
  warm_hz <- tryCatch(
    distance_laplace(y = model$y, X_lambda = model$X_processes[[1]],
                     X_sigma = model$X_processes[[2]], cutpoints = model$cutpoints,
                     key = "hazard", transect = model$transect, mixture = mixture,
                     K_max = K_max, quad_order = model$quad_order,
                     max_iter = as.integer(max_iter), tol = 1e-6, verbose = FALSE),
    error = function(e) NULL)
  eb0 <- if (!is.null(warm_hz) && is.finite(warm_hz$eta_b %||% NA_real_))
    as.numeric(warm_hz$eta_b) else log(2)
  prof <- function(eta_b) { r <- tryCatch(fit_at(eta_b), error = function(e) NULL)
                            if (is.null(r) || !is.finite(r$log_marginal)) -Inf
                            else r$log_marginal }
  opt <- stats::optimize(function(eb) -prof(eb), interval = eb0 + c(-1.2, 1.2))
  eb_hat <- opt$minimum
  ref <- fit_at(eb_hat)
  if (is.null(ref)) return(NULL)
  # Profile SE for the log-shape: FD curvature of the profile log-marginal.
  h <- 0.06; lp0 <- -opt$objective
  curv <- -(prof(eb_hat + h) - 2 * lp0 + prof(eb_hat - h)) / h^2
  var_eb <- if (is.finite(curv) && curv > 0) 1 / curv else NA_real_
  ref$vcov  <- .tobs_distance_insert_shape_vcov(ref$vcov, p_lam, p_sig, is_nb, var_eb)
  ref$eta_b <- eb_hat; ref$shape <- exp(eb_hat)
  ref
}


# Default `build_fun` for .tobs_fit_count_re: assemble the `raw` object
# build_nmix_fit consumes (the (lambda, p[, log_r]) coefficient layout the
# N-mixture and removal families share) and pack it with the appended RE block.
.tobs_count_re_build_nmix <- function(ref, model, design, K_max, mixture) {
  raw <- list(
    mixture     = mixture,
    beta_lambda = ref$beta_lambda,
    beta_p      = ref$beta_p,
    log_r       = ref$log_r,
    r           = ref$r,
    vcov        = ref$vcov,
    log_lik     = ref$log_marginal,
    converged   = ref$converged,
    group_ok    = ref$group_ok,
    groups_failed = ref$groups_failed,
    n_iter      = ref$n_iter,
    K_max       = ref$K_max %||% K_max)
  re_post <- list(arm = ref$arm, design = design,
                  Sigma_list = ref$Sigma_list,
                  b = ref$b, b_var = ref$b_var,
                  n_quad = ref$n_quad, lkj_eta = ref$lkj_eta)
  build_nmix_fit(raw, model, spatial = NULL, re_post = re_post)
}


# Shared grouped-RE fit for a single-species count model (N-mixture / removal /
# distance / ...): warm-start the betas with the no-RE Laplace fit (`warm_fun`),
# seed Sigma at a diagonal 0.25 per term, refine through the family AGHQ helper
# (`aghq_fun`, one of .tobs_nmix_re_aghq / .tobs_removal_re_aghq /
# .tobs_distance_re_aghq), and pack the refined fit through `build_fun`. The
# applicability gate, warm-start, Sigma/BLUP seeding, and AGHQ call are identical
# across families because the per-site marginals share the
# (lambda, <det>[, log_r]) coefficient layout and the REGroupOracle interface;
# `det_arm` names the detection process ("p" for N-mixture / removal, "sigma" for
# distance) and `build_fun(ref, model, design, K_max, mixture)` owns the
# family-specific fit packing. `family_label` flavours the error text.
.tobs_fit_count_re <- function(model, re, warm_fun, aghq_fun, family_label,
                               mixture = "P", K_max = NULL,
                               max_iter = 100L, tol = 1e-6, verbose = TRUE,
                               n_quad = 1L, lkj_eta = 1.5, theta_prior_sd = 100,
                               det_arm = "p",
                               build_fun = .tobs_count_re_build_nmix) {
  if (inherits(re, "tobs_re")) re <- list(re)
  arms <- .tobs_re_split_two_arms(
    re, model, "lambda", det_arm,
    sprintf(paste0("A random effect shared across the abundance (lambda) and ",
                   "detection (%s) arms is not supported on the %s AGHQ path."),
            det_arm, family_label))
  # v1: single arm. Both populated -> reject (cross-arm AGHQ needs a joint
  # two-arm oracle, a separate engine path).
  if (length(arms$lambda) && length(arms[[det_arm]])) {
    stop(sprintf(paste0("Random effects on BOTH the abundance and detection ",
                        "arms in one %s fit are not yet supported; the AGHQ ",
                        "path integrates one arm at a time. Put the RE on ",
                        "lambda OR %s, not both."), family_label, det_arm),
         call. = FALSE)
  }
  design <- if (length(arms$lambda)) arms$lambda else arms[[det_arm]]

  X_lambda <- model$X_processes[[1]]
  X_p      <- model$X_processes[[2]]
  y_long   <- model$y_long

  warm <- tryCatch(
    warm_fun(model, mixture = mixture, K_max = K_max,
             max_iter = as.integer(max_iter), tol = as.numeric(tol)),
    error = function(e) NULL)
  beta_lambda_init <- if (!is.null(warm)) warm$beta_lambda
                      else c(log(max(mean(y_long %||% 1), 0.1)),
                             rep(0, ncol(X_lambda) - 1L))
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

  build_fun(ref, model, design, K_max, mixture)
}
