# =============================================================================
# Occupancy nested-Laplace path. Mirrors `.tobs_laplace` but builds a
# multi-block latent prior (spatial + temporal + iid) and routes the
# occupancy M-step block through tulpa::tulpa_nested_laplace() via the
# generic EM engine's per-block dispatcher. Single-season occupancy only;
# other occupancy variants raise a clear error.
# =============================================================================


#' Build a multi-block latent prior list for tulpa::tulpa_nested_laplace()
#'
#' Converts the tobs-level spatial / temporal / RE specs into the
#' list-of-blocks shape that `tulpa::tulpa_nested_laplace()` expects under
#' its multi-block dispatch (`.is_multi_block_prior`).
#'
#' @param spatial Optional `tobs_spatial` term (from an `icar()` / `bym2()`
#'   formula term). Only BYM2 / ICAR are wired through to the multi-block
#'   engine at present; GP / multiscale_gp / SVC are not yet supported and
#'   raise.
#' @param temporal Optional `tobs_temporal` term (from a `temporal()` formula
#'   term). Types `"ar1"`, `"rw1"`, `"rw2"`, `"iid"` are supported.
#' @param re Optional list of `tobs_re` terms (from `re()` formula terms).
#'   Only `model = "iid"` terms are converted to IID latent blocks here;
#'   correlated structures (ar1 / rw1 / rw2 on RE groups) are passed through
#'   as temporal-like blocks.
#' @param model A `tobs_model` from `.tobs_build_model()`. The structured
#'   terms carry pre-resolved index codes; `model` pins `N = model$n_sites`
#'   for single-season fits.
#'
#' @return `NULL` when no latent block is supplied; a single-block list
#'   when exactly one is supplied; a list-of-blocks otherwise. Each block
#'   is the minimal field set that the tulpa multi-block dispatch fills
#'   defaults around (`.NL_REGISTRY[[type]]$defaults`).
#'
#' @keywords internal
.tobs_to_multi_block_prior <- function(spatial = NULL, temporal = NULL,
                                       re = NULL, model) {
  if (!inherits(model, "tobs_model")) {
    stop("model must be a tobs_model object", call. = FALSE)
  }
  supported <- c("single", "integrated", "dynamic", "count")
  if (!model$model_type %in% supported) {
    stop("Nested Laplace is wired for single-season, integrated, ",
         "dynamic occupancy, and the count / relative-abundance GLMM; got ",
         "model_type = '", model$model_type, "'.",
         call. = FALSE)
  }

  # The latent prior is defined over sites (one GMRF / temporal / iid unit per
  # site). The state ("occ") M-step block has one row per site for single /
  # integrated / dynamic and one row per (site, species) for community, so each
  # state row is mapped to its site via `site_of_row`. The structured terms
  # (spatial / temporal / re) were resolved against the site-level `data`, so
  # their per-site index vectors are expanded through the same map.
  dims        <- .tobs_state_block_dims(model)
  n_sites     <- dims$n_sites
  site_of_row <- dims$site_of_row

  blocks <- list()
  if (!is.null(spatial)) {
    # A spatial bar (or two-term intercept + weighted areal field) carries one
    # OR several spatial specs (the intercept field first, then per-covariate
    # varying-coefficient fields). `.tobs_resolve_occu_spatial_fields()` returns
    # the ordered list; each becomes its own latent block (the weighted ones
    # carry a per-site `svc_weight`). A single plain field stays a length-1 list.
    sp_fields <- .tobs_resolve_occu_spatial_fields(spatial, model)
    for (sf in sp_fields) {
      blocks <- c(blocks,
                  list(.tobs_block_from_spatial(sf, n_sites, site_of_row, model)))
    }
  }
  if (!is.null(temporal)) {
    blocks <- c(blocks,
                list(.tobs_block_from_temporal(temporal, model, n_sites,
                                               site_of_row)))
  }
  if (!is.null(re)) {
    for (r in re) {
      blocks <- c(blocks,
                  list(.tobs_block_from_re(r, model, n_sites, site_of_row)))
    }
  }

  if (length(blocks) == 0L) return(NULL)
  # Always return a list-of-blocks, even for a single block, so the fit routes
  # through tulpa's multi-block dispatch (cpp_nested_laplace_multi). A length-1
  # list is numerically identical to the single-block kernel -- both build a
  # length-1 LatentBlock vector for run_multi_block_nested_laplace -- but the
  # multi-block path is the one that honours `det_prob`, the per-site detection
  # probability for the marginalized `bernoulli` state likelihood whose
  # calibrated predictive variance `predict(type = "state")` reports as
  # psi_lower / psi_upper.
  blocks
}


# Latent-field contribution to the state linear predictor at each row, read
# from the current nested fit's mode. The mode is c(beta (p), latent...), where
# the latent tail concatenates each multi-block's components. A block's eta
# contribution at row i is `field_block[idx(i)]`, with `idx` the block's
# per-row unit index (spatial_idx / temporal_idx / obs_idx) and the eta-entering
# component being the first `max(idx)` entries (BYM2 stores 2*n units: the
# combined effect that enters the predictor comes first, the structured-only
# auxiliary second). Returns a length-`n_rows` numeric offset, or `numeric(0)`
# before the first M-step has produced a mode.
.nested_eta_offset <- function(latent_prior, occ_fit, p_occ, n_rows) {
  if (is.null(occ_fit) || is.null(occ_fit$mode)) return(numeric(0))
  .nested_field_from_mode(latent_prior, occ_fit$mode, p_occ, n_rows)
}

# Grid-marginalised state linear predictor from the engine's per-cell fitted
# eta (`tulpa_nested_laplace(... )$fitted_eta`, [n_grid x N], beta + field at
# each grid cell). Returns `sum_k w_k eta[k, ]` -- exact for every prior,
# including bym2 (whose predictor the engine reconstructs with the right
# per-cell mixing scales). NULL when the engine did not return fitted_eta (older
# tulpa) so callers can fall back to the mode-based reconstruction.
.nested_eta_marginal <- function(occ_fit, n_rows) {
  fe <- occ_fit$fitted_eta
  if (is.null(fe) || !is.matrix(fe) || ncol(fe) != n_rows) return(NULL)
  w <- occ_fit$weights
  if (is.null(w) || length(w) != nrow(fe)) return(NULL)
  as.numeric(crossprod(fe, w / sum(w)))
}

# Field-aware state linear predictor for a nested-Laplace E-step.
#
# The E-step weight P(z = 1 | y) must see the latent block: the field informs
# which undetected units are occupied. Without it the EM converges to the
# fixed-effect-only fixed point, the field cannot track the data, and the inner
# nested-Laplace then fits an unconstrained field to the residual -- which shows
# up as a field invented where the truth is zero, a real field roughly doubled,
# and (through the logistic conditional-vs-marginal factor
# sqrt(1 + 0.346 sigma^2)) a state slope inflated by exactly that factor.
#
# Prefer the engine's exact per-cell fitted eta marginalised over the
# hyperparameter grid -- correct for every prior including bym2. Fall back to the
# grid-weighted mode reconstruction (exact for d_fac = 1 priors; skips bym2) when
# the engine did not return fitted_eta, and to the fixed-effect eta before the
# first M-step has produced a mode.
#
# `eta_fixed` is the state arm's X %*% beta (plus any SPDE offset the caller
# already added); `n_rows` is the state block's row count -- one per site for
# single / integrated / dynamic (psi1).
.nested_state_eta <- function(eta_fixed, latent_prior, occ_fit, p_occ, n_rows) {
  if (is.null(latent_prior)) return(eta_fixed)
  eta_marg <- .nested_eta_marginal(occ_fit, n_rows)
  if (!is.null(eta_marg)) return(eta_marg)
  lat_off <- .nested_eta_offset(latent_prior, occ_fit, p_occ, n_rows)
  if (length(lat_off) == n_rows) eta_fixed + lat_off else eta_fixed
}

# Core: latent-field eta contribution at each state row from one latent vector
# `mode_vec` = c(beta (p_occ), field...). Used both for the EM E-step offset
# (with the grid-weighted mode) and for per-grid-cell marginalised prediction
# (with each `modes[k, ]` row).
.nested_field_from_mode <- function(latent_prior, mode_vec, p_occ, n_rows) {
  if (length(mode_vec) <= p_occ) return(numeric(0))
  field  <- mode_vec[(p_occ + 1L):length(mode_vec)]
  blocks <- if (!is.null(latent_prior$type)) list(latent_prior) else latent_prior

  offset <- numeric(n_rows)
  pos <- 0L
  for (b in blocks) {
    nu  <- .nl_block_field_len(b)
    if (pos + nu > length(field)) break          # guard against layout drift
    block_field <- field[(pos + 1L):(pos + nu)]
    if (identical(b$type, "spde")) {
      # The SPDE field enters eta as the FEM projection (A u), with A re-rowed
      # onto the state rows (b$A is n_rows x n_mesh). d_fac = 1, so the per-row
      # contribution is the many-to-one (A u)_row.
      if (!is.null(b$A) && nrow(b$A) == n_rows && ncol(b$A) == nu) {
        offset <- offset + as.numeric(b$A %*% block_field)
      }
    } else {
      idx <- .nl_block_unit_idx(b)
      # Only blocks whose eta contribution is exactly `x[idx]` (d_fac = 1) are
      # added. A bym2 block's eta mixes its two components with hyperparameter-
      # dependent scales, so a first-n approximation would be wrong-scaled;
      # skipping it (but advancing past its length) keeps the predictor unbiased
      # and leaves multi-block alignment intact. A varying-coefficient (SVC)
      # areal block carries a per-row design weight (svc_weight); its eta
      # contribution is weight[row] * x[idx[row]], so the reconstruction must
      # apply the weight to match the inner Newton's compute_eta.
      if (.nl_block_exact_reconstruct(b) && length(idx) == n_rows) {
        contrib <- block_field[idx]
        if (!is.null(b$svc_weight) && length(b$svc_weight) == n_rows) {
          contrib <- contrib * as.numeric(b$svc_weight)
        }
        offset <- offset + contrib
      }
    }
    pos <- pos + nu
  }
  offset
}


# Whether a latent block's eta contribution is exactly `x[idx]` (d_fac = 1), so
# the per-row linear predictor can be reconstructed in R from the returned modes
# without replicating tulpa's per-prior mixing scales. True for icar /
# car_proper / rw1 / rw2 / ar1 / iid (see tulpa src/nested_laplace.cpp, where
# each sets `block.d_fac = 1.0`). BYM2 mixes a structured (phi) and an
# unstructured (theta) component with hyperparameter-dependent scales
# (sigma*sqrt(rho)*scale and sigma*sqrt(1-rho)), so its eta is not `x[idx]`;
# reconstructing it belongs in the engine, not here.
.nl_block_exact_reconstruct <- function(b) {
  isTRUE(b$type %in% c("icar", "car_proper", "rw1", "rw2", "ar1", "iid"))
}

# Gauss-Hermite nodes / weights (physicists', weight exp(-x^2)) via Golub-Welsch:
# the nodes are the eigenvalues of the symmetric tridiagonal Jacobi matrix
# (zero diagonal, off-diagonal sqrt(k/2)) and the weights are sqrt(pi) times the
# squared first eigenvector components. No hardcoded node tables, no external
# dependency. Used to integrate plogis(eta) over the within-cell Gaussian when
# marginalising psi.
.gauss_hermite <- function(n = 15L) {
  n <- as.integer(n)
  if (n < 2L) return(list(x = 0, w = sqrt(pi)))
  k <- seq_len(n - 1L)
  b <- sqrt(k / 2)
  J <- matrix(0, n, n)
  J[cbind(k, k + 1L)] <- b
  J[cbind(k + 1L, k)] <- b
  e   <- eigen(J, symmetric = TRUE)
  ord <- order(e$values)
  list(x = e$values[ord], w = sqrt(pi) * (e$vectors[1L, ord])^2)
}

# Posterior MEAN of psi = plogis(eta) at each state row, marginalised over BOTH
# the hyperparameter grid (weights `w`) and the within-cell Gaussian
# eta_i ~ N(fe[k, i], fev[k, i]). plogis is integrated over each cell's Gaussian
# by Gauss-Hermite, then weighted by w -- the marginalise-derived-quantities
# rule (no plug-in of the cell mode). `fev = NULL` collapses to the plug-in mean
# over cell modes (older tulpa without per-cell variance).
.nested_psi_mean <- function(fe, fev, w) {
  w <- w / sum(w)
  if (is.null(fev)) return(as.numeric(crossprod(plogis(fe), w)))
  s     <- sqrt(pmax(fev, 0))
  gh    <- .gauss_hermite(15L)
  c0    <- 1 / sqrt(pi)
  ecell <- matrix(0, nrow(fe), ncol(fe))
  for (g in seq_along(gh$x)) {
    ecell <- ecell + (gh$w[g] * c0) * plogis(fe + sqrt(2) * s * gh$x[g])
  }
  as.numeric(crossprod(ecell, w))
}

# Equal-tailed psi credible interval at each state row. The eta posterior is a
# Gaussian mixture over grid cells, F(t) = sum_k w_k Phi((t - fe[k, i]) /
# s[k, i]); plogis is monotone, so psi quantiles are plogis of the eta-mixture
# quantiles, found by bracketing the mixture CDF. A per-cell sd floor keeps the
# root-find well-posed when a cell's within-cell variance is ~0 (the mixture
# then reduces to a weighted quantile of the cell modes). Returns an
# [n_rows x length(probs)] matrix on the psi scale.
.nested_psi_quantiles <- function(fe, fev, w, probs = c(0.025, 0.975)) {
  w      <- w / sum(w)
  s      <- sqrt(pmax(fev, 0))
  n_rows <- ncol(fe)
  out    <- matrix(NA_real_, n_rows, length(probs))
  for (i in seq_len(n_rows)) {
    m   <- fe[, i]
    si  <- pmax(s[, i], 1e-6)
    cdf <- function(t) sum(w * pnorm(t, mean = m, sd = si))
    lo  <- min(m - 8 * si)
    hi  <- max(m + 8 * si)
    for (j in seq_along(probs)) {
      pp  <- probs[j]
      flo <- cdf(lo) - pp
      fhi <- cdf(hi) - pp
      q <- if (flo >= 0) lo
           else if (fhi <= 0) hi
           else stats::uniroot(function(t) cdf(t) - pp, c(lo, hi),
                               tol = 1e-6)$root
      out[i, j] <- plogis(q)
    }
  }
  out
}

# Marginalised state-level posterior of psi for a nested-Laplace fit: the
# posterior MEAN and an equal-tailed 95% credible interval at every state row,
# integrated over the hyperparameter grid (per the marginalise-derived-
# quantities rule -- a weighted summary of the eta posterior, not a plug-in of
# the grid-weighted mode). For single-season occupancy the `heldout` rows are
# the INLA NA-response prediction targets (sites with an all-missing detection
# history, interpolated by the latent field).
#
# Primary path uses the engine's per-cell fitted eta and its predictive variance
# (`occ_fit$fitted_eta` / `$fitted_eta_var`, [n_grid x N]): the per-row eta
# posterior is the Gaussian mixture sum_k w_k N(fe[k, i], fev[k, i]). psi is
# plogis(eta), so the mean is the Gauss-Hermite integral of plogis over each
# cell's Gaussian (weighted by w) and the interval is plogis of the eta-mixture
# quantiles. Exact for every prior including bym2 (the engine reconstructs eta
# and its variance with the right per-cell mixing scales).
#
# When the engine returns `fitted_eta` but not `fitted_eta_var` (older tulpa),
# only the marginalised mean is reported and the interval columns are NA. When
# `fitted_eta` is absent too, eta is reconstructed from the modes -- exact only
# for d_fac = 1 priors (icar / car / temporal / iid); a bym2 fit then returns
# NULL so `predict(type = "state")` reports it is unavailable.
#
# Calibration note. The occupancy M-step is an M-inflated pseudo-binomial
# (n_trials = M), so the per-cell variance is the working-model variance. At a
# held-out node (no likelihood term) it is dominated by the prior conditional
# (kriging) variance, which is M-independent, so held-out intervals are
# calibrated. Observed-site variance is M-deflated (over-confident); those rows
# carry the point estimate, and their interval is conditional on the working
# model. Computed on the (autoscaled) fitting design; psi is on the probability
# scale and so is invariant to the centering / scaling of X.
# Exact-marginal occupancy state pass. The EM converges via the M-inflated
# pseudo-binomial M-step (good for the mode and the detection estimate), but its
# field curvature and grid weights are M-distorted -- the occ M-step weights the
# data ~M times the prior, and the unit-trial Hessian is the EM complete-data
# information, which overstates the marginal information (Louis' identity).
# So, exactly as `.tobs_occu_marginal_refine()` does for the non-spatial fixed
# effects, refine the state field with one nested-Laplace pass on the TRUE
# marginalized state likelihood: each site contributes a Bernoulli on
# D_i = 1{>=1 detection} with mean q_i * sigma(eta), where q_i is the per-site
# probability of detecting at least once given occupancy, read off the converged
# detection estimate. The latent occupancy state is integrated out, so this pass
# carries the calibrated marginal mode, curvature (fitted_eta_var) and grid
# weights with no M-inflation. Held-out sites have no valid visits -> q_i = 0,
# so they drop from the likelihood and are interpolated by the field (the INLA
# NA-response mechanism, now without the n_trials = 0 hack). Detection is plugged
# in at its point estimate (the field's dominant uncertainty is the field).
.tobs_occu_state_marginal_fit <- function(model, em_result, latent_prior,
                                          max_iter = 50L, tol = 1e-6,
                                          n_threads = 1L) {
  X_occ       <- model$X_processes[[1]]
  X_det       <- model$X_processes[[2]]
  X_det_visit <- model$X_det_visit
  y           <- model$y
  n_sites     <- model$n_sites
  max_visits  <- ncol(y)
  p_det       <- ncol(X_det)

  valid_mat <- y >= 0
  n_valid   <- rowSums(valid_mat)
  D_i       <- as.integer(rowSums(valid_mat & (y == 1L), na.rm = TRUE) > 0)

  beta_det <- extract_beta(em_result$fits$det,
                           p_det + (if (is.null(X_det_visit)) 0L
                                    else ncol(X_det_visit)))

  # q_i = P(>=1 detection | occupied) at the converged detection estimate.
  if (is.null(X_det_visit)) {
    p_site <- plogis(as.vector(X_det %*% beta_det[seq_len(p_det)]))
    q_i    <- 1 - (1 - p_site)^n_valid
  } else {
    eta_site  <- as.vector(X_det %*% beta_det[seq_len(p_det)])
    eta_visit <- as.vector(X_det_visit %*% beta_det[(p_det + 1L):length(beta_det)])
    logit_p   <- matrix(eta_site, n_sites, max_visits) +
                 matrix(eta_visit, n_sites, max_visits, byrow = TRUE)
    logit_p   <- .tobs_clamp_eta(logit_p)
    log_1mp   <- -(pmax(logit_p, 0) + log1p(exp(-abs(logit_p))))
    log_1mp[!valid_mat] <- 0
    q_i <- 1 - exp(rowSums(log_1mp))
  }
  q_i[n_valid == 0L] <- 0                     # held-out: no information

  # The marginalized occupancy state likelihood is owned by tulpaObs (a scaled
  # Bernoulli, latent state integrated out); build it as a tulpa LikelihoodSpec
  # over (D_i, q_i) and route it through the nested grid via `likelihood`.
  lik <- occ_make_nested_likelihood(y = as.numeric(D_i), det_prob = q_i)
  # tulpa keeps statistical args top-level; perf/numerical knobs go in `control`.
  tulpa::tulpa_nested_laplace(
    y = D_i, n_trials = rep(1L, n_sites), X = X_occ,
    prior = latent_prior, likelihood = lik,
    control = list(max_iter = as.integer(max_iter),
                   tol = as.numeric(tol),
                   n_threads = as.integer(n_threads))
  )
}

.tobs_nested_state_posterior <- function(model, occ_fit, latent_prior,
                                         heldout = NULL) {
  X_occ  <- model$X_processes[[1]]
  n_rows <- nrow(X_occ)
  w      <- occ_fit$weights

  fe  <- occ_fit$fitted_eta
  fev <- occ_fit$fitted_eta_var
  psi_lower <- rep(NA_real_, n_rows)
  psi_upper <- rep(NA_real_, n_rows)

  if (!is.null(fe) && is.matrix(fe) && ncol(fe) == n_rows &&
      !is.null(w) && length(w) == nrow(fe)) {
    have_var <- !is.null(fev) && is.matrix(fev) && all(dim(fev) == dim(fe))
    psi_mean <- .nested_psi_mean(fe, if (have_var) fev else NULL, w)
    if (have_var) {
      q <- .nested_psi_quantiles(fe, fev, w, probs = c(0.025, 0.975))
      psi_lower <- q[, 1L]
      psi_upper <- q[, 2L]
    }
  } else {
    # Fallback (no engine fitted_eta): reconstruct eta from the modes, exact
    # only for d_fac = 1 priors; bym2 / other mixed-scale priors -> NULL.
    if (is.null(occ_fit$modes) || !is.matrix(occ_fit$modes)) return(NULL)
    blocks <- if (!is.null(latent_prior$type)) list(latent_prior)
              else latent_prior
    if (!all(vapply(blocks, .nl_block_exact_reconstruct, logical(1))))
      return(NULL)
    modes  <- occ_fit$modes
    p_occ  <- ncol(X_occ)
    n_grid <- nrow(modes)
    if (ncol(modes) < p_occ) return(NULL)
    if (is.null(w)) w <- rep(1 / n_grid, n_grid)
    w <- w / sum(w)
    eta <- X_occ %*% t(modes[, seq_len(p_occ), drop = FALSE])
    for (k in seq_len(n_grid)) {
      fld <- .nested_field_from_mode(latent_prior, modes[k, ], p_occ, n_rows)
      if (length(fld) == n_rows) eta[, k] <- eta[, k] + fld
    }
    psi_mean <- as.numeric(plogis(eta) %*% w)
  }

  heldout_flag <- logical(n_rows)
  if (!is.null(heldout) && length(heldout) > 0L) heldout_flag[heldout] <- TRUE
  out <- data.frame(row = seq_len(n_rows), psi = psi_mean,
                    psi_lower = psi_lower, psi_upper = psi_upper,
                    heldout = heldout_flag)
  if (isTRUE(getOption("tobs.nested.debug"))) {
    attr(out, "engine") <- list(fitted_eta = fe, fitted_eta_var = fev,
                                weights = w)
  }
  out
}

# Length of a block's latent vector in the concatenated mode tail.
.nl_block_field_len <- function(b) {
  type <- b$type
  if (type %in% c("icar", "car_proper")) return(as.integer(b$n_spatial_units))
  if (type == "bym2") return(2L * as.integer(b$n_spatial_units))
  if (type %in% c("ar1", "rw1", "rw2")) return(as.integer(b$n_times))
  if (type == "iid") return(as.integer(b$n_units))
  if (type == "spde") return(as.integer(b$n_mesh))
  # Unknown block: assume the per-row index spans the units.
  as.integer(max(.nl_block_unit_idx(b)))
}

# Per-row 1-based unit index a block maps each state row through.
.nl_block_unit_idx <- function(b) {
  if (!is.null(b$spatial_idx))  return(as.integer(b$spatial_idx))
  if (!is.null(b$temporal_idx)) return(as.integer(b$temporal_idx))
  if (!is.null(b$obs_idx))      return(as.integer(b$obs_idx))
  integer(0)
}


# Identify held-out state rows for INLA NA-response prediction. A single-season
# occupancy site with no valid visits (an all-NA / all-missing detection
# history) is a prediction target: dropped from the likelihood (`n_trials = 0`)
# and interpolated by the latent field. Returns integer site indices, or NULL
# when there are none. Other model types are not yet wired for held-out
# prediction (the state row -> "site with a missing response" mapping differs
# for community / dynamic / integrated), so they return NULL.
.tobs_heldout_sites <- function(model) {
  if (!identical(model$model_type, "single")) return(NULL)
  n_valid <- rowSums(model$y >= 0)
  ho <- which(n_valid == 0L)
  if (length(ho) == 0L) NULL else as.integer(ho)
}


# Resolve the state ("occ") M-step block geometry: how many rows it has and
# which site each row belongs to. single / integrated / dynamic (psi1) carry
# one state row per site (identity map); community carries n_sites * n_species
# rows ordered site-major (`rep(site, each = n_species)`), so a site-level
# latent field is shared across the species at a site.
.tobs_state_block_dims <- function(model) {
  n_state_rows <- nrow(model$X_processes[[1]])
  n_sites      <- as.integer(model$n_sites)
  if (!is.null(model$n_species) &&
      n_state_rows == n_sites * as.integer(model$n_species)) {
    site_of_row <- rep(seq_len(n_sites), each = as.integer(model$n_species))
  } else if (n_state_rows == n_sites) {
    site_of_row <- seq_len(n_sites)
  } else {
    stop(sprintf(
      "Cannot map %d state-block rows onto %d sites for the nested latent prior.",
      n_state_rows, n_sites), call. = FALSE)
  }
  list(n_state_rows = n_state_rows, n_sites = n_sites,
       site_of_row = as.integer(site_of_row))
}


# Resolve the spatial term carried on the occupancy formula into the ordered
# list of areal field specs the nested-Laplace path builds blocks from: the
# unweighted intercept field first, then any varying-coefficient (weighted)
# fields. Three shapes flow in:
#   * a plain areal term -- icar()/bym2()/car_proper() -- one field, unweighted.
#   * a weighted areal term -- icar(graph, weight = col, group_var = node) -- a
#     single varying-coefficient field (the SVC slope on `col`); valid here on
#     its own, but typically paired with an intercept term.
#   * an independent varying-coefficient bar -- spatial(~ 1 + x || node,
#     graph = adj) -- desugars to the intercept field plus one weighted field
#     per covariate column, all on the same graph keyed by the bar's node index.
# The correlated bar (single `|`, a free-Sigma MCAR field) is fitted as one
# coupled cross-covariance block on the occu_cover() joint engine; it is not
# wired on the single-arm occu() path, so it errors with a pointer there.
.tobs_resolve_occu_spatial_fields <- function(spatial, model) {
  if (!inherits(spatial, "tobs_spatial")) {
    stop("`spatial` must be a tobs_spatial object", call. = FALSE)
  }
  if (isTRUE(spatial$is_bar)) {
    if (isTRUE(spatial$correlated)) {
      stop("occu(): a correlated spatial bar (`|`, free-Sigma MCAR) is fitted ",
           "on the occu_cover() joint engine, not on the single-arm occu() ",
           "nested-Laplace path. Use the independent bar `||` for separate ",
           "per-coefficient fields, or move the model to occu_cover().",
           call. = FALSE)
    }
    fields <- .tobs_expand_spatial_bar(spatial, model$data)
  } else if (isTRUE(spatial$is_multifield)) {
    # The explicit two-term form: an intercept field plus weighted SVC field(s),
    # already separate `tobs_spatial` specs. Order them intercept-first so the
    # field tables / sigma labels match the bar form.
    fields <- spatial$fields
    weighted <- vapply(fields, function(f) !is.null(f$weight), logical(1))
    if (!any(!weighted)) {
      stop("occu(): a varying-coefficient spatial field needs an unweighted ",
           "intercept field (e.g. icar(graph = adj, group_var = \"cell\")) ",
           "alongside the weighted term(s).", call. = FALSE)
    }
    fields <- c(fields[!weighted], fields[weighted])
  } else {
    fields <- list(spatial)
  }
  # The single-arm path builds one block per field; the weighted (SVC) fields
  # ride the same graph as the intercept field. Order is intercept-first, which
  # the bar expansion already guarantees; a lone weighted term is allowed.
  fields
}


# Build a multi-block-shaped spatial block from a tobs_spatial spec. The
# multi-block C++ entry reads spatial_idx as 1-based per the existing
# .nl_block_spec_for_cpp() contract.
#
# `model` supplies the site -> field-node map: with an areal `group_var` the
# occupancy units (sites, one per data row) map onto fewer field nodes (cells),
# so spatial_idx[row] = data[[group_var]][site_of_row[row]] and the field has
# one node per graph cell. Without group_var the field is one node per site
# (identity), so spatial_idx = site_of_row. A weighted (varying-coefficient)
# term carries a per-site `weight` column that becomes the block's per-row
# `svc_weight`, expanded onto the state rows through the same site_of_row map.
#
# Per-block grids are narrower than `.NL_REGISTRY`'s single-block defaults
# (which target ~20-25 cells per block). Three-block combinations
# (e.g. BYM2 + AR1 + IID) at the single-block defaults exceed the
# multi-block hard cap (2048 cells); the narrower defaults below keep
# typical combos under ~250 cells. Users who need finer integration can
# pass `*_grid` overrides directly via the tobs_* spec attributes (when
# present) -- the helper passes them through if set.
.tobs_block_from_spatial <- function(spatial, n_sites, site_of_row,
                                     model = NULL) {
  if (!inherits(spatial, "tobs_spatial")) {
    stop("`spatial` must be a tobs_spatial object", call. = FALSE)
  }
  type <- spatial$type
  if (type == "spde") {
    # Continuous Matern field on the multi-block nested-Laplace path. The
    # mesh projection A (n_sites x n_mesh) is re-rowed onto the n_state_rows
    # state block via site_of_row (one state row per site for single /
    # integrated / dynamic; community broadcasts the site field across the
    # species at a site), so a state row's field contribution is (A u)_row,
    # the many-to-one FEM projection -- not a one-node spatial_idx. The FEM
    # matrices / n_mesh / Matern priors are unchanged; only A is re-rowed.
    sp <- spatial$tulpa_spec
    if (as.integer(sp$n_mesh) <= 0L) {
      stop("SPDE mesh has 0 nodes; rebuild the mesh (see the tulpaMesh ",
           "zero-triangle note: use cutoff = 0 with the default max_edge).",
           call. = FALSE)
    }
    A_b   <- sp$A[site_of_row, , drop = FALSE]
    A_csc <- methods::as(A_b, "CsparseMatrix")
    n_state_rows <- length(site_of_row)
    out <- list(
      type        = "spde",
      n_mesh      = as.integer(sp$n_mesh),
      n_obs       = as.integer(n_state_rows),
      A           = A_b,                 # broadcast projection (eta offset reader)
      A_x         = as.numeric(A_csc@x),
      A_i         = as.integer(A_csc@i),
      A_p         = as.integer(A_csc@p),
      C0_diag     = as.numeric(sp$C0_diag),
      G1_x        = as.numeric(sp$G1_x),
      G1_i        = as.integer(sp$G1_i),
      G1_p        = as.integer(sp$G1_p),
      nu          = as.numeric(sp$nu),
      prior_range = as.numeric(sp$prior_range),
      prior_sigma = as.numeric(sp$prior_sigma)
    )
    if (!is.null(spatial$range_grid)) out$range_grid <- as.numeric(spatial$range_grid)
    if (!is.null(spatial$sigma_grid)) out$sigma_grid <- as.numeric(spatial$sigma_grid)
    return(out)
  }
  if (type %in% c("gp", "multiscale_gp")) {
    stop("Continuous-field spatial type '", type, "' is not yet wired into the ",
         "multi-block nested-Laplace path (method = 'nested_laplace'); the ",
         "FEM/basis projection for this type is not assembled here. The ",
         "spde() Matern field IS available on this path, and gp()/spde() are ",
         "available on the single-Laplace path (method = 'laplace').",
         call. = FALSE)
  }
  if (!type %in% c("bym2", "icar", "car_proper")) {
    stop("Spatial type '", type, "' is not yet wired into the multi-block ",
         "nested-Laplace path (supported: bym2, icar, car_proper). ",
         "Use `method = 'laplace'` or open an issue if you need this type.",
         call. = FALSE)
  }
  # Site -> field-node map. With an areal `group_var` the field has one node per
  # graph cell and many sites can share a node (e.g. cell-year sites sharing one
  # cell); spatial_idx[site] = data[[group_var]][site]. Without group_var the
  # field has one node per site (identity). Either way the per-state-row index is
  # the per-site index expanded through site_of_row.
  n_nodes <- as.integer(spatial$n_units)
  gv      <- spatial$group_var
  if (!is.null(gv)) {
    if (is.null(model) || is.null(model$data) || !gv %in% names(model$data)) {
      stop(sprintf("spatial group_var '%s' is not a column of the model data.",
                   gv), call. = FALSE)
    }
    site_node <- as.integer(model$data[[gv]])
    if (length(site_node) != n_sites || anyNA(site_node) ||
        min(site_node) < 1L || max(site_node) > n_nodes) {
      stop(sprintf(paste0(
        "spatial group_var '%s' must be an integer cell index in 1..%d, one ",
        "per site (%d sites)."), gv, n_nodes, n_sites), call. = FALSE)
    }
    spatial_idx <- site_node[site_of_row]
  } else {
    if (n_nodes != n_sites) {
      stop(sprintf(paste0(
        "spatial has %d units but the model has %d sites; one spatial unit per ",
        "site is required for the nested-Laplace latent field, or map sites to ",
        "cells with group_var = \"<col>\" on the areal term."),
        n_nodes, n_sites), call. = FALSE)
    }
    spatial_idx <- site_of_row
  }
  out <- list(
    type            = type,
    spatial_idx     = spatial_idx,
    n_spatial_units = n_nodes,
    adj_row_ptr     = as.integer(spatial$adj_row_ptr),
    adj_col_idx     = as.integer(spatial$adj_col_idx),
    n_neighbors     = as.integer(spatial$n_neighbors)
  )
  # A weighted (varying-coefficient) areal term carries a per-site `weight`
  # column; it becomes the block's per-state-row design weight so the field's
  # contribution is weight[i] * z[cell_i] (the areal f(cell, weight, ...)). The
  # weight is per observation (one per site), expanded onto state rows like the
  # node index. Only icar carries the SVC weight on this path.
  if (!is.null(spatial$weight)) {
    if (!identical(type, "icar")) {
      stop(sprintf(paste0(
        "a varying-coefficient (weighted) areal field uses the intrinsic CAR ",
        "(icar) on the occu() nested-Laplace path; model = \"%s\" is not ",
        "supported."), type), call. = FALSE)
    }
    w_site <- as.numeric(spatial$weight)
    if (length(w_site) != n_sites || any(!is.finite(w_site))) {
      stop(sprintf(paste0(
        "spatial weight must be a finite per-site numeric vector of length %d."),
        n_sites), call. = FALSE)
    }
    out$svc_weight <- w_site[site_of_row]
  }
  if (type == "bym2") {
    # A tulpa multi-block prior: the engine multiplies the structured block by
    # this, where the term carries the Riebler constant.
    out$scale_factor <- if (is.null(spatial$scale_factor)) 1.0 else
      .bym2_engine_scale(as.numeric(spatial$scale_factor))
    if (!is.null(spatial$sigma_grid)) out$sigma_grid <- spatial$sigma_grid
    if (!is.null(spatial$rho_grid))   out$rho_grid   <- spatial$rho_grid
    if (is.null(out$sigma_grid) && is.null(out$rho_grid)) {
      sg <- exp(seq(log(0.2), log(2.0), length.out = 3))
      rg <- c(0.3, 0.7)
      gr <- expand.grid(sigma = sg, rho = rg)
      # Ours, not the user's, so the engine may recentre it; `expand.grid()`
      # drops the marker, hence after.
      out$sigma_grid <- tulpa::auto_grid(gr$sigma)
      out$rho_grid   <- tulpa::auto_grid(gr$rho)
    }
  } else if (type == "icar") {
    if (!is.null(spatial$tau_grid)) out$tau_grid <- spatial$tau_grid
  } else if (type == "car_proper") {
    if (!is.null(spatial$tau_grid)) out$tau_grid <- spatial$tau_grid
    if (!is.null(spatial$rho_grid)) out$rho_grid <- spatial$rho_grid
  }
  out
}


# Build a multi-block temporal block from a tobs_temporal spec. Resolves
# the time variable from model$data when given as a string; passes integer
# indices through. `model = "iid"` becomes an iid block (no Q assembly).
.tobs_block_from_temporal <- function(temporal, model, n_sites, site_of_row) {
  if (!inherits(temporal, "tobs_temporal")) {
    stop("`temporal` must be a tobs_temporal object", call. = FALSE)
  }
  type <- temporal$type
  if (!type %in% c("ar1", "rw1", "rw2", "iid")) {
    stop("Temporal type '", type, "' is not supported by the multi-block ",
         "nested-Laplace path (supported: ar1, rw1, rw2, iid).",
         call. = FALSE)
  }
  # Index codes were resolved per-site when the temporal() term was
  # constructed; expand to one entry per state-block row via site_of_row.
  time_idx_site <- as.integer(temporal$time_idx)
  if (length(time_idx_site) != n_sites) {
    stop(sprintf(
      "Resolved temporal index has length %d but the model has %d sites.",
      length(time_idx_site), n_sites), call. = FALSE)
  }
  time_idx <- time_idx_site[site_of_row]
  n_times <- if (!is.null(temporal$n_times)) as.integer(temporal$n_times)
             else max(time_idx, na.rm = TRUE)

  if (type == "iid") {
    out <- list(type = "iid", obs_idx = time_idx, n_units = as.integer(n_times))
    if (!is.null(temporal$sigma_grid)) out$sigma_grid <- temporal$sigma_grid
    return(out)
  }

  out <- list(
    type         = type,
    temporal_idx = time_idx,
    n_times      = as.integer(n_times)
  )
  if (type == "rw1") out$cyclic <- isTRUE(temporal$cyclic)
  if (!is.null(temporal$tau_grid)) out$tau_grid <- temporal$tau_grid
  if (type == "ar1") {
    if (!is.null(temporal$rho_grid)) out$rho_grid <- temporal$rho_grid
    if (is.null(out$tau_grid) && is.null(out$rho_grid)) {
      # 3 x 2 = 6 cells per block keeps a BYM2 + AR1 + IID combo under cap.
      tg <- exp(seq(log(0.5), log(20), length.out = 3))
      rg <- c(0.3, 0.8)
      gr <- expand.grid(tau = tg, rho = rg)
      out$tau_grid <- tulpa::auto_grid(gr$tau)
      out$rho_grid <- tulpa::auto_grid(gr$rho)
    }
  }
  out
}


# Build a multi-block iid (or temporal-on-groups) block from a tobs_re spec.
# A bare `tobs_re(group = "x")` (default model = "iid") becomes an iid block.
# If the user asks for ar1/rw1/rw2 on the group, route through the temporal
# block builder so the same dispatch logic handles both.
.tobs_block_from_re <- function(re, model, n_sites, site_of_row) {
  if (!inherits(re, "tobs_re")) {
    stop("`re` element must be a tobs_re object", call. = FALSE)
  }
  if (identical(re$type, "slope")) {
    stop("Random slopes are not supported by the multi-block nested-Laplace ",
         "path. Use `method = 'laplace'` for uncorrelated slopes (`(x || g)`, ",
         "`(0 + x | g)`) or `method = 'nuts'` for correlated slopes ",
         "(`(1 + x | g)`).", call. = FALSE)
  }

  # Group codes were resolved per-site when the re() term was constructed;
  # expand to one entry per state-block row via site_of_row.
  grp_idx_site <- as.integer(re$group_idx)
  if (length(grp_idx_site) != n_sites) {
    stop(sprintf(
      "Resolved RE group index has length %d but the model has %d sites.",
      length(grp_idx_site), n_sites), call. = FALSE)
  }
  grp_idx <- grp_idx_site[site_of_row]
  n_units <- if (!is.null(re$n_groups)) as.integer(re$n_groups)
             else max(grp_idx, na.rm = TRUE)

  model_name <- re$model %||% "iid"
  if (model_name == "iid") {
    out <- list(type = "iid", obs_idx = grp_idx, n_units = as.integer(n_units))
    if (!is.null(re$sigma_grid)) {
      out$sigma_grid <- re$sigma_grid
    } else {
      # 3-point grid keeps multi-block combos comfortably under the cap.
      out$sigma_grid <- tulpa::auto_grid(exp(seq(log(0.2), log(2.0), length.out = 3)))
    }
    out
  } else if (model_name %in% c("ar1", "rw1", "rw2")) {
    out <- list(type = model_name, temporal_idx = grp_idx,
                n_times = as.integer(n_units))
    if (model_name == "rw1") out$cyclic <- FALSE
    if (!is.null(re$tau_grid)) out$tau_grid <- re$tau_grid
    out
  } else {
    stop("RE model '", model_name, "' is not supported by the multi-block ",
         "nested-Laplace path. Supported RE models: 'iid', 'ar1', 'rw1', 'rw2'.",
         call. = FALSE)
  }
}


# =============================================================================
# Driver — thin wrapper over `.tobs_laplace(..., latent_prior = )`
#
# The multi-block latent prior built from `spatial` / `temporal` / `re` is
# attached to the state ("occ") M-step block, which makes tulpa's generic EM
# engine route that block through `tulpa::tulpa_nested_laplace()` instead of
# `tulpa::tulpa_laplace()`. There is a single set of per-model-type callbacks
# (`build_*_callbacks` in R/laplace.R) shared with the single-Laplace path; the
# attachment and held-out encoding live in `.tobs_laplace_nested()`.
#
# Known approximation. The E-step computes psi[i] = plogis(X[i,]*beta_occ)
# only -- the latent block contribution (spatial / temporal / iid) enters
# through the M-step block log-marginal, not psi. The EM fixed point is the
# fixed-effect MAP under the marginal hyperparameter posterior, with the latent
# block estimated by the inner nested-Laplace fit on the pseudo-binomial
# response.
# =============================================================================

# Decide whether a standalone occu() nested-Laplace fit reroutes through the
# joint direct-grid engine (.tobs_fit_occu_joint) instead of the EM fixed-point
# path (.tobs_em_nested_laplace). The EM path oscillates / does not converge on
# the varying-coefficient (SVC) occupancy bar at EVA scale; the joint engine
# integrates the field hyperparameters on a direct outer grid, so it cannot
# oscillate. The reroute is scoped to exactly the SVC case the joint single-arm
# engine covers:
#   * single-season occupancy (`single`),
#   * an areal ICAR spatial term carrying a varying-coefficient structure -- an
#     independent (`||`) varying-coefficient bar, the explicit intercept +
#     weighted-trend form, or a weighted areal term,
#   * with NO temporal / re block (those need the EM multi-block latent prior;
#     the single-arm joint engine carries the areal field only).
# Everything else -- a plain single intercept field, a correlated (`|`) MCAR
# bar, bym2, temporal / re structure -- stays on the EM path unchanged.
.tobs_occu_reroute_to_joint <- function(model, spatial, temporal, re) {
  if (!identical(model$model_type, "single")) return(FALSE)
  if (!is.null(temporal) || !is.null(re))      return(FALSE)
  if (is.null(spatial) || !inherits(spatial, "tobs_spatial")) return(FALSE)
  if (isTRUE(spatial$is_bar)) {
    # The independent (`||`) bar desugars to intercept + per-coefficient trend
    # fields -- the SVC case. A correlated (`|`) bar is a free-Sigma MCAR field
    # fitted on the occu_cover joint engine, not the single-arm occu() path.
    return(!isTRUE(spatial$correlated))
  }
  if (isTRUE(spatial$is_multifield)) return(TRUE)
  # A lone weighted areal term (a single SVC slope field) is also covered.
  isTRUE(!is.null(spatial$weight) && identical(spatial$type, "icar"))
}


#' Fit an occupancy tobs model via nested-Laplace
#'
#' Internal driver: assembles a multi-block latent prior from `spatial`,
#' `temporal`, and `re`, then routes the state M-step block through
#' `tulpa::tulpa_nested_laplace()` via `.tobs_laplace(latent_prior = )`.
#' Supports single-season, integrated, community, and dynamic occupancy (the
#' state block is `"occ"` in every callback set).
#'
#' @param heldout_state Optional integer indices of state-block rows to treat
#'   as held-out (INLA NA-response prediction targets): they are dropped from
#'   the likelihood (`n_trials = 0`) but kept in the design so their latent
#'   value is informed by the prior. `NULL` (default) fits with no held-out
#'   rows.
#' @keywords internal
.tobs_em_nested_laplace <- function(model, spatial = NULL, temporal = NULL,
                                    re = NULL, priors = NULL,
                                    sigma_beta = 10,
                                    max_iter = 25L, tol = 1e-3,
                                    damping = 0.3,
                                    heldout_state = NULL,
                                    verbose = TRUE) {
  if (!inherits(model, "tobs_model")) {
    stop("model must be a tobs_model object", call. = FALSE)
  }

  multi_prior <- .tobs_to_multi_block_prior(
    spatial = spatial, temporal = temporal, re = re, model = model
  )
  if (is.null(multi_prior)) {
    stop("nested_laplace engine requires at least one latent block ",
         "(spatial, temporal, or re); none were supplied. Use ",
         "`method = 'laplace'` for a fit with no latent structure.",
         call. = FALSE)
  }

  fit <- .tobs_laplace(
    model, spatial = NULL, re = NULL, priors = priors,
    max_iter = max_iter, tol = tol, damping = damping,
    correction = "none",
    latent_prior = multi_prior, heldout_state = heldout_state,
    verbose = verbose
  )
  fit$temporal <- temporal
  fit$re <- re
  # Expose the spatial term so downstream code (predict / diagnostics) and the
  # field-shape check can read the mesh projection `fit$spatial$tulpa_spec$A`,
  # mirroring the single-Laplace fit. The latent realization is in
  # `fit$spatial_field`; for a continuous SPDE block it is the n_mesh field.
  if (!is.null(spatial)) fit$spatial <- spatial
  fit
}
