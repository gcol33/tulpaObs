# =============================================================================
# ccd_outer.R - mode-centred CCD for the outer field-hyperparameter integration
# of the in-package spatial / community fitters.
#
# The areal-BFGS families (dyn_abun, fp_occu), the community N-mixture Newton
# path (ms_abun + areal), and the SPDE community path each integrate their field
# hyperparameters (tau, rho, sigma, range) over a FIXED preset tensor grid. A
# fixed grid can land between nodes and miss a sharp likelihood mode, and a
# weakly-identified axis is swept coarsely or pinned. This is the same failure
# mode the tulpa engine's MCAR / SPDE / re-cov fits removed by mode-centring a
# central-composite design (CCD) at the marginal-likelihood mode.
#
# `.tobs_ccd_outer_grid()` reuses the engine's EXPORTED CCD primitives
# (`tulpa::ccd_grid` / `ccd_to_theta` / `ccd_weights`) and PSIS outer-accuracy
# diagnostic (`tulpa:::.nested_is_pareto_k`) so the node layout, design weights,
# and k-hat are the single source of truth shared with the engine. It owns only
# the per-site mode-find (an inner solve per evaluation, mirrored from
# `tulpa:::fit_spde_nested_ccd`).
#
# Each axis declares its physical box and a transform so the integral is taken in
# an unconstrained coordinate with a flat prior there -- log for a positive scale
# (tau / sigma / range), identity for a bounded correlation (rho), matching the
# implicit prior of the fixed grids those axes replace (flat over log-spaced tau,
# flat over linear rho).
#
# Returns NULL to DECLINE, in which case the caller runs its fixed tensor grid:
# when the outer mode-find fails, lands on the box boundary, or the outer
# curvature is ill-conditioned (a genuinely weakly-identified axis). This is the
# documented net the engine CCD uses -- the fixed grids stay as the fallback.
# =============================================================================

# One outer integration axis.
#   name  : reported hyperparameter name
#   tag   : "log" (positive scale, u = log theta) or "identity" (u = theta)
#   lower / upper : physical box; bounds the mode-find and clamps the CCD nodes
#   start : physical starting value for the mode-find
.tobs_ccd_axis <- function(name, tag, lower, upper, start) {
  if (!tag %in% c("log", "identity"))
    stop("CCD axis tag must be 'log' or 'identity'.", call. = FALSE)
  list(name = name, tag = tag, lower = lower, upper = upper, start = start)
}

# Physical <-> unconstrained (u) maps per axis tag.
.tobs_ccd_fwd <- function(axis, theta) if (axis$tag == "log") log(theta) else theta
.tobs_ccd_inv <- function(axis, u)     if (axis$tag == "log") exp(u)     else u

# Mode-centred CCD over outer field hyperparameters.
#
# `eval_logm(theta_phys)` returns the inner log-marginal (a scalar) at a PHYSICAL
# hyperparameter vector (the same value the fixed-grid loop integrates), or a
# non-finite value at an invalid node. `axes` is a list of `.tobs_ccd_axis()`.
#
# Each `eval_logm` call is one inner Laplace/EM solve, so the mode-find is bounded
# (`max_modefind`) and the outer PSIS k-hat is diagnosed on a small sample. CCD is
# attempted only for k >= 2 axes: a single positive hyperparameter is integrated
# more cheaply and robustly by its 1D grid, which is where CCD offers no node
# saving over the tensor product anyway. The multi-axis cases (proper-CAR
# (tau, rho), BYM2 (sigma, rho), SPDE (range, sigma[, r])) are the ones whose
# fixed tensor either pins or coarsely sweeps a weakly-identified axis.
#
# Returns NULL (decline -> caller uses the fixed grid) or a list:
#   nodes    : [n_node x k] physical hyperparameter values, colnames = axis names
#   dnode    : corrected R-INLA CCD design weights (length n_node)
#   u_hat    : mode in u-space
#   L        : Cholesky factor of the outer posterior covariance in u-space
#   kind     : per-node CCD label ("center"/"axial"/"factorial")
#   n_points : node count
#   pareto_k : outer PSIS k-hat for the Gaussian proposal (NA if not diagnosed)
.tobs_ccd_outer_grid <- function(eval_logm, axes, f0_mult = 1.1,
                                 diagnose_k = TRUE, k_samples = 40L,
                                 max_modefind = 40L,
                                 verbose = getOption("tulpaobs.ccd_verbose", FALSE)) {
  decline <- function(why) { if (isTRUE(verbose)) message("[ccd] decline: ", why); NULL }
  k     <- length(axes)
  if (k < 2L) return(decline("single axis -> the 1D grid is cheaper than a CCD"))
  u0    <- vapply(axes, function(a) .tobs_ccd_fwd(a, a$start), 0.0)
  lower <- vapply(axes, function(a) .tobs_ccd_fwd(a, a$lower), 0.0)
  upper <- vapply(axes, function(a) .tobs_ccd_fwd(a, a$upper), 0.0)

  phys <- function(u) vapply(seq_len(k), function(j) .tobs_ccd_inv(axes[[j]], u[j]), 0.0)
  # Negative outer log-marginal; minimised by the mode-find. Invalid / non-finite
  # nodes return a large finite penalty so L-BFGS-B keeps its box.
  neg <- function(u) {
    lm <- tryCatch(eval_logm(phys(u)), error = function(e) NA_real_)
    if (!is.finite(lm)) 1e10 else -lm
  }

  # Bounded mode-find. `maxit` is kept small because every evaluation is a full
  # inner solve; accepting a maxit-capped interior point (convergence 1) is fine
  # since the CCD design weights correct the Gaussian approximation and the node
  # log-marginals carry the actual integrand.
  op <- tryCatch(
    stats::optim(u0, neg, method = "L-BFGS-B", lower = lower, upper = upper,
                 control = list(factr = 1e7, maxit = as.integer(max_modefind),
                                ndeps = rep(5e-2, k))),
    error = function(e) NULL)
  if (is.null(op) || !op$convergence %in% c(0L, 1L) || any(!is.finite(op$par)))
    return(decline("mode-find failed"))
  # Mode pinned to the box -> the data are uninformative for some axis; decline.
  if (any(abs(op$par - lower) < 1e-3) || any(abs(op$par - upper) < 1e-3))
    return(decline("mode pinned to the box (weakly-identified axis)"))
  u_hat <- op$par

  # Outer curvature. `neg` is the negative log-marginal, so its Hessian at the
  # mode IS the precision of the Gaussian approximation to p(theta | y); its
  # inverse is the posterior covariance whose Cholesky scales the CCD design.
  H <- tryCatch(stats::optimHess(u_hat, neg), error = function(e) NULL)
  if (is.null(H) || any(!is.finite(H))) return(decline("non-finite outer Hessian"))
  H <- 0.5 * (H + t(H))
  ev <- tryCatch(eigen(H, symmetric = TRUE, only.values = TRUE)$values,
                 error = function(e) NA_real_)
  if (any(!is.finite(ev)) || max(abs(ev)) <= 0 || min(ev) <= 1e-6 * max(abs(ev)))
    return(decline("ill-conditioned outer curvature"))
  Sigma <- tryCatch(solve(H), error = function(e) NULL)
  if (is.null(Sigma)) return(decline("outer Hessian not invertible"))
  L <- tryCatch(t(chol(Sigma)), error = function(e) NULL)   # Sigma = L L'
  if (is.null(L)) return(decline("outer covariance not positive-definite"))

  ccd    <- tulpa::ccd_grid(k, f_0 = sqrt(k) * f0_mult)
  dnode  <- tulpa::ccd_weights(ccd)
  u_grid <- tulpa::ccd_to_theta(ccd$z, u_hat, L, log_scale = FALSE)   # nodes in u-space
  for (j in seq_len(k))
    u_grid[, j] <- pmin(pmax(u_grid[, j], lower[j]), upper[j])        # snap to the box
  nodes <- matrix(NA_real_, nrow(u_grid), k)
  for (i in seq_len(nrow(u_grid))) nodes[i, ] <- phys(u_grid[i, ])
  colnames(nodes) <- vapply(axes, function(a) a$name, "")

  pareto_k <- NA_real_
  if (isTRUE(diagnose_k)) {
    # Proposal N(u_hat, L L') in u-space; target is exp(eval_logm(phys(u))) with a
    # flat prior in u (matching the integration weights), so no extra Jacobian.
    lt <- function(U) vapply(seq_len(nrow(U)), function(i) {
      lm <- tryCatch(eval_logm(phys(U[i, ])), error = function(e) NA_real_)
      if (is.finite(lm)) lm else -Inf
    }, 0.0)
    kd <- tryCatch(tulpa:::.nested_is_pareto_k(u_hat, L, lt, as.integer(k_samples)),
                   error = function(e) list(pareto_k = NA_real_))
    pareto_k <- kd$pareto_k
  }

  list(nodes = nodes, dnode = dnode, u_hat = u_hat, L = L,
       kind = ccd$kind, n_points = ccd$n_points, pareto_k = pareto_k)
}
