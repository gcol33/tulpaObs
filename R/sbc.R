# =============================================================================
# Simulation-based calibration for tobs families.
#
# tulpa owns the SBC machinery -- the predictive shapes, the within-atom PIT,
# the CRPS closed forms, the exact simultaneous ECDF bands, and both drivers,
# all behind `tulpa::sbc()`. This file owns the other side of that contract: the
# callbacks that turn a fitted `tobs_fit` into the `model =` list the POSTERIOR
# experiment reads, so calibration is measured on the model a deliverable
# actually ships rather than on a toy fixture.
#
# WHY THE POSTERIOR EXPERIMENT AND NOT THE PRIOR-PREDICTIVE ONE. The nested door
# puts no prior on the fixed effects, so they cannot be drawn from and the
# prior-predictive experiment refuses them. Posterior SBC (Sailynoja, Schmitt,
# Buerkner & Vehtari 2026, Algorithm 2) draws the truth from the posterior at an
# observed data set, which is proper whatever the prior was:
#
#   theta'  ~ pi(theta | y_obs)          truth, from the fit under test
#   y_rep   ~ pi(y | theta')             replicate, on FRESH cells
#   theta'' ~ pi(theta | y_rep, y_obs)   the augmented posterior
#   u       = P(theta'' < theta' | y_rep, y_obs)
#
# Both stages use the algorithm's OWN approximate posterior, which is what turns
# a non-uniform rank into a verdict on the approximation.
#
# THREE CONSTRAINTS THE CONSTRUCTION RESTS ON, none of them optional:
#
#   FRESH CELLS. `occu_cover` couples its occurrence and cover arms through one
#   shared areal field, and theta carries no per-cell field value. A replicate
#   re-observing the cells y_obs already saw is therefore dependent on y_obs
#   given theta and breaks the factorization. Every replicate here is drawn on
#   new cells whose graph is block-diagonal against the observed one, and
#   `group_ids()` is supplied so `tulpa::sbc()` VERIFIES the observable half of
#   that (disjoint labels) instead of taking it on trust.
#
#   PARAMETERS AND DERIVED QUANTITIES, NOT THE LATENT FIELD. The per-cell field
#   is high-dimensional and is integrated out by the fit; what is ranked is the
#   arm fixed effects, the field SD on each arm, the copy scale alpha, and the
#   dispersion where the fit estimates one.
#
#   THE DERIVED SCALES ARE MARGINALIZED, NOT PLUGGED IN. `alpha` and the
#   cover-arm field SD `alpha * sigma` are read PER DRAW off the outer grid
#   through `.tobs_joint_draws()`, whose cells are sampled by their own
#   normalized grid weight. Their reported posterior is therefore the weighted
#   mixture over the joint grid, not a function of component modes.
#
# THE REGISTRY. Everything above is family-agnostic and lives in the shared
# driver; a family is one entry in `.TOBS_SBC_REGISTRY` naming its replicate
# generator, its fitting call, its posterior-draw reader, and optionally the
# joint statistic behind the rank arm. The pooling, the grouping labels, the
# arm construction, the mis-scaled controls and the fixed-quantity guard are
# written once here and shared, and the fitting call is declared as a slot ->
# argument map, so a second family is data entry rather than a second harness.
# The per-simulation seed split is `tulpa::sbc()`'s, not this file's.
# =============================================================================


# ---------------------------------------------------------------------------
# 1. The data shape the shared pooling / refitting operates on
#
# One row per site, plus the visit-level covariate matrices and the areal graph.
# It is deliberately the shape `tobs()` already takes, so `refit()` is a plain
# `tobs()` call rather than a private entry point.
# ---------------------------------------------------------------------------

# Recombine the site-level and visit-level halves a binder split the user's
# arm formula into. The intercept belongs to the site-level half.
.tobs_sbc_recombine <- function(f_site, f_visit) {
  lab <- character(0)
  icpt <- 1L
  if (!is.null(f_site)) {
    tt <- stats::terms(f_site)
    lab <- attr(tt, "term.labels")
    icpt <- attr(tt, "intercept")
  }
  if (!is.null(f_visit)) {
    lab <- c(lab, attr(stats::terms(f_visit), "term.labels"))
  }
  rhs <- if (length(lab)) paste(lab, collapse = " + ") else "1"
  if (!icpt) rhs <- paste(rhs, "- 1")
  stats::as.formula(paste("~", rhs), env = globalenv())
}

# Invert a visit-level design back into the per-visit covariate matrices
# `tobs(visits = )` takes. Only a plain numeric main effect inverts: a factor
# contrast or an interaction column is not a covariate the binder would rebuild,
# so it is refused here rather than silently reconstructed as something else.
.tobs_sbc_visit_matrices <- function(model) {
  n <- model$n_sites; J <- model$max_visits
  out <- list()
  for (slot in c("X_det_visit", "X_pos_visit")) {
    X <- model[[slot]]
    if (is.null(X) || !ncol(X)) next
    for (nm in colnames(X)) {
      if (!is.null(out[[nm]])) next
      if (!grepl("^[.A-Za-z][.A-Za-z0-9_]*$", nm)) {
        stop("SBC rebuilds the replicate design from the fitted one, which ",
             "needs every visit-level design column to be a plain numeric ",
             "covariate. Column ", sQuote(nm), " of ", slot, " is not ",
             "(a factor contrast or an interaction cannot be inverted into a ",
             "`visits` matrix). Refit with numeric visit covariates to run SBC ",
             "on this model.", call. = FALSE)
      }
      out[[nm]] <- matrix(as.numeric(X[, nm]), n, J, byrow = TRUE)
    }
  }
  out
}

# One data set, in the shape every refit, pool and `group_ids()` call reads. A
# family-specific extra (the multiscale family's cell count) rides after the six
# shared slots.
.tobs_sbc_dataset <- function(cells, y, site, y_pos = NULL, visits = NULL,
                              graph = NULL, ...) {
  c(list(cells = cells, y = y, y_pos = y_pos, visits = visits, graph = graph,
         site = site), list(...))
}

# The observed data set, in the shape `refit()` reads.
.tobs_sbc_data_from_fit <- function(fit) {
  m <- fit$model
  y <- m$y
  if (!is.null(m$valid)) y[!m$valid] <- NA
  .tobs_sbc_dataset(m$data, y, seq_len(m$n_sites), y_pos = m$y_pos,
                    visits = .tobs_sbc_visit_matrices(m),
                    graph = fit$spatial$graph)
}

# `rbind()` on two ZERO-COLUMN data frames collapses to 0 rows, 0 columns
# (base R silently drops the row count when there is no column to align on)
# -- an intercept-only fixed-effects fit (a `~ 1` arm, no site covariates)
# hits this on every pool call, not a corner case a fixture happens to avoid.
.tobs_sbc_rbind_cells <- function(a, b) {
  if (ncol(a) == 0L && ncol(b) == 0L) {
    return(data.frame(row.names = seq_len(nrow(a) + nrow(b))))
  }
  rbind(a, b)
}

# Stack the replicate below the observed data on the SITE axis (axis 1),
# whatever shape the response takes; every other axis -- visits, seasons,
# species, distance bins -- is shared design, not something a replicate draws
# fresh. A vector concatenates, a matrix binds rows, a higher-rank array is
# filled along axis 1 on the double scale, and a named list (one leaf per source,
# species or response matrix) stacks leaf by leaf, matched by name, since the
# leaves need not share a visit count. `NULL` on the observed side means the
# family carries no such slot.
.tobs_sbc_stack <- function(a, b) {
  if (is.null(a)) return(NULL)
  if (is.list(a) && !is.data.frame(a)) {
    return(Map(.tobs_sbc_stack, a, b[names(a)]))
  }
  d <- dim(a)
  if (is.null(d)) return(c(a, b))
  if (length(d) == 2L) return(rbind(a, b))
  n_o <- d[1L]; n_r <- dim(b)[1L]
  z <- rbind(matrix(as.double(a), n_o), matrix(as.double(b), n_r))
  dim(z) <- c(n_o + n_r, d[-1L])
  z
}

# Stack two data sets on FRESH sites. The replicate's site labels are offset
# past the observed ones and the two graphs go block-diagonal, so the pooled
# field carries two a-priori independent blocks and `group_ids()` sees
# n_obs + n_rep distinct labels -- which is the observable half of the
# conditional-independence premise `tulpa::sbc()` checks.
.tobs_sbc_pool <- function(obs, rep) {
  n_o <- length(obs$site)
  n_r <- length(rep$site)
  graph <- NULL
  if (!is.null(obs$graph)) {
    g <- matrix(0L, n_o + n_r, n_o + n_r)
    g[seq_len(n_o), seq_len(n_o)] <- as.matrix(obs$graph)
    g[n_o + seq_len(n_r), n_o + seq_len(n_r)] <- as.matrix(rep$graph)
    graph <- g
  }
  .tobs_sbc_dataset(.tobs_sbc_rbind_cells(obs$cells, rep$cells),
                    .tobs_sbc_stack(obs$y, rep$y),
                    c(obs$site, rep$site + n_o),
                    y_pos  = .tobs_sbc_stack(obs$y_pos, rep$y_pos),
                    visits = .tobs_sbc_stack(obs$visits, rep$visits),
                    graph  = graph)
}

.tobs_sbc_groups <- function(data) data$site


# ---------------------------------------------------------------------------
# 2. Rebuilding the fitting call
#
# The spec carries the arm formulas the binder recorded (spatial terms already
# stripped), the family object, the method and the fitting control, so `refit()`
# is the same `tobs()` call the observed fit came from with the pooled data and
# the pooled graph in place.
# ---------------------------------------------------------------------------

.tobs_sbc_default_control <- function() list(verbose = FALSE, progress = TRUE)

# The refit control every spec below starts from (#358). A NUTS fit records
# the sampler knobs it ACTUALLY resolved on `fit$sampler_control` (every NUTS
# fit-assembly writer stamps it -- see `.tobs_sampler_control_snapshot()` in
# engine_defaults.R); defaulting the replicate control from that instead of
# the (engine, family) table's sampler-sized defaults means a fit built with
# a short chain (for speed, or because that is what converges at a given N)
# is refit at that same size across the n.sim replicates, not at a size meant
# for one real fit. The caller's own `fit.control` still wins over both, via
# the same `modifyList()` precedence the old two-argument call had.
.tobs_sbc_control <- function(fit, fit.control) {
  base <- .tobs_sbc_default_control()
  if (identical(fit$method, "nuts") && length(fit$sampler_control))
    base <- utils::modifyList(base, fit$sampler_control)
  utils::modifyList(base, as.list(fit.control))
}

.tobs_sbc_spec <- function(fit, fit.control) {
  m  <- fit$model
  fm <- m$formulas
  fam <- attr(fit, "tobs_family")
  spatial <- fit$spatial
  if (!is.null(spatial) && !identical(spatial$type, "icar")) {
    stop("SBC rebuilds the replicate's areal graph, which is written today for ",
         "an intrinsic `icar()` field; this fit carries a ", spatial$type,
         " field. Refit with icar() to run SBC on it.", call. = FALSE)
  }
  ctl <- .tobs_sbc_control(fit, fit.control)
  if (!is.null(fit$joint_fit) && is.null(ctl$engine)) ctl$engine <- "joint"
  list(model      = m,
       family     = fam,
       method     = fit$method,
       control    = ctl,
       positive   = fit$positive,
       occ_fe     = fm$occ,
       det        = .tobs_sbc_recombine(fm$det, fm$det_visit),
       pos        = .tobs_sbc_recombine(fm$pos, fm$pos_visit),
       has_field  = !is.null(spatial),
       site_cell  = m$site_cell %||% seq_len(m$n_sites))
}

# The occupancy formula with the areal term re-attached on the graph the
# replicate / pooled data set actually carries. The term is written into a fresh
# environment holding that graph, which is how `tobs()` resolves it.
.tobs_sbc_occ_formula <- function(spec, graph) {
  lab <- attr(stats::terms(spec$occ_fe), "term.labels")
  if (is.null(graph)) {
    return(stats::as.formula(paste("~", if (length(lab))
      paste(lab, collapse = " + ") else "1"), env = globalenv()))
  }
  env <- new.env(parent = globalenv())
  env$.tobs_sbc_graph <- graph
  f <- stats::as.formula(
    paste("~", paste(c(lab, "icar(graph = .tobs_sbc_graph)"), collapse = " + ")))
  environment(f) <- env
  f
}

# The copy amplitude the observed fit ASKED for, written back as the `share()`
# argument that asks for it again.
#
# The stored formula does not carry the `share()` term (it is stripped when the
# structured terms are resolved) and the fit does not carry the control the
# formula compiled into, so a rebuilt call reverts to the engine's default axis
# and quietly analyses the replicate under a different model than the observed
# data -- the invariant `test-sbc.R`'s "rebuilds the same call" asserts, and
# already violated for a fit that pins its axis. `.occu_cover_jc_postprocess()`
# records the request for exactly this.
#
# The REQUEST, never the realised axis. `control$adaptive.grid` refines a
# defaulted (`auto`) axis against the data and only densifies a stated one, so
# handing back the nodes a fit realised asks for something different from what
# it asked for, and the rebuilt call stops reproducing it.
.tobs_sbc_copy_arg <- function(spec) {
  req <- spec$alpha_request
  if (is.null(req)) return(NULL)
  # A DEFAULTED axis arrives here as nodes too, carrying the `auto` mark, and
  # writing those back would state them -- turning a defaulted axis into a
  # stated one, which is the same category error as handing back a realised
  # one. A bare `share(spatial())` is what asks for the default, so that is
  # what a defaulted axis rebuilds as.
  if (!is.null(req$alpha_grid)) {
    if (isTRUE(tulpa::is_auto_grid(req$alpha_grid))) return(NULL)
    return(as.call(list(quote(grid), as.numeric(req$alpha_grid))))
  }
  if (!is.null(req$alpha_n))
    return(as.call(list(quote(grid), n = as.integer(req$alpha_n))))
  NULL
}

# The cover formula, carrying the shared field's copy onto the cover arm when
# the observed fit carried one. Without the copy the cover arm never sees the
# field and `alpha` is pinned at zero, which is a different model.
.tobs_sbc_pos_formula <- function(spec) {
  # `control$alpha.grid` / `control$alpha.grid.trend` are the other spelling of
  # the same coupling, and `occu_cover()` refuses either one beside a formula
  # copy, so the copy is written only when the control carries neither. Both are
  # read by exact key: `alpha.grid` is a strict prefix of `alpha.grid.trend`, so
  # a `$` read matches the trend key alone and returns NULL when both are set.
  ctl <- spec$control
  if (!isTRUE(spec$has_copy) || !is.null(ctl[["alpha.grid"]]) ||
      !is.null(ctl[["alpha.grid.trend"]])) {
    return(spec$pos)
  }
  lab <- attr(stats::terms(spec$pos), "term.labels")
  amp <- .tobs_sbc_copy_arg(spec)
  cp  <- if (is.null(amp)) quote(share(spatial())) else
    as.call(list(quote(share), quote(spatial()), alpha = amp))
  rhs <- Reduce(function(a, b) call("+", a, b),
                c(lapply(lab, str2lang), list(cp)))
  f <- stats::as.formula(call("~", rhs))
  environment(f) <- globalenv()
  f
}

.tobs_sbc_refit_occu_cover <- function(spec, data) {
  y_pos <- data$y_pos
  if (!is.null(y_pos)) y_pos[is.na(y_pos)] <- 0
  suppressWarnings(tobs(
    formula   = .tobs_sbc_occ_formula(spec, data$graph),
    data      = data$cells,
    family    = spec$family,
    detection = spec$det,
    positive  = .tobs_sbc_pos_formula(spec),
    y         = data$y,
    y.pos     = y_pos,
    visits    = data$visits,
    method    = spec$method,
    control   = spec$control))
}


# ---------------------------------------------------------------------------
# 3. Posterior draws over the scored quantities
#
# ONE named matrix per family, `[n.draws x n_quantities]`, and both the truth
# (`draw_theta`) and the predictives (`arms`) are read from it -- so a quantity
# can never be simulated under one name and scored under another.
#
# A column with no spread is DROPPED rather than scored: a dispersion the joint
# engine holds fixed in the cell-coupling spec has a point-mass posterior, and
# ranking a truth that equals it produces a degenerate PIT that reads as a
# defect. Which columns those are is decided ONCE, from a probe of the observed
# fit, and then applied to every later draw -- deciding it per call would drop
# every column of the single-draw read `draw_theta()` takes.
# ---------------------------------------------------------------------------

.TOBS_SBC_CONST_TOL <- 1e-8
.TOBS_SBC_PROBE_DRAWS <- 256L

.tobs_sbc_scored <- function(m) {
  spread <- apply(m, 2L, function(v) diff(range(v)))
  colnames(m)[spread > .TOBS_SBC_CONST_TOL * pmax(1, abs(colMeans(m)))]
}

# The joint nested-Laplace families: coefficients and the field scales come off
# the SAME grid-integrated draw, so the derived cover-arm SD and `alpha` are
# marginalized over the outer grid by construction.
.tobs_sbc_draws_joint_occu_cover <- function(fit, n) {
  d <- .tobs_joint_draws(fit, n = n)
  pi_l <- fit$process_info
  cn <- function(k, pre) paste0(pre, "_", pi_l[[k]]$coef_names)
  cols <- list(d$b$occ, d$b$det, d$b$pos)
  nms  <- c(cn(1L, "psi"), cn(2L, "p"), cn(3L, "pos"))
  M <- cbind(cols[[1]], cols[[2]], cols[[3]])
  colnames(M) <- nms
  # Blocks are NAMES, not positions: a dropped constant column must not be able
  # to shift which coefficients an arm is assembled from.
  blocks <- list(occ = cn(1L, "psi"), det = cn(2L, "p"), pos = cn(3L, "pos"))
  # The field scales and the dispersion are ranked on their own continuized
  # marginal (`$hyper`), not on the node values the bundle scores its latent at.
  jf <- .tobs_joint_fit(fit)
  if (length(d$blocks)) {
    sigma <- .tobs_joint_field_sd(jf$theta_grid, d$cells, d$hyper, 1L)
    alpha <- .tobs_joint_amp(jf$theta_grid, d$cells, d$hyper, 1L, "alpha",
                             default = 0)
    M <- cbind(M, sigma = sigma, sigma_pos_field = alpha * sigma,
               alpha = alpha)
  }
  M <- cbind(M, disp = .tobs_joint_disp(jf, d$cells, d$positive, hyper = d$hyper))
  attr(M, "blocks") <- blocks
  M
}

# The sampled route and the field-free Laplace route to the same table. Both
# carry the arm coefficients in `fit$draws` under the names this table uses:
# `nuts` its chain, the plain Laplace fitter (`.tobs_fit_occu_cover()`,
# `R/occu_cover.R`) its MVN pseudo-draws. What differs is only where the field
# scales come from. `nuts` samples them into `$hyper_draws`, so no grid is
# involved and every quantity is continuous; the plain Laplace fitter carries
# no shared field (that composition is nested-Laplace only) and reports none,
# so its field columns are the zero the joint route reports for a fit with no
# copy (#293).
#
# The conventions are the joint route's, not new ones. `hyper_draws[, "sigma"]`
# is the amplitude against the unscaled intrinsic precision, the same quantity
# `b<k>.sigma` carries on the grid (both satisfy `field_sd = sigma *
# sqrt(scale_q)`), and the dispersion column is already on cover()'s own
# surface -- `log_phi` is a log beta precision and `log_sigma_pos` a log
# residual SD -- which is what `.tobs_joint_disp()` converts the grid axis to.
# So the routes' tables are read on one scale and a rank from one is comparable
# with a rank from another.
.tobs_sbc_draws_sampled_occu_cover <- function(fit, n) {
  D <- fit$draws
  if (is.null(D) || !nrow(D))
    stop("This fit carries no posterior draw matrix to score.", call. = FALSE)
  # SBC asks for an arbitrary number of draws from a sample of fixed size, so
  # the rows are resampled. With replacement only when more are asked for than
  # the sample holds, which is a property of the request rather than of the fit.
  S <- nrow(D)
  take <- sample.int(S, n, replace = n > S)

  pi_l <- fit$process_info
  cn <- function(k, pre) paste0(pre, "_", pi_l[[k]]$coef_names)
  nms <- c(cn(1L, "psi"), cn(2L, "p"), cn(3L, "pos"))
  miss <- setdiff(nms, colnames(D))
  if (length(miss))
    stop("This fit's draws carry no column for: ",
         paste(miss, collapse = ", "), ".", call. = FALSE)
  M <- D[take, nms, drop = FALSE]
  colnames(M) <- nms
  blocks <- list(occ = cn(1L, "psi"), det = cn(2L, "p"), pos = cn(3L, "pos"))

  hd <- fit$hyper_draws
  hcol <- function(nm) {
    if (!is.null(hd) && nm %in% colnames(hd)) as.numeric(hd[take, nm])
    else rep(0, n)
  }
  sigma <- hcol("sigma")
  alpha <- hcol("alpha")
  M <- cbind(M, sigma = sigma, sigma_pos_field = alpha * sigma, alpha = alpha)

  disp_col <- .tobs_sbc_disp_name(fit$model)
  if (!disp_col %in% colnames(D))
    stop("This fit's draws carry no `", disp_col, "` column.", call. = FALSE)
  M <- cbind(M, disp = exp(as.numeric(D[take, disp_col])))

  attr(M, "blocks") <- blocks
  M
}

# The cover arm's dispersion coordinate, named for the positive family.
.tobs_sbc_disp_name <- function(m) {
  if (identical(m$positive, "beta")) "log_phi" else "log_sigma_pos"
}

# One entry in the registry, because the family is one family. Which route a
# fit took is a property of the fit, and reading it here is what lets the same
# 11 quantities be scored on the grid-integrated and the sampled posterior and
# compared.
.tobs_sbc_draws_occu_cover <- function(fit, n) {
  if (identical(fit$method, "nuts") || is.null(.tobs_joint_fit(fit)))
    return(.tobs_sbc_draws_sampled_occu_cover(fit, n))
  .tobs_sbc_draws_joint_occu_cover(fit, n)
}


# ---------------------------------------------------------------------------
# 4. The joint-statistic rank arm
#
# ONE scalar per simulation that depends on every scored quantity at once, so a
# posterior getting each marginal right while getting their joint dependence
# wrong is still caught. The rank of the truth's value among the same statistic
# at `n.ref` reference draws from the augmented posterior is uniform on
# 0..n_ref under correct inference for ANY fixed function of theta -- validity
# comes from the exchangeability of truth and references, not from the statistic
# being the exact marginal likelihood -- so what a family supplies here is
# chosen for POWER.
#
# For `occu_cover` that statistic is the data's log-likelihood with the shared
# field integrated out CELL BY CELL against a marginal of the field's own
# geo-mean width, `sigma * sqrt(scale_q)` (the scale the engine's raw-Q block
# actually carries, see section 5), and the copy alpha * sigma onto the cover
# arm. It moves with every coefficient, with both field scales and with the
# dispersion; it is not the model's exact marginal likelihood, which would have
# to integrate the field's cross-cell dependence and its per-cell marginal
# variances too.
# ---------------------------------------------------------------------------

# The dispersion scored when theta carries none: a quantity with no posterior
# spread is dropped from theta, which for the cover dispersion means the fit held
# it, and `cover_pos_disp` is that held value.
.tobs_sbc_held_disp <- function(model) {
  v <- model$cover_pos_disp
  if (is.null(v)) {
    stop("theta carries no `disp` and the fit records no held cover dispersion ",
         "(`model$cover_pos_disp`).", call. = FALSE)
  }
  v
}

.tobs_sbc_loglik_occu_cover <- function(fit, theta, n_quad = 15L) {
  m  <- fit$model
  bl <- attr(theta, "blocks")
  gh <- .gauss_hermite_prob(as.integer(n_quad))
  sigma <- if ("sigma" %in% names(theta)) theta[["sigma"]] else 0
  sigma <- sigma * .tobs_sbc_field_scale(fit$spatial$graph)
  alpha <- if ("alpha" %in% names(theta)) theta[["alpha"]] else 0
  disp  <- if ("disp" %in% names(theta)) theta[["disp"]] else
             .tobs_sbc_held_disp(m)
  # Each site's integral is taken independently, which is what makes this the
  # cell-by-cell surrogate rather than the model's own marginal; one shared
  # quadrature node per pass is therefore correct, and the rows are combined
  # afterwards.
  M <- vapply(seq_along(gh$nodes), function(k) {
    u <- gh$nodes[k]
    e <- .occu_cover_eta_from_par(m, theta[bl$occ], theta[bl$det], theta[bl$pos],
                                  off_occ = sigma * u,
                                  off_pos = alpha * sigma * u)
    .occu_cover_site_ll(m, e$psi, e$p_mat, e$ep_mat, log(disp)) +
      log(gh$weights[k])
  }, numeric(m$n_sites))
  M <- matrix(M, m$n_sites)
  mx <- apply(M, 1L, max)
  sum(mx + log(rowSums(exp(M - mx))))
}


# ---------------------------------------------------------------------------
# 5. Replicate generators
#
# Conditionally independent of the observed data given theta: a FRESH field is
# drawn on a graph isomorphic to the observed one, over new cells. The
# covariate design is the observed design re-instantiated on those cells, so
# the replicate is a second realization of the same survey rather than a second
# observation of the same cells.
#
# The linear predictors come from `.occu_cover_eta_from_par()`, the function the
# fitter and every diagnostic assemble eta with, so a generator cannot drift
# from the convention the likelihood is written in.
#
# THE FIELD'S SCALE IS THE ENGINE'S, NOT THE SIMULATOR'S. The joint
# nested-Laplace path hands the ICAR block the RAW graph precision Q = D - W at
# tau = 1 (`add_icar_prior(..., tau = 1.0, ...)` for a copy block) and carries
# the amplitude in the arm scale, so the field the likelihood integrates is
#   off_occ = sigma * x,   x ~ N(0, Q_aug^-1),   Q_aug = Q + sum_c 1_c 1_c'/J_c
# whose geo-mean marginal SD is sqrt(scale_q), the Sorbye-Rue constant of the
# graph -- NOT 1. `.occu_cover_draw_icar_field()` returns the NORMALISED draw
# (geo-mean marginal SD 1), which is the convention `simulate_occu_cover()` and
# the sampled-hyper NUTS route state their field in, so the replicate multiplies
# that constant back in. Drawing the normalised field instead would generate at
# `sigma` and refit under `sigma * sqrt(scale_q)`: both arm field SDs shift by
# one common factor while `alpha`, their ratio, stays clean.
# ---------------------------------------------------------------------------

# sqrt(scale_q) for a graph, cached on the last graph seen. The rank arm scores
# n.ref + 1 parameter vectors per simulation against the SAME pooled graph, and
# the constant costs an eigendecomposition, so recomputing it per call would
# dominate the run. `identical()` on the graph is the key: a hash could collide,
# and a wrong constant here is a silently mis-scaled experiment.
.TOBS_SBC_SCALE_CACHE <- new.env(parent = emptyenv())

.tobs_sbc_field_scale <- function(graph) {
  if (is.null(graph)) return(1)
  hit <- .TOBS_SBC_SCALE_CACHE$entry
  if (!is.null(hit) && identical(hit$graph, graph)) return(hit$scale)
  s <- sqrt(.occu_cover_icar_scale(as.matrix(graph)))
  .TOBS_SBC_SCALE_CACHE$entry <- list(graph = graph, scale = s)
  s
}

.tobs_sbc_sim_occu_cover <- function(spec, theta, seed) {
  set.seed(seed)
  m  <- spec$model
  bl <- attr(theta, "blocks")
  n <- m$n_sites; J <- m$max_visits
  sigma <- if ("sigma" %in% names(theta)) theta[["sigma"]] else 0
  alpha <- if ("alpha" %in% names(theta)) theta[["alpha"]] else 0
  disp  <- if ("disp" %in% names(theta)) theta[["disp"]] else
             .tobs_sbc_held_disp(m)

  f_cell <- if (spec$has_field)
    (spec$field_scale %||% 1) *
      as.numeric(.occu_cover_draw_icar_field(spec$model_graph, 1L)) else
    rep(0, max(spec$site_cell))
  f <- f_cell[spec$site_cell]

  e <- .occu_cover_eta_from_par(m, theta[bl$occ], theta[bl$det], theta[bl$pos],
                               off_occ = sigma * f, off_pos = alpha * sigma * f)
  z <- stats::rbinom(n, 1L, e$psi)
  det <- matrix(stats::rbinom(n * J, 1L, as.numeric(e$p_mat)), n, J)
  y <- det * z
  y_pos <- matrix(0, n, J)
  hit <- which(y == 1L)
  if (length(hit)) {
    y_pos[hit] <- .tobs_draw_positive_cover(e$ep_mat[hit], disp, spec$positive)
  }
  if (!is.null(m$valid)) { y[!m$valid] <- NA; y_pos[!m$valid] <- NA }

  obs <- spec$data_obs
  list(cells  = obs$cells,
       y      = y,
       y_pos  = y_pos,
       visits = obs$visits,
       graph  = spec$model_graph,
       site   = obs$site)
}


# ---------------------------------------------------------------------------
# 6. Families whose site marginals multiply
#
# `occu`, `abun`, `count`, `removal`, `distance`, `fp_occu`, `royle_nichols`
# and `occu_ttd` each fit ONE data set in which nothing latent is shared between
# two sites: every quantity the fit reports is in theta, and the likelihood is a
# product over sites. Three consequences carry the whole construction, so a
# family here supplies its fitting call and nothing else.
#
#   INDEPENDENCE IS STRUCTURAL, NOT ARRANGED. With no unmodelled quantity shared
#   between two blocks of sites, a replicate generated at theta on a second
#   block is independent of the observed data given theta. The replicate is
#   still given its own site labels and stacked below the observed rows, so
#   `group_ids()` hands `tulpa::sbc()` disjoint labels to CHECK rather than an
#   argument to trust -- the same premise `occu_cover` buys with a
#   block-diagonal graph, here already held by the model. A fit carrying a
#   structured term is refused: its field is exactly the shared unmodelled
#   quantity this rests on not existing, and theta would carry no value for it.
#
#   THE REPLICATE IS THE FAMILY'S OWN POSTERIOR-PREDICTIVE KERNEL, evaluated at
#   a `draws` matrix holding the single row theta. `simulate()` picks a draw and
#   generates from the law the fitter maximised, so there is no second copy of
#   the likelihood for a generator to drift from.
#
#   THE JOINT STATISTIC IS THE FAMILY'S OWN MARGINAL, `.tobs_pointwise_loglik()`
#   summed over sites -- the exact log-likelihood here, where `occu_cover` can
#   only afford a cell-by-cell surrogate. It reads a whole draw matrix, so the
#   truth and every reference are scored in ONE kernel call.
# ---------------------------------------------------------------------------

# The two rejections every registered spec needs, in ONE call so a spec cannot
# carry one and omit the other.
#
# A structured term makes a latent quantity that theta does not carry shared
# across sites, which is the one thing this construction denies. A visit-level
# observation design is rebuilt by no family's simulate() kernel, so the refit
# would score coefficients absent from the data it sees. Either is refused with
# a pointer rather than scored on the coefficients alone, which would report a
# calibration the experiment did not measure. On ms_abun(), for instance, the
# binder stores `X_det_visit` and appends its columns to the detection
# coefficient names while `formulas$det` keeps only the site-level arm, so
# without the visit half theta would be drawn from a posterior carrying
# visit-level coefficients and refit under a model that has none.
.tobs_sbc_reject_unsupported <- function(fit) {
  .tobs_sbc_reject_visit_design(fit)
  st <- fit$model$structured_terms
  present <- character(0)
  if (is.list(st)) present <- names(Filter(Negate(is.null), st))
  if (!is.null(fit$spatial))   present <- unique(c(present, "spatial"))
  if (!is.null(fit$temporal))  present <- unique(c(present, "temporal"))
  # `fit$re` is the term spec every route stamps when a `re()`/`(1|g)` term
  # is in the formula (laplace.R, em_nested_laplace.R, occu_fit.R,
  # nuts_chains.R), the same rule spatial/temporal are caught by above.
  # `fit$re_effects` is only the deterministic-laplace path's computed BLUP
  # summary -- absent on a nested_laplace RE fit, which is what let one pass
  # this gate and reach the nested route's "requires at least one latent
  # block" refusal instead of the intended one below.
  if (!is.null(fit$re) || !is.null(fit$re_effects)) {
    present <- unique(c(present, "re"))
  }
  if (!length(present)) return(invisible(NULL))
  stop("SBC on family ", sQuote(attr(fit, "tobs_family")$name %||% "?"),
       " is registered for a fit whose sites are conditionally independent ",
       "given the scored parameters; this fit carries a structured term (",
       paste(present, collapse = ", "), "), whose latent values theta does not ",
       "hold. Refit without it, or use the coupled `occu_cover` route, which ",
       "draws its replicate on fresh cells.", call. = FALSE)
}

# The -1 a binder writes for a visit not conducted. Every response these
# families take is non-negative, so one rule reads both conventions.
.tobs_sbc_clean_y <- function(y) {
  if (is.null(dim(y))) y <- matrix(y, ncol = 1L)
  storage.mode(y) <- "double"
  y[y < 0] <- NA
  y
}

.tobs_sbc_data_simple <- function(fit, resp) {
  m <- fit$model
  y <- .tobs_sbc_clean_y(m[[resp]])
  cells <- m$data
  if (!is.data.frame(cells) || nrow(cells) != nrow(y)) {
    stop("SBC rebuilds the replicate from the fitted model's own site-level ",
         "data frame, which must have one row per site (", nrow(y), "); this ",
         "fit's has ", if (is.data.frame(cells)) nrow(cells) else "none", ".",
         call. = FALSE)
  }
  .tobs_sbc_dataset(cells, y, seq_len(nrow(y)),
                    visits = .tobs_sbc_visit_matrices(m))
}

# The observed data set of a family whose response sits on the model under one
# shared site axis: `y` (a matrix, an array of any rank, or a named list of
# them) and optionally `y_pos`.
.tobs_sbc_data_model <- function(fit, y = fit$model$y, y_pos = NULL) {
  m <- fit$model
  .tobs_sbc_dataset(m$data, y, seq_len(m$n_sites), y_pos = y_pos)
}

# A visit-level design on the observation arm is REFUSED, not rebuilt. The
# replicate comes from the family's own `simulate()` kernel, and the
# single-season one assembles detection from the SITE-level design alone -- so a
# visit-level column would be scored by the refit and absent from the data the
# refit sees, which is a different model on each side and nothing in the ranks
# would say so.
.tobs_sbc_reject_visit_design <- function(fit) {
  m <- fit$model
  slots <- grep("_visit$", names(m), value = TRUE)
  slots <- slots[vapply(slots, function(s)
    is.matrix(m[[s]]) && ncol(m[[s]]) > 0L, logical(1))]
  if (!length(slots)) return(invisible(NULL))
  stop("SBC on family ", sQuote(attr(fit, "tobs_family")$name %||% "?"),
       " generates its replicate from the family's own simulate() kernel, ",
       "which assembles the observation arm from the site-level design; this ",
       "fit carries a visit-level one (", paste(slots, collapse = ", "),
       "). Refit with site-level observation covariates to run SBC on it.",
       call. = FALSE)
}

# A registration covers a declared scope of its family. A fit outside it is
# refused, naming what it carries, rather than scored under a model the
# replicate generator or the pooled refit does not reproduce.
.tobs_sbc_scope <- function(family, ok, scope, got, why = NULL) {
  if (isTRUE(ok)) return(invisible(NULL))
  stop("SBC on ", family, "() is registered for ", scope, " only (got ", got,
       ")", if (!is.null(why)) paste0("; ", why), ".", call. = FALSE)
}

# Constant-rate scope of a multi-season family: none of the named arms carries
# a season-varying design (`model$<arm>_season_varying`), which the replicate
# and the pooled refit do not reproduce.
.tobs_sbc_constant_rate <- function(family, arms) {
  function(fit) {
    sv <- vapply(arms, function(a)
      isTRUE(fit$model[[paste0(a, "_season_varying")]]), logical(1))
    .tobs_sbc_scope(family, !any(sv), "constant-rate fits",
                    paste0("a season-varying (", paste(arms[sv], collapse = ", "),
                           ") design"))
  }
}

# Full-overlap scope of a multi-source family: every source observes every
# site, which is the only case the exact marginal the rank arm reuses is wired
# for. `maps` are the per-source row -> site maps, 1-based.
.tobs_sbc_full_overlap <- function(family, maps, n_sites) {
  full <- vapply(maps, function(mp)
    identical(as.integer(mp), seq_len(n_sites)), logical(1))
  .tobs_sbc_scope(family, all(full),
                  "full-overlap fits (every source observes every site)",
                  "partial overlap")
}

# The fitting call a registered family is refit with, declared as data.
#
# `formulas` maps each `tobs()` formula argument to the model's own arm-formula
# slot -- recombined from the site-level half the binder recorded -- or to a
# function of the fit building it; the first entry is the state formula and is
# passed as `formula`. `response(data, spec)` returns the response arguments for
# a data set. `family` / `method` / `control` override what the observed fit
# records where a family's front door needs it rebuilt, and `fields(fit)` adds
# whatever the family's replicate generator reads off the spec. Every spec goes
# through the one rejection gate first and the family's `scope` check next.
.tobs_sbc_spec_fn <- function(formulas, response = .tobs_sbc_response_y,
                              scope = NULL, family = NULL, method = NULL,
                              control = .tobs_sbc_control, fields = NULL) {
  function(fit, fit.control) {
    .tobs_sbc_reject_unsupported(fit)
    if (!is.null(scope)) scope(fit)
    m <- fit$model
    args <- lapply(formulas, function(s)
      if (is.function(s)) s(fit) else .tobs_sbc_recombine(m$formulas[[s]], NULL))
    c(list(model    = m,
           fit_obs  = fit,
           family   = if (is.null(family)) attr(fit, "tobs_family") else family(fit),
           method   = method %||% fit$method %||% "laplace",
           control  = control(fit, fit.control),
           args     = args,
           response = response),
      if (!is.null(fields)) fields(fit))
  }
}

.tobs_sbc_response_y <- function(data, spec) list(y = data$y)

.tobs_sbc_response_species <- function(data, spec) {
  list(y = data$y, species = spec$model$species_names)
}

# The same `tobs()` call the observed fit came from, on the data set handed in.
.tobs_sbc_refit <- function(spec, data) {
  a <- spec$args
  args <- c(a[1L], list(data = data$cells, family = spec$family), a[-1L],
            spec$response(data, spec),
            list(method = spec$method, control = spec$control), spec$extra)
  names(args)[1L] <- "formula"
  if (length(data$visits)) args$visits <- data$visits
  suppressWarnings(do.call(tobs, args))
}

# The replicate at theta: the family's own posterior-predictive kernel on a
# draw matrix holding the single row theta.
.tobs_sbc_sim_simple <- function(spec, theta, seed) {
  set.seed(seed)
  f <- spec$fit_obs
  D <- matrix(theta, nrow = 1L)
  colnames(D) <- names(theta)
  f$draws <- D
  obs <- spec$data_obs
  .tobs_sbc_dataset(obs$cells,
                    .tobs_sbc_clean_y(stats::simulate(f, nsim = 1L)[[1L]]),
                    obs$site, visits = obs$visits)
}

# The fit's posterior-draw path, read rather than re-sampled: a row of the
# matrix the fitter itself drew from N(mode, V). At the driver's `n.draws` the
# selection is a permutation of that matrix, so the reported predictive is the
# fit's own; at `draw_theta`'s single row it is one draw from it.
.tobs_sbc_draws_fit <- function(fit, n) {
  D <- fit$draws
  if (is.null(D) || !is.matrix(D)) {
    stop("SBC reads the posterior from the fit's own draw matrix; ",
         "`fit$draws` is missing or not a matrix.", call. = FALSE)
  }
  if (is.null(colnames(D))) colnames(D) <- fit$col_names %||% names(fit$means)
  idx <- sample.int(nrow(D), n, replace = n > nrow(D))
  D[idx, , drop = FALSE]
}

# The exact marginal log-likelihood at each row of `Theta`, summed over sites.
# One kernel call for the whole reference set.
.tobs_sbc_loglik_many_simple <- function(fit, Theta) {
  f <- fit
  f$draws <- Theta
  rowSums(.tobs_pointwise_loglik(f, n.draws = nrow(Theta)))
}

# One registry row over the shared pieces: the posterior read off `fit$draws`,
# the family's own `simulate()` kernel as the replicate, the generic refit, the
# family's `.tobs_pointwise_loglik()` marginal as the rank statistic, and the
# shared site-axis pool. A family replaces whichever it cannot share.
.tobs_sbc_entry <- function(spec, data = .tobs_sbc_data_model,
                            simulate = .tobs_sbc_sim_simple,
                            draws = .tobs_sbc_draws_fit,
                            refit = .tobs_sbc_refit,
                            loglik_many = .tobs_sbc_loglik_many_simple,
                            pool = NULL) {
  e <- list(spec = spec, data = data, draws = draws, simulate = simulate,
            refit = refit, loglik_many = loglik_many)
  if (!is.null(pool)) e$pool <- pool
  e
}

# The row of a family whose site marginals multiply. `state` / `det` name the
# arm formulas in the fitted model, `resp` its response slot, `extra` any
# further arguments the family's front door takes (a list, or a function of the
# model), `y_vector` a front door taking the response as a plain vector.
.tobs_sbc_simple_entry <- function(state, det = NULL, resp = "y",
                                   extra = NULL, y_vector = FALSE) {
  .tobs_sbc_entry(
    spec = .tobs_sbc_spec_fn(
      formulas = c(formula = state, detection = det),
      response = if (y_vector) function(data, spec)
        list(y = as.numeric(data$y)) else .tobs_sbc_response_y,
      scope = function(fit) {
        if (is.null(fit$model$formulas[[state]])) {
          stop("SBC rebuilds the fitting call from the model's own arm ",
               "formulas; slot ", sQuote(state), " is absent.", call. = FALSE)
        }
      },
      fields = function(fit)
        list(extra = if (is.function(extra)) extra(fit$model) else extra)),
    data = function(fit) .tobs_sbc_data_simple(fit, resp = resp))
}


# ---------------------------------------------------------------------------
# 6a-bis. The community group's shared draws / simulate / loglik_many
#
# Seven families rank a fixed species set to one design: theta carries one
# block of length P per species, species-major, with theta_s = mu + b_s,
# followed by any SHARED coordinates that carry no per-species deviation. The
# three functions over that layout were written out seven times, differing only
# in per-arm bookkeeping. The cost of that is on the record: the #226 fix (mu
# and b_s are not independent in the posterior, so they draw jointly) is a few
# lines inside `.tobs_sbc_community_b_draws()`, but landing it meant editing
# all seven copies, with nothing to guarantee the seventh got it.
#
# An ARM is one per-species coefficient block: `name` (what the family's kernel
# calls it), `prefix` (its theta-column prefix), `pi` (its slot in
# `model$process_info`, where the coefficient names live), `blup` and `coef`
# (its `fit$ms_community` slots -- the BLUP the draw centres on, and the
# coefficient matrix `simulate()` overwrites). Prefix and slot are independent:
# `ms_occu_cover`'s first arm is written `psi_` in theta and `occ` in
# `ms_community`.
#
# A GLOBAL is a shared block with no per-species deviation: `ms_dyn_occu`'s
# gamma / eps transition coefficients, `ms_occu_cover`'s cover-arm dispersion.
# Its columns sit after every species block, are copied straight off the
# community draw, and are written back onto `fit$means` by name. The dispersion
# is the one whose column name is not a prefix + coefficient-name product, so a
# global names its columns with a function of the model rather than a slot.
# ---------------------------------------------------------------------------

.tobs_sbc_arm <- function(name, prefix, pi, blup = paste0("blup_", name),
                          coef = paste0("coef_", name)) {
  list(name = name, prefix = prefix, pi = pi, blup = blup, coef = coef)
}

# A shared block whose columns are `prefix` + the coefficient names at
# `process_info[[pi]]`.
.tobs_sbc_global <- function(name, prefix, pi) {
  list(name = name,
       cols = function(m) paste0(prefix, m$process_info[[pi]]$coef_names))
}

# A shared block naming its own columns (the cover-arm dispersion, whose name
# depends on the positive family rather than on a design).
.tobs_sbc_global_named <- function(name, cols) list(name = name, cols = cols)

# An arm / global list may be a plain list or a function of the model, since
# `ms_int_occu`'s detection arms are one per source.
.tobs_sbc_resolve_blocks <- function(x, m) if (is.function(x)) x(m) else x

# Theta column names + the per-arm coefficient names every community generic
# below reads. Species-major over the arms, then the globals.
.tobs_sbc_community_names <- function(m, arms, globals) {
  arms  <- .tobs_sbc_resolve_blocks(arms, m)
  coefs <- lapply(arms, function(a) m$process_info[[a$pi]]$coef_names)
  par   <- unlist(lapply(seq_along(arms),
                         function(k) paste0(arms[[k]]$prefix, coefs[[k]])),
                  use.names = FALSE)
  sp_cols <- as.vector(outer(par, m$species_names,
                             function(p, sp) paste0(sp, "_", p)))
  globals <- .tobs_sbc_resolve_blocks(globals, m)
  gcols   <- lapply(globals, function(g) g$cols(m))
  names(gcols) <- vapply(globals, `[[`, character(1), "name")
  list(arms = arms, coefs = coefs, par = par, sp_cols = sp_cols,
       globals = globals, gcols = gcols,
       global_cols = unlist(gcols, use.names = FALSE) %||% character(0),
       cols = c(sp_cols, unlist(gcols, use.names = FALSE)))
}

# Species s's coefficients out of one theta row, as a list keyed by arm name.
.tobs_sbc_community_beta <- function(th, nm, m, s) {
  sp <- m$species_names[s]
  out <- lapply(seq_along(nm$arms), function(k)
    th[paste0(sp, "_", nm$arms[[k]]$prefix, nm$coefs[[k]])])
  names(out) <- vapply(nm$arms, `[[`, character(1), "name")
  out
}

# The shared blocks out of one theta row, keyed by global name. Unnamed: a
# kernel multiplies these into a design, and a stray name would ride out on the
# scalar `loglik_many` returns.
.tobs_sbc_community_globals <- function(th, nm) {
  lapply(nm$gcols, function(cols) unname(th[cols]))
}

# The community posterior draw: each species' block is the community mean draw
# plus that species' deviation, drawn JOINTLY with it ( -- see
# .tobs_sbc_community_b_draws), and each shared block is copied straight off
# the same draw.
.tobs_sbc_draws_community <- function(fit, n, nm) {
  m  <- fit$model
  S  <- m$n_species; P <- length(nm$par)
  cm <- fit$ms_community
  family <- attr(fit, "tobs_family")$name %||% m$model_type %||% "?"
  .tobs_sbc_require_laplace_covariance(fit, family, "fit$vcov", fit$vcov)
  mu_draws <- .tobs_sbc_mvn_draws(fit$means, fit$vcov, n)  # n x (P + globals)
  # The draw inherits its column names from `fit$vcov`'s dimnames; a shared
  # block is read off it by NAME, so name them here rather than rely on that.
  if (is.null(colnames(mu_draws))) colnames(mu_draws) <- names(fit$means)
  mu_arm <- mu_draws[, seq_len(P), drop = FALSE]
  M <- matrix(NA_real_, n, length(nm$cols)); colnames(M) <- nm$cols
  for (s in seq_len(S)) {
    b_hat_s <- unlist(lapply(nm$arms, function(a) cm[[a$blup]][s, ]),
                      use.names = FALSE)
    .tobs_sbc_require_laplace_covariance(
      fit, family, sprintf("fit$ms_community$Bf/Cinv[[%d]]", s),
      if (is.null(cm$Bf[[s]]) || is.null(cm$Cinv[[s]])) NULL else TRUE)
    # `Bf[[s]]` spans the FULL (mu, globals) vector, not just the RE arms: the
    # community EM's Newton solve treats a shared block as jointly informative
    # about each species' deviation even where it carries no RE of its own, so
    # the draw conditions on all of `mu_draws`.
    b_draws <- .tobs_sbc_community_b_draws(mu_draws, fit$means, b_hat_s,
                                           cm$Bf[[s]], cm$Cinv[[s]], n)
    M[, (s - 1L) * P + seq_len(P)] <- mu_arm + b_draws
  }
  if (length(nm$global_cols)) {
    M[, nm$global_cols] <- mu_draws[, nm$global_cols]
  }
  M
}

# The replicate at one theta row: write each arm's [S x p] coefficient matrix
# and each shared block back onto the fit, then let the family's own
# `simulate()` handler draw the response.
.tobs_sbc_sim_community <- function(spec, theta, seed, nm, y_pos = FALSE) {
  set.seed(seed)
  f <- spec$fit_obs
  m <- f$model
  S <- m$n_species
  for (k in seq_along(nm$arms)) {
    a <- nm$arms[[k]]; cn <- nm$coefs[[k]]
    # do.call(rbind, lapply(...)), not t(vapply(...)) -- vapply collapses a
    # length-1 per-species result (a single detection coefficient, e.g.
    # `detection = ~ 1`) to a plain vector, and t() of that transposes the
    # wrong axis (1 x S instead of S x 1).
    B <- do.call(rbind, lapply(seq_len(S), function(s)
      theta[paste0(m$species_names[s], "_", a$prefix, cn)]))
    dimnames(B) <- list(m$species_names, cn)
    f$ms_community[[a$coef]] <- B
  }
  if (length(nm$global_cols)) {
    f$means[nm$global_cols] <- theta[nm$global_cols]
  }
  rep <- stats::simulate(f, nsim = 1L)[[1L]]
  obs <- spec$data_obs
  .tobs_sbc_dataset(obs$cells, if (y_pos) rep$y else rep, obs$site,
                    y_pos = if (y_pos) rep$y_pos else NULL)
}

# The exact marginal log-likelihood at each row of `Theta`, summed over
# species. `ll` is the family's own per-species kernel -- the one thing that
# genuinely differs; `setup` builds whatever it needs once per call (species
# views, a distance engine, a long-form response), and `row` whatever depends
# on the shared blocks alone, once per row rather than once per species.
.tobs_sbc_loglik_many_community <- function(fit, Theta, nm, ll,
                                            setup = NULL, row = NULL) {
  m   <- fit$model
  S   <- m$n_species
  aux <- if (is.null(setup)) NULL else setup(fit)
  vapply(seq_len(nrow(Theta)), function(i) {
    th <- Theta[i, ]
    gl <- .tobs_sbc_community_globals(th, nm)
    rw <- if (is.null(row)) NULL else row(m, gl, aux)
    total <- 0
    for (s in seq_len(S)) {
      total <- total + ll(m, s, .tobs_sbc_community_beta(th, nm, m, s),
                          gl, aux, rw)
    }
    total
  }, numeric(1))
}

# One community registry row. `arms` / `globals` declare the layout; `ll`
# (+ optional `setup` / `row`) the family's kernel; `y_pos` says the family's
# `simulate()` returns both response arrays. `spec` is the family's declared
# fitting call; `data` / `refit` / `pool` default to the shared ones.
.tobs_sbc_community_entry <- function(arms, ll, spec,
                                      data = .tobs_sbc_data_model,
                                      refit = .tobs_sbc_refit,
                                      globals = list(), setup = NULL,
                                      row = NULL, pool = NULL, y_pos = FALSE) {
  nmf <- function(m) .tobs_sbc_community_names(m, arms, globals)
  entry <- list(
    spec     = spec,
    data     = data,
    # The resolved theta layout, exposed so a consumer can read the family's
    # ranked quantity off `ms_community` without restating its arm list.
    names    = nmf,
    draws    = function(fit, n) .tobs_sbc_draws_community(fit, n, nmf(fit$model)),
    simulate = function(spec, theta, seed)
      .tobs_sbc_sim_community(spec, theta, seed, nmf(spec$model), y_pos = y_pos),
    refit    = refit,
    loglik_many = function(fit, Theta)
      .tobs_sbc_loglik_many_community(fit, Theta, nmf(fit$model),
                                      ll = ll, setup = setup, row = row))
  if (!is.null(pool)) entry$pool <- pool
  entry
}


# ---------------------------------------------------------------------------
# 6b. Multi-season families: dyn_occu(), dyn_abun(), distsamp_open(),
# t_occu(). `model$y` is a 3D array whose first axis is the site; the shared
# pool stacks that axis and leaves the season axis (and its transition
# structure) alone, so only the fitting call and, where the family has no
# `simulate()` handler at a given theta, the replicate are their own.
# Constant-rate only: a season-varying detection / colonization / extinction
# (dyn_occu) or omega / gamma (dyn_abun) design is not seen by
# `.tobs_sbc_reject_visit_design()` and is refused explicitly.
# ---------------------------------------------------------------------------

.tobs_sbc_spec_dyn_occu <- .tobs_sbc_spec_fn(
  c(formula = "occ", detection = "det", colonization = "col",
    extinction = "ext"),
  scope = .tobs_sbc_constant_rate("dyn_occu", c("det", "col", "ext")))

# Forward-simulate the HMM: z_1 ~ Bern(psi1), z_t | z_{t-1} ~ colonization /
# extinction, y_{i,j,t} | z_{i,t} ~ Bern(p_i) over every visit. Missing visits
# in the fitted model's OWN design (its own `NA` pattern, not the pooled
# data's) are masked back to `NA` so the replicate's visit design matches the
# one it will be pooled against.
.tobs_sbc_sim_dyn_occu <- function(spec, theta, seed) {
  set.seed(seed)
  m <- spec$model
  pi_list <- m$process_info
  off <- 0L
  get_beta <- function(k) {
    b <- theta[off + seq_len(pi_list[[k]]$p)]
    off <<- off + pi_list[[k]]$p
    b
  }
  beta_psi1 <- get_beta(1L); beta_p <- get_beta(2L)
  beta_gam  <- get_beta(3L); beta_eps <- get_beta(4L)

  psi1 <- plogis(as.vector(m$X_processes[[1L]] %*% beta_psi1))
  p    <- plogis(as.vector(m$X_processes[[2L]] %*% beta_p))
  gam  <- plogis(as.vector(m$X_processes[[3L]] %*% beta_gam))
  eps  <- plogis(as.vector(m$X_processes[[4L]] %*% beta_eps))

  n <- m$n_sites; Tn <- m$n_seasons; J <- dim(m$y)[2L]
  z <- matrix(0L, n, Tn)
  z[, 1L] <- stats::rbinom(n, 1L, psi1)
  for (t in 2:Tn) {
    z[, t] <- ifelse(z[, t - 1L] == 1L,
                     stats::rbinom(n, 1L, 1 - eps),
                     stats::rbinom(n, 1L, gam))
  }
  y <- array(NA_real_, dim = c(n, J, Tn))
  for (t in seq_len(Tn)) {
    y[, , t] <- z[, t] * matrix(stats::rbinom(n * J, 1L, rep(p, J)), n, J)
  }
  y[is.na(m$y)] <- NA_real_

  .tobs_sbc_dataset(m$data, y, seq_len(n))
}

# dyn_abun(): the Dail-Madsen open N-mixture has a working `simulate()` handler,
# so its replicate is the family's own kernel.
.tobs_sbc_spec_dyn_abun <- .tobs_sbc_spec_fn(
  c(formula = "lambda", detection = "p", omega = "omega", gamma = "gamma"),
  scope = .tobs_sbc_constant_rate("dyn_abun", c("omega", "gamma")))

# distsamp_open(): the "season" axis is the open-population primary period
# rather than a revisit season, but the shape and the pooling rule are the
# same; `fit$draws` and the `.tobs_pointwise_loglik()` dispatch are the
# standard ones. Constant-dynamics, Poisson only (the alternative dynamics /
# negbin / zero-inflated layers are not reproduced by the rebuilt front door).
.tobs_sbc_spec_distsamp_open <- .tobs_sbc_spec_fn(
  c(formula = "lambda", detection = "sigma", omega = "omega", gamma = "gamma"),
  family = function(fit)
    distsamp_open(transect = fit$model$transect,
                  cutpoints = fit$model$cutpoints,
                  mixture = "poisson", dynamics = "constant"),
  scope = function(fit) {
    dyn <- fit$model$dynamics %||% "constant"
    mix <- fit$model$mixture %||% "poisson"
    .tobs_sbc_scope("distsamp_open",
                    identical(dyn, "constant") && identical(mix, "poisson"),
                    "constant-dynamics, Poisson fits",
                    sprintf("dynamics = \"%s\", mixture = \"%s\"", dyn, mix),
                    "the pooled refit does not yet reproduce it")
  })


# ---------------------------------------------------------------------------
# 6c. Multi-source and multi-response families: int_occu(), gdistremoval(),
# occu_multi(), dyn_int_occu(). `y` is a NAMED list of matrices (or, for
# dyn_int_occu, of 3D arrays) sharing one site axis -- one per source, per
# response (yDist / yRem) or per species -- and the shared pool stacks each
# leaf on that axis, matched by name.
#
# The two source-list families store their site maps on different bases, and
# both checks below are right for their own: int_occu()'s binder writes the
# 0-based maps its C++ marginal reads (`R/occu.R`, `site_maps <- ... m - 1L`),
# ms_int_occu()'s writes 1-based ones it indexes R arrays with
# (`R/ms_int_occu.R`, `yi[md, , s]`). `.tobs_sbc_full_overlap()` takes 1-based
# maps, so int_occu's are shifted by one before the check.
# ---------------------------------------------------------------------------

.tobs_sbc_int_occu_src_names <- function(m) {
  vapply(seq_len(m$n_sources), function(s) m$process_info[[1L + s]]$name,
        character(1))
}

# Full-overlap fits only: `.tobs_ploglik_integrated`, which `loglik_many`
# reuses unchanged, is itself only wired for that case.
.tobs_sbc_spec_int_occu <- .tobs_sbc_spec_fn(
  list(formula   = "occ",
       detection = function(fit)
         lapply(.tobs_sbc_int_occu_src_names(fit$model), function(nm)
           .tobs_sbc_recombine(fit$model$formulas[[nm]], NULL))),
  response = function(data, spec)
    list(y = data$y[.tobs_sbc_int_occu_src_names(spec$model)]),
  scope = function(fit)
    .tobs_sbc_full_overlap("int_occu", lapply(fit$model$site_maps, `+`, 1L),
                           fit$model$n_sites))

.tobs_sbc_data_int_occu <- function(fit) {
  .tobs_sbc_data_model(fit, stats::setNames(
    fit$model$y_sources, .tobs_sbc_int_occu_src_names(fit$model)))
}

# Shared latent occupancy z ~ Bern(psi), independently observed through each
# source's own detection process (per-visit Bernoulli(p_s)); the ordered
# detected/all-zero event per source stays consistent with the same z draw,
# matching what the exact marginal (`.tobs_ploglik_integrated`) integrates
# over.
.tobs_sbc_sim_int_occu <- function(spec, theta, seed) {
  set.seed(seed)
  m <- spec$model
  pi_list <- m$process_info
  n_src <- m$n_sources

  beta_psi <- theta[seq_len(pi_list[[1L]]$p)]
  psi <- plogis(as.vector(m$X_processes[[1L]] %*% beta_psi))
  z <- stats::rbinom(m$n_sites, 1L, psi)

  off <- pi_list[[1L]]$p
  y <- vector("list", n_src)
  for (s in seq_len(n_src)) {
    pp <- pi_list[[1L + s]]
    beta_s <- theta[off + seq_len(pp$p)]
    off <- off + pp$p
    p_s <- plogis(as.vector(m$X_processes[[1L + s]] %*% beta_s))
    J_s <- ncol(m$y_sources[[s]])
    ys <- z * matrix(stats::rbinom(m$n_sites * J_s, 1L, rep(p_s, J_s)),
                     m$n_sites, J_s)
    ys[m$y_sources[[s]] < 0] <- -1
    y[[s]] <- ys
  }
  names(y) <- .tobs_sbc_int_occu_src_names(m)
  .tobs_sbc_dataset(m$data, y, seq_len(m$n_sites))
}

# gdistremoval(): two response matrices (yDist band counts, yRem period counts)
# that must stay row-consistent; the family's `simulate()` handler returning
# `list(yDist =, yRem =)` is the replicate generator.
.tobs_sbc_spec_gdistremoval <- .tobs_sbc_spec_fn(
  c(formula = "lambda", detection = "sigma", removal = "r"),
  response = function(data, spec)
    list(y = data$y$yDist, y.rem = data$y$yRem))

.tobs_sbc_data_gdistremoval <- function(fit) {
  .tobs_sbc_data_model(fit, list(yDist = fit$model$y, yRem = fit$model$y_rem))
}

.tobs_sbc_sim_gdistremoval <- function(spec, theta, seed) {
  set.seed(seed)
  f <- spec$fit_obs
  D <- matrix(theta, nrow = 1L)
  colnames(D) <- names(theta)
  f$draws <- D
  rep <- stats::simulate(f, nsim = 1L)[[1L]]
  list(cells = spec$model$data, y = list(yDist = rep$yDist, yRem = rep$yRem),
       y_pos = NULL, visits = NULL, graph = NULL,
       site = seq_len(spec$model$n_sites))
}

# occu_multi(): one detection matrix per species. Species are NOT independent
# given a shared z -- the joint state is one draw from the log-linear
# multi-species model, then each species is independently detected given its
# own z_k -- so the replicate wraps the family's own
# `.tobs_simulate_occu_multi()` handler. Species names are attached explicitly,
# since the binder does not always name `model$y`.
.tobs_sbc_spec_occu_multi <- .tobs_sbc_spec_fn(
  c(formula = "state", detection = "det"),
  response = function(data, spec)
    list(y = data$y[spec$model$species], species = spec$model$species))

.tobs_sbc_data_occu_multi <- function(fit) {
  .tobs_sbc_data_model(fit, stats::setNames(fit$model$y, fit$model$species))
}

.tobs_sbc_sim_occu_multi <- function(spec, theta, seed) {
  set.seed(seed)
  f <- spec$fit_obs
  D <- matrix(theta, nrow = 1L)
  colnames(D) <- names(theta)
  f$draws <- D
  rep <- stats::simulate(f, nsim = 1L)[[1L]]
  obs <- spec$data_obs
  list(cells = obs$cells, y = rep, y_pos = NULL, visits = NULL, graph = NULL,
       site = obs$site)
}

# dyn_int_occu(): the product of the multi-season and multi-source shapes, a
# named list of per-source 3D [n_sites x max_visits_s x T] arrays. The family's
# own `.tobs_simulate_dyn_int_occu()` already reproduces each source's own NA
# pattern, so partial season overlap needs no special-casing. The binder only
# accepts full SITE overlap across sources and has no season-varying option,
# so no scope check is needed.
.tobs_sbc_spec_dyn_int_occu <- .tobs_sbc_spec_fn(
  c(formula = "psi", colonization = "gamma", extinction = "eps",
    detection = "p"),
  response = function(data, spec)
    list(y = data$y[spec$model$sources], sources = spec$model$sources))

.tobs_sbc_data_dyn_int_occu <- function(fit) {
  .tobs_sbc_data_model(fit, stats::setNames(fit$model$y, fit$model$sources))
}

.tobs_sbc_sim_dyn_int_occu <- function(spec, theta, seed) {
  set.seed(seed)
  f <- spec$fit_obs
  D <- matrix(theta, nrow = 1L)
  colnames(D) <- names(theta)
  f$draws <- D
  rep <- stats::simulate(f, nsim = 1L)[[1L]]
  obs <- spec$data_obs
  list(cells = obs$cells, y = rep, y_pos = NULL, visits = NULL, graph = NULL,
       site = obs$site)
}


# ---------------------------------------------------------------------------
# 6d. Two-block multiarm fits: occu_categorical() and cover(). A
# `tobs_multiarm_fit` carries NO `fit$draws` -- its two arms (presence + class,
# presence + positive) are independent Laplace-Gaussian blocks fit as two
# separate GLMs, not a joint MVN pseudo-draw matrix -- so `draws()` samples the
# two blocks independently (they ARE independent by construction) and stacks
# them into one theta row, and `loglik_many()` scores the two-arm likelihood
# directly, since neither family has a `.tobs_pointwise_loglik` dispatch. The
# response is a plain length-N vector (one observation per unit, no visit or
# season axis). Non-spatial by construction; neither takes a `detection`
# formula, so the visit-design half of the rejection gate never fires on them.
# ---------------------------------------------------------------------------

# The covariance-based draw (`mu ~ N(means, vcov)`, or an arm's `beta ~
# N(beta, V)`) exists only on the fitter that actually computes it -- the
# deterministic Laplace-EM / two-Laplace routes. A sampled fit (nuts,
# pg_gibbs) or a nested-engine fit this family's registration was never
# extended to structurally never populates that field; that is a property of
# `fit$method`, not a symptom of a bad fit. Refuse HERE, by name, before the
# NULL/missing value reaches `.tobs_sbc_mvn_draws()` and is misreported as
# "non-converged or rank-deficient".
.tobs_sbc_require_laplace_covariance <- function(fit, family, label, value) {
  if (!is.null(value) && !anyNA(value)) return(invisible(NULL))
  stop(sprintf(paste0(
    "SBC on %s() draws its posterior from the fit's own Laplace covariance ",
    "(%s), which method = '%s' does not provide here -- this is a scope ",
    "limit of the registered draw path for %s(), not a sign this fit failed ",
    "to converge. Refit under the Laplace-EM method this registration reads, ",
    "or wait for a chain-draw registration for %s()."),
    family, label, fit$method %||% "?", family, family), call. = FALSE)
}

# SBC draws the posterior from the fit's own covariance, so a missing one is an
# error here rather than the point mass `.rmvn()` returns. A refit's per-species
# Cinv/vcov can be numerically singular or genuinely indefinite on a
# data-sparse replicate (ms_distance's acceptance run at n.sim = 100: "the
# leading minor of order 3 is not positive" survived a 1e-4-relative ridge
# retry); `.rmvn()` eigen-clips that to the nearest valid covariance.
.tobs_sbc_mvn_draws <- function(mean, V, n) {
  if (is.null(V) || anyNA(V)) {
    stop("SBC draws the posterior from this fit's own covariance, which is ",
         "missing or non-finite here (a non-converged or rank-deficient fit).",
         call. = FALSE)
  }
  .rmvn(n, mean, V)
}

# Independent Gaussian blocks, drawn in order and bound column-wise. Each block
# is `list(mean, V, names)`.
.tobs_sbc_draws_blocks <- function(n, blocks) {
  M <- do.call(cbind, lapply(blocks, function(b)
    .tobs_sbc_mvn_draws(b$mean, b$V, n)))
  colnames(M) <- unlist(lapply(blocks, `[[`, "names"), use.names = FALSE)
  M
}

# The response vector and the unit count sit on the fit's encoding.
.tobs_sbc_data_vector <- function(fit) {
  enc <- fit$encoding
  .tobs_sbc_dataset(enc$data, enc$y, seq_len(enc$N))
}

# occu_categorical(): the encoding's own formula, unrecombined, rebuilt through
# the family's front door on the classes the fit recorded, under the plain
# Laplace route and the caller's control alone.
.tobs_sbc_spec_occu_categorical <- .tobs_sbc_spec_fn(
  list(formula = function(fit) fit$encoding$formula),
  family  = function(fit) occu_categorical(classes = fit$class_labels),
  method  = "laplace",
  control = function(fit, fit.control)
    utils::modifyList(list(), as.list(fit.control)))

# Columns named `occ_<coef>` / `class_<label>_<coef>`, matching `beta_class`'s
# column-major `vec()` order so a theta row unpacks back into the same matrix
# shape it was drawn from.
.tobs_sbc_draws_occu_categorical <- function(fit, n) {
  .tobs_sbc_draws_blocks(n, list(
    list(mean = fit$beta_occ, V = fit$vcov_occ,
         names = paste0("occ_", names(fit$beta_occ))),
    list(mean = as.numeric(fit$beta_class), V = fit$vcov_class,
         names = as.vector(outer(rownames(fit$beta_class),
                                 colnames(fit$beta_class),
                                 function(r, cl) paste0("class_", cl, "_", r))))))
}

.tobs_sbc_theta_occu_categorical <- function(fit, theta) {
  occ_names <- names(fit$beta_occ)
  beta_occ <- theta[paste0("occ_", occ_names)]
  names(beta_occ) <- occ_names
  Beta <- matrix(theta[grepl("^class_", names(theta))],
                 nrow = nrow(fit$beta_class), ncol = ncol(fit$beta_class),
                 dimnames = dimnames(fit$beta_class))
  list(beta_occ = beta_occ, Beta = Beta)
}

.tobs_sbc_sim_occu_categorical <- function(spec, theta, seed) {
  set.seed(seed)
  f <- spec$fit_obs
  th <- .tobs_sbc_theta_occu_categorical(f, theta)
  obs <- spec$data_obs
  X <- stats::model.matrix(spec$args$formula, obs$cells)
  psi <- stats::plogis(as.numeric(X %*% th$beta_occ))
  present <- stats::rbinom(nrow(X), 1L, psi)
  eta <- X %*% th$Beta
  E <- exp(pmin(eta, 700)); denom <- 1 + rowSums(E)
  P <- cbind(E / denom, 1 / denom)
  K <- f$K
  cls <- apply(P, 1L, function(pr) sample.int(K, 1L, prob = pr))
  y <- ifelse(present == 1L, cls, 0L)
  .tobs_sbc_dataset(obs$cells, as.integer(y), obs$site)
}

# The two-arm log-likelihood at each reference theta row, read straight off
# the encoding the refit already carries -- the same closed forms
# `fit_occu_categorical()` fits, evaluated rather than optimized.
.tobs_sbc_loglik_many_occu_categorical <- function(fit, Theta) {
  enc <- fit$encoding
  vapply(seq_len(nrow(Theta)), function(i) {
    th <- .tobs_sbc_theta_occu_categorical(fit, Theta[i, ])
    eta_occ <- enc$X_occ %*% th$beta_occ
    ll_occ <- sum(stats::dbinom(enc$present, 1L, stats::plogis(eta_occ),
                                log = TRUE))
    eta_cls <- enc$X_class %*% th$Beta
    E <- exp(pmin(eta_cls, 700)); denom <- 1 + rowSums(E)
    pmat <- cbind(E / denom, 1 / denom)
    ll_occ + sum(log(pmat[cbind(seq_along(enc$cls), enc$cls)]))
  }, numeric(1))
}

# cover(): `encoding$data` / `encoding$y` are the raw fitting inputs the refit
# call is rebuilt from, and `V_occ` / `V_pos` the full per-arm coefficient
# covariance beside the diagonal `se_occ` / `se_pos` -- needed to draw an arm's
# coefficients JOINTLY rather than assume independence across them, which the
# ms_occu near-miss (#226) is a direct warning against.
#
# Scope: the two-Laplace `.dispatch_cover()` routes (`laplace`, `laplace_sla`),
# the only ones populating `V_occ` / `V_pos`, and `positive = "lognormal"`
# (`simulate_cover()`'s own generator covers lognormal/gaussian, and
# lognormal's `log(y_pos) ~ N(eta_pos, sigma_pos)` is simpler to get right than
# beta's mu/phi -> shape1/shape2 reparameterization). `sigma_pos` is NOT part
# of theta -- the fitted point estimate carries no SE anywhere in the package
# (`sigma_pos_sd` is hardcoded NA), so it is held fixed at the observed fit's
# value for both the replicate and the reference score. `X_occ` / `X_pos` are
# rebuilt fresh via `model.matrix()` from the arm formula + natural-scale data,
# NEVER read off `encoding$occ_data$X` / `pos_data$X` -- those are autoscaled
# design matrices paired with SCALED coefficients internally, while
# `beta_occ` / `beta_pos` on the fit are already unscaled to natural units
# (`.unscale_beta_vec`); mixing the two would reproduce the exact class of bug
# #225 found in `int_occu()`.
.tobs_sbc_spec_cover <- .tobs_sbc_spec_fn(
  list(formula  = function(fit) .tobs_sbc_recombine(fit$encoding$fe_occ, NULL),
       presence = function(fit) .tobs_sbc_recombine(fit$encoding$fe_occ, NULL),
       positive = function(fit) .tobs_sbc_recombine(fit$encoding$fe_pos, NULL)),
  scope = function(fit) {
    .tobs_sbc_scope("cover", identical(fit$positive, "lognormal"),
                    "positive = \"lognormal\"",
                    paste0("\"", fit$positive, "\""),
                    "the beta/beta_oi/lognormal_trunc/ordinal/gaussian arms are follow-ups")
    .tobs_sbc_scope("cover",
                    (fit$method %||% "laplace") %in% c("laplace", "laplace_sla"),
                    "method = \"laplace\"/\"laplace_sla\"",
                    paste0("\"", fit$method, "\""),
                    paste0("it draws the posterior from fit$V_occ/V_pos, the ",
                           "two-Laplace dispatcher's own per-arm covariance, ",
                           "which the nuts / nested_laplace(_sla) routes do not ",
                           "populate -- a scope limit of the registered draw ",
                           "path, not a sign this fit failed to converge"))
  },
  fields = function(fit) list(sigma_pos = fit$sigma_pos))

.tobs_sbc_draws_cover <- function(fit, n) {
  .tobs_sbc_draws_blocks(n, list(
    list(mean = fit$beta_occ, V = fit$V_occ,
         names = paste0("occ_", names(fit$beta_occ))),
    list(mean = fit$beta_pos, V = fit$V_pos,
         names = paste0("pos_", names(fit$beta_pos)))))
}

.tobs_sbc_sim_cover <- function(spec, theta, seed) {
  set.seed(seed)
  obs <- spec$data_obs
  X_occ <- stats::model.matrix(spec$args$presence, obs$cells)
  X_pos <- stats::model.matrix(spec$args$positive, obs$cells)
  beta_occ <- theta[paste0("occ_", colnames(X_occ))]
  beta_pos <- theta[paste0("pos_", colnames(X_pos))]
  eta_occ <- as.vector(X_occ %*% beta_occ)
  eta_pos <- as.vector(X_pos %*% beta_pos)
  n <- nrow(X_occ)
  occur <- stats::rbinom(n, 1L, stats::plogis(eta_occ))
  log_cover <- stats::rnorm(n, eta_pos, spec$sigma_pos)
  cover <- ifelse(occur == 1L, pmin(exp(log_cover), 1), 0)
  .tobs_sbc_dataset(obs$cells, cover, obs$site)
}

.tobs_sbc_loglik_many_cover <- function(fit, Theta) {
  enc <- fit$encoding
  X_occ <- stats::model.matrix(enc$fe_occ, enc$data)
  X_pos <- stats::model.matrix(enc$fe_pos, enc$data)
  occur <- as.integer(enc$y > 0)
  is_pos <- occur == 1L
  y_pos <- enc$y[is_pos]
  X_pos_obs <- X_pos[is_pos, , drop = FALSE]
  sigma_pos <- fit$sigma_pos
  occ_names <- colnames(X_occ); pos_names <- colnames(X_pos)
  vapply(seq_len(nrow(Theta)), function(i) {
    th <- Theta[i, ]
    beta_occ <- th[paste0("occ_", occ_names)]
    beta_pos <- th[paste0("pos_", pos_names)]
    eta_occ <- as.vector(X_occ %*% beta_occ)
    ll_occ <- sum(stats::dbinom(occur, 1L, stats::plogis(eta_occ), log = TRUE))
    eta_pos <- as.vector(X_pos_obs %*% beta_pos)
    ll_occ + sum(stats::dnorm(log(y_pos), eta_pos, sigma_pos, log = TRUE))
  }, numeric(1))
}


# ---------------------------------------------------------------------------
# 6e. t_occu(): a Polya-Gamma Gibbs family, but `.tobs_pg_finalize_fit()`
# already reports the pooled cross-chain draws as `fit$draws` in the standard
# shape (real posterior samples, not a Gaussian pseudo-draw matrix) -- so
# `.tobs_sbc_draws_fit` reads it unchanged, no Gibbs-aware draws() needed.
# `model$y` is `[n_sites x n_seasons x max_visits]` (season BEFORE visits --
# the one family with this axis order), but the shared pool only ever stacks
# axis 1 and leaves the others alone, so it is still exactly right. What
# genuinely IS custom: the family has no `simulate()` handler at all
# (`simulate_t_occu()` draws a FRESH truth, not a replicate at a given theta)
# and no `.tobs_pointwise_loglik` dispatch, so `simulate` is hand-written
# (drawing a fresh AR1 year-effect sequence at the theta's own (sigma, rho) --
# the same generative model `simulate_t_occu()` uses, parameterized by theta
# instead of drawing its own truth). It is also fully-observed only (no
# NA-visit masking, matching `simulate_t_occu()`'s own output).
#
# The log_lik rank arm needs the AR1 year effect eta MARGINALIZED (theta
# carries only its (sigma, rho) hyperparameters, not eta itself), which is a
# real T-dimensional integral -- but a Laplace one, using the SAME Louis
# (1982) mixture identity already established elsewhere in this codebase for
# an additive logit offset into a 2-state mixture (e.g. the ms_dyn_occu psi1
# areal-field oracle: score = w - psi, curv = psi(1-psi) - w(1-w), w = the
# posterior P(z=1|data)). Per season t, eta_t enters every site's occupancy
# logit as a shared additive offset, so summing that identity over sites
# gives the score/curvature of the data term in eta_t; the AR1 prior
# contributes its own tridiagonal precision Q(rho)/sigma^2 (`.t_occu_ar1_Q`,
# the SAME precision the Gibbs sampler's GMRF update already uses). Newton on
# the T x T system (diag(curv) + Q/sigma^2) finds the mode; the Laplace
# log-marginal is the joint log-density at the mode plus the Gaussian
# normalizing correction. FD-validated against the joint-density gradient in
# eta and cross-checked against a brute-force dense-grid sum at T = 2.
# ---------------------------------------------------------------------------

.tobs_sbc_spec_t_occu <- .tobs_sbc_spec_fn(c(formula = "occ", detection = "det"))

.tobs_sbc_sim_t_occu <- function(spec, theta, seed) {
  set.seed(seed)
  m <- spec$model
  occ_names <- m$process_info[[1L]]$coef_names
  det_names <- m$process_info[[2L]]$coef_names
  beta_psi <- theta[paste0("psi_", occ_names)]
  beta_p   <- theta[paste0("p_",   det_names)]
  sigma <- exp(theta[["log_sigma_ar1"]]); rho <- theta[["rho_ar1"]]

  T_s <- m$n_seasons; n <- m$n_sites; J <- m$max_visits
  eta <- numeric(T_s)
  eta[1L] <- stats::rnorm(1, 0, sigma / sqrt(pmax(1 - rho^2, 1e-6)))
  for (t in 2:T_s) eta[t] <- rho * eta[t - 1L] + stats::rnorm(1, 0, sigma)
  eta <- eta - mean(eta)

  lin  <- as.vector(m$X_occ %*% beta_psi)
  p_it <- stats::plogis(as.vector(m$X_det %*% beta_p))
  z <- matrix(0L, n, T_s)
  for (t in seq_len(T_s)) z[, t] <- stats::rbinom(n, 1L, stats::plogis(lin + eta[t]))

  obs <- spec$data_obs
  y <- array(0L, dim = c(n, T_s, J))
  for (t in seq_len(T_s)) for (j in seq_len(J))
    y[, t, j] <- ifelse(z[, t] == 1L, stats::rbinom(n, 1L, p_it), 0L)
  .tobs_sbc_dataset(obs$cells, y, obs$site)
}

# Newton mode-finder for the T-dim AR1 year effect at fixed (lin, p_i, rho,
# sigma). `nvis`/`kdet` are the model's own [n_sites x T] detection
# suff-stats. Returns the mode plus the pieces the Laplace formula needs
# (H = the Gaussian precision at the mode, Pmix = the per-(site,season)
# mixture probability at the mode).
.tobs_t_occu_eta_laplace <- function(lin, p_i, nvis, kdet, T_s, rho, sigma,
                                     max_iter = 50L, tol = 1e-9) {
  n <- length(lin)
  Pmat <- matrix(p_i, n, T_s)
  Bmat <- Pmat^kdet * (1 - Pmat)^(nvis - kdet)
  Amat <- (kdet == 0) * 1
  prec <- .t_occu_ar1_Q(T_s, rho) / sigma^2
  eta <- numeric(T_s)
  H <- prec; psi <- matrix(stats::plogis(lin), n, T_s); Pmix <- Amat
  for (it in seq_len(max_iter)) {
    psi  <- stats::plogis(outer(lin, eta, "+"))
    Pmix <- pmax((1 - psi) * Amat + psi * Bmat, 1e-300)
    w    <- psi * Bmat / Pmix
    score <- colSums(w - psi)
    curv  <- pmax(colSums(psi * (1 - psi) - w * (1 - w)), 1e-8)
    H <- diag(curv, T_s) + prec
    grad <- score - as.vector(prec %*% eta)
    step <- tryCatch(solve(H, grad), error = function(e) NULL)
    if (is.null(step) || anyNA(step)) break
    eta <- eta + step
    if (max(abs(step)) < tol) break
  }
  list(eta = eta, H = H, Pmix = Pmix)
}

# Laplace approximation to log integral_eta P(data|eta,theta) p(eta|theta)
# d(eta) -- see the section header for the derivation. `|Q(rho)| = 1 - rho^2`
# for the standard unit-innovation AR1 precision (`.t_occu_ar1_Q`), and the
# -T/2 log(2 pi) / +T/2 log(2 pi) terms from the prior normalizer and the
# Laplace formula cancel exactly.
.tobs_t_occu_laplace_ll <- function(theta, model) {
  occ_names <- model$process_info[[1L]]$coef_names
  det_names <- model$process_info[[2L]]$coef_names
  beta_psi <- theta[paste0("psi_", occ_names)]
  beta_p   <- theta[paste0("p_",   det_names)]
  sigma <- exp(theta[["log_sigma_ar1"]]); rho <- theta[["rho_ar1"]]
  T_s <- model$n_seasons

  lin <- as.vector(model$X_occ %*% beta_psi)
  p_i <- stats::plogis(as.vector(model$X_det %*% beta_p))
  fitm <- .tobs_t_occu_eta_laplace(lin, p_i, model$nvis, model$kdet, T_s, rho, sigma)
  eta <- fitm$eta
  prec <- .t_occu_ar1_Q(T_s, rho) / sigma^2

  ll_data  <- sum(log(fitm$Pmix))
  ll_prior <- -T_s * log(sigma) + 0.5 * log(pmax(1 - rho^2, 1e-300)) -
    0.5 * as.numeric(crossprod(eta, prec %*% eta))
  logdetH <- as.numeric(determinant(fitm$H, logarithm = TRUE)$modulus)
  ll_data + ll_prior - 0.5 * logdetH
}

.tobs_sbc_loglik_many_t_occu <- function(fit, Theta) {
  m <- fit$model
  vapply(seq_len(nrow(Theta)), function(i)
    .tobs_t_occu_laplace_ll(Theta[i, ], m), numeric(1))
}


# ---------------------------------------------------------------------------
# 6f. occu_multiscale_cover(): a THREE-level hurdle -- cell occupancy psi
# gates plot availability theta, which gates per-visit detection p and the
# cover hurdle pos -- fit under `method = "laplace"` (non-spatial: iid cells,
# field fixed at 0). Unlike `cover()`, this is the STANDARD single-block fit
# shape: `fit$means`/`fit$draws` are a real joint MVN (dispersion
# `log_sigma_pos` optimized jointly with the coefficients via BFGS +
# `optimHess`, so it carries a proper SE, unlike `cover()`'s post-hoc MOM
# estimate) and `.tobs_pointwise_loglik` already dispatches for this
# model_type -- so `draws`/`loglik_many` are the shared generic ones.
#
# What IS custom, and why: the exchangeable UNIT for this family's SBC
# premise is the CELL, not the plot -- z_c (cell occupancy) is the top
# hierarchical draw, and every plot within a cell shares it, so plots are
# not conditionally independent given theta the way one row is in every
# other registered family. `site` in the pooled data is therefore the
# per-plot CELL INDEX (`model$plot_cell`), not a running plot counter, so
# the driver's fresh-groups premise (unique `site` labels) counts unique
# CELLS. Pooling offsets the replicate's cell indices past the observed
# cell COUNT (`model$n_cells`), not the plot count.
#
# `formula` is REQUIRED to carry an `icar(graph =, group.var = "<col>")`
# term declaring the plot -> cell map, even under the non-spatial engine
# (the graph itself is ignored there, but the term's presence and its
# group_var column are not optional syntax). Rather than track and
# reproduce the OBSERVED fit's own group_var column name through pooling, a
# fixed synthetic column (`.sbc_cell`) is written into `cells` by both
# `data()` and `simulate()`, with a same-size dummy graph built at refit
# time -- the group_var mapping the refit reads is entirely the SBC
# adapter's own, decoupled from whatever the original fit used.
#
# `simulate()` has no family handler to inject theta into, so it is a
# hand-written three-level generator: z_c ~ Bernoulli(psi_c) per cell,
# a_j | z_{cell(j)} ~ Bernoulli(theta_j) per plot, y_jv | a_j ~
# Bernoulli(p_j) per visit, cover | y_jv=1 ~ Lognormal(eta_pos_j,
# sigma_pos). `positive = "lognormal"` only, matching `cover()`'s scope
# (no beta-arm generator to mirror, and dispersion carries a real SE here
# so there is no dispersion-fixing complication to work around). No
# visit-level covariates (the rejection gate refuses them).
# ---------------------------------------------------------------------------

.tobs_sbc_spec_occu_mscale_cover <- .tobs_sbc_spec_fn(
  c(formula = "psi", detection = "p", availability = "theta",
    positive = "pos"),
  response = function(data, spec) list(y = data$y, y.pos = data$y_pos),
  scope = function(fit) {
    pos <- fit$model$positive %||% "lognormal"
    .tobs_sbc_scope("occu_multiscale_cover", identical(pos, "lognormal"),
                    "positive = \"lognormal\"", paste0("\"", pos, "\""))
  })

.tobs_sbc_data_occu_mscale_cover <- function(fit) {
  m <- fit$model
  cells <- m$data
  cells$.sbc_cell <- m$plot_cell
  .tobs_sbc_dataset(cells, m$y, m$plot_cell, y_pos = m$y_pos,
                    n_cells = m$n_cells)
}

.tobs_sbc_pool_occu_mscale_cover <- function(obs, rep) {
  n_c_o <- obs$n_cells
  rep$cells$.sbc_cell <- rep$cells$.sbc_cell + n_c_o
  out <- .tobs_sbc_pool(obs, rep)
  out$site <- c(obs$site, rep$site + n_c_o)
  out$n_cells <- n_c_o + rep$n_cells
  out
}

# The state formula with the plot -> cell term re-attached on a dummy graph of
# the pooled cell count; the term is resolved in this frame, which holds it.
.tobs_sbc_refit_occu_mscale_cover <- function(spec, data) {
  n_cells <- max(data$cells$.sbc_cell)
  dummy_graph <- matrix(0, n_cells, n_cells)
  tt <- stats::terms(spec$args$formula)
  rhs <- paste(c(attr(tt, "term.labels"),
                 'icar(graph = dummy_graph, group.var = ".sbc_cell")'),
               collapse = " + ")
  if (!attr(tt, "intercept")) rhs <- paste(rhs, "- 1")
  spec$args$formula <- stats::as.formula(paste("~", rhs), env = environment())
  .tobs_sbc_refit(spec, data)
}

.tobs_sbc_sim_occu_mscale_cover <- function(spec, theta, seed) {
  set.seed(seed)
  m <- spec$model
  psi_names   <- m$process_info[[1L]]$coef_names
  theta_names <- m$process_info[[2L]]$coef_names
  det_names   <- m$process_info[[3L]]$coef_names
  pos_names   <- m$process_info[[4L]]$coef_names
  beta_psi   <- theta[paste0("psi_",   psi_names)]
  beta_theta <- theta[paste0("theta_", theta_names)]
  beta_p     <- theta[paste0("p_",     det_names)]
  beta_pos   <- theta[paste0("pos_",   pos_names)]
  sigma_pos  <- exp(theta[["log_sigma_pos"]])

  n_cells <- m$n_cells; n_plots <- m$n_plots; J <- m$max_visits
  psi_c   <- stats::plogis(as.vector(m$X_psi %*% beta_psi))
  z_c     <- stats::rbinom(n_cells, 1L, psi_c)

  theta_j <- stats::plogis(as.vector(m$X_theta %*% beta_theta))
  a_j     <- stats::rbinom(n_plots, 1L, theta_j) * z_c[m$plot_cell]

  p_j <- stats::plogis(as.vector(m$X_p_site %*% beta_p))
  eta_pos_j <- as.vector(m$X_pos_site %*% beta_pos)

  y <- matrix(0L, n_plots, J)
  y_pos <- matrix(0, n_plots, J)
  for (v in seq_len(J)) {
    y[, v] <- ifelse(a_j == 1L, stats::rbinom(n_plots, 1L, p_j), 0L)
    cov_draw <- exp(stats::rnorm(n_plots, eta_pos_j, sigma_pos))
    y_pos[, v] <- ifelse(y[, v] == 1L, pmin(cov_draw, 1), 0)
  }

  obs <- spec$data_obs
  .tobs_sbc_dataset(obs$cells, y, obs$site, y_pos = y_pos,
                    n_cells = obs$n_cells)
}


# ---------------------------------------------------------------------------
# 6g. The community families: ms_occu(), ms_occu_cover(), ms_int_occu(),
# ms_count(), jsdm(), ms_distance(), ms_dyn_occu(), ms_abun().
#
# Every one ranks a FIXED species set's own coefficients (mu + b_s per
# species), not the community means with fresh species drawn per replicate:
# theta is the per-species REALIZED coefficient vector for every species
# (column names `<species>_<arm>_<coef>`), followed by any shared globals. The
# draw is the joint mu/b_s draw of `.tobs_sbc_community_b_draws()` (#226: mu
# and b_s are NOT independent in the posterior), the replicate is the family's
# own `simulate()` handler with the candidate theta written onto
# `ms_community`'s coefficient matrices (section 6a-bis), and the rank
# statistic sums each species' own exact marginal -- the kernel the fitter
# optimizes. The response is a site-axis-first array or a named list of them,
# pooled by the shared stacker.
#
# SPECIES-COUNT SCOPED. At a small species count the Laplace-Gaussian
# approximation to a genuinely non-Gaussian community posterior breaks down
# and posterior SBC fails hard; at a community size typical of real
# ecological data the same plain Laplace-EM calibrates cleanly. The fixtures
# (`.SBC_REG_FIXTURES`) are sized on that evidence, per family -- do not shrink
# them without re-running the investigation behind the registry comment.
# ---------------------------------------------------------------------------

# mu and b_s are NOT independent in the posterior -- conditional on a draw
# of mu, b_s's mean shifts by -Cinv_s %*% t(Bf_s) %*% (mu_draw - mu_hat)
# (Bf_s the mu-b_s cross-Hessian block from the community EM's own Newton
# solve, R/community_em.R), while its covariance is UNCHANGED at Cinv_s.
# Derived by block-inverting the joint (mu, b_s) arrowhead precision [Sf +
# Bf_s Cinv_s Bf_s', Bf_s; Bf_s', Cinv_s^-1] and validated to machine
# precision against a direct construction of that matrix on random SPD inputs.
# Drawing mu and each b_s independently under-covers on whichever species'
# deviation most strongly trades off against the community mean.
.tobs_sbc_community_b_draws <- function(mu_draws, mu_hat, b_hat_s, Bf_s, Cinv_s, n) {
  eps <- .tobs_sbc_mvn_draws(rep(0, length(b_hat_s)), Cinv_s, n)
  delta <- sweep(mu_draws, 2L, mu_hat, "-")
  shift <- (delta %*% Bf_s) %*% Cinv_s
  sweep(eps - shift, 2L, b_hat_s, "+")
}

# ms_occu(): the exact per-species two-state marginal the fitter itself
# optimizes, reached through the shared multi-source form with a single
# detection arm.
.tobs_sbc_spec_ms_occu <- .tobs_sbc_spec_fn(
  c(formula = "occ", detection = "det"),
  response = .tobs_sbc_response_species)

.tobs_sbc_ll_ms_occu <- function(m, s, b, g, aux, row) {
  .ms_int_occu_sp_ll(as.numeric(m$X_occ %*% b$psi),
                     list(as.numeric(m$X_det %*% b$p)), m$summaries[[s]])
}

# ms_occu_cover(): the occ + p + pos analogue of ms_occu, with ONE extra
# coordinate in `fit$means` / `vcov` beyond the three arms: the shared
# cover-arm dispersion (`log_sigma_pos` for lognormal, `log_phi` for beta), a
# single community-level scalar with no per-species deviation. `Bf_s` is
# `(P+1) x P` (mu+disp -> b_s), so the conditional shift already reads the
# disp coordinate, while the disp value itself is copied straight off the
# (mu, disp) draw. The replicate handler reads `ms_community$coef_occ` /
# `coef_p` / `coef_pos` (note "occ", not the "psi" prefix of `fit$means`) and
# the last element of `fit$means`. `y` and `y_pos` are both
# `[site x visit x species]`. Non-spatial `laplace`, `positive = "lognormal"`
# only (beta's mu/phi reparameterization is a separate correctness question).
.tobs_sbc_spec_ms_occu_cover <- .tobs_sbc_spec_fn(
  c(formula = "occ", detection = "det", positive = "pos"),
  response = function(data, spec)
    list(y = data$y, y.pos = data$y_pos, species = spec$model$species_names),
  scope = function(fit)
    .tobs_sbc_scope("ms_occu_cover", identical(fit$model$positive, "lognormal"),
                    "positive = \"lognormal\"",
                    paste0("\"", fit$model$positive, "\"")))

.tobs_sbc_data_ms_occu_cover <- function(fit) {
  .tobs_sbc_data_model(fit, y_pos = fit$model$y_pos)
}

.tobs_sbc_setup_ms_occu_cover <- function(fit) {
  lapply(seq_len(fit$model$n_species),
         function(s) .ms_occu_cover_species_view(fit$model, s))
}

.tobs_sbc_ll_ms_occu_cover <- function(m, s, b, g, aux, row) {
  .occu_cover_sp_ll(aux[[s]], b$occ, b$p, b$pos, g$disp)
}

# ms_int_occu(): the community analogue of int_occu(), a list of per-source 3D
# [n_sites x visits_d x n_species] arrays. Named by `process_names`
# ("p1".."pD"), NOT `source_names`: the family's own
# `.tobs_simulate_ms_int_occu()` names its replicate list that way, and
# `source_names` defaults to a DIFFERENT label ("src1".."srcD") whenever `y`
# came in unnamed, which would silently mismatch the two lists' keys at
# pooling. Full site-overlap only, as for int_occu().
.tobs_sbc_spec_ms_int_occu <- .tobs_sbc_spec_fn(
  list(formula   = "occ",
       detection = function(fit)
         lapply(fit$model$formulas$det, function(f) .tobs_sbc_recombine(f, NULL))),
  response = function(data, spec)
    list(y = data$y[spec$model$process_names],
         species = spec$model$species_names),
  scope = function(fit)
    .tobs_sbc_full_overlap("ms_int_occu", fit$model$site_maps, fit$model$n_sites))

.tobs_sbc_data_ms_int_occu <- function(fit) {
  .tobs_sbc_data_model(fit, stats::setNames(fit$model$y, fit$model$process_names))
}

# One occupancy arm plus one detection arm per source. `p<d>` matches
# `model$process_names`, the SAME labels `ms_community`'s `coef_p<d>` /
# `blup_p<d>` fields are keyed by (build_ms_int_occu_fit()).
.tobs_sbc_arms_ms_int_occu <- function(m) {
  c(list(.tobs_sbc_arm("psi", "psi_", 1L)),
    lapply(seq_len(m$n_sources), function(d) {
      nm <- paste0("p", d)
      .tobs_sbc_arm(nm, paste0(nm, "_"), 1L + d)
    }))
}

# The multi-source per-species two-state marginal, the same kernel `ms_occu`
# reaches with a single detection arm.
.tobs_sbc_ll_ms_int_occu <- function(m, s, b, g, aux, row) {
  .ms_int_occu_sp_ll(
    as.numeric(m$X_psi %*% b$psi),
    lapply(seq_len(m$n_sources),
           function(d) as.numeric(m$X_p[[d]] %*% b[[paste0("p", d)]])),
    m$summaries[[s]])
}

# ms_count(): community relative-abundance GLMM, no detection or latent state;
# `model$y` is a plain [n_sites x n_species] matrix. `response = "poisson"`
# only (negbin/gaussian/bernoulli/binomial carry an extra per-species
# dispersion arm). jsdm() is `ms_count(response = "bernoulli")` under the hood
# -- the same `fit$model` shape and the same community Laplace-EM -- so its
# spec is ms_count()'s minus the response check, rebuilt through `jsdm()`, and
# only the kernel differs.
.TOBS_SBC_MS_COUNT_FORMULAS <- c(formula = "mu")

.tobs_sbc_spec_ms_count <- .tobs_sbc_spec_fn(
  .TOBS_SBC_MS_COUNT_FORMULAS, response = .tobs_sbc_response_species,
  scope = function(fit)
    .tobs_sbc_scope("ms_count", identical(fit$model$response, "poisson"),
                    "response = \"poisson\"",
                    paste0("\"", fit$model$response, "\""),
                    "negbin/gaussian/bernoulli/binomial are follow-ups"))

.tobs_sbc_spec_jsdm <- .tobs_sbc_spec_fn(
  .TOBS_SBC_MS_COUNT_FORMULAS, response = .tobs_sbc_response_species,
  family = function(fit) jsdm())

.tobs_sbc_ll_ms_count_pois <- function(m, s, b, g, aux, row) {
  .ms_count_ll_pois(m$summaries[[s]], b$mu)
}

.tobs_sbc_ll_ms_count_bern <- function(m, s, b, g, aux, row) {
  .ms_count_ll_bern(m$summaries[[s]], b$mu)
}

# ms_distance(): community binned distance sampling (msDS), the two-arm
# (lambda/sigma) sibling of ms_occu; `model$y` is `[n_sites x n_bins x
# n_species]`. The replicate handler draws through `cpp_simulate_distance`,
# the SAME kernel the likelihood integrates against, and the rank statistic
# sums each species' binned-distance marginal off the family's own sweep.
# `key = "halfnorm"`, `mixture = "poisson"` only, matching
# `.tobs_simulate_ms_distance()`'s own restriction (the hazard key's log-shape
# rides the community EM's `global` slot, not `ms_community`, so neither
# `simulate()` nor this registration can reach it -- #227).
.tobs_sbc_spec_ms_distance <- .tobs_sbc_spec_fn(
  c(formula = "lambda", detection = "sigma"),
  response = .tobs_sbc_response_species,
  family = function(fit)
    ms_distance(key = fit$model$key, transect = fit$model$transect,
                cutpoints = fit$model$cutpoints, mixture = fit$model$mixture),
  scope = function(fit) {
    m <- fit$model
    .tobs_sbc_scope("ms_distance",
                    identical(m$key, "halfnorm") && identical(m$mixture, "poisson"),
                    "key = \"halfnorm\", mixture = \"poisson\"",
                    sprintf("key = \"%s\", mixture = \"%s\"", m$key, m$mixture),
                    paste0("the hazard key's log-shape and negbin's dispersion ",
                           "are not carried by simulate()"))
  })

.tobs_sbc_setup_ms_distance <- function(fit) {
  m <- fit$model
  list(eng = .tobs_ms_distance_engine(m),
       X_lambda = m$X_processes[[1L]], X_sigma = m$X_processes[[2L]])
}

.tobs_sbc_ll_ms_distance <- function(m, s, b, g, aux, row) {
  sw <- aux$eng$sweep(s, as.numeric(aux$X_lambda %*% b$lambda),
                      as.numeric(aux$X_sigma %*% b$sigma), value_only = TRUE)
  sum(sw$log_lik)
}

# ms_dyn_occu(): community DYNAMIC occupancy, `model$y` a 4D `[n_sites x
# max_visits x n_seasons x n_species]` array. Two arms carry a per-species
# deviation (psi1, p); the season-to-season colonization / extinction
# coefficients (gamma, eps) are SHARED community globals with no per-species
# term. `fit$means` / `fit$vcov` carry all four jointly, and
# `ms_community$Bf[[s]]` is `(P + G) x P` -- the cross-Hessian block against
# the WHOLE (mu, global) vector, since the community EM's Newton solve finds
# gamma/eps jointly informative about each species' deviation even though they
# carry no random effect of their own -- so each b_s conditions on the full
# draw and the globals copy straight out of it. The rank statistic is each
# species' exact HMM forward at its psi1/p slice plus the shared gamma/eps.
.tobs_sbc_spec_ms_dyn_occu <- .tobs_sbc_spec_fn(
  c(formula = "occ", detection = "det", colonization = "col",
    extinction = "ext"),
  response = .tobs_sbc_response_species)

# The transition probabilities depend only on the SHARED gamma / eps blocks, so
# they are built once per theta row rather than once per species.
.tobs_sbc_row_ms_dyn_occu <- function(m, g, aux) {
  list(gamma = as.numeric(stats::plogis(m$X_gamma %*% g$gamma)),
       eps   = as.numeric(stats::plogis(m$X_eps   %*% g$eps)))
}

.tobs_sbc_setup_ms_dyn_occu <- function(fit) {
  m <- fit$model
  lapply(seq_len(m$n_species), function(s)
    .ms_dyn_occu_emit_stats(m$y[, , , s], m$valid[, , , s],
                            m$n_sites, m$n_seasons))
}

.tobs_sbc_ll_ms_dyn_occu <- function(m, s, b, g, aux, row) {
  p <- as.numeric(stats::plogis(m$X_p %*% b$p))
  .ms_dyn_occu_fwd_ll_vec(
    as.numeric(stats::plogis(m$X_psi1 %*% b$psi1)),
    row$gamma, row$eps,
    .ms_dyn_occu_emissions(p, aux[[s]]$nvalid, aux[[s]]$ndet),
    m$n_sites, m$n_seasons)
}

# ms_abun(): community N-mixture (msNMix). Unlike every other community family
# here, its engine is `tulpa::tulpa_re_aghq()`, which exposes the per-species
# posterior covariance and the lambda/p cross-covariance (`blup_cross` /
# `blup_cov_g`) that N-mixture's own lambda/p identifiability ridge makes
# genuinely nonzero. They are only present when the fit ran the AGHQ/joint path
# or the default Laplace-EM (`tulpa::tulpa_re_em()`), so the scope check
# requires them explicitly rather than falling back to an independent
# (mu, b_s) draw. The refit runs under the ORIGINAL fit's optimizer and
# quadrature, so every replicate's posterior is the one the observed fit
# reported. Poisson only (negbin/zip/zinb add a further per-species RE arm
# this theta layout does not carry).
.tobs_sbc_reject_ms_abun_scope <- function(fit) {
  .tobs_sbc_scope("ms_abun", identical(fit$mixture, "poisson"),
                  "mixture = \"poisson\"", paste0("\"", fit$mixture, "\""),
                  paste0("negbin/zip/zinb add a further per-species RE arm ",
                         "this registration does not cover"))
  if (is.null(fit$ms_community$Cinv) || is.null(fit$ms_community$Bf)) {
    stop("SBC on ms_abun() needs the per-species posterior covariance and ",
         "mode/theta cross-Hessian (fit$ms_community$Cinv/Bf), which this fit ",
         "does not carry; a joint (mu, b_s) draw cannot be built without them. ",
         "Refit with the default Laplace-EM or control = list(optimizer = ",
         "\"joint_fd\") rather than drawing mu and b_s independently.",
         call. = FALSE)
  }
}

.tobs_sbc_spec_ms_abun <- .tobs_sbc_spec_fn(
  c(formula = "lambda", detection = "det"),
  response = .tobs_sbc_response_species,
  family = function(fit) ms_abun(mixture = fit$mixture),
  scope = .tobs_sbc_reject_ms_abun_scope,
  control = function(fit, fit.control) {
    ctl <- .tobs_sbc_control(fit, fit.control)
    ctl$optimizer <- fit$ms_community$optimizer
    ctl[["n.quad"]] <- fit$ms_community$n_quad
    ctl
  })

# The community N-mixture kernel: each species' own exact closed-form Royle
# marginal, the SAME `nmix_site_marginal()` the fitter optimizes.
.tobs_sbc_setup_ms_abun <- function(fit) {
  m  <- fit$model
  lf <- .tobs_ms_nmix_longform(m)
  list(lf = lf, X_lambda = m$X_processes[[1L]],
       K_max = fit$K_max %||% (max(lf$y) + 100L))
}

.tobs_sbc_ll_ms_abun <- function(m, s, b, g, aux, row) {
  k  <- aux$lf$species_idx == s
  Xp <- aux$lf$X_p[k, , drop = FALSE]
  nmix_site_marginal(aux$lf$y[k], aux$lf$site_idx[k], aux$X_lambda, Xp,
                     mixture = "P", K_max = aux$K_max)$eval(
                       as.numeric(aux$X_lambda %*% b$lambda),
                       as.numeric(Xp %*% b$p))$log_lik
}


# ---------------------------------------------------------------------------
# 7. The registry
#
# A family is one row. Everything not named here -- pooling, the grouping
# labels, the arm construction, the mis-scaled controls -- is the shared
# driver's; the per-simulation seed split is `tulpa::sbc()`'s.
# ---------------------------------------------------------------------------

.TOBS_SBC_REGISTRY <- list(
  occu_cover = list(
    draws    = .tobs_sbc_draws_occu_cover,
    simulate = .tobs_sbc_sim_occu_cover,
    refit    = .tobs_sbc_refit_occu_cover,
    loglik   = .tobs_sbc_loglik_occu_cover),

  occu          = .tobs_sbc_simple_entry("occ",    "det"),
  abun          = .tobs_sbc_simple_entry("lambda", "det"),
  removal       = .tobs_sbc_simple_entry("lambda", "det"),
  distance      = .tobs_sbc_simple_entry("lambda", "sigma"),
  royle_nichols = .tobs_sbc_simple_entry("lambda", "r"),
  occu_ttd      = .tobs_sbc_simple_entry("psi",    "rate"),
  fp_occu       = .tobs_sbc_simple_entry(
    "psi", "p11",
    extra = function(m) list(p10 = m$formulas$p10, certainty = m$formulas$b)),
  # `.tobs_simulate_count()` (R/count_methods.R) is registered as count()'s
  # `simulate()` handler, so the shared replicate already draws at the sampled
  # theta row the SBC driver writes into `f$draws`.
  count         = .tobs_sbc_simple_entry("occ", resp = "y_count", y_vector = TRUE),
  # Independent-observer double_observer(): both p1/p2 read ONE `detection`
  # formula slot (`fit$model$formulas$p`). The dependent (`type =
  # "dependent"`) protocol additionally carries a fixed per-site `primary`
  # assignment the front door takes as `...`, threaded through as `extra`.
  double_observer = .tobs_sbc_simple_entry(
    "lambda", "p",
    extra = function(m) if (!is.null(m$primary)) list(primary = m$primary)
                        else NULL),

  dyn_occu      = .tobs_sbc_entry(.tobs_sbc_spec_dyn_occu,
                                  simulate = .tobs_sbc_sim_dyn_occu),
  dyn_abun      = .tobs_sbc_entry(.tobs_sbc_spec_dyn_abun),
  distsamp_open = .tobs_sbc_entry(.tobs_sbc_spec_distsamp_open),
  t_occu        = .tobs_sbc_entry(.tobs_sbc_spec_t_occu,
                                  simulate = .tobs_sbc_sim_t_occu,
                                  loglik_many = .tobs_sbc_loglik_many_t_occu),

  int_occu      = .tobs_sbc_entry(.tobs_sbc_spec_int_occu,
                                  data = .tobs_sbc_data_int_occu,
                                  simulate = .tobs_sbc_sim_int_occu),
  gdistremoval  = .tobs_sbc_entry(.tobs_sbc_spec_gdistremoval,
                                  data = .tobs_sbc_data_gdistremoval,
                                  simulate = .tobs_sbc_sim_gdistremoval),
  occu_multi    = .tobs_sbc_entry(.tobs_sbc_spec_occu_multi,
                                  data = .tobs_sbc_data_occu_multi,
                                  simulate = .tobs_sbc_sim_occu_multi),
  dyn_int_occu  = .tobs_sbc_entry(.tobs_sbc_spec_dyn_int_occu,
                                  data = .tobs_sbc_data_dyn_int_occu,
                                  simulate = .tobs_sbc_sim_dyn_int_occu),

  occu_categorical = .tobs_sbc_entry(
    .tobs_sbc_spec_occu_categorical,
    data        = .tobs_sbc_data_vector,
    draws       = .tobs_sbc_draws_occu_categorical,
    simulate    = .tobs_sbc_sim_occu_categorical,
    loglik_many = .tobs_sbc_loglik_many_occu_categorical),
  cover = .tobs_sbc_entry(
    .tobs_sbc_spec_cover,
    data        = .tobs_sbc_data_vector,
    draws       = .tobs_sbc_draws_cover,
    simulate    = .tobs_sbc_sim_cover,
    loglik_many = .tobs_sbc_loglik_many_cover),

  occu_multiscale_cover = .tobs_sbc_entry(
    .tobs_sbc_spec_occu_mscale_cover,
    data     = .tobs_sbc_data_occu_mscale_cover,
    pool     = .tobs_sbc_pool_occu_mscale_cover,
    simulate = .tobs_sbc_sim_occu_mscale_cover,
    refit    = .tobs_sbc_refit_occu_mscale_cover),

  # ms_occu(): at S=5 posterior SBC collapsed (p_unif as low as 0, several
  # quantities below 1e-3). AGHQ variance-component debiasing and a
  # from-scratch Vf/Cinv/Bf consistency fix were BOTH ruled out as the cause; a
  # direct comparison against `method = "nuts"`'s exact joint posterior (Rhat
  # 1.011, ESS 523, 0 divergences) showed Laplace's Vf and Cinv 3-15x too
  # narrow AND its point estimates measurably off at S=5 -- the
  # Laplace-Gaussian approximation breaking down, not a subtle attenuation. At
  # S=20 the same plain Laplace-EM calibrates cleanly: 5 seeds, min p_unif
  # 0.0017-0.032, 0 of 81 quantities below 1e-3, no reproducible failing
  # coefficient. The fixture is deliberately S=20.
  ms_occu = .tobs_sbc_community_entry(
    arms  = list(.tobs_sbc_arm("psi", "psi_", 1L),
                 .tobs_sbc_arm("p",   "p_",   2L)),
    ll    = .tobs_sbc_ll_ms_occu,
    spec  = .tobs_sbc_spec_ms_occu),
  # ms_occu_cover(): safe to attempt because its Cinv is kept consistent with
  # Sigma even under AGHQ debiasing (commit `03b87ad`).
  ms_occu_cover = .tobs_sbc_community_entry(
    arms    = list(.tobs_sbc_arm("occ", "psi_", 1L),
                   .tobs_sbc_arm("p",   "p_",   2L),
                   .tobs_sbc_arm("pos", "pos_", 3L)),
    globals = list(.tobs_sbc_global_named("disp", .tobs_sbc_disp_name)),
    ll      = .tobs_sbc_ll_ms_occu_cover,
    setup   = .tobs_sbc_setup_ms_occu_cover,
    y_pos   = TRUE,
    spec    = .tobs_sbc_spec_ms_occu_cover,
    data    = .tobs_sbc_data_ms_occu_cover),
  # ms_int_occu(): ms_occu's failure mode at a small species count (a
  # three-seed probe pins `sp3_p2_(Intercept)` at p_unif ~0.003). At S=14
  # (this family's own recovery-test fixture, N=140, 2 sources) the plain
  # Laplace-EM calibrates cleanly: 5 seeds, min p_unif 0.0013-0.052, 0 of 43
  # quantities below 1e-3.
  ms_int_occu = .tobs_sbc_community_entry(
    arms  = .tobs_sbc_arms_ms_int_occu,
    ll    = .tobs_sbc_ll_ms_int_occu,
    spec  = .tobs_sbc_spec_ms_int_occu,
    data  = .tobs_sbc_data_ms_int_occu),
  # ms_count(): a small fixture pins `sp3_mu_(Intercept)` at p_unif ~1e-6 on
  # seeds 0, 1, 2 -- worse than ms_occu / ms_int_occu, plausibly because there
  # is no detection arm to dilute it. At S=20 it calibrates cleanly: 5 seeds,
  # min p_unif 0.0016-0.086, 0 of 41 quantities below 1e-3. jsdm() shares the
  # engine and the S=20 fixture.
  ms_count = .tobs_sbc_community_entry(
    arms  = list(.tobs_sbc_arm("mu", "mu_", 1L)),
    ll    = .tobs_sbc_ll_ms_count_pois,
    spec  = .tobs_sbc_spec_ms_count),
  jsdm = .tobs_sbc_community_entry(
    arms  = list(.tobs_sbc_arm("mu", "mu_", 1L)),
    ll    = .tobs_sbc_ll_ms_count_bern,
    spec  = .tobs_sbc_spec_jsdm),
  ms_distance = .tobs_sbc_community_entry(
    arms  = list(.tobs_sbc_arm("lambda", "lambda_", 1L),
                 .tobs_sbc_arm("sigma",  "sigma_",  2L)),
    ll    = .tobs_sbc_ll_ms_distance,
    setup = .tobs_sbc_setup_ms_distance,
    spec  = .tobs_sbc_spec_ms_distance),
  ms_dyn_occu = .tobs_sbc_community_entry(
    arms    = list(.tobs_sbc_arm("psi1", "psi1_", 1L),
                   .tobs_sbc_arm("p",    "p_",    2L)),
    globals = list(.tobs_sbc_global("gamma", "gamma_", 3L),
                   .tobs_sbc_global("eps",   "eps_",   4L)),
    ll      = .tobs_sbc_ll_ms_dyn_occu,
    setup   = .tobs_sbc_setup_ms_dyn_occu,
    row     = .tobs_sbc_row_ms_dyn_occu,
    spec    = .tobs_sbc_spec_ms_dyn_occu),
  ms_abun = .tobs_sbc_community_entry(
    arms  = list(.tobs_sbc_arm("lambda", "lambda_", 1L),
                 .tobs_sbc_arm("p",      "p_",      2L)),
    ll    = .tobs_sbc_ll_ms_abun,
    setup = .tobs_sbc_setup_ms_abun,
    spec  = .tobs_sbc_spec_ms_abun)
)

# Every entry supplies the callbacks the driver reads. `loglik` / `loglik_many`
# are the optional rank arm; `spec` / `data` / `pool` fall back to the shared
# ones.
.TOBS_SBC_REQUIRED <- c("draws", "simulate", "refit")

.tobs_sbc_registered <- function() sort(names(.TOBS_SBC_REGISTRY))


# ---------------------------------------------------------------------------
# 8. The shared driver
# ---------------------------------------------------------------------------

# The predictives one fit reports, plus the deliberately mis-scaled controls.
#
# A control arm is the SAME draw set rescaled about its own mean, so it differs
# from the reported posterior in width alone and in nothing else. A calibration
# band that nothing can fail is not evidence, so a run that declares controls is
# asserting that the instrument still reacts.
.tobs_sbc_arms <- function(fit, data, reg, n.draws, n.ref, controls,
                           bad.factor, scored) {
  M <- reg$draws(fit, n.draws)
  blk <- attr(M, "blocks")
  qs <- intersect(scored, colnames(M))
  post <- lapply(qs, function(q) tulpa::sbc_draws(M[, q]))
  names(post) <- qs

  if ((!is.null(reg$loglik) || !is.null(reg$loglik_many)) && n.ref > 0L) {
    idx <- if (nrow(M) <= n.ref) seq_len(nrow(M)) else
      round(seq(1, nrow(M), length.out = n.ref))
    # The statistic is evaluated on the FULL parameter vector, references and
    # truth alike; restricting to the scored columns would leave a nuisance the
    # two sides then carry at different values.
    mk <- function(i) {
      th <- M[i, ]; names(th) <- colnames(M)
      attr(th, "blocks") <- blk
      th
    }
    th0 <- data$theta
    if (!is.null(reg$loglik_many)) {
      # A vectorized statistic scores the whole reference set in one call; the
      # truth goes through the SAME call on a one-row matrix, so the two sides
      # cannot reach it by different routes. It is ordered by NAME against the
      # reference columns, since the statistic reads its layout positionally.
      T0 <- matrix(th0[colnames(M)], nrow = 1L,
                   dimnames = list(NULL, colnames(M)))
      ll_ref <- reg$loglik_many(fit, M[idx, , drop = FALSE])
      ll0    <- reg$loglik_many(fit, T0)
    } else {
      ll_ref <- vapply(idx, function(i) reg$loglik(fit, mk(i)), numeric(1))
      ll0 <- reg$loglik(fit, th0)
    }
    post$log_lik <- tulpa::sbc_rank(sum(ll_ref < ll0), length(ll_ref))
  }

  arms <- list(posterior = post)
  scale_arm <- function(mult) {
    a <- lapply(qs, function(q) {
      v <- M[, q]
      tulpa::sbc_draws(mean(v) + mult * (v - mean(v)))
    })
    names(a) <- qs
    a
  }
  if ("wide" %in% controls)   arms$wide   <- scale_arm(bad.factor)
  if ("narrow" %in% controls) arms$narrow <- scale_arm(1 / bad.factor)
  arms
}

# A dispersion the joint engine holds FIXED is set per data set, before the
# outer grid, so the observed fit and every augmented refit hold it at different
# values. The replicate is then generated at one dispersion and scored under
# another, which biases the cover-arm ranks without anything in the result
# saying so. Putting it on the outer grid estimates it, removes the mismatch,
# and adds it to what is scored.
.tobs_sbc_check_fixed_dispersion <- function(fit, fixed) {
  if (!"disp" %in% fixed) return(invisible(NULL))
  if (is.null(fit$model$cover_pos_disp)) return(invisible(NULL))
  warning("this fit holds the cover dispersion fixed at ",
          signif(fit$model$cover_pos_disp, 4),
          ", a value set per data set rather than estimated, so the replicate ",
          "is generated at the observed fit's value while each refit uses its ",
          "own. Put it on the outer grid -- fit.control = list(phi.grid.pos = ",
          "<grid>) -- to estimate and score it instead.", call. = FALSE)
  invisible(NULL)
}

.tobs_sbc_build_model <- function(object, n.draws, n.ref, controls, bad.factor,
                                  fit.control) {
  fam <- attr(object, "tobs_family")$name
  reg <- .TOBS_SBC_REGISTRY[[fam]]
  if (is.null(reg)) {
    stop("SBC is not registered for family ", sQuote(fam), ". Registered: ",
         paste(.tobs_sbc_registered(), collapse = ", "),
         ". Add an entry to `.TOBS_SBC_REGISTRY` (a simulate() and, for the ",
         "rank arm, a log-likelihood) to cover it.", call. = FALSE)
  }
  spec <- (reg$spec %||% .tobs_sbc_spec)(object, fit.control)
  spec$model_graph <- object$spatial$graph
  # The replicate's field is drawn at the scale the engine's own block carries,
  # resolved ONCE here so every generator call reads one constant (section 5).
  spec$field_scale <- .tobs_sbc_field_scale(spec$model_graph)
  spec$data_obs <- (reg$data %||% .tobs_sbc_data_from_fit)(object)

  probe <- reg$draws(object, .TOBS_SBC_PROBE_DRAWS)
  scored <- .tobs_sbc_scored(probe)
  fixed <- setdiff(colnames(probe), scored)
  # Whether the shared field is copied onto the cover arm is read off the fit's
  # own draws: the copy is what gives `alpha` a posterior at all, and without it
  # the cover arm never sees the field.
  spec$has_copy <- "alpha" %in% scored
  spec$alpha_request <- object$spatial$alpha_request
  .tobs_sbc_check_fixed_dispersion(object, fixed)

  structure(list(
    data_obs = spec$data_obs,
    fit = function(data) reg$refit(spec, data),
    draw_theta = function(fit, seed) {
      set.seed(seed)
      M <- reg$draws(fit, 1L)
      th <- as.numeric(M[1L, ])
      names(th) <- colnames(M)
      attr(th, "blocks") <- attr(M, "blocks")
      th
    },
    simulate = function(theta, seed) reg$simulate(spec, theta, seed),
    pool = reg$pool %||% .tobs_sbc_pool,
    arms = function(fit, data) .tobs_sbc_arms(fit, data, reg, n.draws, n.ref,
                                              controls, bad.factor, scored),
    group_ids = .tobs_sbc_groups),
    family = fam, fixed = fixed, quantities = scored)
}


#' Simulation-based calibration for a fitted tobs model
#'
#' Runs the POSTERIOR simulation-based calibration experiment (Talts et al.
#' 2018; Sailynoja, Schmitt, Buerkner and Vehtari 2026, Algorithm 2) on a
#' fitted `tobs_fit`, through the engine's `tulpa::sbc()` front door. The truth
#' is drawn from the fit's own posterior at the observed data, a replicate is
#' simulated at that truth on FRESH cells, the model is refitted on the two data
#' sets together, and the truth's rank under that augmented posterior is
#' recorded. Under exact inference those ranks are Uniform(0, 1); the departure
#' from uniformity is a verdict on the approximation.
#'
#' @section What is scored:
#' The arm fixed effects, the areal field SD on each arm, the copy scale
#' `alpha`, and the dispersion when the fit estimates one. The per-cell field
#' itself is not scored -- it is integrated out by the fit and carries no truth
#' theta could hold. The cover-arm field SD and `alpha` are read per draw off
#' the outer hyperparameter grid, whose cells are sampled by their own
#' normalized weight, so what is ranked is the grid-marginalized posterior of
#' each rather than a function of component modes.
#'
#' A quantity whose posterior has no spread -- a dispersion the joint engine
#' holds fixed in its cell-coupling spec, for instance -- is dropped rather than
#' scored, and named in `attr(x, "fixed")`.
#'
#' An extra `log_lik` quantity ranks a joint statistic of the whole parameter
#' vector, which catches a posterior getting each marginal right while getting
#' the dependence between them wrong.
#'
#' @section Fresh cells:
#' The replicate is drawn on new cells whose areal graph is block-diagonal
#' against the observed one. `occu_cover` couples its occurrence and cover arms
#' through one shared field, and theta carries no per-cell field value, so a
#' replicate re-observing the observed cells would depend on the observed data
#' given theta and break the factorization posterior SBC rests on. `group_ids`
#' is supplied to `tulpa::sbc()`, which verifies the observable half of that
#' (disjoint group labels) rather than assuming it.
#'
#' @section Controls:
#' `controls = c("wide", "narrow")` adds arms reporting the same draws rescaled
#' about their mean by `bad.factor` and its reciprocal. They are deliberately
#' mis-scaled posteriors and are expected to leave the band: a calibration read
#' that nothing can fail is not evidence.
#'
#' @param object A fitted `tobs_fit`. Its family must be registered; see
#'   Details for the roster.
#' @param n.sim Number of simulations. Each costs one refit on the pooled data.
#' @param n.draws Posterior draws per fit backing the reported predictives.
#' @param n.ref Reference draws behind the `log_lik` rank. `0` drops that arm.
#' @param quantities Optional character vector restricting what is scored.
#' @param controls Mis-scaled control arms to add: any of `"wide"`, `"narrow"`.
#' @param bad.factor The control arms' SD multiplier.
#' @param level Simultaneous band level.
#' @param seed Seed offset. Simulation `i` draws its truth at `seed + i`.
#' @param model.only Return the callback list instead of running the
#'   experiment, for inspection or for passing to `tulpa::sbc()` directly.
#' @param fit.control Merged into the `control` list of every refit.
#' @param control Passed to `tulpa::sbc()` (`progress`, `rand_seed`).
#' @param ... Unused.
#'
#' @details
#' Registered families: `occu_cover` (the coupled occupancy + cover hurdle on
#' the joint nested-Laplace engine), and `occu`, `abun`, `count`, `removal`,
#' `distance`, `fp_occu`, `royle_nichols`, `occu_ttd` and `double_observer`,
#' whose site marginals multiply. For those the replicate is the family's own
#' `simulate()` kernel at the drawn theta, the rank arm is the family's exact marginal
#' log-likelihood, and two kinds of fit are refused rather than approximated: a
#' structured term (its field is a latent quantity shared across sites that
#' theta does not hold, which is what the coupled `occu_cover` route handles by
#' drawing on fresh cells) and a visit-level observation design (the replicate
#' kernel assembles the observation arm from the site-level design, so the two
#' sides would carry different models).
#'
#' `dyn_occu` (constant-rate only) pools on the SITE axis and leaves the
#' season axis alone, since its response is 3D
#' `[n_sites x max_visits x n_seasons]`: a bespoke replicate generator
#' forward-simulates the colonization/extinction HMM, but the posterior
#' draws and the rank statistic are the same shared machinery every other
#' registered family uses.
#'
#' `int_occu` (full source-site overlap only) pools on the SITE axis and
#' leaves the per-source axis alone, since its response is a list of one
#' detection matrix per source; the rank statistic is the shared machinery
#' unchanged.
#'
#' `gdistremoval` pools two response matrices (`yDist` band counts, `yRem`
#' period counts) together on the SITE axis, reusing the same
#' list-of-matrices `pool` `int_occu` uses; its replicate generator is its
#' own `simulate()` handler.
#'
#' `dyn_abun` (constant-rate only) shares `dyn_occu`'s 3D
#' `[n_sites x max_visits x n_seasons]` response and site-axis pooling, but
#' unlike `dyn_occu` it already has a working `simulate()` handler (the
#' Dail-Madsen open N-mixture forward is not a two-state HMM that needs a
#' bespoke generator), so its replicate is the same family-`simulate()` route
#' as `occu`/`abun`/etc.
#'
#' `occu_categorical` is a `tobs_multiarm_fit` with two independent
#' Laplace-Gaussian blocks (presence, class) rather than a joint MVN draw
#' matrix; `draws()` samples the two blocks independently and `loglik_many()`
#' scores the two-arm likelihood directly, since this family has no
#' `.tobs_pointwise_loglik` dispatch to reuse.
#'
#' `distsamp_open` (constant-dynamics, Poisson only) shares `dyn_abun`'s 3D
#' response and site-axis pooling, and its `fit$means`/`fit$draws` are the
#' standard `.tobs_bfgs_marginal_fit()` shape, so `draws`/`simulate`/
#' `loglik_many` are all the shared generic ones.
#'
#' `occu_multi`'s response is a list of S per-species matrices, the same
#' shape `int_occu()` pools; its `simulate()` is custom, since species share
#' one joint multi-species state rather than being independently observed.
#'
#' `dyn_int_occu` is the product of the multi-season and multi-source shapes:
#' a named list of S per-source 3D arrays, each stacked on its site axis.
#'
#' `t_occu` is a `pg_gibbs` family whose `fit$draws` is already the real
#' pooled posterior sample; it has no `simulate()` handler and no
#' `.tobs_pointwise_loglik` dispatch, so both are custom -- `simulate` draws
#' a fresh AR1 year effect at the theta's own `(sigma, rho)`, and
#' `loglik_many` is a Laplace approximation to that year effect's marginal
#' (FD-validated, brute-force cross-checked at T = 2).
#'
#' `cover` is a `tobs_multiarm_fit` with the same two-independent-block shape
#' as `occu_categorical`; `positive = "lognormal"` only for v1, with
#' dispersion held fixed (no SE is reported for it anywhere in the package).
#'
#' `occu_multiscale_cover` is a standard single-block fit (unlike `cover`'s
#' two-block one); the exchangeable unit is the CELL rather than the plot, so
#' pooling and site labeling track cell indices.
#'
#' The community families `ms_occu`, `ms_occu_cover` (`positive =
#' "lognormal"`), `ms_int_occu` (full source-site overlap), `ms_count`
#' (`response = "poisson"`), `jsdm`, `ms_distance` (`key = "halfnorm"`,
#' `mixture = "poisson"`), `ms_dyn_occu` and `ms_abun` (`mixture =
#' "poisson"`) rank a FIXED species set's own realized coefficients (community
#' mean plus each species' deviation, drawn jointly), plus any shared
#' coefficients such as `ms_dyn_occu`'s colonization / extinction and
#' `ms_occu_cover`'s cover dispersion. The replicate is the family's own
#' `simulate()` at the drawn coefficients and the rank statistic sums each
#' species' exact marginal. Calibration of these Laplace-EM posteriors depends
#' on the species count: a small community is measurably non-Gaussian, so they
#' are checked at a community size typical of real data.
#'
#' A family is registered by adding one entry to the internal registry -- its
#' fitting call as a map from arm-formula slots to `tobs()` arguments, a
#' replicate generator where the family's own `simulate()` cannot serve, and
#' optionally a joint statistic; the pooling, the grouping labels, the arm
#' construction and the controls are shared.
#'
#' @return An object of class `tulpa_sbc` from \pkg{tulpa}, carrying the per-quantity
#'   rank ECDFs, the exact simultaneous band, and the uniformity tests. When
#'   `model.only = TRUE`, the callback list instead.
#'
#' @references
#' Talts, S., Betancourt, M., Simpson, D., Vehtari, A. and Gelman, A. (2018).
#' Validating Bayesian inference algorithms with simulation-based calibration.
#' arXiv:1804.06788.
#'
#' Sailynoja, T., Schmitt, M., Buerkner, P.-C. and Vehtari, A. (2026).
#' Posterior SBC: simulation-based calibration checking conditional on data.
#' Statistics and Computing 36:78.
#'
#' @seealso [tulpa::sbc()], [waic()]
#' @examples
#' \donttest{
#' N <- 20L; J <- 3L
#' adj <- matrix(0L, N, N)
#' for (s in seq_len(N)) {
#'   if (s > 1L) adj[s, s - 1L] <- 1L
#'   if (s < N)  adj[s, s + 1L] <- 1L
#' }
#' sim <- simulate_occu_cover(N = N, J = J, positive = "lognormal",
#'                            adj = adj, sigma = 0.7, alpha = 1, seed = 1L)
#' long <- data.frame(site_id = rep(seq_len(N), each = J),
#'                    visit = rep(seq_len(J), times = N),
#'                    y = as.vector(t(sim$y)),
#'                    det_cov1 = sim$visit_data$det_cov1,
#'                    pos_cov1 = sim$visit_data$pos_cov1)
#' od <- tobs_data(long, y = "y", site = "site_id", visit = "visit",
#'                 det.covs = c("det_cov1", "pos_cov1"))
#' y_pos <- sim$y_pos; y_pos[is.na(y_pos)] <- 0
#' # The dispersion goes on the outer grid so it is estimated and scored, and
#' # the field-SD grid is pinned so both stages integrate the same support.
#' ctl <- list(engine = "joint", verbose = FALSE,
#'             sigma.grid   = exp(seq(log(0.15), log(2.0), length.out = 5)),
#'             phi.grid.pos = exp(seq(log(0.20), log(0.90), length.out = 4)))
#' fit <- tobs(~ occ_cov1 + icar(graph = adj),
#'             data = cbind(data.frame(site_id = seq_len(N)), sim$data),
#'             family = occu_cover("lognormal"),
#'             detection = ~ det_cov1,
#'             positive = ~ pos_cov1 + share(spatial()),
#'             y = od$y, y.pos = y_pos, visits = od$det.covs,
#'             method = "nested_laplace", control = ctl)
#' sbc(fit, n.sim = 3L, controls = "narrow", fit.control = ctl)
#' }
#' @export
sbc.tobs_fit <- function(object, n.sim = 100L, n.draws = 1000L, n.ref = 200L,
                         quantities = NULL, controls = character(),
                         bad.factor = 1.25, level = 0.95, seed = 0L,
                         model.only = FALSE, fit.control = list(),
                         control = list(), ...) {
  controls <- if (length(controls))
    match.arg(controls, c("wide", "narrow"), several.ok = TRUE) else character(0)
  model <- .tobs_sbc_build_model(object, as.integer(n.draws), as.integer(n.ref),
                                 controls, bad.factor, fit.control)
  if (isTRUE(model.only)) return(model)
  # tulpa::sbc()'s S3-dispatch parameter is named `object` (the experiment
  # name "posterior"/"prior_predictive"), not `experiment` -- a stale call
  # site name silently fell through to `...`, leaving `object` at its default
  # ("prior_predictive") and erroring there for every family, not just a
  # newly registered one (found while adding double_observer to the
  # registry).
  res <- tulpa::sbc(object = "posterior", model = model,
                    n_sim = as.integer(n.sim), quantities = quantities,
                    level = level, seed = as.integer(seed), control = control)
  res$tobs_family <- attr(model, "family")
  res$fixed_quantities <- attr(model, "fixed")
  res
}
