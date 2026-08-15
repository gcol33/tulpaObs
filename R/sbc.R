# =============================================================================
# Simulation-based calibration for tobs families (gcol33/tulpaObs#207).
#
# tulpa owns the SBC machinery -- the predictive shapes, the within-atom PIT,
# the CRPS closed forms, the exact simultaneous ECDF bands, and both drivers,
# all behind `tulpa::sbc()` (gcol33/tulpa#380). This file owns the other side of
# that contract: the callbacks that turn a fitted `tobs_fit` into the `model =`
# list the POSTERIOR experiment reads, so calibration is measured on the model a
# deliverable actually ships rather than on a toy fixture.
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
# arm construction, the mis-scaled controls, the fixed-quantity guard and the
# seed split are written once and shared, so a second family is data entry
# rather than a second harness.
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

# The observed data set, in the shape `refit()` reads.
.tobs_sbc_data_from_fit <- function(fit) {
  m <- fit$model
  y <- m$y
  if (!is.null(m$valid)) y[!m$valid] <- NA
  list(cells  = m$data,
       y      = y,
       y_pos  = m$y_pos,
       visits = .tobs_sbc_visit_matrices(m),
       graph  = fit$spatial$graph,
       site   = seq_len(m$n_sites))
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
  vis <- lapply(names(obs$visits),
                function(nm) rbind(obs$visits[[nm]], rep$visits[[nm]]))
  names(vis) <- names(obs$visits)
  list(cells  = .tobs_sbc_rbind_cells(obs$cells, rep$cells),
       y      = rbind(obs$y, rep$y),
       y_pos  = if (is.null(obs$y_pos)) NULL else rbind(obs$y_pos, rep$y_pos),
       visits = vis,
       graph  = graph,
       site   = c(obs$site, rep$site + n_o))
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
  ctl <- utils::modifyList(.tobs_sbc_default_control(),
                           as.list(fit.control))
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

# The cover formula, carrying the shared field's copy onto the cover arm when
# the observed fit carried one. Without the copy the cover arm never sees the
# field and `alpha` is pinned at zero, which is a different model.
.tobs_sbc_pos_formula <- function(spec) {
  # `control$alpha.grid` is the other spelling of the same coupling, and
  # `occu_cover()` refuses both at once, so the formula copy is written only
  # when the control does not already carry the scale.
  if (!isTRUE(spec$has_copy) || !is.null(spec$control$alpha.grid)) {
    return(spec$pos)
  }
  lab <- attr(stats::terms(spec$pos), "term.labels")
  f <- stats::as.formula(
    paste("~", paste(c(lab, "copy(spatial())"), collapse = " + ")))
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
    y_pos     = y_pos,
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
  if (length(d$blocks)) {
    b1 <- d$blocks[[1L]]
    M <- cbind(M, sigma = b1$amp_occ, sigma_pos_field = b1$amp_pos,
               alpha = ifelse(b1$amp_occ > 0, b1$amp_pos / b1$amp_occ, 0))
  }
  M <- cbind(M, disp = d$disp)
  attr(M, "blocks") <- blocks
  M
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

.tobs_sbc_loglik_occu_cover <- function(fit, theta, n_quad = 15L) {
  m  <- fit$model
  bl <- attr(theta, "blocks")
  gh <- .gauss_hermite_prob(as.integer(n_quad))
  sigma <- if ("sigma" %in% names(theta)) theta[["sigma"]] else 0
  sigma <- sigma * .tobs_sbc_field_scale(fit$spatial$graph)
  alpha <- if ("alpha" %in% names(theta)) theta[["alpha"]] else 0
  disp  <- if ("disp" %in% names(theta)) theta[["disp"]] else
             (m$cover_pos_disp %||% 1)
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
# one common factor while `alpha`, their ratio, stays clean
# (gcol33/tulpaObs#213).
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

.tobs_sbc_draw_positive <- function(eta, disp, positive) {
  n <- length(eta)
  if (identical(positive, "beta")) {
    mu <- stats::plogis(eta)
    pmin(pmax(stats::rbeta(n, mu * disp, (1 - mu) * disp), 1e-6), 1 - 1e-6)
  } else if (identical(positive, "gaussian")) {
    stats::rnorm(n, eta, disp)
  } else {
    exp(stats::rnorm(n, eta, disp))
  }
}

.tobs_sbc_sim_occu_cover <- function(spec, theta, seed) {
  set.seed(seed)
  m  <- spec$model
  bl <- attr(theta, "blocks")
  n <- m$n_sites; J <- m$max_visits
  sigma <- if ("sigma" %in% names(theta)) theta[["sigma"]] else 0
  alpha <- if ("alpha" %in% names(theta)) theta[["alpha"]] else 0
  disp  <- if ("disp" %in% names(theta)) theta[["disp"]] else
             (m$cover_pos_disp %||% 1)

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
    y_pos[hit] <- .tobs_sbc_draw_positive(e$ep_mat[hit], disp, spec$positive)
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

# A structured term makes a latent quantity that theta does not carry shared
# across sites, which is the one thing this construction denies. Refused with a
# pointer rather than scored on the coefficients alone, which would report a
# calibration the experiment did not measure.
.tobs_sbc_reject_structure <- function(fit) {
  st <- fit$model$structured_terms
  present <- character(0)
  if (is.list(st)) present <- names(Filter(Negate(is.null), st))
  if (!is.null(fit$spatial))   present <- unique(c(present, "spatial"))
  if (!is.null(fit$temporal))  present <- unique(c(present, "temporal"))
  if (!is.null(fit$re_effects)) present <- unique(c(present, "re"))
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
  list(cells  = cells,
       y      = y,
       y_pos  = NULL,
       visits = .tobs_sbc_visit_matrices(m),
       graph  = NULL,
       site   = seq_len(nrow(y)))
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

.tobs_sbc_spec_simple <- function(fit, fit.control, state, det) {
  .tobs_sbc_reject_structure(fit)
  .tobs_sbc_reject_visit_design(fit)
  m <- fit$model
  if (is.null(m$formulas[[state]])) {
    stop("SBC rebuilds the fitting call from the model's own arm formulas; ",
         "slot ", sQuote(state), " is absent.", call. = FALSE)
  }
  list(model   = m,
       fit_obs = fit,
       family  = attr(fit, "tobs_family"),
       method  = fit$method,
       control = utils::modifyList(.tobs_sbc_default_control(),
                                   as.list(fit.control)),
       state   = .tobs_sbc_recombine(m$formulas[[state]], NULL),
       det     = if (is.null(det)) NULL else
                   .tobs_sbc_recombine(m$formulas[[det]], NULL))
}

.tobs_sbc_refit_simple <- function(spec, data) {
  args <- list(formula = spec$state,
               data    = data$cells,
               family  = spec$family,
               y       = if (isTRUE(spec$y_vector)) as.numeric(data$y) else data$y,
               method  = spec$method,
               control = spec$control)
  if (!is.null(spec$det)) args$detection <- spec$det
  if (length(data$visits)) args$visits <- data$visits
  suppressWarnings(do.call(tobs, c(args, spec$extra)))
}

.tobs_sbc_sim_simple <- function(spec, theta, seed) {
  set.seed(seed)
  f <- spec$fit_obs
  D <- matrix(theta, nrow = 1L)
  colnames(D) <- names(theta)
  f$draws <- D
  obs <- spec$data_obs
  list(cells  = obs$cells,
       y      = .tobs_sbc_clean_y(spec$replicate(f)),
       y_pos  = NULL,
       visits = obs$visits,
       graph  = NULL,
       site   = obs$site)
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

# The Poisson / negative-binomial count GLMM has no latent state and so no
# `simulate()` handler to borrow; its replicate is the response drawn at the
# fitted mean, which is that handler's whole content for this family.
.tobs_sbc_replicate_count <- function(f) {
  m <- f$model
  p <- m$process_info[[1L]]$p
  mu <- exp(as.vector(m$X_occ %*% f$draws[1L, seq_len(p)]))
  switch(m$response,
         poisson = stats::rpois(length(mu), mu),
         negbin  = stats::rnbinom(length(mu), mu = mu,
                                  size = exp(f$draws[1L, p + 1L])),
         stop("SBC on count() is registered for the poisson and negbin ",
              "responses; this fit is ", sQuote(m$response), ".",
              call. = FALSE))
}

# One registry row. `state` / `det` name the arm formulas in the fitted model,
# `resp` its response slot, `extra` any further arguments the family's front
# door takes, `replicate` a generator when the family has no `simulate()`
# handler.
.tobs_sbc_simple_entry <- function(state, det = NULL, resp = "y",
                                   extra = NULL, replicate = NULL,
                                   y_vector = FALSE) {
  list(
    spec = function(fit, fit.control) {
      sp <- .tobs_sbc_spec_simple(fit, fit.control, state = state, det = det)
      sp$extra <- if (is.function(extra)) extra(fit$model) else extra
      sp$replicate <- replicate %||%
        function(f) stats::simulate(f, nsim = 1L)
      sp$y_vector <- y_vector
      sp
    },
    data        = function(fit) .tobs_sbc_data_simple(fit, resp = resp),
    draws       = .tobs_sbc_draws_fit,
    simulate    = .tobs_sbc_sim_simple,
    refit       = .tobs_sbc_refit_simple,
    loglik_many = .tobs_sbc_loglik_many_simple)
}


# ---------------------------------------------------------------------------
# 6b. dyn_occu() (gcol33/tulpaObs#220, multi-season group): pooling on the
# SITE axis, leaving the season axis alone. `model$y` is a 3D
# [n_sites x max_visits x n_seasons] array, so the shared simple-family
# `data`/`pool` (2D `rbind`) cannot be reused; every other piece (draws off
# `fit$draws`, the `dynamic` pointwise log-likelihood already dispatched by
# `.tobs_pointwise_loglik`) is. Constant-rate only (season-varying detection /
# colonization / extinction is a follow-up -- `.tobs_sbc_reject_visit_design`
# alone does not see those designs, so they are refused explicitly).
# ---------------------------------------------------------------------------

.tobs_sbc_reject_season_varying <- function(fit) {
  m <- fit$model
  sv <- c(det = isTRUE(m$det_season_varying),
         col = isTRUE(m$col_season_varying),
         ext = isTRUE(m$ext_season_varying))
  if (!any(sv)) return(invisible(NULL))
  stop("SBC on dyn_occu() is registered for constant-rate fits; this fit ",
       "carries a season-varying (", paste(names(sv)[sv], collapse = ", "),
       ") design, which the replicate generator does not yet reproduce.",
       call. = FALSE)
}

.tobs_sbc_spec_dyn_occu <- function(fit, fit.control) {
  .tobs_sbc_reject_structure(fit)
  .tobs_sbc_reject_visit_design(fit)
  .tobs_sbc_reject_season_varying(fit)
  m <- fit$model
  list(model   = m,
       fit_obs = fit,
       family  = attr(fit, "tobs_family"),
       method  = fit$method,
       control = utils::modifyList(.tobs_sbc_default_control(),
                                   as.list(fit.control)),
       occ = .tobs_sbc_recombine(m$formulas$occ, NULL),
       det = .tobs_sbc_recombine(m$formulas$det, NULL),
       col = .tobs_sbc_recombine(m$formulas$col, NULL),
       ext = .tobs_sbc_recombine(m$formulas$ext, NULL))
}

# Shared by every multi-season family (dyn_occu, dyn_abun): `model$y` is a 3D
# [n_sites x max_visits x n_seasons] array under both, and neither `data` nor
# `pool` reads anything family-specific off it.
.tobs_sbc_data_3d_season <- function(fit) {
  m <- fit$model
  list(cells  = m$data,
       y      = m$y,
       y_pos  = NULL,
       visits = NULL,
       graph  = NULL,
       site   = seq_len(m$n_sites))
}

# Concatenate two 3D [n_sites x max_visits x n_seasons] response arrays on the
# SITE axis; the season axis (and its transition structure) is shared design,
# not something a replicate draws fresh.
.tobs_sbc_pool_3d_season <- function(obs, rep) {
  n_o <- dim(obs$y)[1L]; n_r <- dim(rep$y)[1L]
  y <- array(NA_real_, dim = c(n_o + n_r, dim(obs$y)[2L], dim(obs$y)[3L]))
  y[seq_len(n_o), , ] <- obs$y
  y[n_o + seq_len(n_r), , ] <- rep$y
  list(cells = .tobs_sbc_rbind_cells(obs$cells, rep$cells), y = y, y_pos = NULL, visits = NULL,
       graph = NULL, site = c(obs$site, rep$site + n_o))
}

.tobs_sbc_refit_dyn_occu <- function(spec, data) {
  suppressWarnings(do.call(tobs, list(
    formula = spec$occ, data = data$cells, family = spec$family,
    detection = spec$det, colonization = spec$col, extinction = spec$ext,
    y = data$y, method = spec$method, control = spec$control)))
}

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

  list(cells = m$data, y = y, y_pos = NULL, visits = NULL, graph = NULL,
       site = seq_len(n))
}


# ---------------------------------------------------------------------------
# 6b2. dyn_abun() (gcol33/tulpaObs#220, multi-season group): the abundance
# sibling of dyn_occu -- same 3D [n_sites x max_visits x n_seasons] response,
# same site-axis pooling (`.tobs_sbc_data_3d_season` / `.tobs_sbc_pool_3d_season`
# above, shared verbatim). UNLIKE dyn_occu it already has a working
# `simulate()` handler (the Dail-Madsen forward is not a two-state HMM to hand
# -write), so the replicate is the family's own kernel via the shared
# `.tobs_sbc_sim_simple`, exactly the `occu()`-style pattern -- no bespoke
# forward simulator needed. Constant-rate only (season-varying omega/gamma is
# a follow-up, same reasoning as dyn_occu's colonization/extinction).
# ---------------------------------------------------------------------------

.tobs_sbc_reject_dyn_abun_season_varying <- function(fit) {
  m <- fit$model
  sv <- c(omega = isTRUE(m$omega_season_varying),
         gamma = isTRUE(m$gamma_season_varying))
  if (!any(sv)) return(invisible(NULL))
  stop("SBC on dyn_abun() is registered for constant-rate fits; this fit ",
       "carries a season-varying (", paste(names(sv)[sv], collapse = ", "),
       ") design, which the pooled refit does not yet reproduce.",
       call. = FALSE)
}

.tobs_sbc_spec_dyn_abun <- function(fit, fit.control) {
  .tobs_sbc_reject_structure(fit)
  .tobs_sbc_reject_visit_design(fit)
  .tobs_sbc_reject_dyn_abun_season_varying(fit)
  m <- fit$model
  list(model     = m,
       fit_obs   = fit,
       family    = attr(fit, "tobs_family"),
       method    = fit$method,
       control   = utils::modifyList(.tobs_sbc_default_control(),
                                     as.list(fit.control)),
       lambda    = .tobs_sbc_recombine(m$formulas$lambda, NULL),
       p         = .tobs_sbc_recombine(m$formulas$p,      NULL),
       omega     = .tobs_sbc_recombine(m$formulas$omega,  NULL),
       gamma     = .tobs_sbc_recombine(m$formulas$gamma,  NULL),
       replicate = function(f) stats::simulate(f, nsim = 1L))
}

.tobs_sbc_refit_dyn_abun <- function(spec, data) {
  suppressWarnings(do.call(tobs, list(
    formula = spec$lambda, data = data$cells, family = spec$family,
    detection = spec$p, omega = spec$omega, gamma = spec$gamma,
    y = data$y, method = spec$method, control = spec$control)))
}


# ---------------------------------------------------------------------------
# 6c. int_occu() (gcol33/tulpaObs#220, multi-source group): pooling on the
# SITE axis, leaving the per-source axis alone. `model$y_sources` is a list
# of one detection matrix per source (`model$site_maps` the row -> site map),
# so this is a second response shape the shared simple-family `data`/`pool`
# cannot read. Full-overlap fits only (every source observes every site) --
# `.tobs_ploglik_integrated`, which `loglik_many` reuses unchanged, is itself
# only wired for that case.
# ---------------------------------------------------------------------------

.tobs_sbc_int_occu_src_names <- function(m) {
  vapply(seq_len(m$n_sources), function(s) m$process_info[[1L + s]]$name,
        character(1))
}

.tobs_sbc_spec_int_occu <- function(fit, fit.control) {
  .tobs_sbc_reject_structure(fit)
  .tobs_sbc_reject_visit_design(fit)
  m <- fit$model
  full <- vapply(seq_len(m$n_sources), function(s)
    identical(as.integer(m$site_maps[[s]]), seq_len(m$n_sites) - 1L), logical(1))
  if (!all(full)) {
    stop("SBC on int_occu() is registered for full-overlap fits (every ",
         "source observes every site); this fit has partial overlap.",
         call. = FALSE)
  }
  src_nm <- .tobs_sbc_int_occu_src_names(m)
  list(model     = m,
       fit_obs   = fit,
       family    = attr(fit, "tobs_family"),
       method    = fit$method,
       control   = utils::modifyList(.tobs_sbc_default_control(),
                                     as.list(fit.control)),
       occ       = .tobs_sbc_recombine(m$formulas$occ, NULL),
       det_list  = lapply(src_nm, function(nm) .tobs_sbc_recombine(m$formulas[[nm]], NULL)),
       src_names = src_nm)
}

.tobs_sbc_data_int_occu <- function(fit) {
  m <- fit$model
  list(cells  = m$data,
       y      = stats::setNames(m$y_sources, .tobs_sbc_int_occu_src_names(m)),
       y_pos  = NULL,
       visits = NULL,
       graph  = NULL,
       site   = seq_len(m$n_sites))
}

# Concatenate each source's own matrix on the SITE axis; sources stay
# separate lists (they need not share a visit count).
# Generic over any family whose `y` is a NAMED LIST of matrices sharing one
# site axis -- int_occu()'s per-source detection matrices, gdistremoval()'s
# yDist/yRem -- so both registry entries below share this one `pool`.
.tobs_sbc_pool_named_matrices <- function(obs, rep) {
  n_o <- length(obs$site); n_r <- length(rep$site)
  y <- Map(function(a, b) rbind(a, b), obs$y, rep$y[names(obs$y)])
  list(cells = .tobs_sbc_rbind_cells(obs$cells, rep$cells), y = y, y_pos = NULL, visits = NULL,
       graph = NULL, site = c(obs$site, rep$site + n_o))
}

.tobs_sbc_refit_int_occu <- function(spec, data) {
  suppressWarnings(do.call(tobs, list(
    formula = spec$occ, data = data$cells, family = spec$family,
    detection = spec$det_list, y = data$y[spec$src_names],
    method = spec$method, control = spec$control)))
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
  list(cells = m$data, y = y, y_pos = NULL, visits = NULL, graph = NULL,
       site = seq_len(m$n_sites))
}


# ---------------------------------------------------------------------------
# 6d. gdistremoval() (gcol33/tulpaObs#220, multi-response group): a single
# time-step family whose response is TWO matrices (yDist band counts, yRem
# period counts) that must stay row-consistent, pooled together on the site
# axis via the same `.tobs_sbc_pool_named_matrices` int_occu()'s per-source
# matrices use. Already has a `simulate()` handler returning
# `list(yDist=, yRem=)`, reused here as the replicate generator the way the
# occu()-family entries reuse their own `simulate()` handler.
# ---------------------------------------------------------------------------

.tobs_sbc_spec_gdistremoval <- function(fit, fit.control) {
  .tobs_sbc_reject_structure(fit)
  .tobs_sbc_reject_visit_design(fit)
  m <- fit$model
  list(model   = m,
       fit_obs = fit,
       family  = attr(fit, "tobs_family"),
       method  = fit$method,
       control = utils::modifyList(.tobs_sbc_default_control(),
                                   as.list(fit.control)),
       lambda    = .tobs_sbc_recombine(m$formulas$lambda, NULL),
       detection = .tobs_sbc_recombine(m$formulas$sigma,  NULL),
       removal   = .tobs_sbc_recombine(m$formulas$r,      NULL))
}

.tobs_sbc_data_gdistremoval <- function(fit) {
  m <- fit$model
  list(cells  = m$data,
       y      = list(yDist = m$y, yRem = m$y_rem),
       y_pos  = NULL,
       visits = NULL,
       graph  = NULL,
       site   = seq_len(m$n_sites))
}

.tobs_sbc_refit_gdistremoval <- function(spec, data) {
  suppressWarnings(do.call(tobs, list(
    formula = spec$lambda, data = data$cells, family = spec$family,
    detection = spec$detection, removal = spec$removal,
    y = data$y$yDist, y_rem = data$y$yRem,
    method = spec$method, control = spec$control)))
}

.tobs_sbc_sim_gdistremoval <- function(spec, theta, seed) {
  set.seed(seed)
  f <- spec$fit_obs
  D <- matrix(theta, nrow = 1L)
  colnames(D) <- names(theta)
  f$draws <- D
  rep <- stats::simulate(f, nsim = 1L)
  list(cells = spec$model$data, y = list(yDist = rep$yDist, yRem = rep$yRem),
       y_pos = NULL, visits = NULL, graph = NULL,
       site = seq_len(spec$model$n_sites))
}


# ---------------------------------------------------------------------------
# 6e. occu_categorical() (gcol33/tulpaObs#220, multiarm-S3 group): presence +
# nominal-class hurdle. A `tobs_multiarm_fit` with NO `fit$model` / `fit$draws`
# -- the two arms are independent Laplace-Gaussian blocks (`beta_occ`/
# `vcov_occ`, `beta_class`/`vcov_class`), fit as two separate GLMs, not a
# joint MVN pseudo-draw matrix off a shared `tulpa_laplace()` postprocessor.
# `draws()` therefore samples the two blocks independently (they ARE
# independent by construction) and stacks them into one theta row;
# `loglik_many()` scores the two-arm likelihood directly off the encoding,
# since there is no `.tobs_pointwise_loglik` dispatch for this family to
# reuse. Non-spatial by construction (no spatial/NUTS route exists yet) and
# takes no `detection` formula, so neither `.tobs_sbc_reject_structure` nor
# `.tobs_sbc_reject_visit_design` apply here.
# ---------------------------------------------------------------------------

.tobs_sbc_mvn_draws <- function(mean, V, n) {
  p <- length(mean)
  if (is.null(V) || anyNA(V)) {
    stop("SBC draws the posterior from this fit's own covariance, which is ",
         "missing or non-finite here (a non-converged or rank-deficient fit).",
         call. = FALSE)
  }
  Vs <- (V + t(V)) / 2
  # A refit's per-species Cinv/vcov can occasionally be numerically singular
  # OR genuinely indefinite, not merely ill-conditioned (a community family
  # fit on a data-sparse simulated replicate -- observed on ms_distance's own
  # acceptance run at n.sim=100, "the leading minor of order 3 is not
  # positive" survived even a 1e-4-relative ridge retry), which a bare
  # chol() has no recourse for. Eigen-clip to the nearest valid covariance
  # instead of guessing at a ridge scale: a small floor on the eigenvalues
  # leaves a well-conditioned Vs untouched (floor << its smallest eigenvalue)
  # and guarantees chol() succeeds regardless of how indefinite the input is.
  L <- tryCatch(chol(Vs), error = function(e) {
    ee <- eigen(Vs, symmetric = TRUE)
    floor_val <- max(1e-8, 1e-6 * max(ee$values))
    vals <- pmax(ee$values, floor_val)
    Vpd <- ee$vectors %*% (vals * t(ee$vectors))
    chol((Vpd + t(Vpd)) / 2)
  })
  Z <- matrix(stats::rnorm(n * p), n, p)
  sweep(Z %*% L, 2L, mean, `+`)
}

# Shared by every family whose response is a plain length-N vector (one
# observation per unit, no visit/season axis) and whose fit stores its
# encoding at `fit$encoding` with `data`/`y`/`N` fields -- occu_categorical,
# cover.
.tobs_sbc_data_vector <- function(fit) {
  enc <- fit$encoding
  list(cells = enc$data, y = enc$y, y_pos = NULL, visits = NULL,
       graph = NULL, site = seq_len(enc$N))
}

.tobs_sbc_pool_vector <- function(obs, rep) {
  list(cells = .tobs_sbc_rbind_cells(obs$cells, rep$cells), y = c(obs$y, rep$y),
       y_pos = NULL, visits = NULL, graph = NULL,
       site = c(obs$site, rep$site + length(obs$site)))
}

.tobs_sbc_spec_occu_categorical <- function(fit, fit.control) {
  list(formula = fit$encoding$formula, class_labels = fit$class_labels,
       control = utils::modifyList(list(), as.list(fit.control)),
       fit_obs = fit)
}

# Columns named `occ_<coef>` / `class_<label>_<coef>`, matching `beta_class`'s
# column-major `vec()` order so a theta row unpacks back into the same matrix
# shape it was drawn from.
.tobs_sbc_draws_occu_categorical <- function(fit, n) {
  occ_names   <- names(fit$beta_occ)
  class_names <- as.vector(outer(rownames(fit$beta_class),
                                 colnames(fit$beta_class),
                                 function(r, cl) paste0("class_", cl, "_", r)))
  D_occ <- .tobs_sbc_mvn_draws(fit$beta_occ, fit$vcov_occ, n)
  D_cls <- .tobs_sbc_mvn_draws(as.numeric(fit$beta_class), fit$vcov_class, n)
  M <- cbind(D_occ, D_cls)
  colnames(M) <- c(paste0("occ_", occ_names), class_names)
  M
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
  X <- stats::model.matrix(spec$formula, obs$cells)
  psi <- stats::plogis(as.numeric(X %*% th$beta_occ))
  present <- stats::rbinom(nrow(X), 1L, psi)
  eta <- X %*% th$Beta
  E <- exp(pmin(eta, 700)); denom <- 1 + rowSums(E)
  P <- cbind(E / denom, 1 / denom)
  K <- f$K
  cls <- apply(P, 1L, function(pr) sample.int(K, 1L, prob = pr))
  y <- ifelse(present == 1L, cls, 0L)
  list(cells = obs$cells, y = as.integer(y), y_pos = NULL, visits = NULL,
       graph = NULL, site = obs$site)
}

.tobs_sbc_refit_occu_categorical <- function(spec, data) {
  suppressWarnings(do.call(tobs, list(
    formula = spec$formula, data = data$cells,
    family = occu_categorical(classes = spec$class_labels),
    y = data$y, method = "laplace", control = spec$control)))
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


# ---------------------------------------------------------------------------
# 6f. distsamp_open() (gcol33/tulpaObs#220, multi-season group): the same 3D
# [n_sites x n_bins x n_seasons] response shape and site-axis pooling as
# dyn_abun (section 6b2) -- the "season" axis here is the open-population
# primary period rather than a revisit season, but the shape and the pooling
# rule (stack on sites, leave the season axis alone) are identical, so
# `.tobs_sbc_data_3d_season`/`.tobs_sbc_pool_3d_season` are reused unchanged.
# `fit$means`/`fit$draws`/`fit$model` are the standard shape too
# (`.tobs_bfgs_marginal_fit()`), and `.tobs_pointwise_loglik` already
# dispatches for this model_type, so `draws`/`loglik_many` are also the
# shared generic ones -- only `spec`/`refit` are custom, to route the
# family's four arm formulas (`lambda`, `sigma`, `omega`, `gamma`) and rebuild
# `distsamp_open(cutpoints=, transect=)`. Constant-dynamics, Poisson only for
# v1 (the alternative dynamics / negbin / zero-inflated layers are follow-ups,
# same reasoning as dyn_abun's season-varying-rate exclusion).
# ---------------------------------------------------------------------------

.tobs_sbc_reject_distsamp_open_scope <- function(fit) {
  m <- fit$model
  dyn <- m$dynamics %||% "constant"
  mix <- m$mixture %||% "poisson"
  if (identical(dyn, "constant") && identical(mix, "poisson")) return(invisible(NULL))
  stop("SBC on distsamp_open() is registered for constant-dynamics, Poisson ",
       "fits; this fit carries dynamics = ", sQuote(dyn), ", mixture = ",
       sQuote(mix), ", which the pooled refit does not yet reproduce.",
       call. = FALSE)
}

.tobs_sbc_spec_distsamp_open <- function(fit, fit.control) {
  .tobs_sbc_reject_structure(fit)
  .tobs_sbc_reject_visit_design(fit)
  .tobs_sbc_reject_distsamp_open_scope(fit)
  m <- fit$model
  list(model     = m,
       fit_obs   = fit,
       family    = attr(fit, "tobs_family"),
       method    = fit$method,
       control   = utils::modifyList(.tobs_sbc_default_control(),
                                     as.list(fit.control)),
       lambda    = .tobs_sbc_recombine(m$formulas$lambda, NULL),
       sigma     = .tobs_sbc_recombine(m$formulas$sigma,  NULL),
       omega     = .tobs_sbc_recombine(m$formulas$omega,  NULL),
       gamma     = .tobs_sbc_recombine(m$formulas$gamma,  NULL),
       cutpoints = m$cutpoints, transect = m$transect,
       replicate = function(f) stats::simulate(f, nsim = 1L))
}

.tobs_sbc_refit_distsamp_open <- function(spec, data) {
  suppressWarnings(do.call(tobs, list(
    formula = spec$lambda, data = data$cells,
    family = distsamp_open(transect = spec$transect, cutpoints = spec$cutpoints,
                           mixture = "poisson", dynamics = "constant"),
    detection = spec$sigma, omega = spec$omega, gamma = spec$gamma,
    y = data$y, method = spec$method, control = spec$control)))
}


# ---------------------------------------------------------------------------
# 6g. occu_multi() (gcol33/tulpaObs#220, multi-response group): `model$y` is a
# list of S per-species detection matrices sharing one site axis -- the same
# shape as int_occu()'s per-source list, so `.tobs_sbc_pool_named_matrices` is
# reused unchanged for pooling (species names attached explicitly in `data()`,
# since the binder does not always name `model$y`). `fit$means`/`fit$draws`
# are the standard `.tobs_bfgs_marginal_fit()` shape and
# `.tobs_pointwise_loglik` already dispatches for this model_type, so
# `draws`/`loglik_many` are the shared generic ones. `simulate` is custom
# (unlike int_occu, species are NOT independent given a shared z -- the joint
# state is one draw from the log-linear multi-species model, then each
# species is independently detected given its own z_k), wrapping the
# family's own `.tobs_simulate_occu_multi()` handler exactly as
# `.tobs_sbc_sim_gdistremoval()` wraps `simulate_gdistremoval()`'s handler --
# `.tobs_sbc_sim_simple`'s cleaning step assumes a single matrix/array
# response, not a list, so it cannot be reused here either.
# ---------------------------------------------------------------------------

.tobs_sbc_data_occu_multi <- function(fit) {
  m <- fit$model
  list(cells = m$data, y = stats::setNames(m$y, m$species), y_pos = NULL,
       visits = NULL, graph = NULL, site = seq_len(m$n_sites))
}

.tobs_sbc_spec_occu_multi <- function(fit, fit.control) {
  .tobs_sbc_reject_structure(fit)
  .tobs_sbc_reject_visit_design(fit)
  m <- fit$model
  list(model     = m,
       fit_obs   = fit,
       family    = attr(fit, "tobs_family"),
       method    = fit$method,
       control   = utils::modifyList(.tobs_sbc_default_control(),
                                     as.list(fit.control)),
       state     = .tobs_sbc_recombine(m$formulas$state, NULL),
       det       = .tobs_sbc_recombine(m$formulas$det,   NULL),
       species   = m$species)
}

.tobs_sbc_refit_occu_multi <- function(spec, data) {
  suppressWarnings(do.call(tobs, list(
    formula = spec$state, data = data$cells, family = spec$family,
    detection = spec$det, y = data$y[spec$species], species = spec$species,
    method = spec$method, control = spec$control)))
}

.tobs_sbc_sim_occu_multi <- function(spec, theta, seed) {
  set.seed(seed)
  f <- spec$fit_obs
  D <- matrix(theta, nrow = 1L)
  colnames(D) <- names(theta)
  f$draws <- D
  rep <- stats::simulate(f, nsim = 1L)
  obs <- spec$data_obs
  list(cells = obs$cells, y = rep, y_pos = NULL, visits = NULL, graph = NULL,
       site = obs$site)
}


# ---------------------------------------------------------------------------
# 6h. dyn_int_occu() (gcol33/tulpaObs#220): the product of the multi-season
# and multi-source shapes -- `model$y` is a NAMED list of S per-source 3D
# [n_sites x max_visits_s x T] arrays sharing the site axis. Neither the
# named-matrices pool (2D `rbind`, int_occu()/gdistremoval()) nor the
# 3D-season pool (a single array, dyn_occu()/dyn_abun()/distsamp_open()) fits
# alone, so a new generic `.tobs_sbc_pool_named_3d` composes both: site-axis
# stack WITHIN each source's own 3D array. `simulate()` is custom, wrapping
# the family's own `.tobs_simulate_dyn_int_occu()` handler -- which already
# reproduces each source's own NA/missingness pattern, so partial season
# overlap needs no special-casing here, it falls out of the family's own
# simulate(). `fit$means`/`fit$draws` are the standard shape and
# `.tobs_pointwise_loglik` already dispatches for this model_type, so
# `draws`/`loglik_many` are the shared generic ones. The binder itself only
# accepts full SITE overlap across sources (every source's array shares
# `n_sites`), so no extra reject is needed there (unlike int_occu()); v1 has
# no season-varying-rate option to reject either (the family's only mode).
# ---------------------------------------------------------------------------

.tobs_sbc_pool_named_3d <- function(obs, rep) {
  n_o <- length(obs$site); n_r <- length(rep$site)
  y <- Map(function(a, b) {
    d <- dim(a)
    z <- array(NA_real_, dim = c(n_o + n_r, d[2L], d[3L]))
    z[seq_len(n_o), , ] <- a
    z[n_o + seq_len(n_r), , ] <- b
    z
  }, obs$y, rep$y[names(obs$y)])
  list(cells = .tobs_sbc_rbind_cells(obs$cells, rep$cells), y = y, y_pos = NULL, visits = NULL,
       graph = NULL, site = c(obs$site, rep$site + n_o))
}

.tobs_sbc_data_dyn_int_occu <- function(fit) {
  m <- fit$model
  list(cells = m$data, y = stats::setNames(m$y, m$sources), y_pos = NULL,
       visits = NULL, graph = NULL, site = seq_len(m$n_sites))
}

.tobs_sbc_spec_dyn_int_occu <- function(fit, fit.control) {
  .tobs_sbc_reject_structure(fit)
  .tobs_sbc_reject_visit_design(fit)
  m <- fit$model
  list(model     = m,
       fit_obs   = fit,
       family    = attr(fit, "tobs_family"),
       method    = fit$method,
       control   = utils::modifyList(.tobs_sbc_default_control(),
                                     as.list(fit.control)),
       state     = .tobs_sbc_recombine(m$formulas$psi,   NULL),
       col       = .tobs_sbc_recombine(m$formulas$gamma, NULL),
       ext       = .tobs_sbc_recombine(m$formulas$eps,   NULL),
       det       = .tobs_sbc_recombine(m$formulas$p,     NULL),
       sources   = m$sources)
}

.tobs_sbc_refit_dyn_int_occu <- function(spec, data) {
  suppressWarnings(do.call(tobs, list(
    formula = spec$state, data = data$cells, family = spec$family,
    colonization = spec$col, extinction = spec$ext, detection = spec$det,
    y = data$y[spec$sources], sources = spec$sources,
    method = spec$method, control = spec$control)))
}

.tobs_sbc_sim_dyn_int_occu <- function(spec, theta, seed) {
  set.seed(seed)
  f <- spec$fit_obs
  D <- matrix(theta, nrow = 1L)
  colnames(D) <- names(theta)
  f$draws <- D
  rep <- stats::simulate(f, nsim = 1L)
  obs <- spec$data_obs
  list(cells = obs$cells, y = rep, y_pos = NULL, visits = NULL, graph = NULL,
       site = obs$site)
}


# ---------------------------------------------------------------------------
# 6i. t_occu() (gcol33/tulpaObs#220): a Polya-Gamma Gibbs family, but
# `.tobs_pg_finalize_fit()` already reports the pooled cross-chain draws as
# `fit$draws` in the standard shape (real posterior samples, not a Gaussian
# pseudo-draw matrix) -- so `.tobs_sbc_draws_fit` reads it unchanged, no
# Gibbs-aware draws() needed. `model$y` is `[n_sites x n_seasons x
# max_visits]` (season BEFORE visits -- the one family with this axis order),
# but the 3D-season pool only ever stacks axis 1 and leaves the other two
# alone, so `.tobs_sbc_data_3d_season`/`.tobs_sbc_pool_3d_season` are still
# exactly right. What genuinely IS custom: the family has no `simulate()`
# handler at all (`simulate_t_occu()` draws a FRESH truth, not a replicate at
# a given theta) and no `.tobs_pointwise_loglik` dispatch, so `simulate` is
# hand-written (drawing a fresh AR1 year-effect sequence at the theta's own
# (sigma, rho) -- the same generative model `simulate_t_occu()` uses,
# parameterized by theta instead of drawing its own truth). v1 is also
# fully-observed only (no NA-visit masking, matching `simulate_t_occu()`'s
# own output).
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

.tobs_sbc_spec_t_occu <- function(fit, fit.control) {
  .tobs_sbc_reject_structure(fit)
  .tobs_sbc_reject_visit_design(fit)
  m <- fit$model
  list(model   = m,
       fit_obs = fit,
       family  = attr(fit, "tobs_family"),
       method  = fit$method,
       control = utils::modifyList(.tobs_sbc_default_control(),
                                   as.list(fit.control)),
       occ     = .tobs_sbc_recombine(m$formulas$occ, NULL),
       det     = .tobs_sbc_recombine(m$formulas$det, NULL))
}

.tobs_sbc_refit_t_occu <- function(spec, data) {
  suppressWarnings(do.call(tobs, list(
    formula = spec$occ, data = data$cells, family = spec$family,
    detection = spec$det, y = data$y, method = spec$method,
    control = spec$control)))
}

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
  list(cells = obs$cells, y = y, y_pos = NULL, visits = NULL, graph = NULL,
       site = obs$site)
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
# 6j. ms_occu() (gcol33/tulpaObs#220, community group): the design decision
# the issue asks to state rather than default -- rank a FIXED species set's
# own coefficients (mu + b_s per species), not the community means with fresh
# species drawn per replicate. theta is the per-species REALIZED coefficient
# vector for every species (S x P, column names `<species>_psi_<coef>` /
# `<species>_p_<coef>`), not the P-length community mean `fit$means` alone.
#
# `.tobs_community_em()` (R/community_em.R) computes `Cinv` (per-species
# posterior covariance Cov(b_s|y), Louis 1982) and `Bf` (the mu-b_s
# cross-Hessian block from the same Newton solve), both exposed as
# `fit$ms_community$Cinv`/`Bf` (R/ms_occu.R). `draws()` samples the community
# mean once per row (`mu ~ N(means, vcov)`) then draws each species'
# deviation CONDITIONAL on that mu draw via `.tobs_sbc_community_b_draws`
# (`b_s | mu ~ N(blup_s - Cinv_s t(Bf_s)(mu_draw-means), Cinv_s)`), theta_s =
# mu + b_s -- NOT independently, which was #226 (drawing mu and b_s as
# independent ignores their posterior cross-covariance and biases Var(theta_s)
# on whichever species trades off most against the community mean). See the
# registry entry below for the species-count scope this fix needed to actually
# register cleanly (S=20, not the S=5 fixture the bug was originally found on).
#
# `simulate()` reuses the family's OWN `.tobs_simulate_ms_occu()` handler
# (the validated `cpp_simulate_ms_occu` kernel, which already respects the
# observed model's own visit-validity pattern) rather than a hand-written
# generator -- injecting the candidate theta by overwriting
# `f$ms_community$coef_psi`/`coef_p` (the per-species coefficient matrices
# `.tobs_simulate_ms_occu` actually reads), not `f$draws` (which that handler
# does not consult). `loglik_many()` sums each species' own two-state
# marginal (`.ms_int_occu_sp_ll`, the exact kernel `.tobs_fit_ms_occu()`
# already optimizes) at that species' theta slice -- the natural per-species
# analogue of every other family's marginal rank arm. Pooling/data reuse the
# generic 3D-season helpers unchanged: `model$y` is `[site x visit x
# species]`, and stacking axis 1 (sites) while leaving axes 2-3 alone is
# exactly what they already do, regardless of what those axes mean.
# ---------------------------------------------------------------------------

.tobs_sbc_spec_ms_occu <- function(fit, fit.control) {
  .tobs_sbc_reject_structure(fit)
  .tobs_sbc_reject_visit_design(fit)
  m <- fit$model
  list(model   = m,
       fit_obs = fit,
       family  = attr(fit, "tobs_family"),
       method  = fit$method,
       control = utils::modifyList(.tobs_sbc_default_control(),
                                   as.list(fit.control)),
       occ     = .tobs_sbc_recombine(m$formulas$occ, NULL),
       det     = .tobs_sbc_recombine(m$formulas$det, NULL),
       species = m$species_names)
}

.tobs_sbc_refit_ms_occu <- function(spec, data) {
  suppressWarnings(do.call(tobs, list(
    formula = spec$occ, data = data$cells, family = spec$family,
    detection = spec$det, y = data$y, species = spec$species,
    method = spec$method, control = spec$control)))
}

# theta column names: `<species>_psi_<coef>` / `<species>_p_<coef>`, one
# block of length P per species, species-major -- the layout `draws()`,
# `simulate()` and `loglik_many()` all share.
.tobs_sbc_ms_occu_names <- function(m) {
  occ_names <- m$process_info[[1L]]$coef_names
  det_names <- m$process_info[[2L]]$coef_names
  par <- c(paste0("psi_", occ_names), paste0("p_", det_names))
  list(occ = occ_names, det = det_names, par = par,
       cols = as.vector(outer(par, m$species_names,
                              function(p, sp) paste0(sp, "_", p))))
}

# gcol33/tulpaObs#226: mu and b_s are NOT independent in the posterior --
# conditional on a draw of mu, b_s's mean shifts by
# -Cinv_s %*% t(Bf_s) %*% (mu_draw - mu_hat) (Bf_s the mu-b_s cross-Hessian
# block from the community EM's own Newton solve, R/community_em.R), while
# its covariance is UNCHANGED at Cinv_s. Derived by block-inverting the joint
# (mu, b_s) arrowhead precision [Sf + Bf_s Cinv_s Bf_s', Bf_s; Bf_s',
# Cinv_s^-1] and validated to machine precision against a direct construction
# of that matrix on random SPD inputs (dev_notes probe, 2026-08-12). Drawing
# mu and each b_s independently -- what every affected family did before this
# -- is exactly the #226 bug: it under-covers on whichever species' deviation
# most strongly trades off against the community mean.
.tobs_sbc_community_b_draws <- function(mu_draws, mu_hat, b_hat_s, Bf_s, Cinv_s, n) {
  eps <- .tobs_sbc_mvn_draws(rep(0, length(b_hat_s)), Cinv_s, n)
  delta <- sweep(mu_draws, 2L, mu_hat, "-")
  shift <- (delta %*% Bf_s) %*% Cinv_s
  sweep(eps - shift, 2L, b_hat_s, "+")
}

.tobs_sbc_draws_ms_occu <- function(fit, n) {
  m <- fit$model
  nm <- .tobs_sbc_ms_occu_names(m)
  S <- m$n_species; P <- length(nm$par)
  cm <- fit$ms_community
  mu_draws <- .tobs_sbc_mvn_draws(fit$means, fit$vcov, n)   # n x P
  M <- matrix(NA_real_, n, S * P)
  for (s in seq_len(S)) {
    b_hat_s <- c(cm$blup_psi[s, ], cm$blup_p[s, ])
    b_draws <- .tobs_sbc_community_b_draws(mu_draws, fit$means, b_hat_s,
                                           cm$Bf[[s]], cm$Cinv[[s]], n)
    M[, (s - 1L) * P + seq_len(P)] <- mu_draws + b_draws
  }
  colnames(M) <- nm$cols
  M
}

.tobs_sbc_sim_ms_occu <- function(spec, theta, seed) {
  set.seed(seed)
  f <- spec$fit_obs
  m <- f$model
  nm <- .tobs_sbc_ms_occu_names(m)
  S <- m$n_species
  # do.call(rbind, lapply(...)), not t(vapply(...)) -- vapply collapses a
  # length-1 per-species result (a single detection coefficient, e.g.
  # `detection = ~ 1`) to a plain vector, and t() of that transposes the
  # wrong axis (1 x S instead of S x 1).
  coef_psi <- do.call(rbind, lapply(seq_len(S), function(s)
    theta[paste0(m$species_names[s], "_psi_", nm$occ)]))
  coef_p <- do.call(rbind, lapply(seq_len(S), function(s)
    theta[paste0(m$species_names[s], "_p_", nm$det)]))
  dimnames(coef_psi) <- list(m$species_names, nm$occ)
  dimnames(coef_p)   <- list(m$species_names, nm$det)
  f$ms_community$coef_psi <- coef_psi
  f$ms_community$coef_p   <- coef_p
  rep <- stats::simulate(f, nsim = 1L)
  obs <- spec$data_obs
  list(cells = obs$cells, y = rep, y_pos = NULL, visits = NULL, graph = NULL,
       site = obs$site)
}

.tobs_sbc_loglik_many_ms_occu <- function(fit, Theta) {
  m <- fit$model
  nm <- .tobs_sbc_ms_occu_names(m)
  S <- m$n_species
  vapply(seq_len(nrow(Theta)), function(i) {
    th <- Theta[i, ]
    ll <- 0
    for (s in seq_len(S)) {
      beta_psi <- th[paste0(m$species_names[s], "_psi_", nm$occ)]
      beta_p   <- th[paste0(m$species_names[s], "_p_",   nm$det)]
      eta_psi <- as.numeric(m$X_occ %*% beta_psi)
      eta_p   <- as.numeric(m$X_det %*% beta_p)
      ll <- ll + .ms_int_occu_sp_ll(eta_psi, list(eta_p), m$summaries[[s]])
    }
    ll
  }, numeric(1))
}


# ---------------------------------------------------------------------------
# 6j-bis. ms_occu_cover() (gcol33/tulpaObs#220, community group): the
# occ+p+pos analogue of ms_occu -- same "rank a fixed species set" design
# (theta_s = mu_arm + b_s per species, one block per arm, species-major),
# same joint mu/b_s draw (gcol33/tulpaObs#226 part 1, `.tobs_sbc_community_b_draws`),
# now safe to attempt because `Cinv` is kept consistent with `Sigma` even
# under AGHQ debiasing (#226 part 2 fix specific to this family, commit
# `03b87ad` -- `ms_occu`/`ms_int_occu`/`ms_count` have no AGHQ path at all and
# stay deferred).
#
# `fit$means`/`vcov` carry ONE extra coordinate beyond the three arms' P
# dimensions: the shared cover-arm dispersion (`log_sigma_pos` for lognormal,
# `log_phi` for beta), a SINGLE community-level scalar with no per-species
# deviation (`b_list[[s]]` in `R/ms_occu_cover.R` is P-long, arm coefs only).
# `.tobs_sbc_community_b_draws(mu_draws, mu_hat, ...)` needs no change for
# this -- `Bf_s` is `(P+1) x P` (mu+disp -> b_s), so the shift term already
# reads the disp coordinate; a species' theta gets the disp-informed shift on
# its arm coefs, while the disp value itself is copied straight off the
# u=(mu,disp) draw with no per-species addition.
#
# Reuses `.occu_cover_sp_ll` (the exact per-species two-state marginal
# `.tobs_fit_ms_occu_cover()` already optimizes, via
# `.ms_occu_cover_species_view()`) for `loglik_many` -- same "family's own
# kernel" pattern as `ms_occu`/`ms_int_occu`. `simulate()` reuses
# `.tobs_simulate_ms_occu_cover()`, overwriting `f$ms_community$coef_occ`/
# `coef_p`/`coef_pos` (the field names THAT handler reads -- note these are
# "occ"/"p"/"pos", not the "psi"/"p"/"pos" prefix `fit$means`'s own
# coefficient names use) and the last element of `f$means` (the shared
# log-dispersion `.tobs_simulate_ms_occu_cover` reads via
# `object$means[[length(object$means)]]`).
#
# `y`/`y_pos` are both `[site x visit x species]` -- the 3D-season pool
# generalizes to carry both response arrays, not just one.
#
# v1 scope: non-spatial `laplace`, `positive = "lognormal"` only (matches
# `occu_cover`/`cover`'s own v1 scope -- beta's mu/phi reparameterization is a
# separate correctness question, not attempted here).
# ---------------------------------------------------------------------------

.tobs_sbc_spec_ms_occu_cover <- function(fit, fit.control) {
  .tobs_sbc_reject_structure(fit)
  .tobs_sbc_reject_visit_design(fit)
  m <- fit$model
  if (!identical(m$positive, "lognormal")) {
    stop("SBC on ms_occu_cover() is registered for positive = \"lognormal\" ",
         "only (gcol33/tulpaObs#220).", call. = FALSE)
  }
  list(model   = m,
       fit_obs = fit,
       family  = attr(fit, "tobs_family"),
       method  = fit$method,
       control = utils::modifyList(.tobs_sbc_default_control(),
                                   as.list(fit.control)),
       occ     = .tobs_sbc_recombine(m$formulas$occ, NULL),
       det     = .tobs_sbc_recombine(m$formulas$det, NULL),
       pos     = .tobs_sbc_recombine(m$formulas$pos, NULL),
       species = m$species_names)
}

.tobs_sbc_refit_ms_occu_cover <- function(spec, data) {
  suppressWarnings(do.call(tobs, list(
    formula = spec$occ, data = data$cells, family = spec$family,
    detection = spec$det, positive = spec$pos,
    y = data$y, y_pos = data$y_pos, species = spec$species,
    method = spec$method, control = spec$control)))
}

# theta columns: S blocks of `psi_<coef>`/`p_<coef>`/`pos_<coef>` (species-
# major, matching ms_occu's layout), plus ONE trailing shared dispersion
# column (`log_sigma_pos`/`log_phi`) -- not per-species.
.tobs_sbc_ms_occu_cover_names <- function(m) {
  occ_names <- m$process_info[[1L]]$coef_names
  det_names <- m$process_info[[2L]]$coef_names
  pos_names <- m$process_info[[3L]]$coef_names
  par <- c(paste0("psi_", occ_names), paste0("p_", det_names),
          paste0("pos_", pos_names))
  disp_name <- if (identical(m$positive, "beta")) "log_phi" else "log_sigma_pos"
  list(occ = occ_names, det = det_names, pos = pos_names, par = par,
       disp_name = disp_name,
       cols = c(as.vector(outer(par, m$species_names,
                                function(p, sp) paste0(sp, "_", p))),
               disp_name))
}

.tobs_sbc_data_ms_occu_cover <- function(fit) {
  m <- fit$model
  list(cells = m$data, y = m$y, y_pos = m$y_pos, visits = NULL, graph = NULL,
       site = seq_len(m$n_sites))
}

.tobs_sbc_pool_ms_occu_cover <- function(obs, rep) {
  n_o <- dim(obs$y)[1L]; n_r <- dim(rep$y)[1L]
  y <- array(NA_real_, dim = c(n_o + n_r, dim(obs$y)[2L], dim(obs$y)[3L]))
  y[seq_len(n_o), , ] <- obs$y
  y[n_o + seq_len(n_r), , ] <- rep$y
  y_pos <- array(NA_real_, dim = c(n_o + n_r, dim(obs$y_pos)[2L], dim(obs$y_pos)[3L]))
  y_pos[seq_len(n_o), , ] <- obs$y_pos
  y_pos[n_o + seq_len(n_r), , ] <- rep$y_pos
  list(cells = .tobs_sbc_rbind_cells(obs$cells, rep$cells), y = y, y_pos = y_pos,
       visits = NULL, graph = NULL, site = c(obs$site, rep$site + n_o))
}

.tobs_sbc_draws_ms_occu_cover <- function(fit, n) {
  m <- fit$model
  nm <- .tobs_sbc_ms_occu_cover_names(m)
  S <- m$n_species; P <- length(nm$par)
  cm <- fit$ms_community
  mu_draws <- .tobs_sbc_mvn_draws(fit$means, fit$vcov, n)   # n x (P+1)
  mu_arm <- mu_draws[, seq_len(P), drop = FALSE]
  M <- matrix(NA_real_, n, S * P + 1L)
  for (s in seq_len(S)) {
    b_hat_s <- c(cm$blup_occ[s, ], cm$blup_p[s, ], cm$blup_pos[s, ])
    b_draws <- .tobs_sbc_community_b_draws(mu_draws, fit$means, b_hat_s,
                                           cm$Bf[[s]], cm$Cinv[[s]], n)
    M[, (s - 1L) * P + seq_len(P)] <- mu_arm + b_draws
  }
  M[, S * P + 1L] <- mu_draws[, P + 1L]
  colnames(M) <- nm$cols
  M
}

.tobs_sbc_sim_ms_occu_cover <- function(spec, theta, seed) {
  set.seed(seed)
  f <- spec$fit_obs
  m <- f$model
  nm <- .tobs_sbc_ms_occu_cover_names(m)
  S <- m$n_species
  coef_occ <- do.call(rbind, lapply(seq_len(S), function(s)
    theta[paste0(m$species_names[s], "_psi_", nm$occ)]))
  coef_p <- do.call(rbind, lapply(seq_len(S), function(s)
    theta[paste0(m$species_names[s], "_p_", nm$det)]))
  coef_pos <- do.call(rbind, lapply(seq_len(S), function(s)
    theta[paste0(m$species_names[s], "_pos_", nm$pos)]))
  dimnames(coef_occ) <- list(m$species_names, nm$occ)
  dimnames(coef_p)   <- list(m$species_names, nm$det)
  dimnames(coef_pos) <- list(m$species_names, nm$pos)
  f$ms_community$coef_occ <- coef_occ
  f$ms_community$coef_p   <- coef_p
  f$ms_community$coef_pos <- coef_pos
  f$means[[length(f$means)]] <- theta[[nm$disp_name]]
  rep <- stats::simulate(f, nsim = 1L)
  obs <- spec$data_obs
  list(cells = obs$cells, y = rep$y, y_pos = rep$y_pos, visits = NULL,
       graph = NULL, site = obs$site)
}

.tobs_sbc_loglik_many_ms_occu_cover <- function(fit, Theta) {
  m <- fit$model
  nm <- .tobs_sbc_ms_occu_cover_names(m)
  S <- m$n_species
  views <- lapply(seq_len(S), function(s) .ms_occu_cover_species_view(m, s))
  vapply(seq_len(nrow(Theta)), function(i) {
    th <- Theta[i, ]
    ld <- th[[nm$disp_name]]
    ll <- 0
    for (s in seq_len(S)) {
      bo   <- th[paste0(m$species_names[s], "_psi_", nm$occ)]
      bp   <- th[paste0(m$species_names[s], "_p_",   nm$det)]
      bpos <- th[paste0(m$species_names[s], "_pos_", nm$pos)]
      ll <- ll + .occu_cover_sp_ll(views[[s]], bo, bp, bpos, ld)
    }
    ll
  }, numeric(1))
}


# ---------------------------------------------------------------------------
# 6k. cover() (gcol33/tulpaObs#220, multiarm-S3 group): a `tobs_multiarm_fit`
# with the same two-independent-Laplace-Gaussian-block shape as
# occu_categorical (presence, positive), which is why `.tobs_sbc_data_vector`/
# `.tobs_sbc_pool_vector` (renamed from the occu_categorical-only spelling,
# section 6e) are reused unchanged here too. Two things this family did NOT
# already carry, both now added (small, independently useful, kept regardless
# of SBC): `encoding$data`/`encoding$y` (the raw fitting inputs, needed to
# rebuild a refit call -- occu_categorical needed the identical addition) and
# `V_occ`/`V_pos` (the full per-arm coefficient covariance, previously only
# the diagonal `se_occ`/`se_pos` was stored -- needed to draw a JOINT sample
# of an arm's coefficients rather than assume independence across them,
# which the ms_occu near-miss (#226) is a direct warning against).
#
# v1 scope: non-spatial `laplace`, `positive = "lognormal"` only
# (`simulate_cover()`'s own generator only covers lognormal/gaussian, and
# lognormal's `log(y_pos) ~ N(eta_pos, sigma_pos)` is simpler to get right
# than beta's mu/phi -> shape1/shape2 reparameterization). `sigma_pos` is
# NOT part of theta -- the fitted point estimate carries no SE anywhere in
# the package (`sigma_pos_sd` is hardcoded NA), so it is held fixed at the
# observed fit's value for both the replicate and the reference score,
# matching the established `.tobs_sbc_check_fixed_dispersion` convention
# elsewhere in this codebase. `X_occ`/`X_pos` are rebuilt fresh via
# `model.matrix()` from the arm formula + natural-scale data, NEVER read off
# `encoding$occ_data$X`/`pos_data$X` -- those are autoscaled design matrices
# (gcol33/tulpaObs#9) paired with SCALED coefficients internally, while
# `beta_occ`/`beta_pos` on the fit are already unscaled to natural units
# (`.unscale_beta_vec`); mixing the two would reproduce the exact class of
# bug #225 found in `int_occu()`.
# ---------------------------------------------------------------------------

.tobs_sbc_spec_cover <- function(fit, fit.control) {
  .tobs_sbc_reject_structure(fit)
  if (!identical(fit$positive, "lognormal")) {
    stop("SBC on cover() is registered for positive = \"lognormal\" only ",
         "(got \"", fit$positive, "\"); the beta/beta_oi/lognormal_trunc/",
         "ordinal/gaussian arms are follow-ups.", call. = FALSE)
  }
  enc <- fit$encoding
  list(presence  = .tobs_sbc_recombine(enc$fe_occ, NULL),
       positive  = .tobs_sbc_recombine(enc$fe_pos, NULL),
       family    = attr(fit, "tobs_family"),
       method    = fit$method %||% "laplace",
       control   = utils::modifyList(.tobs_sbc_default_control(),
                                     as.list(fit.control)),
       sigma_pos = fit$sigma_pos,
       fit_obs   = fit)
}

.tobs_sbc_refit_cover <- function(spec, data) {
  suppressWarnings(do.call(tobs, list(
    formula = spec$presence, data = data$cells, family = spec$family,
    presence = spec$presence, positive = spec$positive, y = data$y,
    method = spec$method, control = spec$control)))
}

.tobs_sbc_draws_cover <- function(fit, n) {
  occ_names <- names(fit$beta_occ)
  pos_names <- names(fit$beta_pos)
  D_occ <- .tobs_sbc_mvn_draws(fit$beta_occ, fit$V_occ, n)
  D_pos <- .tobs_sbc_mvn_draws(fit$beta_pos, fit$V_pos, n)
  M <- cbind(D_occ, D_pos)
  colnames(M) <- c(paste0("occ_", occ_names), paste0("pos_", pos_names))
  M
}

.tobs_sbc_sim_cover <- function(spec, theta, seed) {
  set.seed(seed)
  obs <- spec$data_obs
  X_occ <- stats::model.matrix(spec$presence, obs$cells)
  X_pos <- stats::model.matrix(spec$positive, obs$cells)
  beta_occ <- theta[paste0("occ_", colnames(X_occ))]
  beta_pos <- theta[paste0("pos_", colnames(X_pos))]
  eta_occ <- as.vector(X_occ %*% beta_occ)
  eta_pos <- as.vector(X_pos %*% beta_pos)
  n <- nrow(X_occ)
  occur <- stats::rbinom(n, 1L, stats::plogis(eta_occ))
  log_cover <- stats::rnorm(n, eta_pos, spec$sigma_pos)
  cover <- ifelse(occur == 1L, pmin(exp(log_cover), 1), 0)
  list(cells = obs$cells, y = cover, y_pos = NULL, visits = NULL, graph = NULL,
       site = obs$site)
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
# 6l. ms_int_occu() (gcol33/tulpaObs#220, community group): REGISTERED --
# the community analogue of int_occu(): `model$y` is a list of D per-source
# 3D [n_sites x visits_d x n_species] arrays sharing the site axis
# (`.tobs_sbc_pool_named_3d` from dyn_int_occu() reused unchanged), ranking
# the fixed species set's own realized coefficients via the joint mu/b_s
# draw `.tobs_sbc_community_b_draws` (#226 part 1). Originally found to share
# ms_occu's exact failure mode (a direct probe at three seeds found
# `sp3_p2_(Intercept)` stuck at p_unif ~0.0029-0.0030, reproducible not
# noise) -- and, like ms_occu, the actual cause turned out to be species
# count, not a bug: at S=14 (matching this family's own recovery-test
# fixture) the plain Laplace-EM calibrates cleanly, no debiasing needed (see
# the registry entry below for the numbers). `simulate()` injects
# `ms_community$coef_psi`/`coef_p<d>` directly rather than the `f$draws`
# trick since `.tobs_simulate_ms_int_occu()` doesn't consult `object$draws`;
# `loglik_many` sums each species' own exact two-state marginal via
# `.ms_int_occu_sp_ll`. Full site-overlap only (matching int_occu()'s
# existing gate).
# ---------------------------------------------------------------------------

.tobs_sbc_data_ms_int_occu <- function(fit) {
  m <- fit$model
  # Named by `process_names` ("p1".."pD"), NOT `source_names` -- the family's
  # own `.tobs_simulate_ms_int_occu()` names its replicate list that way
  # (`names(srcs) <- model$process_names`), and `source_names` defaults to a
  # DIFFERENT label ("src1".."srcD") whenever `y` came in unnamed, which would
  # silently mismatch the two lists' keys at pooling.
  list(cells = m$data, y = stats::setNames(m$y, m$process_names), y_pos = NULL,
       visits = NULL, graph = NULL, site = seq_len(m$n_sites))
}

.tobs_sbc_spec_ms_int_occu <- function(fit, fit.control) {
  .tobs_sbc_reject_structure(fit)
  m <- fit$model
  full <- vapply(m$site_maps, function(mp)
    identical(as.integer(mp), seq_len(m$n_sites)), logical(1))
  if (!all(full)) {
    stop("SBC on ms_int_occu() is registered for full-overlap fits (every ",
         "source observes every site); this fit has partial overlap.",
         call. = FALSE)
  }
  list(model     = m,
       fit_obs   = fit,
       family    = attr(fit, "tobs_family"),
       method    = fit$method,
       control   = utils::modifyList(.tobs_sbc_default_control(),
                                     as.list(fit.control)),
       occ       = .tobs_sbc_recombine(m$formulas$occ, NULL),
       det_list  = lapply(m$formulas$det, function(f) .tobs_sbc_recombine(f, NULL)),
       species   = m$species_names,
       sources   = m$process_names)
}

.tobs_sbc_refit_ms_int_occu <- function(spec, data) {
  suppressWarnings(do.call(tobs, list(
    formula = spec$occ, data = data$cells, family = spec$family,
    detection = spec$det_list, y = data$y[spec$sources], species = spec$species,
    method = spec$method, control = spec$control)))
}

# theta column names: `<species>_psi_<coef>` / `<species>_p<d>_<coef>`, one
# block of length P per species, species-major -- the layout `draws()`,
# `simulate()` and `loglik_many()` all share. `arm_names` c("p1", ..., "pD")
# matches `model$process_names`, the SAME labels `ms_community`'s
# `coef_p<d>` / `blup_p<d>` fields are keyed by (build_ms_int_occu_fit()).
.tobs_sbc_ms_int_occu_names <- function(m) {
  D <- m$n_sources
  occ_names <- m$process_info[[1L]]$coef_names
  det_names_list <- lapply(seq_len(D), function(d) m$process_info[[1L + d]]$coef_names)
  par <- c(paste0("psi_", occ_names),
          unlist(lapply(seq_len(D), function(d)
            paste0("p", d, "_", det_names_list[[d]]))))
  list(occ_names = occ_names, det_names_list = det_names_list, par = par,
       cols = as.vector(outer(par, m$species_names,
                              function(p, sp) paste0(sp, "_", p))))
}

.tobs_sbc_draws_ms_int_occu <- function(fit, n) {
  m <- fit$model
  nm <- .tobs_sbc_ms_int_occu_names(m)
  D <- m$n_sources; S <- m$n_species; P <- length(nm$par)
  cm <- fit$ms_community
  mu_draws <- .tobs_sbc_mvn_draws(fit$means, fit$vcov, n)   # n x P
  M <- matrix(NA_real_, n, S * P)
  for (s in seq_len(S)) {
    b_hat_s <- c(cm$blup_psi[s, ],
                unlist(lapply(seq_len(D), function(d) cm[[paste0("blup_p", d)]][s, ])))
    # gcol33/tulpaObs#226 (see .tobs_sbc_community_b_draws): mu and b_s draw
    # jointly, not independently.
    b_draws <- .tobs_sbc_community_b_draws(mu_draws, fit$means, b_hat_s,
                                           cm$Bf[[s]], cm$Cinv[[s]], n)
    M[, (s - 1L) * P + seq_len(P)] <- mu_draws + b_draws
  }
  colnames(M) <- nm$cols
  M
}

.tobs_sbc_sim_ms_int_occu <- function(spec, theta, seed) {
  set.seed(seed)
  f <- spec$fit_obs
  m <- f$model
  nm <- .tobs_sbc_ms_int_occu_names(m)
  D <- m$n_sources; S <- m$n_species
  coef_psi <- do.call(rbind, lapply(seq_len(S), function(s)
    theta[paste0(m$species_names[s], "_psi_", nm$occ_names)]))
  dimnames(coef_psi) <- list(m$species_names, nm$occ_names)
  f$ms_community$coef_psi <- coef_psi
  for (d in seq_len(D)) {
    arm <- paste0("p", d); cn <- nm$det_names_list[[d]]
    coef_d <- do.call(rbind, lapply(seq_len(S), function(s)
      theta[paste0(m$species_names[s], "_", arm, "_", cn)]))
    dimnames(coef_d) <- list(m$species_names, cn)
    f$ms_community[[paste0("coef_", arm)]] <- coef_d
  }
  rep <- stats::simulate(f, nsim = 1L)
  obs <- spec$data_obs
  list(cells = obs$cells, y = rep, y_pos = NULL, visits = NULL, graph = NULL,
       site = obs$site)
}

.tobs_sbc_loglik_many_ms_int_occu <- function(fit, Theta) {
  m <- fit$model
  nm <- .tobs_sbc_ms_int_occu_names(m)
  D <- m$n_sources; S <- m$n_species
  vapply(seq_len(nrow(Theta)), function(i) {
    th <- Theta[i, ]
    ll <- 0
    for (s in seq_len(S)) {
      beta_psi <- th[paste0(m$species_names[s], "_psi_", nm$occ_names)]
      eta_psi <- as.numeric(m$X_psi %*% beta_psi)
      eta_p <- lapply(seq_len(D), function(d) {
        arm <- paste0("p", d)
        beta_d <- th[paste0(m$species_names[s], "_", arm, "_", nm$det_names_list[[d]])]
        as.numeric(m$X_p[[d]] %*% beta_d)
      })
      ll <- ll + .ms_int_occu_sp_ll(eta_psi, eta_p, m$summaries[[s]])
    }
    ll
  }, numeric(1))
}


# ---------------------------------------------------------------------------
# 6m. occu_multiscale_cover() (gcol33/tulpaObs#220, multiarm-S3 group): a
# THREE-level hurdle -- cell occupancy psi gates plot availability theta,
# which gates per-visit detection p and the cover hurdle pos -- fit under
# `method = "laplace"` (non-spatial: iid cells, field fixed at 0). Unlike
# `cover()`, this is the STANDARD single-block fit shape: `fit$means`/
# `fit$draws` are a real joint MVN (dispersion `log_sigma_pos` optimized
# jointly with the coefficients via BFGS + `optimHess`, so it carries a
# proper SE, unlike `cover()`'s post-hoc MOM estimate) and
# `.tobs_pointwise_loglik` already dispatches for this model_type -- so
# `draws`/`loglik_many` are the shared generic ones, unlike every other
# family in this section.
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
# `formula` is REQUIRED to carry an `icar(graph =, group_var = "<col>")`
# term declaring the plot -> cell map, even under the non-spatial engine
# (the graph itself is ignored there, but the term's presence and its
# group_var column are not optional syntax). Rather than track and
# reproduce the OBSERVED fit's own group_var column name through pooling, a
# fixed synthetic column (`.sbc_cell`) is written into `cells` by both
# `data()` and `simulate()`, with a same-size dummy graph built at refit
# time -- the group_var mapping the refit reads is entirely the SBC
# adapter's own, decoupled from whatever the original fit used.
#
# `simulate()` has no family handler to inject theta into (unlike `cover()`
# reusing `.tobs_simulate_ms_occu()`-style machinery), so it is a
# hand-written three-level generator: z_c ~ Bernoulli(psi_c) per cell,
# a_j | z_{cell(j)} ~ Bernoulli(theta_j) per plot, y_jv | a_j ~
# Bernoulli(p_j) per visit, cover | y_jv=1 ~ Lognormal(eta_pos_j,
# sigma_pos). `positive = "lognormal"` only, matching `cover()`'s v1 scope
# and reasoning (no beta-arm generator to mirror, and dispersion carries a
# real SE here so there is no dispersion-fixing complication to work
# around). No visit-level covariates for v1 (`.tobs_sbc_reject_visit_design`
# rejects them).
# ---------------------------------------------------------------------------

.tobs_sbc_reject_occu_mscale_cover_scope <- function(fit) {
  m <- fit$model
  if (!identical(m$positive %||% "lognormal", "lognormal")) {
    stop("SBC on occu_multiscale_cover() is registered for positive = ",
         "\"lognormal\" only (got \"", m$positive, "\").", call. = FALSE)
  }
}

.tobs_sbc_spec_occu_mscale_cover <- function(fit, fit.control) {
  .tobs_sbc_reject_structure(fit)
  .tobs_sbc_reject_visit_design(fit)
  .tobs_sbc_reject_occu_mscale_cover_scope(fit)
  m <- fit$model
  list(model   = m,
       fit_obs = fit,
       family  = attr(fit, "tobs_family"),
       method  = fit$method %||% "laplace",
       control = utils::modifyList(.tobs_sbc_default_control(),
                                   as.list(fit.control)),
       psi     = .tobs_sbc_recombine(m$formulas$psi,   NULL),
       theta   = .tobs_sbc_recombine(m$formulas$theta, NULL),
       det     = .tobs_sbc_recombine(m$formulas$p,     NULL),
       pos     = .tobs_sbc_recombine(m$formulas$pos,   NULL))
}

.tobs_sbc_data_occu_mscale_cover <- function(fit) {
  m <- fit$model
  cells <- m$data
  cells$.sbc_cell <- m$plot_cell
  list(cells = cells, y = m$y, y_pos = m$y_pos, visits = NULL, graph = NULL,
       site = m$plot_cell, n_cells = m$n_cells)
}

.tobs_sbc_pool_occu_mscale_cover <- function(obs, rep) {
  n_c_o <- obs$n_cells
  rep_cells <- rep$cells
  rep_cells$.sbc_cell <- rep_cells$.sbc_cell + n_c_o
  list(cells   = .tobs_sbc_rbind_cells(obs$cells, rep_cells),
       y       = rbind(obs$y, rep$y),
       y_pos   = rbind(obs$y_pos, rep$y_pos),
       visits  = NULL, graph = NULL,
       site    = c(obs$site, rep$site + n_c_o),
       n_cells = n_c_o + rep$n_cells)
}

.tobs_sbc_refit_occu_mscale_cover <- function(spec, data) {
  n_cells <- max(data$cells$.sbc_cell)
  dummy_graph <- matrix(0, n_cells, n_cells)
  fe_labels <- attr(stats::terms(spec$psi), "term.labels")
  icpt <- attr(stats::terms(spec$psi), "intercept")
  rhs <- paste(c(fe_labels,
                'icar(graph = dummy_graph, group_var = ".sbc_cell")'),
              collapse = " + ")
  if (!icpt) rhs <- paste(rhs, "- 1")
  full_formula <- stats::as.formula(paste("~", rhs), env = environment())
  suppressWarnings(do.call(tobs, list(
    formula = full_formula, data = data$cells, family = spec$family,
    detection = spec$det, availability = spec$theta, positive = spec$pos,
    y = data$y, y_pos = data$y_pos,
    method = spec$method, control = spec$control)))
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
  list(cells = obs$cells, y = y, y_pos = y_pos, visits = NULL, graph = NULL,
       site = obs$site, n_cells = obs$n_cells)
}


# ---------------------------------------------------------------------------
# 6n. ms_count() (gcol33/tulpaObs#220, community group): community relative-
# abundance GLMM, no detection/latent state (the abundance analogue of
# ms_occu). `model$y` is a plain [n_sites x n_species] matrix (no visit
# axis), so the SHARED generic `.tobs_sbc_pool` (2D `rbind`) applies
# unchanged. Ranks the FIXED species set's own realized coefficients (mu +
# b_s per species), the same design as ms_occu/ms_int_occu -- MULTI-SEED
# tested (0, 1, 2) before registering, per the lesson #226 and the
# ms_int_occu near-miss both taught: a single seed clearing the 1e-3 bar is
# not evidence of calibration. `simulate()` is custom, overwriting
# `ms_community$coef_mu` directly (`.tobs_simulate_ms_count()` reads that,
# not `object$draws`, the same mechanism cover()/ms_occu()/ms_int_occu()
# use for the identical reason). `loglik_many` sums each species' own exact
# Poisson log-likelihood (`.ms_count_ll_pois`, the same kernel the fitter
# optimizes) at that species' theta slice. `response = "poisson"` only for
# v1 (negbin/gaussian/bernoulli/binomial carry an extra per-species
# dispersion arm, a follow-up); `jsdm()` is `ms_count(response =
# "bernoulli")` under the hood and is NOT covered by this registration --
# a different response family with its own calibration to check. Originally
# deferred on a small fixture's failure (matching ms_occu/ms_int_occu); like
# both of those, the cause was species count, not a bug, and the registry
# entry below uses a fixture at S=20 -- see gcol33/tulpaObs#226.
# ---------------------------------------------------------------------------

.tobs_sbc_data_ms_count <- function(fit) {
  m <- fit$model
  list(cells = m$data, y = m$y, y_pos = NULL, visits = NULL, graph = NULL,
       site = seq_len(m$n_sites))
}

.tobs_sbc_spec_ms_count <- function(fit, fit.control) {
  .tobs_sbc_reject_structure(fit)
  m <- fit$model
  if (!identical(m$response, "poisson")) {
    stop("SBC on ms_count() is registered for response = \"poisson\" only ",
         "(got \"", m$response, "\"); negbin/gaussian/bernoulli/binomial are ",
         "follow-ups.", call. = FALSE)
  }
  list(model   = m,
       fit_obs = fit,
       family  = attr(fit, "tobs_family"),
       method  = fit$method,
       control = utils::modifyList(.tobs_sbc_default_control(),
                                   as.list(fit.control)),
       formula = .tobs_sbc_recombine(m$formulas$mu, NULL),
       species = m$species_names)
}

.tobs_sbc_refit_ms_count <- function(spec, data) {
  suppressWarnings(do.call(tobs, list(
    formula = spec$formula, data = data$cells, family = spec$family,
    y = data$y, species = spec$species,
    method = spec$method, control = spec$control)))
}

.tobs_sbc_ms_count_names <- function(m) {
  mu_names <- m$process_info[[1L]]$coef_names
  list(mu_names = mu_names,
       cols = as.vector(outer(mu_names, m$species_names,
                              function(p, sp) paste0(sp, "_mu_", p))))
}

.tobs_sbc_draws_ms_count <- function(fit, n) {
  m <- fit$model
  nm <- .tobs_sbc_ms_count_names(m)
  S <- m$n_species; P <- length(nm$mu_names)
  cm <- fit$ms_community
  mu_draws <- .tobs_sbc_mvn_draws(fit$means, fit$vcov, n)
  M <- matrix(NA_real_, n, S * P)
  for (s in seq_len(S)) {
    # gcol33/tulpaObs#226 (see .tobs_sbc_community_b_draws): mu and b_s draw
    # jointly, not independently.
    b_draws <- .tobs_sbc_community_b_draws(mu_draws, fit$means, cm$blup_mu[s, ],
                                           cm$Bf[[s]], cm$Cinv[[s]], n)
    M[, (s - 1L) * P + seq_len(P)] <- mu_draws + b_draws
  }
  colnames(M) <- nm$cols
  M
}

.tobs_sbc_sim_ms_count <- function(spec, theta, seed) {
  set.seed(seed)
  f <- spec$fit_obs
  m <- f$model
  nm <- .tobs_sbc_ms_count_names(m)
  S <- m$n_species
  coef_mu <- do.call(rbind, lapply(seq_len(S), function(s)
    theta[paste0(m$species_names[s], "_mu_", nm$mu_names)]))
  dimnames(coef_mu) <- list(m$species_names, nm$mu_names)
  f$ms_community$coef_mu <- coef_mu
  rep <- stats::simulate(f, nsim = 1L)
  obs <- spec$data_obs
  list(cells = obs$cells, y = rep, y_pos = NULL, visits = NULL, graph = NULL,
       site = obs$site)
}

.tobs_sbc_loglik_many_ms_count <- function(fit, Theta) {
  m <- fit$model
  nm <- .tobs_sbc_ms_count_names(m)
  S <- m$n_species
  vapply(seq_len(nrow(Theta)), function(i) {
    th <- Theta[i, ]
    ll <- 0
    for (s in seq_len(S)) {
      beta <- th[paste0(m$species_names[s], "_mu_", nm$mu_names)]
      ll <- ll + .ms_count_ll_pois(m$summaries[[s]], beta)
    }
    ll
  }, numeric(1))
}


# ---------------------------------------------------------------------------
# 6o. jsdm() (gcol33/tulpaObs#220, community group): `jsdm()` is
# `ms_count(response = "bernoulli")` under the hood (`.dispatch_jsdm` calls
# `.tobs_build_ms_count(..., response = "bernoulli")` and routes through the
# same community Laplace-EM / latent driver / NUTS target as ms_count() --
# gcol33/tulpaObs#121). The fit object shares ms_count()'s `fit$model`
# shape byte for byte, so `data`/`draws`/`simulate` are the exact same
# generic helpers ms_count() already registers (`.tobs_sbc_data_ms_count`,
# `.tobs_sbc_draws_ms_count`, `.tobs_sbc_sim_ms_count` -- the last reads
# `f$ms_community$coef_mu` and calls `stats::simulate()`, neither of which
# is response-specific). `attr(fit, "tobs_family")` is "jsdm", not
# "ms_count", though (the family ctor's own name), so the SBC dispatch
# needs its own registry row regardless. Only `spec` (constructs `jsdm()`,
# not `ms_count(response=)`) and `loglik_many` (the Bernoulli kernel
# `.ms_count_ll_bern`, not `.ms_count_ll_pois`) are new.
# ---------------------------------------------------------------------------

.tobs_sbc_spec_jsdm <- function(fit, fit.control) {
  .tobs_sbc_reject_structure(fit)
  m <- fit$model
  list(model   = m,
       fit_obs = fit,
       family  = attr(fit, "tobs_family"),
       method  = fit$method,
       control = utils::modifyList(.tobs_sbc_default_control(),
                                   as.list(fit.control)),
       formula = .tobs_sbc_recombine(m$formulas$mu, NULL),
       species = m$species_names)
}

.tobs_sbc_refit_jsdm <- function(spec, data) {
  suppressWarnings(do.call(tobs, list(
    formula = spec$formula, data = data$cells, family = jsdm(),
    y = data$y, species = spec$species,
    method = spec$method, control = spec$control)))
}

.tobs_sbc_loglik_many_jsdm <- function(fit, Theta) {
  m <- fit$model
  nm <- .tobs_sbc_ms_count_names(m)
  S <- m$n_species
  vapply(seq_len(nrow(Theta)), function(i) {
    th <- Theta[i, ]
    ll <- 0
    for (s in seq_len(S)) {
      beta <- th[paste0(m$species_names[s], "_mu_", nm$mu_names)]
      ll <- ll + .ms_count_ll_bern(m$summaries[[s]], beta)
    }
    ll
  }, numeric(1))
}


# ---------------------------------------------------------------------------
# 6p. ms_distance() (gcol33/tulpaObs#220, community group): community binned
# distance sampling (the spAbundance msDS analogue), the two-arm
# (lambda/sigma) sibling of ms_occu -- same "rank a fixed species set" design
# (theta_s = mu_arm + b_s per species, species-major), same joint mu/b_s draw
# (`.tobs_sbc_community_b_draws`, #226 part 1). `.tobs_community_em()`
# already exposes `Cinv`/`Bf` on this family's fit (R/ms_distance.R,
# `build_ms_distance_fit()` -- added alongside the #226 fix, not something
# this registration needed to add). `model$y` is a plain 3D
# `[n_sites x n_bins x n_species]` array, structurally identical to
# ms_occu()'s `[n_sites x max_visits x n_species]` -- the generic
# `.tobs_sbc_data_3d_season`/`.tobs_sbc_pool_3d_season` pair (stacking axis 1
# while leaving axes 2-3 alone) applies unchanged, species standing in for
# season. `simulate()` reuses the family's own `.tobs_simulate_ms_distance()`
# handler (draws through `cpp_simulate_distance`, the SAME kernel the
# likelihood integrates against), overwriting `f$ms_community$coef_lambda`/
# `coef_sigma`. `loglik_many()` sums each species' own exact binned-distance
# marginal via `.tobs_ms_distance_engine(m)$sweep(s, eta_lam, eta_sig,
# value_only = TRUE)$log_lik` (the same kernel the fitter and
# `simulate_ms_distance()` optimize against/simulate through), NOT a
# hand-rolled likelihood.
#
# Registered `key = "halfnorm"`, `mixture = "poisson"` only, matching
# `.tobs_simulate_ms_distance()`'s own restriction (the hazard key's
# log-shape rides the community EM's `global` slot, not `ms_community`, so
# neither `simulate()` nor this registration can reach it yet -- #227);
# negbin is a second per-species dispersion arm this registration does not
# cover either.
# ---------------------------------------------------------------------------

.tobs_sbc_reject_ms_distance_scope <- function(fit) {
  m <- fit$model
  if (!identical(m$key, "halfnorm") || !identical(m$mixture, "poisson")) {
    stop("SBC on ms_distance() is registered for key = \"halfnorm\", ",
         "mixture = \"poisson\" only (the hazard key's log-shape and negbin's ",
         "dispersion are not carried by simulate(); gcol33/tulpaObs#227). Got ",
         "key = \"", m$key, "\", mixture = \"", m$mixture, "\".", call. = FALSE)
  }
}

# theta column names: `<species>_lambda_<coef>` / `<species>_sigma_<coef>`,
# one block of length P per species, species-major.
.tobs_sbc_ms_distance_names <- function(m) {
  lam_nm <- m$process_info[[1L]]$coef_names
  sig_nm <- m$process_info[[2L]]$coef_names
  par <- c(paste0("lambda_", lam_nm), paste0("sigma_", sig_nm))
  list(lam = lam_nm, sig = sig_nm, par = par,
       cols = as.vector(outer(par, m$species_names,
                              function(p, sp) paste0(sp, "_", p))))
}

.tobs_sbc_spec_ms_distance <- function(fit, fit.control) {
  .tobs_sbc_reject_structure(fit)
  .tobs_sbc_reject_ms_distance_scope(fit)
  m <- fit$model
  list(model     = m,
       fit_obs   = fit,
       family    = attr(fit, "tobs_family"),
       method    = fit$method,
       control   = utils::modifyList(.tobs_sbc_default_control(),
                                     as.list(fit.control)),
       lambda    = .tobs_sbc_recombine(m$formulas$lambda, NULL),
       sigma     = .tobs_sbc_recombine(m$formulas$sigma, NULL),
       cutpoints = m$cutpoints, key = m$key, transect = m$transect,
       mixture   = m$mixture,
       species   = m$species_names)
}

.tobs_sbc_refit_ms_distance <- function(spec, data) {
  suppressWarnings(do.call(tobs, list(
    formula = spec$lambda, data = data$cells,
    family = ms_distance(key = spec$key, transect = spec$transect,
                         cutpoints = spec$cutpoints, mixture = spec$mixture),
    detection = spec$sigma, y = data$y, species = spec$species,
    method = spec$method, control = spec$control)))
}

.tobs_sbc_draws_ms_distance <- function(fit, n) {
  m <- fit$model
  nm <- .tobs_sbc_ms_distance_names(m)
  S <- m$n_species; P <- length(nm$par)
  cm <- fit$ms_community
  mu_draws <- .tobs_sbc_mvn_draws(fit$means, fit$vcov, n)
  M <- matrix(NA_real_, n, S * P)
  for (s in seq_len(S)) {
    b_hat_s <- c(cm$blup_lambda[s, ], cm$blup_sigma[s, ])
    b_draws <- .tobs_sbc_community_b_draws(mu_draws, fit$means, b_hat_s,
                                           cm$Bf[[s]], cm$Cinv[[s]], n)
    M[, (s - 1L) * P + seq_len(P)] <- mu_draws + b_draws
  }
  colnames(M) <- nm$cols
  M
}

.tobs_sbc_sim_ms_distance <- function(spec, theta, seed) {
  set.seed(seed)
  f <- spec$fit_obs
  m <- f$model
  nm <- .tobs_sbc_ms_distance_names(m)
  S <- m$n_species
  coef_lambda <- do.call(rbind, lapply(seq_len(S), function(s)
    theta[paste0(m$species_names[s], "_lambda_", nm$lam)]))
  coef_sigma <- do.call(rbind, lapply(seq_len(S), function(s)
    theta[paste0(m$species_names[s], "_sigma_", nm$sig)]))
  dimnames(coef_lambda) <- list(m$species_names, nm$lam)
  dimnames(coef_sigma)  <- list(m$species_names, nm$sig)
  f$ms_community$coef_lambda <- coef_lambda
  f$ms_community$coef_sigma  <- coef_sigma
  rep <- stats::simulate(f, nsim = 1L)
  obs <- spec$data_obs
  list(cells = obs$cells, y = rep, y_pos = NULL, visits = NULL, graph = NULL,
       site = obs$site)
}

.tobs_sbc_loglik_many_ms_distance <- function(fit, Theta) {
  m <- fit$model
  nm <- .tobs_sbc_ms_distance_names(m)
  S <- m$n_species
  eng <- .tobs_ms_distance_engine(m)
  X_lambda <- m$X_processes[[1L]]; X_sigma <- m$X_processes[[2L]]
  vapply(seq_len(nrow(Theta)), function(i) {
    th <- Theta[i, ]
    ll <- 0
    for (s in seq_len(S)) {
      beta_lam <- th[paste0(m$species_names[s], "_lambda_", nm$lam)]
      beta_sig <- th[paste0(m$species_names[s], "_sigma_", nm$sig)]
      eta_lam <- as.numeric(X_lambda %*% beta_lam)
      eta_sig <- as.numeric(X_sigma %*% beta_sig)
      sw <- eng$sweep(s, eta_lam, eta_sig, value_only = TRUE)
      ll <- ll + sum(sw$log_lik)
    }
    ll
  }, numeric(1))
}


# ---------------------------------------------------------------------------
# 6q. ms_dyn_occu() (gcol33/tulpaObs#220, community group): community DYNAMIC
# (multi-season) occupancy, the HMM-forward analogue of ms_occu. `model$y` is
# a 4D `[n_sites x max_visits x n_seasons x n_species]` array -- ONE axis
# past every existing pool helper, so this section adds the generic
# `.tobs_sbc_data_4d_species`/`.tobs_sbc_pool_4d_species` pair (stack axis 1,
# leave axes 2-4 alone; the same recipe as `_3d_season`, one more dimension).
#
# Two of the four arms carry a per-species deviation (psi1, p -- the "rank a
# fixed species set" design, same as every other community family); the
# other two (gamma, eps, season-to-season colonization/extinction) are
# SHARED community globals with no per-species term at all (`R/ms_dyn_occu.R`
# section comment: "colonization/extinction coefficients carry no per-species
# random effect... shared community globals"). `fit$means`/`fit$vcov` already
# carry all four arms jointly (psi1, p, gamma, eps concatenated, ONE MVN), so
# `draws()` draws that whole vector once per row and conditions each
# species' b_s on the FULL draw via `.tobs_sbc_community_b_draws` (#226 part
# 1) -- `fit$ms_community$Bf[[s]]` is `(P + G) x P` (P = the two RE arms'
# width, G = gamma + eps), the cross-Hessian block against the WHOLE (mu,
# global) vector, not just the RE arms: the community EM's Newton solve
# finds gamma/eps jointly informative about each species' deviation even
# though they carry no random effect of their own, verified by inspecting
# `dim(Bf[[1]])` directly rather than assumed from the P x P shape every
# other community family's `Bf` block has. The shared gamma/eps columns
# themselves copy straight out of that same draw, UNCHANGED per species
# (there is nothing to condition -- they are the same shared draw for every
# species in that row, exactly as they are the same shared FIT for every
# species).
#
# `simulate()` reuses the family's own `.tobs_simulate_ms_dyn_occu()` handler
# (drives `cpp_simulate_ms_dyn_occu`, the SAME kernel the likelihood
# integrates against for the RE arms' HMM forward), overwriting
# `f$ms_community$coef_psi1`/`coef_p` (per-species) AND
# `f$means[gamma_names]`/`f$means[eps_names]` (the shared globals that
# handler reads directly off `object$means`, not `ms_community`).
# `loglik_many()` sums each species' own `.ms_dyn_occu_forward_ll()` (the
# exact HMM-forward marginal the fitter's `sp_ll` closure optimizes) at that
# species' psi1/p slice plus the SAME shared gamma/eps every species uses.
# ---------------------------------------------------------------------------

# Concatenate two 4D [n_sites x max_visits x n_seasons x n_species] response
# arrays on the SITE axis; visits/seasons/species are shared design, not
# something a replicate draws fresh. The `_3d_season` recipe one axis over.
.tobs_sbc_data_4d_species <- function(fit) {
  m <- fit$model
  list(cells  = m$data,
       y      = m$y,
       y_pos  = NULL,
       visits = NULL,
       graph  = NULL,
       site   = seq_len(m$n_sites))
}

.tobs_sbc_pool_4d_species <- function(obs, rep) {
  n_o <- dim(obs$y)[1L]; n_r <- dim(rep$y)[1L]
  y <- array(NA_real_, dim = c(n_o + n_r, dim(obs$y)[2L], dim(obs$y)[3L],
                               dim(obs$y)[4L]))
  y[seq_len(n_o), , , ] <- obs$y
  y[n_o + seq_len(n_r), , , ] <- rep$y
  list(cells = .tobs_sbc_rbind_cells(obs$cells, rep$cells), y = y, y_pos = NULL,
       visits = NULL, graph = NULL, site = c(obs$site, rep$site + n_o))
}

.tobs_sbc_spec_ms_dyn_occu <- function(fit, fit.control) {
  .tobs_sbc_reject_structure(fit)
  m <- fit$model
  list(model   = m,
       fit_obs = fit,
       family  = attr(fit, "tobs_family"),
       method  = fit$method,
       control = utils::modifyList(.tobs_sbc_default_control(),
                                   as.list(fit.control)),
       occ     = .tobs_sbc_recombine(m$formulas$occ, NULL),
       det     = .tobs_sbc_recombine(m$formulas$det, NULL),
       col     = .tobs_sbc_recombine(m$formulas$col, NULL),
       ext     = .tobs_sbc_recombine(m$formulas$ext, NULL),
       species = m$species_names)
}

.tobs_sbc_refit_ms_dyn_occu <- function(spec, data) {
  suppressWarnings(do.call(tobs, list(
    formula = spec$occ, data = data$cells, family = spec$family,
    detection = spec$det, colonization = spec$col, extinction = spec$ext,
    y = data$y, species = spec$species,
    method = spec$method, control = spec$control)))
}

# theta column names: `<species>_psi1_<coef>` / `<species>_p_<coef>` (one
# block of length P per species, species-major) followed by the SHARED
# `gamma_<coef>` / `eps_<coef>` columns (no species prefix -- one value per
# row, not per species).
.tobs_sbc_ms_dyn_occu_names <- function(m) {
  pi_list <- m$process_info
  psi1_nm <- pi_list[[1L]]$coef_names
  p_nm    <- pi_list[[2L]]$coef_names
  gam_nm  <- pi_list[[3L]]$coef_names
  eps_nm  <- pi_list[[4L]]$coef_names
  par <- c(paste0("psi1_", psi1_nm), paste0("p_", p_nm))
  sp_cols <- as.vector(outer(par, m$species_names,
                             function(p, sp) paste0(sp, "_", p)))
  global_cols <- c(paste0("gamma_", gam_nm), paste0("eps_", eps_nm))
  list(psi1 = psi1_nm, p = p_nm, gam = gam_nm, eps = eps_nm,
       par = par, sp_cols = sp_cols, global_cols = global_cols,
       cols = c(sp_cols, global_cols))
}

.tobs_sbc_draws_ms_dyn_occu <- function(fit, n) {
  m <- fit$model
  nm <- .tobs_sbc_ms_dyn_occu_names(m)
  S <- m$n_species; P <- length(nm$par)
  cm <- fit$ms_community
  # `Bf[[s]]` is (P + G) x P -- the cross-Hessian block between b_s and the
  # FULL (mu, global) vector, not just the RE arms (verified: dim(Bf[[1]]) =
  # 4x2 for P_psi1=P_p=P_gam=P_eps=1, i.e. (P+G) x P). The community EM's
  # Newton solve treats gamma/eps as jointly informative about each species'
  # deviation even though they carry no RE of their own, so the draw must
  # condition on the FULL mu_draws, not a psi1/p-only slice.
  re_idx <- seq_len(P)
  mu_draws <- .tobs_sbc_mvn_draws(fit$means, fit$vcov, n)   # n x (P + G)
  M <- matrix(NA_real_, n, length(nm$cols))
  colnames(M) <- nm$cols
  for (s in seq_len(S)) {
    b_hat_s <- c(cm$blup_psi1[s, ], cm$blup_p[s, ])
    b_draws <- .tobs_sbc_community_b_draws(mu_draws, fit$means, b_hat_s,
                                           cm$Bf[[s]], cm$Cinv[[s]], n)
    M[, (s - 1L) * P + seq_len(P)] <- mu_draws[, re_idx, drop = FALSE] + b_draws
  }
  M[, nm$global_cols] <- mu_draws[, nm$global_cols]
  M
}

.tobs_sbc_sim_ms_dyn_occu <- function(spec, theta, seed) {
  set.seed(seed)
  f <- spec$fit_obs
  m <- f$model
  nm <- .tobs_sbc_ms_dyn_occu_names(m)
  S <- m$n_species
  coef_psi1 <- do.call(rbind, lapply(seq_len(S), function(s)
    theta[paste0(m$species_names[s], "_psi1_", nm$psi1)]))
  coef_p <- do.call(rbind, lapply(seq_len(S), function(s)
    theta[paste0(m$species_names[s], "_p_", nm$p)]))
  dimnames(coef_psi1) <- list(m$species_names, nm$psi1)
  dimnames(coef_p)    <- list(m$species_names, nm$p)
  f$ms_community$coef_psi1 <- coef_psi1
  f$ms_community$coef_p    <- coef_p
  f$means[paste0("gamma_", nm$gam)] <- theta[paste0("gamma_", nm$gam)]
  f$means[paste0("eps_",   nm$eps)] <- theta[paste0("eps_",   nm$eps)]
  rep <- stats::simulate(f, nsim = 1L)
  obs <- spec$data_obs
  list(cells = obs$cells, y = rep, y_pos = NULL, visits = NULL, graph = NULL,
       site = obs$site)
}

.tobs_sbc_loglik_many_ms_dyn_occu <- function(fit, Theta) {
  m <- fit$model
  nm <- .tobs_sbc_ms_dyn_occu_names(m)
  S <- m$n_species
  n_sites <- m$n_sites; n_seasons <- m$n_seasons
  vapply(seq_len(nrow(Theta)), function(i) {
    th <- Theta[i, ]
    beta_gam <- th[paste0("gamma_", nm$gam)]
    beta_eps <- th[paste0("eps_",   nm$eps)]
    gamma <- as.numeric(stats::plogis(m$X_gamma %*% beta_gam))
    eps   <- as.numeric(stats::plogis(m$X_eps   %*% beta_eps))
    ll <- 0
    for (s in seq_len(S)) {
      beta_psi1 <- th[paste0(m$species_names[s], "_psi1_", nm$psi1)]
      beta_p    <- th[paste0(m$species_names[s], "_p_",    nm$p)]
      psi1 <- as.numeric(stats::plogis(m$X_psi1 %*% beta_psi1))
      p    <- as.numeric(stats::plogis(m$X_p    %*% beta_p))
      ll <- ll + .ms_dyn_occu_forward_ll(psi1, p, gamma, eps,
                                         m$y[, , , s], m$valid[, , , s],
                                         n_sites, n_seasons)
    }
    ll
  }, numeric(1))
}


# ---------------------------------------------------------------------------
# 7. The registry
#
# A family is one row. Everything not named here -- pooling, the grouping
# labels, the arm construction, the mis-scaled controls, the seed split -- is
# the shared driver's.
# ---------------------------------------------------------------------------

.TOBS_SBC_REGISTRY <- list(
  occu_cover = list(
    draws    = .tobs_sbc_draws_joint_occu_cover,
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
  count         = .tobs_sbc_simple_entry(
    "occ", resp = "y_count", y_vector = TRUE,
    replicate = .tobs_sbc_replicate_count),
  # Independent-observer double_observer(): both p1/p2 read ONE `detection`
  # formula slot (`fit$model$formulas$p`), same shape as every other simple
  # entry's site-level (state, det) pair (gcol33/tulpaObs#220). The dependent
  # (`type = "dependent"`) protocol additionally carries a fixed per-site
  # `primary` assignment the front door takes as `...`, which
  # `.tobs_sbc_refit_simple`'s `spec$extra` threads through unchanged.
  double_observer = .tobs_sbc_simple_entry(
    "lambda", "p",
    extra = function(m) if (!is.null(m$primary)) list(primary = m$primary)
                        else NULL),
  # Multi-season group (gcol33/tulpaObs#220): pools on the site axis, leaves
  # the season axis alone (custom `spec`/`data`/`pool`/`simulate`/`refit` --
  # see section 6b -- `draws` and `loglik_many` are the shared simple-family
  # ones, since `fit$draws` and `.tobs_pointwise_loglik` already work for a
  # `dynamic` fit unchanged). Constant-rate only.
  dyn_occu = list(
    spec        = .tobs_sbc_spec_dyn_occu,
    data        = .tobs_sbc_data_3d_season,
    pool        = .tobs_sbc_pool_3d_season,
    draws       = .tobs_sbc_draws_fit,
    simulate    = .tobs_sbc_sim_dyn_occu,
    refit       = .tobs_sbc_refit_dyn_occu,
    loglik_many = .tobs_sbc_loglik_many_simple),
  # dyn_abun() (gcol33/tulpaObs#220): same 3D response + site-axis pooling as
  # dyn_occu (section 6b2), but the replicate is the family's OWN simulate()
  # kernel via the shared simple-family `simulate` (no hand-written forward
  # needed -- the Dail-Madsen open N-mixture already has a working handler).
  dyn_abun = list(
    spec        = .tobs_sbc_spec_dyn_abun,
    data        = .tobs_sbc_data_3d_season,
    pool        = .tobs_sbc_pool_3d_season,
    draws       = .tobs_sbc_draws_fit,
    simulate    = .tobs_sbc_sim_simple,
    refit       = .tobs_sbc_refit_dyn_abun,
    loglik_many = .tobs_sbc_loglik_many_simple),
  # Multi-source group (gcol33/tulpaObs#220): pools on the site axis, leaves
  # the per-source axis alone (custom `spec`/`data`/`pool`/`simulate`/`refit`
  # -- see section 6c). Full-overlap only.
  int_occu = list(
    spec        = .tobs_sbc_spec_int_occu,
    data        = .tobs_sbc_data_int_occu,
    pool        = .tobs_sbc_pool_named_matrices,
    draws       = .tobs_sbc_draws_fit,
    simulate    = .tobs_sbc_sim_int_occu,
    refit       = .tobs_sbc_refit_int_occu,
    loglik_many = .tobs_sbc_loglik_many_simple),
  # Multi-response group (gcol33/tulpaObs#220): y = (yDist, yRem), two
  # response matrices that must stay row-consistent -- gdistremoval() already
  # has a simulate() handler returning list(yDist=, yRem=), reused here as the
  # replicate generator (see section 6d).
  gdistremoval = list(
    spec        = .tobs_sbc_spec_gdistremoval,
    data        = .tobs_sbc_data_gdistremoval,
    pool        = .tobs_sbc_pool_named_matrices,
    draws       = .tobs_sbc_draws_fit,
    simulate    = .tobs_sbc_sim_gdistremoval,
    refit       = .tobs_sbc_refit_gdistremoval,
    loglik_many = .tobs_sbc_loglik_many_simple),
  # Multiarm-S3 group (gcol33/tulpaObs#220): a `tobs_multiarm_fit` with two
  # independent Laplace-Gaussian blocks instead of a joint MVN draw matrix
  # (see section 6e) -- custom draws()/loglik_many() unpack/repack the two
  # blocks, everything else is a plain length-N vector response (no visits,
  # no detection formula).
  occu_categorical = list(
    spec        = .tobs_sbc_spec_occu_categorical,
    data        = .tobs_sbc_data_vector,
    pool        = .tobs_sbc_pool_vector,
    draws       = .tobs_sbc_draws_occu_categorical,
    simulate    = .tobs_sbc_sim_occu_categorical,
    refit       = .tobs_sbc_refit_occu_categorical,
    loglik_many = .tobs_sbc_loglik_many_occu_categorical),
  # distsamp_open() (gcol33/tulpaObs#220): shares dyn_abun's 3D response and
  # site-axis pooling (section 6f) -- constant-dynamics, Poisson only for v1.
  distsamp_open = list(
    spec        = .tobs_sbc_spec_distsamp_open,
    data        = .tobs_sbc_data_3d_season,
    pool        = .tobs_sbc_pool_3d_season,
    draws       = .tobs_sbc_draws_fit,
    simulate    = .tobs_sbc_sim_simple,
    refit       = .tobs_sbc_refit_distsamp_open,
    loglik_many = .tobs_sbc_loglik_many_simple),
  # occu_multi() (gcol33/tulpaObs#220): a list-of-species response, same shape
  # as int_occu()'s list-of-source one (section 6g) -- pooling is the shared
  # named-matrices route; simulate() is custom (joint multi-species state,
  # not independent per-source arms).
  occu_multi = list(
    spec        = .tobs_sbc_spec_occu_multi,
    data        = .tobs_sbc_data_occu_multi,
    pool        = .tobs_sbc_pool_named_matrices,
    draws       = .tobs_sbc_draws_fit,
    simulate    = .tobs_sbc_sim_occu_multi,
    refit       = .tobs_sbc_refit_occu_multi,
    loglik_many = .tobs_sbc_loglik_many_simple),
  # dyn_int_occu() (gcol33/tulpaObs#220): the product shape (section 6h) --
  # a named list of S per-source 3D arrays, pooled with the new
  # `.tobs_sbc_pool_named_3d`; simulate() is custom (family's own handler).
  dyn_int_occu = list(
    spec        = .tobs_sbc_spec_dyn_int_occu,
    data        = .tobs_sbc_data_dyn_int_occu,
    pool        = .tobs_sbc_pool_named_3d,
    draws       = .tobs_sbc_draws_fit,
    simulate    = .tobs_sbc_sim_dyn_int_occu,
    refit       = .tobs_sbc_refit_dyn_int_occu,
    loglik_many = .tobs_sbc_loglik_many_simple),
  # t_occu() (gcol33/tulpaObs#220, section 6i): a pg_gibbs family whose
  # fit$draws is already the real pooled posterior sample (no Gibbs-aware
  # draws() needed); shares the 3D-season pool despite its different axis
  # order (season before visits). loglik_many is a Laplace approximation to
  # the AR1 year effect's marginal (see section 6i).
  t_occu = list(
    spec        = .tobs_sbc_spec_t_occu,
    data        = .tobs_sbc_data_3d_season,
    pool        = .tobs_sbc_pool_3d_season,
    draws       = .tobs_sbc_draws_fit,
    simulate    = .tobs_sbc_sim_t_occu,
    refit       = .tobs_sbc_refit_t_occu,
    loglik_many = .tobs_sbc_loglik_many_t_occu),
  # ms_occu() (gcol33/tulpaObs#220, community group, section 6j): the "rank a
  # fixed species set" design (theta_s = mu + b_s), joint mu/b_s draw via
  # `.tobs_sbc_community_b_draws` (#226 part 1: mu and b_s are NOT independent
  # in the posterior; Cov(mu, b_s) != 0, derivation validated to machine
  # precision by block-inverting the joint (mu, b_s) arrowhead precision).
  # SPECIES-COUNT SCOPED: at S=5 (a first, small fixture) this failed
  # ACCEPTANCE hard (posterior SBC collapsed to p_unif as low as 0, several
  # quantities below 1e-3). Tried AGHQ variance-component debiasing (the fix
  # that works for `ms_occu_cover()`) and a from-scratch Vf/Cinv/Bf
  # consistency fix -- BOTH ruled out as the cause; a direct comparison
  # against `method="nuts"`'s exact joint posterior (Rhat 1.011, ESS 523, 0
  # divergences -- a trustworthy reference) showed Laplace's Vf and Cinv were
  # 3-15x too narrow AND its point estimates measurably off at S=5, not a
  # subtle attenuation but the Laplace-Gaussian approximation breaking down
  # on a genuinely non-Gaussian posterior at that few species. At S=20 (a
  # community size typical of real ecological data), the SAME plain
  # Laplace-EM -- no debiasing at all -- calibrates cleanly: 5 seeds, every
  # posterior min p_unif comfortably above 1e-3 (range 0.0017-0.032), 0
  # quantities below 1e-3 out of 81 possible across all 5 runs, no
  # reproducible failing coefficient. Registered on THAT evidence; the
  # `.SBC_REG_FIXTURES$ms_occu` fixture below is deliberately S=20, not S=5 --
  # do not shrink it without re-running this investigation. `ms_int_occu` and
  # `ms_count` (same `.tobs_community_em()` engine, same original S=5-style
  # failure) are likely fixable the identical way but each needs its OWN
  # species-count check before registering, not an assumption this transfers.
  ms_occu = list(
    spec        = .tobs_sbc_spec_ms_occu,
    data        = .tobs_sbc_data_3d_season,
    pool        = .tobs_sbc_pool_3d_season,
    draws       = .tobs_sbc_draws_ms_occu,
    simulate    = .tobs_sbc_sim_ms_occu,
    refit       = .tobs_sbc_refit_ms_occu,
    loglik_many = .tobs_sbc_loglik_many_ms_occu),
  #
  # ms_occu_cover() (gcol33/tulpaObs#220, community group, section 6j-bis):
  # the occ+p+pos analogue of ms_occu, safe to attempt because its Cinv is
  # kept consistent with Sigma even under AGHQ debiasing (#226 part 2 fix
  # specific to this family, commit `03b87ad`) -- CONTRACT-verified (refit on
  # 0-pooled data reproduces means exactly, pool/draws/simulate/loglik_many
  # all finite and correctly shaped); acceptance not yet run.
  ms_occu_cover = list(
    spec        = .tobs_sbc_spec_ms_occu_cover,
    data        = .tobs_sbc_data_ms_occu_cover,
    pool        = .tobs_sbc_pool_ms_occu_cover,
    draws       = .tobs_sbc_draws_ms_occu_cover,
    simulate    = .tobs_sbc_sim_ms_occu_cover,
    refit       = .tobs_sbc_refit_ms_occu_cover,
    loglik_many = .tobs_sbc_loglik_many_ms_occu_cover),
  #
  # cover() (gcol33/tulpaObs#220, multiarm-S3 group, section 6k): the same
  # two-independent-block shape as occu_categorical, reusing the generalized
  # vector data/pool helpers; positive = "lognormal" only for v1, dispersion
  # held fixed (no SE reported anywhere in the package for it).
  cover = list(
    spec        = .tobs_sbc_spec_cover,
    data        = .tobs_sbc_data_vector,
    pool        = .tobs_sbc_pool_vector,
    draws       = .tobs_sbc_draws_cover,
    simulate    = .tobs_sbc_sim_cover,
    refit       = .tobs_sbc_refit_cover,
    loglik_many = .tobs_sbc_loglik_many_cover),
  # ms_int_occu() (gcol33/tulpaObs#220, community group, section 6l): shared
  # ms_occu's exact failure mode -- both are `.tobs_community_em()` consumers
  # hitting a genuinely non-Gaussian posterior at a small species count, not
  # a fixable bug (verified for ms_occu against method="nuts"'s exact
  # reference; see gcol33/tulpaObs#226). SPECIES-COUNT SCOPED the same way:
  # at S=14 (matching this family's own existing recovery-test fixture, N=140,
  # 2 sources), the plain Laplace-EM (no debiasing) calibrates cleanly -- 5
  # seeds, min p_unif range 0.0013-0.052, 0 quantities below 1e-3 out of 43
  # possible across all 5 runs, no reproducible failing coefficient. Do NOT
  # shrink `.SBC_REG_FIXTURES$ms_int_occu`'s species count for speed.
  ms_int_occu = list(
    spec        = .tobs_sbc_spec_ms_int_occu,
    data        = .tobs_sbc_data_ms_int_occu,
    pool        = .tobs_sbc_pool_named_3d,
    draws       = .tobs_sbc_draws_ms_int_occu,
    simulate    = .tobs_sbc_sim_ms_int_occu,
    refit       = .tobs_sbc_refit_ms_int_occu,
    loglik_many = .tobs_sbc_loglik_many_ms_int_occu),
  #
  # occu_multiscale_cover() (gcol33/tulpaObs#220, multiarm-S3 group, section
  # 6m): the standard single-block fit shape (unlike cover()), so
  # draws/loglik_many are the shared generic ones; site is the per-plot CELL
  # index (the exchangeable unit for this family), not a plot counter.
  occu_multiscale_cover = list(
    spec        = .tobs_sbc_spec_occu_mscale_cover,
    data        = .tobs_sbc_data_occu_mscale_cover,
    pool        = .tobs_sbc_pool_occu_mscale_cover,
    draws       = .tobs_sbc_draws_fit,
    simulate    = .tobs_sbc_sim_occu_mscale_cover,
    refit       = .tobs_sbc_refit_occu_mscale_cover,
    loglik_many = .tobs_sbc_loglik_many_simple),
  # ms_count() (gcol33/tulpaObs#220, community group, section 6n): shared
  # ms_occu's exact failure mode and, like ms_occu/ms_int_occu, the actual
  # cause was species count, not a bug (see gcol33/tulpaObs#226). Multi-seed
  # (0, 1, 2) posterior SBC on a small fixture originally found
  # `sp3_mu_(Intercept)` pinned at p_unif ~9.6e-7-9.9e-7 every time (several
  # other coefficients also suspiciously low, e.g. `sp3_mu_x` ~5e-6) -- worse
  # than ms_occu/ms_int_occu's own small-fixture failures, plausibly because
  # this family has no detection arm to dilute it. At S=20 (matching
  # ms_occu's own resolved scale), the plain Laplace-EM calibrates cleanly:
  # 5 seeds, min p_unif range 0.0016-0.086, 0 quantities below 1e-3 out of 41
  # possible across all 5 runs, no reproducible failing coefficient. No
  # explicit `pool` slot needed -- `y` is a plain 2D `[site x species]`
  # matrix, so the generic default pooling applies unchanged.
  ms_count = list(
    spec        = .tobs_sbc_spec_ms_count,
    data        = .tobs_sbc_data_ms_count,
    draws       = .tobs_sbc_draws_ms_count,
    simulate    = .tobs_sbc_sim_ms_count,
    refit       = .tobs_sbc_refit_ms_count,
    loglik_many = .tobs_sbc_loglik_many_ms_count),
  # jsdm() (gcol33/tulpaObs#220, community group, section 6o): ms_count()'s
  # generic data/draws/simulate helpers reused unchanged (see the section
  # comment); only spec (family = jsdm()) and loglik_many (Bernoulli kernel)
  # are jsdm-specific.
  jsdm = list(
    spec        = .tobs_sbc_spec_jsdm,
    data        = .tobs_sbc_data_ms_count,
    draws       = .tobs_sbc_draws_ms_count,
    simulate    = .tobs_sbc_sim_ms_count,
    refit       = .tobs_sbc_refit_jsdm,
    loglik_many = .tobs_sbc_loglik_many_jsdm),
  # ms_distance() (gcol33/tulpaObs#220, community group, section 6p): the
  # two-arm (lambda/sigma) sibling of ms_occu; `data`/`pool` reuse the
  # generic 3D-season helpers unchanged (model$y is [site x bin x species],
  # structurally identical to ms_occu's [site x visit x species]).
  ms_distance = list(
    spec        = .tobs_sbc_spec_ms_distance,
    data        = .tobs_sbc_data_3d_season,
    pool        = .tobs_sbc_pool_3d_season,
    draws       = .tobs_sbc_draws_ms_distance,
    simulate    = .tobs_sbc_sim_ms_distance,
    refit       = .tobs_sbc_refit_ms_distance,
    loglik_many = .tobs_sbc_loglik_many_ms_distance),
  # ms_dyn_occu() (gcol33/tulpaObs#220, community group, section 6q): the
  # HMM-forward dynamic analogue of ms_occu, with two SHARED (non-species)
  # transition globals (gamma, eps) alongside the two per-species RE arms
  # (psi1, p). New 4D pool (model$y is [site x visit x season x species]).
  ms_dyn_occu = list(
    spec        = .tobs_sbc_spec_ms_dyn_occu,
    data        = .tobs_sbc_data_4d_species,
    pool        = .tobs_sbc_pool_4d_species,
    draws       = .tobs_sbc_draws_ms_dyn_occu,
    simulate    = .tobs_sbc_sim_ms_dyn_occu,
    refit       = .tobs_sbc_refit_ms_dyn_occu,
    loglik_many = .tobs_sbc_loglik_many_ms_dyn_occu)
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
#' unchanged. Registering it surfaced a real bug (`gcol33/tulpaObs#225`,
#' fixed): the detection intercept was silently left in standardized-covariate
#' units on every fit with an autoscaled detection covariate, independent of
#' SBC.
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
#' a named list of S per-source 3D arrays, pooled by a new
#' `.tobs_sbc_pool_named_3d` that composes both existing pooling rules.
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
#' A family is registered by adding one entry to the internal registry -- a
#' replicate generator, a refit call and optionally a joint statistic; the
#' pooling, the grouping labels, the arm construction and the controls are
#' shared.
#'
#' @return An object of class `sbc` from \pkg{tulpa}, carrying the per-quantity
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
#'             sigma.grid   = exp(seq(log(0.15), log(2.0), length.out = 9)),
#'             phi.grid.pos = exp(seq(log(0.20), log(0.90), length.out = 7)))
#' fit <- tobs(~ occ_cov1 + icar(graph = adj),
#'             data = cbind(data.frame(site_id = seq_len(N)), sim$data),
#'             family = occu_cover("lognormal"),
#'             detection = ~ det_cov1,
#'             positive = ~ pos_cov1 + copy(spatial()),
#'             y = od$y, y_pos = y_pos, visits = od$det.covs,
#'             method = "nested_laplace", control = ctl)
#' sbc(fit, n.sim = 20L, controls = "narrow", fit.control = ctl)
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
  # newly registered one (gcol33/tulpaObs#220, found while adding
  # double_observer to the registry).
  res <- tulpa::sbc(object = "posterior", model = model,
                    n_sim = as.integer(n.sim), quantities = quantities,
                    level = level, seed = as.integer(seed), control = control)
  res$tobs_family <- attr(model, "family")
  res$fixed_quantities <- attr(model, "fixed")
  res
}
