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
  list(cells  = rbind(obs$cells, rep$cells),
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
  ctl <- utils::modifyList(list(verbose = FALSE, progress = FALSE),
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
       control = utils::modifyList(list(verbose = FALSE, progress = FALSE),
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
       control = utils::modifyList(list(verbose = FALSE, progress = FALSE),
                                   as.list(fit.control)),
       occ = .tobs_sbc_recombine(m$formulas$occ, NULL),
       det = .tobs_sbc_recombine(m$formulas$det, NULL),
       col = .tobs_sbc_recombine(m$formulas$col, NULL),
       ext = .tobs_sbc_recombine(m$formulas$ext, NULL))
}

.tobs_sbc_data_dyn_occu <- function(fit) {
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
.tobs_sbc_pool_dyn_occu <- function(obs, rep) {
  n_o <- dim(obs$y)[1L]; n_r <- dim(rep$y)[1L]
  y <- array(NA_real_, dim = c(n_o + n_r, dim(obs$y)[2L], dim(obs$y)[3L]))
  y[seq_len(n_o), , ] <- obs$y
  y[n_o + seq_len(n_r), , ] <- rep$y
  list(cells = rbind(obs$cells, rep$cells), y = y, y_pos = NULL, visits = NULL,
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
       control   = utils::modifyList(list(verbose = FALSE, progress = FALSE),
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
  list(cells = rbind(obs$cells, rep$cells), y = y, y_pos = NULL, visits = NULL,
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
       control = utils::modifyList(list(verbose = FALSE, progress = FALSE),
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
    data        = .tobs_sbc_data_dyn_occu,
    pool        = .tobs_sbc_pool_dyn_occu,
    draws       = .tobs_sbc_draws_fit,
    simulate    = .tobs_sbc_sim_dyn_occu,
    refit       = .tobs_sbc_refit_dyn_occu,
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
    loglik_many = .tobs_sbc_loglik_many_simple)
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
