# =============================================================================
# occu_cover_spatial.R - v2 nested-Laplace spatial path for occu_cover().
#
# Adds a cell-level latent ICAR field z[1..n_cells] (sum-to-zero) shared
# across the psi and cover arms, mirroring Michael Glaser's mod.joint
# (`copy = "cell.occ"` linking the besag field on cell.occ to cell.ab).
#
# Augmented linear predictors:
#
#   eta_psi_i   = X_psi[i, ] %*% beta_psi + z[i]
#   eta_p_ij    = X_p[i, ] %*% beta_p_site + X_p_visit[ij, ] %*% beta_p_visit
#   eta_pos_ij  = X_pos[i, ] %*% beta_pos_site + X_pos_visit[ij, ] %*% beta_pos_visit
#                 + alpha * z[i]
#
# This is `cover()`'s joint parameterisation, not "sigma * z + alpha * sigma *
# z": the latter confounds `alpha` with `sigma` when detection is weak.
# Here `z` carries its own marginal variance sigma^2 via the prior precision,
# and `alpha` is the cover-arm scaling alone (matches INLA `copy=` semantics).
#
# z's prior is the ICAR (intrinsic CAR / Besag) penalty (1/(2 sigma^2)) * z' Q z with
# Q = D - W (degree minus adjacency, the unscaled graph Laplacian). The
# field is rank-deficient (kernel = constant vector), so a soft sum-to-zero
# penalty 0.5 * kappa * (sum(z))^2 is added to identify the psi intercept
# separately from the field mean.
#
# Fit:
#   joint Laplace MAP via optim(BFGS) on
#     par = c(beta_psi, beta_p, beta_pos, log_disp, z, alpha, log_sigma)
#   SEs from the inverse observed-Fisher Hessian (full joint covariance).
#
# v2.0 limitations (intentional, will be addressed in v3):
#   - ICAR only (rho = 1). BYM2 with free rho mixing requires another
#     latent block; for now bym2() in the formula is read as ICAR with a
#     warning the first time. Most of Michael's `mod.joint` runs use
#     besag (= ICAR), so the head-to-head still holds.
#   - Joint optim grows linearly in n_cells. For n_cells > ~2000 the
#     dense Hessian becomes the bottleneck; sparse-Q + sparse linear
#     algebra is v3.
#   - No outer-grid integration of (sigma, alpha); they're point-estimated
#     at the joint MAP with their SE read off the Hessian. The full
#     nested-Laplace integration of these hyperpriors is v3.
# =============================================================================


# ---------------------------------------------------------------------------
# Build the ICAR precision matrix Q = D - W from an adjacency.
# ---------------------------------------------------------------------------
.occu_cover_icar_Q <- function(adj) {
  n <- nrow(adj)
  D <- rowSums(adj)
  if (any(D == 0L)) {
    isolated <- which(D == 0L)
    stop(sprintf("ICAR graph has %d isolated node(s) (no neighbours): %s. ",
                 length(isolated),
                 paste(utils::head(isolated, 5L), collapse = ", ")),
         "Drop them or connect them before fitting.", call. = FALSE)
  }
  Q <- -as.matrix(adj)
  diag(Q) <- D
  Q
}


# ---------------------------------------------------------------------------
# Sorbye-Rue scaling factor: geometric mean of diag(Q^+) under the sum-to-zero
# constraint. Multiplying Q by this factor gives a precision whose
# generalized-inverse diagonal has geometric mean 1, so `sigma` in the linear
# predictor `sigma * z` is the geo-mean marginal SD of the field.
# Matches INLA's `scale.model = TRUE` (Sorbye & Rue 2014), the convention
# Michael's `mod.joint` uses on the besag field.
# ---------------------------------------------------------------------------
.occu_cover_icar_scale <- function(adj) {
  n <- nrow(adj)
  Q <- .occu_cover_icar_Q(adj)
  eig <- eigen(Q, symmetric = TRUE)
  pos <- eig$values > 1e-10
  if (!any(pos)) {
    stop("ICAR graph has no positive eigenvalues; check connectivity.",
         call. = FALSE)
  }
  V_p <- eig$vectors[, pos, drop = FALSE]
  Qinv_diag <- rowSums(V_p^2 /
                          matrix(eig$values[pos], n, sum(pos), byrow = TRUE))
  exp(mean(log(Qinv_diag)))
}


# ---------------------------------------------------------------------------
# Pull the coupled spatial field(s) from the psi formula. Returns NULL when no
# spatial term is present, otherwise a list with `fe` (the fixed-effects psi
# formula) and `fields` (the ordered field specs: the one unweighted intercept
# field first, then any weighted spatially-varying-coefficient fields).
#
# A weighted areal term (`icar(graph = adj, weight = col)`) is a second coupled
# field -- a spatially-varying coefficient on `col` sharing the same areal
# graph -- the formula-DSL spelling of the trend field that `control$trend`
# also produces (gcol33/tulpaObs#15). The joint_coupled engine couples N such
# fields; the legacy single-field v2/v3 engines take only the intercept field.
# ---------------------------------------------------------------------------
.occu_cover_spatial_fields <- function(formula, data, arm_fields = list()) {
  bind <- .tobs_bind_formulas(list(psi = formula), data)
  if (length(bind$terms) == 0L && length(arm_fields) == 0L) return(NULL)
  # `.tobs_bind_formulas` returns terms wrapped in `list(spec = ..., process = ...)`.
  spatial <- Filter(function(t) inherits(t$spec, "tobs_spatial"), bind$terms)
  # Placed arm fields (detection / positive): evaluate each lifted field call and
  # tag its spec with the arm, in the occurrence formula's environment (where the
  # occurrence field's own graph symbol resolves), then route alongside the
  # occurrence formula's own fields.
  if (length(arm_fields)) {
    extra <- lapply(arm_fields, function(af) {
      spec <- .tobs_eval_arm_field(af$call, af$arm, data, environment(formula))
      # Only the bar form carries an arm-specific field; a plain areal
      # constructor placed in an arm formula has no such spelling.
      if (!isTRUE(spec$is_bar)) {
        stop(sprintf(paste0(
          "occu_cover(): a field placed in the %s formula must use the bar form ",
          "spatial(~ 1 + w || cell, graph = adj); got %s()."),
          af$arm, spec$type %||% class(spec)[1L]), call. = FALSE)
      }
      list(spec = spec)
    })
    spatial <- c(spatial, extra)
  }
  if (length(spatial) == 0L) return(NULL)
  # A varying-coefficient bar (`spatial(~ 1 + w || node, graph = adj)`,
  # gcol33/tulpaObs#61) desugars in place to the intercept field + per-covariate
  # trend fields, the same pair the two-term form produces -- INDEPENDENT areal
  # fields, each with its own precision. A correlated bar (single `|`, the
  # free-Sigma separable MCAR of gcol33/tulpaObs#64) declares the SAME per-field
  # design, but the intercept + coefficient fields share a free cross-covariance
  # Sigma on the occupancy arm and are copied onto the cover arm together with one
  # amplitude alpha (gcol33/tulpaObs#63). It is fitted as one coupled MCAR block
  # rather than independent fields, so it must be the only spatial term (one field
  # structure per fit); `correlated` flags it for the joint-coupled fitter.
  specs <- list()
  correlated <- FALSE
  # Arm-specific INDEPENDENT fields, keyed by internal arm slot: "pos" (cover) and
  # "p" (detection). A single-arm `to =` bar (placement in the positive / detection
  # formula, lifted here) becomes a separate latent field on that arm ALONE, with
  # its own precision and NO cross-arm copy of the occupancy field's alpha. This is
  # the opt-in for arm-structured trends the alpha copy cannot express (the copy
  # collapses to a global slope when the shapes differ, e.g. flattening
  # delta_cover_cond). Distinct from the shared `||` desugar (which copies a
  # psi-anchored field onto the cover arm with alpha).
  armspec  <- list()
  arm_slot <- c(positive = "pos", detection = "p")
  for (t in spatial) {
    if (isTRUE(t$spec$is_bar)) {
      to <- t$spec$to %||% .tobs_cover_arms
      if (isTRUE(t$spec$correlated)) {
        if (correlated || length(specs) > 0L || length(spatial) > 1L) {
          stop(paste0(
            "occu_cover(): a correlated spatial bar (`|`) must be the only ",
            "spatial term in the psi formula (one MCAR field structure per ",
            "fit). Drop the other areal term(s), or use the INDEPENDENT ",
            "spelling `||` to combine separate per-coefficient fields."),
            call. = FALSE)
        }
        specs <- .tobs_expand_spatial_bar(t$spec, data)
        correlated <- TRUE
      } else if (length(to) == 1L) {
        if (identical(to, "presence")) {
          stop(paste0(
            "occu_cover(): there is no separate presence arm; the occupancy ",
            "field is the occurrence spatial() term. Write the field in the ",
            "`occurrence` formula, or place it in `positive` (cover) / ",
            "`detection`."), call. = FALSE)
        }
        if (!to %in% names(arm_slot)) {
          stop(sprintf(paste0(
            "occu_cover(): an arm-specific field belongs in the positive (cover) ",
            "or detection formula; got arm \"%s\"."), to), call. = FALSE)
        }
        slot <- arm_slot[[to]]
        if (!is.null(armspec[[slot]])) {
          stop(sprintf(paste0(
            "occu_cover(): at most one arm-specific field per arm (%s); ",
            "combine the coefficient fields into one bar in that arm's formula, ",
            "e.g. %s = ~ x + spatial(~ 1 + w || cell, graph = adj)."),
            to, to), call. = FALSE)
        }
        armspec[[slot]] <- .tobs_armspecific_bar_fields(t$spec, data)
      } else {
        specs <- c(specs, .cover_desugar_spatial_bar(t$spec, data))
      }
    } else {
      specs[[length(specs) + 1L]] <- t$spec
    }
  }

  # Soft guard (gcol33/tulpaObs#62): a bare `| / ||` RE bar whose grouping factor
  # is also an areal term's graph-node group_var is fitted as an IID random effect,
  # not a spatial field (the engine-bar-idiom papercut). Informs (message) without
  # rejecting -- RE bars are legitimate -- and is silent when the bar's factor is
  # unrelated to any spatial term.
  .tobs_cover_bar_re_guard(formula, specs)

  bad <- Filter(function(s) !s$type %in% c("icar", "bym2"), specs)
  if (length(bad)) {
    stop(sprintf(
      "occu_cover() spatial path supports icar() or bym2() in the psi formula; got %s().",
      bad[[1L]]$type), call. = FALSE)
  }
  if (any(vapply(specs, function(s) identical(s$type, "bym2"), logical(1)))) {
    warning("occu_cover() reads bym2() as ICAR (rho fixed to 1); ",
            "BYM2 with free rho mixing is the v3 escape hatch.", call. = FALSE)
  }

  # Exactly one unweighted (intercept) field is the shared base field; the rest
  # are weighted SVC fields. The cell-coupling model has one node per cell, so
  # every coupled field shares the base graph (only the weight column differs).
  weighted <- vapply(specs, function(s) !is.null(s$weight), logical(1))
  base <- specs[!weighted]
  if (length(base) == 0L) {
    stop("occu_cover() spatial requires one unweighted intercept field ",
         "(e.g. icar(graph = adj)); only weighted SVC field(s) were given.",
         call. = FALSE)
  }
  if (length(base) > 1L) {
    stop("occu_cover() spatial supports exactly one unweighted intercept ",
         "field; got ", length(base), ". Additional coupled fields must be ",
         "weighted SVC terms, e.g. icar(graph = adj, weight = year).",
         call. = FALSE)
  }
  base_graph <- base[[1L]]$graph
  for (s in specs[weighted]) {
    if (!identical(dim(s$graph), dim(base_graph)) ||
        !all(s$graph == base_graph)) {
      stop("occu_cover() coupled fields must share the same areal graph as ",
           "the intercept field (same nodes / adjacency).", call. = FALSE)
    }
  }

  # Optional group_var maps each site (occupancy unit) to a field node, so the
  # site count can exceed the node count (e.g. site = cell-year sharing one cell
  # field). All coupled fields must name the same group_var (or none).
  gvs <- unique(Filter(Negate(is.null), lapply(specs, function(s) s$group_var)))
  if (length(gvs) > 1L) {
    stop("occu_cover() coupled fields must share a single group_var (or none).",
         call. = FALSE)
  }
  group_var <- if (length(gvs) == 1L) gvs[[1L]] else NULL

  # An optional per-group random INTERCEPT on the occupancy arm, layered on the
  # shared field (gcol33/tulpaObs#56, the consumer of tulpa#86's field + per-group
  # RE composition in the joint cell-coupling engine). It joins the joint fit as a
  # single `iid` prior block whose per-cell offset rides the occupancy predictor;
  # its variance integrates on the outer grid alongside the field sigma / alpha.
  # Scope: one random-intercept term -- a scalar per group -- maps onto the one
  # iid block. A random slope or a correlated multi-coefficient block has no
  # scalar-per-group form here and errors (the non-spatial cover-hurdle EM is the
  # route for richer RE).
  re_terms <- Filter(function(t) inherits(t$spec, "tobs_re"), bind$terms)
  re_spec  <- NULL
  if (length(re_terms) > 0L) {
    if (length(re_terms) > 1L) {
      stop("occu_cover() spatial + RE supports a single random-intercept term ",
           "on the occupancy formula; got ", length(re_terms), ".", call. = FALSE)
    }
    rs <- re_terms[[1L]]$spec
    if (!identical(rs$type, "intercept")) {
      stop("occu_cover() spatial + RE supports a random INTERCEPT only ",
           "(e.g. (1 | group) / re(group)); a random slope or correlated block ",
           "is not wired on the shared-field joint engine.", call. = FALSE)
    }
    re_spec <- list(group_idx = as.integer(rs$group_idx),
                    n_groups  = as.integer(rs$n_groups))
  }

  # Correlated (`|`) MCAR field requirements (gcol33/tulpaObs#63): at least one
  # coefficient beyond the intercept (a single field has no cross-covariance),
  # intrinsic CAR only, and no per-group RE (the MCAR block already spans the
  # full coupled field structure; a layered iid block is not wired with it yet).
  if (correlated) {
    if (length(specs) < 2L) {
      stop(paste0(
        "occu_cover(): a correlated spatial bar (`|`) needs at least one ",
        "coefficient beyond the intercept (e.g. spatial(~ 1 + x | cell, ",
        "graph = adj)); a single field has no cross-covariance to estimate. ",
        "Use icar()/`||` for an uncorrelated field."), call. = FALSE)
    }
    if (!all(vapply(specs, function(s) identical(s$type, "icar"), logical(1)))) {
      stop(paste0(
        "occu_cover(): a correlated spatial bar (`|`) uses the intrinsic CAR ",
        "(icar); bym2 is not supported for the MCAR field."), call. = FALSE)
    }
    if (!is.null(re_spec)) {
      stop(paste0(
        "occu_cover(): a correlated spatial bar (`|`) does not compose with a ",
        "per-group occupancy random effect; fit one or the other."),
        call. = FALSE)
    }
  }

  # Arm-specific fields (gcol33/tulpaObs#110, extended to the detection arm): each
  # is a separate, non-copied block on ONE arm (cover or detection); it composes
  # with the shared occupancy field (which still drives psi and, via the alpha
  # copy, delta_cover_exp) but NOT with the correlated `|` MCAR field (that already
  # spans the whole coupled structure with its own copy). Like the occu_cover
  # shared field, an arm-specific field is fitted as ICAR (rho fixed to 1); bym2 /
  # car on the bar is read as ICAR. Every arm-specific field shares the occupancy
  # field's areal graph (one node set).
  arm_label <- c(pos = "positive", p = "detection")
  for (slot in names(armspec)) {
    af <- armspec[[slot]]
    if (correlated) {
      stop(sprintf(paste0(
        "occu_cover(): an arm-specific field (to = \"%s\") does not compose with ",
        "a correlated `|` MCAR field; use one spatial structure."),
        arm_label[[slot]]), call. = FALSE)
    }
    if (!identical(af$type, "icar")) {
      warning(sprintf(paste0(
        "occu_cover() reads the arm-specific %s field as ICAR (rho fixed to 1); ",
        "bym2/car with free mixing on that arm is not yet wired."),
        arm_label[[slot]]), call. = FALSE)
      af$type <- "icar"
    }
    if (!identical(dim(af$graph), dim(base_graph)) ||
        !all(af$graph == base_graph)) {
      stop(sprintf(paste0(
        "occu_cover(): the arm-specific %s field (to = \"%s\") must share the ",
        "same areal graph as the occupancy field (same nodes / adjacency)."),
        arm_label[[slot]], arm_label[[slot]]), call. = FALSE)
    }
    armspec[[slot]] <- af
  }

  list(fe = bind$fe$psi, fields = c(base, specs[weighted]),
       armspec = armspec, pos_armspec = armspec[["pos"]],
       group_var = group_var, re = re_spec, correlated = correlated)
}


# ---------------------------------------------------------------------------
# NLP for the spatial v2 fit.
#
# par layout (offsets recorded once in the fitter):
#   beta_psi      [1, p_psi]
#   beta_p        [p_psi + 1, p_psi + p_p]
#   beta_pos      [p_psi + p_p + 1, p_psi + p_p + p_pos]
#   log_disp      [.. + 1]
#   z             [.. + (1 .. n_cells)]
#   alpha         [.. + 1]
#   log_sigma     [.. + 1]
# ---------------------------------------------------------------------------
.tobs_occu_cover_spatial_nlp <- function(par, model, Q, scale_q, pmean, pprec,
                                          kappa_sum = 1e4) {
  cl <- .tobs_clamp_eta

  pi_list <- model$process_info
  p_psi   <- pi_list[[1L]]$p
  p_p     <- pi_list[[2L]]$p
  p_pos   <- pi_list[[3L]]$p
  n_cells <- model$n_sites
  max_visits <- model$max_visits

  off <- 0L
  beta_psi <- par[off + seq_len(p_psi)]; off <- off + p_psi
  beta_p   <- par[off + seq_len(p_p)];   off <- off + p_p
  beta_pos <- par[off + seq_len(p_pos)]; off <- off + p_pos
  log_disp <- par[off + 1L];             off <- off + 1L
  z        <- par[off + seq_len(n_cells)]; off <- off + n_cells
  alpha    <- par[off + 1L];             off <- off + 1L
  log_sigma <- par[off + 1L]
  sigma    <- exp(log_sigma)

  # Psi linear predictor: cell-level + field (z carries its own marginal SD
  # sigma via the prior precision; no sigma multiplier here).
  eta_psi <- as.numeric(model$X_occ %*% beta_psi) + z
  psi     <- stats::plogis(cl(eta_psi))

  # Detection: site-level + visit-level (no field — detection is observation-
  # process, the shared field belongs to the latent state and the cover arm).
  bp_site  <- beta_p[seq_len(ncol(model$X_det_site))]
  bp_visit <- if (!is.null(model$X_det_visit)) {
    beta_p[ncol(model$X_det_site) + seq_len(ncol(model$X_det_visit))]
  } else {
    numeric(0)
  }
  eta_p_site <- as.numeric(model$X_det_site %*% bp_site)
  p_mat <- matrix(eta_p_site, n_cells, max_visits)
  if (length(bp_visit)) {
    eta_p_visit <- as.numeric(model$X_det_visit %*% bp_visit)
    p_mat <- p_mat + matrix(eta_p_visit, n_cells, max_visits, byrow = TRUE)
  }
  p_mat <- stats::plogis(cl(p_mat))

  # Cover-arm linear predictor (site-level + visit-level + alpha * sigma * z).
  bpos_site  <- beta_pos[seq_len(ncol(model$X_pos_site))]
  bpos_visit <- if (!is.null(model$X_pos_visit)) {
    beta_pos[ncol(model$X_pos_site) + seq_len(ncol(model$X_pos_visit))]
  } else {
    numeric(0)
  }
  eta_pos_site <- as.numeric(model$X_pos_site %*% bpos_site)
  ep_mat <- matrix(eta_pos_site, n_cells, max_visits)
  if (length(bpos_visit)) {
    eta_pos_visit <- as.numeric(model$X_pos_visit %*% bpos_visit)
    ep_mat <- ep_mat + matrix(eta_pos_visit, n_cells, max_visits, byrow = TRUE)
  }
  ep_mat <- ep_mat + matrix(alpha * z, n_cells, max_visits)

  valid <- model$valid
  y     <- model$y
  y_pos <- model$y_pos

  log_p   <- ifelse(valid, log(p_mat),     0)
  log_1mp <- ifelse(valid, log(1 - p_mat), 0)

  pos_mask <- valid & (y == 1L)
  log_f_pos <- matrix(0, n_cells, max_visits)
  if (identical(model$positive, "beta")) {
    phi_d <- exp(log_disp)
    mu_pos <- stats::plogis(cl(ep_mat))
    a <- mu_pos * phi_d
    b <- (1 - mu_pos) * phi_d
    dens <- lgamma(phi_d) - lgamma(a) - lgamma(b) +
            (a - 1) * log(y_pos) + (b - 1) * log(1 - y_pos)
    log_f_pos[pos_mask] <- dens[pos_mask]
  } else {
    sigma_pos <- exp(log_disp)
    dens <- -log(y_pos) - log(sigma_pos) - 0.5 * log(2 * pi) -
            0.5 * ((log(y_pos) - ep_mat) / sigma_pos)^2
    log_f_pos[pos_mask] <- dens[pos_mask]
  }

  log_h <- ifelse(valid,
                  ifelse(y == 1L, log_p + log_f_pos, log_1mp),
                  0)

  any_det <- rowSums(y * valid, na.rm = FALSE) > 0
  log_psi   <- log(pmax(psi, 1e-300))
  log_1mpsi <- log(pmax(1 - psi, 1e-300))

  det_ll <- log_psi + rowSums(log_h)
  ln_a <- log_psi   + rowSums(log_1mp)
  ln_b <- log_1mpsi
  nodet_ll <- .tobs_logsumexp2(ln_a, ln_b)

  ll <- sum(ifelse(any_det, det_ll, nodet_ll))

  # Priors and penalties.
  # - Gaussian prior on betas (and on alpha, log_sigma via pprec).
  # - ICAR prior on z: 0.5 * z' Q z.
  # - Soft sum-to-zero on z: 0.5 * kappa_sum * sum(z)^2.
  beta_penalty <- 0.5 * sum(pprec * (par - pmean)^2)
  # Sorbye-Rue scaled ICAR with marginal variance sigma^2:
  #   z ~ N(0, sigma^2 * (scale_q * Q)^{-1})
  #   -log p(z) = 0.5 * (scale_q / sigma^2) * z' Q z + (n_eff/2) * log(sigma^2)
  # The log-determinant term is constant in z but needed for sigma's score.
  # n_eff = rank of Q under sum-to-zero = n_cells - 1.
  inv_sig2 <- exp(-2 * log_sigma)
  z_prior  <- 0.5 * inv_sig2 * scale_q *
              as.numeric(crossprod(z, Q %*% z)) +
              (n_cells - 1) * log_sigma
  z_sumzero <- 0.5 * kappa_sum * sum(z)^2

  -ll + beta_penalty + z_prior + z_sumzero
}


# ---------------------------------------------------------------------------
# Spatial fitter.
# ---------------------------------------------------------------------------
.tobs_fit_occu_cover_spatial <- function(model, adj,
                                          priors    = NULL,
                                          max.iter  = 300L,
                                          tol       = 1e-6,
                                          verbose   = TRUE,
                                          sigma.beta = 5,
                                          ...) {
  pi_list <- model$process_info
  p_psi   <- pi_list[[1L]]$p
  p_p     <- pi_list[[2L]]$p
  p_pos   <- pi_list[[3L]]$p
  n_cells <- model$n_sites

  if (nrow(adj) != n_cells) {
    stop(sprintf("Spatial graph has %d nodes but data has %d cells. ",
                 nrow(adj), n_cells),
         "Pass an adjacency matching the cell grid (one node per cell).",
         call. = FALSE)
  }

  Q       <- .occu_cover_icar_Q(adj)
  scale_q <- .occu_cover_icar_scale(adj)
  n_par <- p_psi + p_p + p_pos + 1L + n_cells + 2L

  par_names <- c(
    paste0("psi_", pi_list[[1L]]$coef_names),
    paste0("p_",   pi_list[[2L]]$coef_names),
    paste0("pos_", pi_list[[3L]]$coef_names),
    if (identical(model$positive, "beta")) "log_phi" else "log_sigma_pos",
    sprintf("z[%d]", seq_len(n_cells)),
    "alpha", "log_sigma"
  )

  # Warm starts.
  start <- numeric(n_par)
  any_det <- rowSums(model$y * model$valid) > 0
  det_rate <- max(mean(any_det), 1e-3)
  start[1L] <- stats::qlogis(min(max(det_rate, 1e-3), 1 - 1e-3))

  pos_vals <- model$y_pos[model$valid & model$y == 1L]
  disp_idx <- p_psi + p_p + p_pos + 1L
  pos_int_idx <- p_psi + p_p + 1L
  if (length(pos_vals) > 0L) {
    if (identical(model$positive, "beta")) {
      start[pos_int_idx] <- stats::qlogis(min(max(mean(pos_vals), 1e-3), 1 - 1e-3))
      start[disp_idx]    <- log(10)
    } else {
      start[pos_int_idx] <- mean(log(pos_vals))
      start[disp_idx]    <- log(stats::sd(log(pos_vals)) + 0.1)
    }
  } else {
    start[disp_idx] <- if (identical(model$positive, "beta")) log(10) else log(0.4)
  }

  # Warm-start z from each cell's per-cell mean cover residual.
  # If cell i has positive observations with cover values, log(mean(cover))
  # minus the population mean gives a rough estimate of the cell deviation
  # on the cover-arm scale; scale to roughly the prior amplitude. Starting
  # z at 0 traps the optim in (small z, small sigma, large alpha) ridge.
  z_idx_init <- p_psi + p_p + p_pos + 1L + seq_len(n_cells)
  cell_pos_mean <- rep(NA_real_, n_cells)
  for (i in seq_len(n_cells)) {
    vals <- model$y_pos[i, model$valid[i, ] & model$y[i, ] == 1L]
    if (length(vals) > 0L) {
      cell_pos_mean[i] <- if (identical(model$positive, "beta"))
                            stats::qlogis(min(max(mean(vals), 1e-3), 1 - 1e-3))
                          else
                            mean(log(vals))
    }
  }
  z_init <- cell_pos_mean - mean(cell_pos_mean, na.rm = TRUE)
  z_init[is.na(z_init)] <- 0
  # Soft-shrink toward zero and rescale to a modest amplitude (sigma ~ 0.5).
  if (stats::sd(z_init) > 0) z_init <- 0.5 * z_init / stats::sd(z_init)
  start[z_idx_init] <- z_init

  start[n_par - 1L] <- 1.0   # alpha — same sign, same scale anchor
  start[n_par]      <- 0     # log_sigma — start at sigma = 1

  # Gaussian-prior precision aligned with par. Penalize only fixed-effect
  # betas (not z, not alpha, not log_sigma — those carry their own priors).
  pmean <- numeric(n_par)
  pprec <- numeric(n_par)
  if (isTRUE(is.null(priors)) || !isFALSE(priors)) {
    beta_idx <- c(seq_len(p_psi),
                  p_psi + seq_len(p_p),
                  p_psi + p_p + seq_len(p_pos))
    pprec[beta_idx] <- 1 / (sigma.beta^2)
  }
  # Anchors on the hyperparameters to break the (z, alpha, sigma) ridge.
  # alpha ~ N(0, 1^2) — symmetric "no-effect" prior; lets alpha be either
  # sign / scale without nudging toward 1. log_sigma ~ N(0, 0.6^2) — sigma
  # in roughly [0.3, 3.3], wide enough to not dominate when data has signal.
  # The joint Laplace still slides this ridge when the cover arm dominates
  # in info; v3 (proper nested-Laplace with z profiled out) is the structural
  # fix.
  pmean[n_par - 1L] <- 0
  pprec[n_par - 1L] <- 1 / 1
  pmean[n_par]      <- 0
  pprec[n_par]      <- 1 / (0.6^2)

  # max_calls ceiling for the v2 outer BFGS: n_par + 3 calls per iter
  # (FD gradient + line search), * max.iter.
  calls_per_iter_est <- n_par + 3L
  max_calls_est <- as.integer(max.iter) * calls_per_iter_est
  report <- if (isTRUE(verbose) || is.null(verbose) || verbose != FALSE)
              .tobs_progress_reporter("occu_cover v2", throttle = 5,
                                       max_calls = max_calls_est)
            else
              function(...) invisible(NULL)

  wrapped_nlp <- function(par) {
    val <- .tobs_occu_cover_spatial_nlp(par, model, Q, scale_q, pmean, pprec)
    report(val)
    val
  }

  opt <- stats::optim(start, wrapped_nlp,
                       method = "BFGS", hessian = TRUE,
                       control = list(maxit = max.iter, reltol = tol,
                                      trace = if (isTRUE(verbose)) 1L else 0L))

  # Post-process: re-center z to sum to zero (the soft penalty makes
  # mean(z) ~ 0 but not exactly), absorbing the residual into the two
  # intercepts. Linear predictors are psi = beta + z and pos = beta + alpha * z,
  # so the offsets are +z_mean on psi[1] and +alpha * z_mean on pos[1].
  z_idx     <- p_psi + p_p + p_pos + 1L + seq_len(n_cells)
  alpha_idx <- n_par - 1L
  z_mean    <- mean(opt$par[z_idx])
  alpha_hat <- opt$par[alpha_idx]
  opt$par[z_idx] <- opt$par[z_idx] - z_mean
  opt$par[1L]                  <- opt$par[1L]                  + z_mean
  opt$par[p_psi + p_p + 1L]    <- opt$par[p_psi + p_p + 1L]    + alpha_hat * z_mean

  V <- tryCatch(solve(opt$hessian), error = function(e) NULL)
  if (is.null(V)) {
    warning("occu_cover spatial: Hessian not invertible; SEs unreliable.",
            call. = FALSE)
    V <- matrix(NA_real_, n_par, n_par)
  }
  se <- sqrt(pmax(diag(V), 0))

  means <- opt$par
  names(means) <- par_names
  names(se)    <- par_names
  dimnames(V)  <- list(par_names, par_names)

  # Split betas from field for cleaner downstream summaries.
  beta_names <- par_names[-z_idx]
  field_names <- par_names[z_idx]
  beta_idx_all <- setdiff(seq_len(n_par), z_idx)

  # Pseudo-draws on the FIXED-EFFECT / hyperparameter block only (not z;
  # the field draws are large and the user typically wants the cell-level
  # marginal summaries computed below).
  n_draws <- 1000L
  draws <- .occu_cover_rmvn(n_draws, means[beta_idx_all], V[beta_idx_all, beta_idx_all, drop = FALSE])
  colnames(draws) <- beta_names

  # Per-cell field summary (posterior mean + 95% Wald CI).
  z_mean_post <- means[z_idx]
  z_sd_post   <- se[z_idx]
  field_table <- data.frame(
    cell      = seq_len(n_cells),
    z_mean    = z_mean_post,
    z_sd      = z_sd_post,
    z_lower   = z_mean_post - 1.96 * z_sd_post,
    z_upper   = z_mean_post + 1.96 * z_sd_post
  )

  structure(c(list(
    draws        = draws,
    means        = means,
    sds          = se,
    vcov         = V,
    n_samples    = n_draws,
    n_params     = n_par,
    log_prob     = rep(-opt$value, n_draws),
    log_lik      = -opt$value,
    N            = sum(model$valid)),
    .tobs_na_nuts_diagnostics(n_draws),
    list(
    col_names    = beta_names,
    param_names  = beta_names,
    process_info = pi_list,
    model        = model,
    spatial      = list(type = "icar", graph = adj,
                        sigma_mean = exp(means["log_sigma"]),
                        alpha_mean = means["alpha"]),
    spatial_field = z_mean_post,
    field_table  = field_table,
    method       = "nested_laplace",
    positive     = model$positive,
    convergence  = list(converged = opt$convergence == 0L,
                        n_iter    = opt$counts[1L])
  )), class = c("tobs_fit", "tulpa_fit"))
}
