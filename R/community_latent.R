# =============================================================================
# community_latent.R - term resolution for the community latent structure
#
# A community model carrying a shared areal field and/or latent factors on one
# arm has the linear predictor
#
#   eta_{s,i} = X_i . (mu + b_s) + sum_k W[i,k] F[u(i),k]
#                                + sum_q lambda_{s,q} zeta_{q,i}
#
# with b_s ~ N(0, Sigma) the per-species coefficient deviations, F the areal
# field(s) (an all-ones weight column for the intercept field, covariate values
# for a varying-coefficient field), and (zeta, lambda) the per-site latent
# factors and their per-species loadings.
#
# The fit is tulpa's block-coordinate ascent `tulpa::tulpa_latent_ascent()`:
# the family's community Laplace-EM with the latent structure as a per-cell
# offset, alternated with a field Newton solve and the marginal-likelihood
# loading EM. Every family reaches it through ONE oracle,
#
#   working(eta_mat) -> list(score = [n_sites x n_species],
#                            curv  = [n_sites x n_species])
#
# the score and curvature of that family's marginal log-likelihood with respect
# to an additive offset on the structured arm's linear predictor, plus
# `data_ll(eta)`, `ll_cell(eta, idx)` (factors) and optionally `cov_curv(eta)`.
# Poisson gives (y - mu, mu); the occupancy two-state marginal gives its own
# pair; a Bernoulli JSDM gives (y - psi, psi(1-psi)). Adding a family to every
# latent route is therefore one callback, not a new fitter.
#
# This file resolves the formula terms into what the engine takes -- the field
# geometry (precision builder, constraints, site -> node map, weights) and the
# factor start -- and writes the engine's result onto the fit.
# =============================================================================


# ---------------------------------------------------------------------------
# Field geometry
# ---------------------------------------------------------------------------

# Candidate kappa grid from the observed-coordinate extent. The Matern range is
# rho = sqrt(8 nu) / kappa, so a grid of ranges spanning a tenth to most of the
# domain becomes a kappa grid.
.tobs_latent_spde_kappa_grid <- function(spec) {
  co  <- spec$obs_coords
  ext <- max(apply(co, 2L, function(z) diff(range(z))))
  nu  <- spec$nu %||% 1
  sqrt(8 * nu) / (ext * c(0.1, 0.2, 0.35, 0.5, 0.8))
}

# Continuous Matern (spde) field geometry. The mesh nodes carry the field and a
# barycentric projector A [n_sites x n_mesh] maps them onto sites, so the driver
# reads Mt = A (node -> site) and M = t(A) (the adjoint that scatters site scores
# onto nodes) -- the same linear map slot the areal group_var incidence uses.
.tobs_latent_spde_setup <- function(fields, model, what) {
  if (length(fields) > 1L || !is.null(fields[[1L]]$weight)) {
    stop(what, " spde() is the single shared intercept field only (no ",
         "varying-coefficient bar); use icar()/car_proper() for those.",
         call. = FALSE)
  }
  spec <- fields[[1L]]$tulpa_spec
  if (is.null(spec) || is.null(spec$A) || is.null(spec$G) ||
      is.null(spec$C0_diag)) {
    stop(what, " spde() needs the tulpa mesh spec (the projector + FEM ",
         "matrices).", call. = FALSE)
  }
  alpha <- (spec$nu %||% 1) + 1
  if (abs(alpha - round(alpha)) > 1e-8) {
    stop(sprintf(paste0("%s spde(nu = %g): the community latent field builds the ",
                        "integer-alpha SPDE precision (alpha = nu + 1 in 2-D), so ",
                        "nu must be a whole number (nu = 1 is the default)."),
                 what, spec$nu), call. = FALSE)
  }
  Ns <- model$n_sites
  if (nrow(spec$A) != Ns) {
    stop(sprintf(paste0("%s spde() projector has %d rows but the model has %d ",
                        "sites; the mesh must be built on the model's ",
                        "coordinates."), what, nrow(spec$A), Ns), call. = FALSE)
  }
  list(graph = NULL, n_nodes = as.integer(spec$n_mesh),
       M = Matrix::t(spec$A), type = "spde", K = 1L,
       W = matrix(1, Ns, 1L), field_labels = "intercept", spec = spec)
}


# Resolve the spatial term(s) on a community formula into the field geometry:
# the per-field per-site weights W [n_sites x K], the field labels, the field
# kind, the adjacency graph, and the node <- site linear map M -- a group_var
# incidence for an areal field with sites > cells, or the barycentric projector
# transpose for a continuous spde mesh. M is NULL when an areal field has one
# node per site. Shared by every community latent fitter.
.tobs_latent_field_setup <- function(spatial, model, what,
                                     allow = c("icar", "car_proper", "bym2")) {
  Ns     <- model$n_sites
  fields <- .tobs_resolve_occu_spatial_fields(spatial, model)
  ptype  <- fields[[1L]]$type %||% "icar"
  if (!all(vapply(fields, function(fd) identical(fd$type %||% "icar", ptype),
                  logical(1))) || !ptype %in% allow) {
    stop(what, " community field supports ",
         paste(sprintf("%s()", allow), collapse = " / "),
         " (one field kind per formula).", call. = FALSE)
  }
  # Continuous Matern field: mesh nodes + projector, no adjacency graph.
  if (identical(ptype, "spde")) {
    return(.tobs_latent_spde_setup(fields, model, what))
  }
  A <- fields[[1L]]$graph
  if (is.null(A)) {
    stop(what, " areal field needs the adjacency graph on the icar() term.",
         call. = FALSE)
  }
  n_nodes <- nrow(A)

  # group_var maps several sites to one field cell; the field has one node per
  # graph cell and the incidence M aggregates the sites onto it.
  Mmap <- NULL
  gv   <- fields[[1L]]$group_var
  if (!is.null(gv)) {
    if (is.null(model$data) || !gv %in% names(model$data)) {
      stop(sprintf("spatial group.var '%s' is not a column of the data.", gv),
           call. = FALSE)
    }
    node_of_site <- as.integer(model$data[[gv]])
    if (length(node_of_site) != Ns || anyNA(node_of_site) ||
        min(node_of_site) < 1L || max(node_of_site) > n_nodes) {
      stop(sprintf(paste0("spatial group.var '%s' must be an integer cell index in ",
                          "1..%d, one per site (%d sites)."), gv, n_nodes, Ns),
           call. = FALSE)
    }
    if (n_nodes != Ns) {
      Mmap <- Matrix::sparseMatrix(i = node_of_site, j = seq_len(Ns), x = 1,
                                   dims = c(n_nodes, Ns))
    }
  } else if (n_nodes != Ns) {
    stop(sprintf(paste0("icar graph has %d nodes but the model has %d sites; add ",
                        "group.var = \"<cell>\" to map sites to cells, or use one ",
                        "node per site."), n_nodes, Ns), call. = FALSE)
  }

  if (identical(ptype, "bym2") &&
      (length(fields) > 1L || !is.null(fields[[1L]]$weight) || !is.null(Mmap))) {
    stop(what, " bym2() field is the single shared intercept field only (no ",
         "varying-coefficient bar / group.var); use icar()/car_proper() for ",
         "those.", call. = FALSE)
  }

  K <- length(fields)
  W <- matrix(1, Ns, K)
  field_labels <- character(K)
  for (k in seq_len(K)) {
    wk <- fields[[k]]$weight
    if (!is.null(wk)) {
      if (length(wk) != Ns) {
        stop(what, " varying-coefficient field weight must be one value per ",
             "site.", call. = FALSE)
      }
      W[, k] <- as.numeric(wk)
      field_labels[k] <- fields[[k]]$weight_label %||% paste0("trend", k - 1L)
    } else {
      field_labels[k] <- "intercept"
    }
  }

  # An intrinsic field (icar, bym2's structured part) has one constant null
  # direction per connected component of the graph, each identified by its own
  # sum-to-zero; a proper CAR is constrained over the whole graph.
  label <- .tobs_graph_components(as.matrix(A))
  components <- if (ptype %in% c("icar", "bym2"))
    unname(split(seq_len(n_nodes), label)) else list(seq_len(n_nodes))
  list(graph = A, n_nodes = n_nodes, M = Mmap, type = ptype, K = K, W = W,
       field_labels = field_labels, components = components)
}

# The geometry `tulpa::tulpa_latent_ascent()` takes, built from a resolved field
# setup. The precision is a closure over the term's own graph or mesh, so the
# engine never sees the field kind:
#
#   icar        Q = D - W (tulpa_car_precision_Q at rho = 1), one sum-to-zero
#               per connected component, no hyperparameter.
#   car_proper  Q(rho) = D - rho W, rho profiled over a small grid by the field
#               marginal; constrained over the whole graph.
#   spde        Q(kappa) the integer-alpha Matern precision on the mesh at unit
#               scale (tulpa_spde_precision_Q, alpha = nu + 1 in 2-D; nu = 1 is
#               kappa^4 C0 + 2 kappa^2 G1 + G1 C0^-1 G1), the driver's per-field
#               tau carrying the overall scale (the SPDE tau^2); proper (full
#               rank), so unconstrained; kappa profiled over a range grid.
#   bym2        structured precision diag(node_prec) (D - W), scaled per
#               connected component and loaded at the reference component's
#               constant (.bym2_scaling); (sigma, rho) profiled.
#
# These grids are this fitter's own field-marginal profile: the ascent never
# hands a prior block to an outer-grid engine, so neither the outer-grid
# registry nor the `auto_grid()` provenance marker applies to them. The CAR
# nodes coincide with the engine's `joint_car_rho` because both place a
# correlation axis where an areal likelihood turns, not because either reads
# the other.
.tobs_latent_geometry <- function(geom) {
  base <- list(type = geom$type, M = geom$M, W = geom$W, n_nodes = geom$n_nodes,
               components = geom$components)
  switch(geom$type,
    spde = {
      spec <- geom$spec
      kappa_grid <- .tobs_latent_spde_kappa_grid(spec)
      c(base, list(
        Q = function(h) tulpa::tulpa_spde_precision_Q(spec, h[["kappa"]],
                                                      tau_spde = 1),
        hyper = c(kappa = stats::median(kappa_grid)),
        hyper_grid = list(kappa = kappa_grid)))
    },
    car_proper = {
      graph <- geom$graph
      c(base, list(
        Q = function(h) tulpa::tulpa_car_precision_Q(graph, rho = h[["rho"]]),
        hyper = c(rho = 0.95),
        hyper_grid = list(rho = c(0.5, 0.8, 0.95, 0.99))))
    },
    bym2 = {
      bsc <- .bym2_scaling(geom$graph)
      Qic <- tulpa::tulpa_car_precision_Q(geom$graph)
      if (!is.null(bsc$node_prec)) Qic <- Matrix::Diagonal(x = bsc$node_prec) %*% Qic
      c(base, list(
        Q = function(h) Qic,
        hyper = c(sigma = 0.6, rho = 0.5),
        # rho varies fastest, so the profile visits the cells sigma-major.
        hyper_grid = list(rho = c(0.25, 0.5, 0.75, 0.9),
                          sigma = c(0.3, 0.6, 1.0, 1.6)),
        scale = bsc$scale))
    },
    {
      Q0 <- tulpa::tulpa_car_precision_Q(geom$graph)
      c(base, list(Q = function(h) Q0, hyper = numeric(0), hyper_grid = NULL))
    })
}



# ---------------------------------------------------------------------------
# Latent factors
# ---------------------------------------------------------------------------

# Initial factors / loadings. The factors start at deterministic, mutually
# distinct smooth columns (a fixed seed-free start keeps a fit reproducible); the
# loadings at a small constant, which the first Newton pass spreads.
.tobs_latent_factor_init <- function(latent, Ns, S) {
  Qk <- as.integer(latent$n_factors %||% 1L)
  if (Qk < 1L) stop("latent(): n.factors must be >= 1.", call. = FALSE)
  if (Qk > S - 1L) {
    stop(sprintf("latent(): n.factors (%d) must be < n_species (%d).", Qk, S),
         call. = FALSE)
  }
  zeta <- matrix(0, Ns, Qk)
  for (q in seq_len(Qk)) zeta[, q] <- scale(cos(seq_len(Ns) * q))[, 1]
  list(n_factors = Qk, zeta = zeta, lambda = matrix(0.1, S, Qk))
}


# ---------------------------------------------------------------------------
# Block coordinate driver
# ---------------------------------------------------------------------------

# Fit a community model with a shared latent structure: a shared areal field
# (`spatial`), latent factors (`latent`), or BOTH (the spatial-factor case), by
# `tulpa::tulpa_latent_ascent()`. The single entry for every community latent
# route across every family: it resolves the terms, hands the engine the
# geometry and the factor start, and returns the engine's result in the
# spellings the attach writers and the fit's `latent_control` read.
#
# `make_oracle(em)` is a family callback returning the working oracle;
# `em_fit(site_off, fac_off, em_prev)` runs the family's community EM at the
# given offsets; `offset_of(em)` returns the arm's [n_sites x n_species]
# coefficient predictor.
#
# The outer budget is resolved per model. A field run stops at `tol` (tulpa's
# field budget, 400, is a cap); the factor block alternates with the
# coefficient block along a slow mode (the per-cell offset change decays
# about 2% per pass, so `tol` 1e-4 would need roughly 300 passes), so a
# factor model takes `factor.outer`, which
# each family sets from its OWN measurement (ms_count's curve: community
# intercept bias +0.0613 / +0.0237 / +0.0001 / -0.0038 at 25 / 60 / 150 / 400
# passes), and an explicit `control$max.outer` wins over both. `factor.starts`
# (the candidate starting directions of the first factor pass, each a full
# loading-EM run) is likewise set per family from the cost of one oracle
# evaluation: a closed-form GLMM density is cheap, the N-mixture oracle sums
# over the latent N per species-site. `n.quad` is the Gauss-Hermite node count
# per factor dimension; 5 is enough for what the estimator reads off the
# integral (test-community-latent-quad.R).
.tobs_community_latent_ascent <- function(spatial, latent, model, what,
                                          make_oracle, em_fit, offset_of,
                                          allow = c("icar", "car_proper", "bym2"),
                                          tol = 1e-4, max.outer = NULL,
                                          factor.outer = 25L,
                                          factor.starts = 8L,
                                          n.quad = NULL, verbose = FALSE) {
  S  <- model$n_species
  Ns <- model$n_sites
  geom <- if (is.null(spatial)) NULL else
    .tobs_latent_field_setup(spatial, model, what = what, allow = allow)
  fac <- if (is.null(latent)) NULL else .tobs_latent_factor_init(latent, Ns, S)

  res <- tulpa::tulpa_latent_ascent(
    make_oracle = make_oracle, em_fit = em_fit, offset_of = offset_of,
    n_sites = Ns, n_species = S,
    geometry = if (is.null(geom)) NULL else .tobs_latent_geometry(geom),
    factors = if (is.null(fac)) NULL else fac[c("zeta", "lambda")],
    control = list(tol = tol, max_outer = max.outer, factor_outer = factor.outer,
                   factor_starts = factor.starts, n_quad = n.quad,
                   verbose = verbose, label = what))

  st <- res$settings
  fd <- res$field
  list(em = res$em, geom = geom,
       # The settings the fit actually RAN under, not the ones it was called
       # with; the factor knobs read NA on a field-only fit.
       settings = list(max.outer = st$max_outer, tol = st$tol,
                       factor.outer = st$factor_outer,
                       factor.starts = st$factor_starts, n.quad = st$n_quad),
       field = if (is.null(fd)) NULL else list(
         Fm = fd$Fm, tau = fd$tau, type = fd$type,
         rho = if (identical(fd$type, "car_proper")) fd$hyper[["rho"]] else NA_real_,
         kappa = if (identical(fd$type, "spde")) fd$hyper[["kappa"]] else NA_real_,
         site_off = fd$site_off, bym2 = fd$bym2),
       factor = res$factor)
}


# ---------------------------------------------------------------------------
# Result assembly
# ---------------------------------------------------------------------------

# Record the block-coordinate settings the driver resolved. Idempotent, and
# called from both attach paths so a field-only, factor-only or spatial-factor
# fit all carry it.
.tobs_latent_attach_settings <- function(fit, res) {
  if (!is.null(res$settings)) fit$latent_control <- res$settings
  fit
}

# Attach the shared-field summary of a block-coordinate latent fit onto `fit`.
# `offset_slot` names the model slot holding the per-site eta contribution that
# fitted() / WAIC add back (e.g. "count_field_offset").
.tobs_latent_attach_field <- function(fit, res, spatial, offset_slot) {
  fit <- .tobs_latent_attach_settings(fit, res)
  fd <- res$field
  if (is.null(fd)) return(fit)
  fit$spatial <- spatial
  fit$method  <- "nested_laplace"
  if (!is.null(fd$bym2)) {
    # The reported field is the combined BYM2 field phi = a v + b u.
    fit$spatial_field <- as.numeric(fd$bym2$phi)
    fit$model[[offset_slot]] <- as.numeric(fd$bym2$phi)
    fit$spatial_hyper <- list(type = "bym2", sigma = fd$bym2$sigma,
                              rho = fd$bym2$rho, scale = fd$bym2$scale,
                              field_labels = "intercept")
    fit$means <- c(fit$means, sigma_field_intercept = fd$bym2$sigma)
    return(fit)
  }
  labels <- res$geom$field_labels
  # spatial_field is the per-node intercept field; the offset slot is the
  # per-site eta contribution (the node field mapped through the site -> cell map).
  fit$spatial_field <- as.numeric(fd$Fm[, 1L])
  fit$model[[offset_slot]] <- fd$site_off
  sigma_field <- 1 / sqrt(fd$tau)
  fit$spatial_hyper <- list(type = fd$type, tau = fd$tau, sigma = sigma_field,
                            rho = fd$rho, field_labels = labels)
  # A continuous Matern (spde) field additionally reports its range: the mesh
  # field lives on n_mesh nodes, and rho = sqrt(8 nu) / kappa is the practical
  # correlation range in the coordinate units.
  if (identical(fd$type, "spde") && is.finite(fd$kappa %||% NA_real_)) {
    nu <- res$geom$spec$nu %||% 1
    fit$spatial_hyper$kappa <- fd$kappa
    fit$spatial_hyper$range <- sqrt(8 * nu) / fd$kappa
    fit$means <- c(fit$means, spde_range = sqrt(8 * nu) / fd$kappa)
  }
  fit$means <- c(fit$means,
                 stats::setNames(sigma_field, paste0("sigma_field_", labels)))
  if (res$geom$K > 1L) {
    fit$trend_fields <- lapply(2:res$geom$K, function(k) as.numeric(fd$Fm[, k]))
    names(fit$trend_fields) <- labels[-1L]
    fit$trend_field <- fit$trend_fields[[1L]]
  }
  fit
}

# Attach the latent-factor summary. The loadings and factors are identified only
# up to rotation, so the recoverable target is the residual species covariance
# Sigma_res = lambda lambda' (reported alongside the implied correlation).
.tobs_latent_attach_factor <- function(fit, res, latent, model, offset_slot) {
  fit <- .tobs_latent_attach_settings(fit, res)
  fc <- res$factor
  if (is.null(fc)) return(fit)
  fit$latent <- latent
  lambda <- fc$lambda
  Sigma_res <- tcrossprod(lambda)
  dimnames(Sigma_res) <- list(model$species_names, model$species_names)
  rownames(lambda) <- model$species_names
  colnames(lambda) <- paste0("factor", seq_len(fc$n_factors))
  S <- nrow(lambda)
  fit$ms_factor <- list(
    n_factors = fc$n_factors, loadings = lambda, factors = fc$zeta,
    residual_cov = Sigma_res,
    residual_cor = stats::cov2cor(Sigma_res + diag(1e-10, S)),
    # Scale diagnostics. residual_cor is row-normalised and so cannot show a
    # magnitude problem; these can.
    magnitude_saturated = isTRUE(fc$magnitude_saturated),
    loading_em_converged = isTRUE(fc$loading_em_converged),
    marginal_loglik = fc$marginal_loglik)
  # The score-matched offset the fit conditioned on, not zeta t(lambda): fitted()
  # and WAIC read this slot, and they should see the predictor the coefficients
  # were estimated against.
  fit$model[[offset_slot]] <- fc$offset
  fit
}