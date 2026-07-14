# =============================================================================
# ms_count_spatial.R - community-spatial relative-abundance GLMM
# (ms_count() + a shared areal field; the spAbundance sfMsAbund analogue,
# gcol33/tulpaObs#117), Poisson.
#
#   log mu_{s,i} = X_i . (mu_beta + b_s) + f_{u(i)}
#   b_s ~ N(0, Sigma_beta),  f ~ ICAR(tau)   (one shared field across species)
#
# Fit by block coordinate ascent, reusing the pure-R community Laplace-EM
# (.tobs_community_em) WITHOUT modifying it: the shared field enters each
# species' log-likelihood as a fixed per-site OFFSET (captured in the sp_ll /
# sp_grad closures), so a coefficient update is an ordinary community EM; the
# field update given the coefficients is a self-contained Poisson-ICAR Laplace
# (an analytic sparse Newton with a closed-form tau M-step). Alternating the two
# converges to the joint mode. This sidesteps the community EM's
# finite-difference Hessian (which does not scale to an O(n_sites) field) and
# needs no C++.
# =============================================================================


# Poisson-ICAR field update for K covariate-weighted areal fields (the community
# SVC generalisation). Given the per-(species, site) coefficient offsets
# `offset_mat` [n_sites x n_species] and the per-field per-site weights `W`
# [n_sites x K] (an all-ones column for the intercept field, covariate values for
# a varying-coefficient field), refine the fields `F` [n_sites x K] and their
# precisions `tau` [K]. The per-site field contribution to eta is
# sum_k W[i,k] F[i,k]; per field the score aggregates the working residual over
# species, weighted by W[,k]. Joint Newton over the K*n field vector (a K x K
# block system, each (k,l) block diag(W[,k] W[,l] w) + (k==l) tau_k Q), then a
# per-field closed-form tau M-step (ICAR rank = n - 1). Each field is demeaned
# (the ICAR null space is the constant; the fixed effects absorb the level).
.ms_count_field_solve <- function(offset_mat, Q, F, tau, W, y_rowsum,
                                  max_iter = 50L, tol = 1e-8,
                                  constrain_mean = TRUE, rankdef = 1L, M = NULL) {
  Ns <- nrow(offset_mat); K <- ncol(W); Nn <- nrow(F)
  # `M` (n_nodes x n_sites incidence) maps sites to field nodes when several
  # sites share a node (a group_var / sites > cells design). Identity when the
  # field has one node per site.
  if (is.null(M)) M <- Matrix::Diagonal(Ns)
  Mt <- Matrix::t(M)
  build_H <- function(w) {
    blocks <- vector("list", K * K)
    for (k in seq_len(K)) for (l in seq_len(K)) {
      Dw <- M %*% Matrix::Diagonal(x = W[, k] * W[, l] * w) %*% Mt   # n_nodes^2
      blocks[[(k - 1L) * K + l]] <- if (k == l) Dw + tau[k] * Q else Dw
    }
    do.call(rbind, lapply(seq_len(K), function(k)
      do.call(cbind, blocks[((k - 1L) * K + 1L):(k * K)])))
  }
  for (it in seq_len(max_iter)) {
    F_site    <- as.matrix(Mt %*% F)                     # n_sites x K (node -> site)
    field_off <- rowSums(W * F_site)                     # per-site field eta
    w   <- rowSums(exp(pmin(offset_mat + field_off, 700)))
    g   <- unlist(lapply(seq_len(K), function(k)
      as.numeric(M %*% (W[, k] * (y_rowsum - w))) - tau[k] * as.numeric(Q %*% F[, k])))
    H   <- build_H(w)
    step <- as.numeric(Matrix::solve(H, g))
    F   <- F + matrix(step, Nn, K)
    # icar's null space is the constant (demean; the fixed effects absorb the
    # level). car_proper keeps the same sum-to-zero convention (rho only sets the
    # dependence strength of the deviations).
    if (isTRUE(constrain_mean)) for (k in seq_len(K)) F[, k] <- F[, k] - mean(F[, k])
    if (max(abs(step)) < tol) break
  }
  F_site    <- as.matrix(Mt %*% F)
  field_off <- rowSums(W * F_site)
  w    <- rowSums(exp(pmin(offset_mat + field_off, 700)))
  Cov  <- Matrix::solve(build_H(w))
  tau_new <- numeric(K)
  df <- Nn - as.integer(rankdef)                         # icar: nodes-1, car: nodes
  for (k in seq_len(K)) {
    idx  <- (k - 1L) * Nn + seq_len(Nn)
    quad <- as.numeric(t(F[, k]) %*% (Q %*% F[, k])) +
            sum(Matrix::diag(Q %*% Cov[idx, idx, drop = FALSE]))
    tau_new[k] <- df / max(quad, 1e-8)
  }
  list(F = F, tau = tau_new)
}


# log pseudo-determinant of the CAR precision on the sum-to-zero subspace: the
# sum of log of the n-1 largest eigenvalues (the smallest eigenvalue's
# eigenvector is ~ the constant, which the sum-to-zero constraint removes). Used
# to compare the field marginal across a proper-CAR rho grid.
.ms_count_car_logdet <- function(Qr) {
  ev <- sort(eigen(as.matrix(Qr), symmetric = TRUE, only.values = TRUE)$values,
             decreasing = TRUE)
  sum(log(pmax(ev[-length(ev)], 1e-10)))
}

# Marginal (Laplace) objective of the field solve at its converged mode, used to
# pick the proper-CAR rho over a small grid: the Poisson data log-likelihood plus
# the field GMRF prior (0.5 (df log tau + log|Q|) - 0.5 tau f'Q f) minus the
# Laplace normaliser 0.5 log|H|. `logdetQ` is log|Q(rho)| (per unit tau).
.ms_count_field_marginal <- function(offset_mat, Q, F, tau, W, y_rowsum, logdetQ,
                                     rankdef = 0L, M = NULL) {
  Ns <- nrow(offset_mat); K <- ncol(W); Nn <- nrow(F)
  if (is.null(M)) M <- Matrix::Diagonal(Ns)
  Mt <- Matrix::t(M)
  field_off <- rowSums(W * as.matrix(Mt %*% F))
  eta <- offset_mat + field_off
  ll  <- sum(y_rowsum * field_off) - sum(exp(pmin(eta, 700))) # data terms in f
  w   <- rowSums(exp(pmin(eta, 700)))
  build_H <- function() {
    blocks <- vector("list", K * K)
    for (k in seq_len(K)) for (l in seq_len(K)) {
      Dw <- M %*% Matrix::Diagonal(x = W[, k] * W[, l] * w) %*% Mt
      blocks[[(k - 1L) * K + l]] <- if (k == l) Dw + tau[k] * Q else Dw
    }
    do.call(rbind, lapply(seq_len(K), function(k)
      do.call(cbind, blocks[((k - 1L) * K + 1L):(k * K)])))
  }
  H <- build_H()
  df <- Nn - as.integer(rankdef)
  prior <- 0; for (k in seq_len(K))
    prior <- prior + 0.5 * (df * log(tau[k]) + logdetQ) -
             0.5 * tau[k] * as.numeric(t(F[, k]) %*% (Q %*% F[, k]))
  ll + prior - 0.5 * as.numeric(Matrix::determinant(H, logarithm = TRUE)$modulus)
}


# Fit a community count model with a shared latent structure: a shared areal
# field (`spatial`, the sfMsAbund / svcMsAbund case), latent factors (`latent`,
# the lfMsAbund case), or BOTH (the spatial-factor case). One block coordinate
# ascent over: (a) the community EM with the combined latent as a per-species
# offset; (b) the multi-field Poisson-ICAR field update (if `spatial`); (c) the
# Poisson factor update (if `latent`, with centred loadings when a field is also
# present so the field owns the shared spatial mean). Single source of truth for
# every latent-count route; the field-only / factor-only fitters are the
# spatial-only / latent-only special cases.
.tobs_fit_ms_count_latent <- function(model, spatial = NULL, latent = NULL,
                                      max.iter = 200L, tol = 1e-4,
                                      sigma.beta = 5, priors = NULL,
                                      max.outer = 25L, verbose = FALSE, ...) {
  if (!identical(model$response %||% "poisson", "poisson")) {
    stop("Community-spatial / latent-factor count (ms_count + a shared field or ",
         "latent()) is Poisson-only in this release (gcol33/tulpaObs#117).",
         call. = FALSE)
  }
  X <- model$X; P <- ncol(X); S <- model$n_species; Ns <- model$n_sites
  su <- model$summaries
  if (any(!model$valid)) {
    stop("ms_count() shared field / latent factors need a complete y (no NA ",
         "species-site cells) (gcol33/tulpaObs#117).", call. = FALSE)
  }
  y_mat    <- matrix(as.numeric(model$y), Ns, S)
  y_rowsum <- rowSums(y_mat)
  has_field  <- !is.null(spatial)
  has_factor <- !is.null(latent)

  # ---- field setup (if a shared areal field) ----
  W <- Q <- NULL; K <- 0L; field_labels <- character(0); Ffield <- NULL; tau <- NULL
  Mmap <- NULL; n_nodes <- Ns
  if (has_field) {
    fields <- .tobs_resolve_occu_spatial_fields(spatial, model)
    A <- fields[[1L]]$graph
    if (is.null(A)) {
      stop("ms_count() areal field needs the adjacency graph on the icar() term.",
           call. = FALSE)
    }
    # group_var maps several sites to one field cell (sites > cells); the field
    # has one node per graph cell and the site -> cell incidence M aggregates.
    # Without group_var the field is one node per site (identity).
    n_nodes <- nrow(A)
    gv <- fields[[1L]]$group_var
    if (!is.null(gv)) {
      if (is.null(model$data) || !gv %in% names(model$data))
        stop(sprintf("spatial group_var '%s' is not a column of the data.", gv),
             call. = FALSE)
      node_of_site <- as.integer(model$data[[gv]])
      if (length(node_of_site) != Ns || anyNA(node_of_site) ||
          min(node_of_site) < 1L || max(node_of_site) > n_nodes)
        stop(sprintf("spatial group_var '%s' must be an integer cell index in ",
                     "1..%d, one per site (%d sites).", gv, n_nodes, Ns),
             call. = FALSE)
      if (n_nodes != Ns)                                 # only build M when it aggregates
        Mmap <- Matrix::sparseMatrix(i = node_of_site, j = seq_len(Ns),
                                     x = 1, dims = c(n_nodes, Ns))
    } else if (n_nodes != Ns) {
      stop(sprintf("icar graph has %d nodes but the model has %d sites; add ",
                   "group_var = \"<cell>\" to map sites to cells, or use one ",
                   "node per site.", n_nodes, Ns), call. = FALSE)
    }
    ptype <- fields[[1L]]$type %||% "icar"
    if (!all(vapply(fields, function(fd) identical(fd$type %||% "icar", ptype),
                    logical(1))) || !ptype %in% c("icar", "car_proper")) {
      stop("ms_count() community field supports icar() or car_proper() (one ",
           "field kind per formula; bym2() is a follow-up, gcol33/tulpaObs#117).",
           call. = FALSE)
    }
    K <- length(fields); W <- matrix(1, Ns, K); field_labels <- character(K)
    for (k in seq_len(K)) {
      wk <- fields[[k]]$weight
      if (!is.null(wk)) {
        if (length(wk) != Ns)
          stop("ms_count() varying-coefficient field weight must be one value ",
               "per site.", call. = FALSE)
        W[, k] <- as.numeric(wk)
        field_labels[k] <- fields[[k]]$weight_label %||% paste0("trend", k - 1L)
      } else field_labels[k] <- "intercept"
    }
    # The field is sum-to-zero (it captures spatial DEVIATIONS; the intercept in X
    # owns the level -- so both are constrained/demeaned, rank n-1). icar fixes
    # the dependence at the intrinsic limit (Q = D - W); car_proper estimates the
    # dependence strength rho over a small grid (Q(rho) = D - rho W) by the field
    # marginal likelihood.
    Dg   <- Matrix::Diagonal(x = rowSums(A))
    Wadj <- methods::as(A, "CsparseMatrix")
    is_car <- identical(ptype, "car_proper")
    constrain_mean <- TRUE
    rankdef  <- 1L
    rho_grid <- if (is_car) c(0.5, 0.8, 0.95, 0.99) else NA_real_
    rho <- if (is_car) 0.95 else NA_real_
    Q <- if (is_car) Dg - rho * Wadj else Dg - Wadj
    Ffield <- matrix(0, n_nodes, K); tau <- rep(1, K)
  }
  Mt_field <- if (!is.null(Mmap)) Matrix::t(Mmap) else NULL
  # Per-site field offset from the (possibly aggregated) node field.
  site_field_off <- function(Ff) {
    Fs <- if (is.null(Mt_field)) Ff else as.matrix(Mt_field %*% Ff)
    rowSums(W * Fs)
  }

  # ---- factor setup (if latent factors) ----
  Qk <- 0L; eta <- lambda <- NULL
  if (has_factor) {
    Qk <- as.integer(latent$n_factors %||% 1L)
    if (Qk < 1L) stop("latent(): n_factors must be >= 1.", call. = FALSE)
    if (Qk > S - 1L)
      stop(sprintf("latent(): n_factors (%d) must be < n_species (%d).", Qk, S),
           call. = FALSE)
    eta <- matrix(0, Ns, Qk)
    for (q in seq_len(Qk)) eta[, q] <- scale(cos(seq_len(Ns) * q))[, 1]
    lambda <- matrix(0.1, S, Qk)
  }

  arm_idx <- list(mu = seq_len(P))
  mu0 <- numeric(P); mu0[1L] <- log(max(mean(y_rowsum / S), 0.1))
  em  <- NULL

  for (outer in seq_len(max.outer)) {
    site_off <- if (has_field) site_field_off(Ffield) else numeric(Ns)
    fac_off  <- if (has_factor) tcrossprod(eta, lambda) else matrix(0, Ns, S)
    # (a) community EM given the combined latent offset.
    sp_ll <- function(s, theta, global) {
      e <- as.numeric(su[[s]]$X %*% theta) + site_off + fac_off[, s]
      sum(stats::dpois(su[[s]]$y, exp(pmin(e, 700)), log = TRUE))
    }
    sp_grad <- function(s, theta, global) {
      mu_s <- exp(pmin(as.numeric(su[[s]]$X %*% theta) + site_off + fac_off[, s], 700))
      as.numeric(crossprod(su[[s]]$X, su[[s]]$y - mu_s))
    }
    em <- .tobs_community_em(
      S = S, P = P, arm_idx = arm_idx, sp_ll = sp_ll, sp_grad = sp_grad,
      init_mu = if (is.null(em)) mu0 else em$mu, init_global = numeric(0),
      penalize_global = FALSE, sigma_beta = sigma.beta, priors = priors,
      sigma_init = 0.3, max_iter = min(as.integer(max.iter), 60L),
      tol = as.numeric(tol), newton_max = 30L, verbose = FALSE)

    offset_mat <- vapply(seq_len(S),
                         function(s) as.numeric(X %*% (em$mu + em$b_list[[s]])),
                         numeric(Ns))
    delta <- 0
    # (b) field update (its "offset" holds the coefficients + the factor part).
    if (has_field) {
      fu <- .ms_count_field_solve(offset_mat + fac_off, Q, Ffield, tau, W, y_rowsum,
                                  constrain_mean = constrain_mean, rankdef = rankdef,
                                  M = Mmap)
      delta <- max(delta, max(abs(fu$F - Ffield))); Ffield <- fu$F; tau <- fu$tau
    }
    # (c) factor update (its "offset" holds the coefficients + the field part).
    if (has_factor) {
      site_off2 <- if (has_field) site_field_off(Ffield) else numeric(Ns)
      gu <- .ms_count_factor_update(offset_mat + matrix(site_off2, Ns, S), y_mat,
                                    eta, lambda, center = has_field)
      delta <- max(delta, max(abs(tcrossprod(gu$eta, gu$lambda) - fac_off)))
      eta <- gu$eta; lambda <- gu$lambda
    }
    if (isTRUE(verbose))
      message(sprintf("[ms_count latent %d] delta=%.2e", outer, delta))
    if (outer > 2L && delta < tol) break
  }

  # Proper-CAR rho selection: at the converged coefficients, re-solve the field
  # over a small rho grid and keep the rho maximising the field marginal.
  if (has_field && is_car) {
    offset_mat <- vapply(seq_len(S),
                         function(s) as.numeric(X %*% (em$mu + em$b_list[[s]])),
                         numeric(Ns))
    off_f <- offset_mat + (if (has_factor) tcrossprod(eta, lambda) else 0)
    best <- list(m = -Inf)
    for (rr in rho_grid) {
      Qr  <- Dg - rr * Wadj
      # log|Q(rho)| on the sum-to-zero subspace: drop the smallest eigenvalue
      # (the constrained constant direction) via the pseudo-determinant.
      ld  <- .ms_count_car_logdet(Qr)
      fr  <- .ms_count_field_solve(off_f, Qr, Ffield, tau, W, y_rowsum,
                                   constrain_mean = TRUE, rankdef = 1L, M = Mmap)
      mm  <- .ms_count_field_marginal(off_f, Qr, fr$F, fr$tau, W, y_rowsum, ld,
                                      rankdef = 1L, M = Mmap)
      if (is.finite(mm) && mm > best$m)
        best <- list(m = mm, rho = rr, F = fr$F, tau = fr$tau)
    }
    if (is.finite(best$m)) { rho <- best$rho; Ffield <- best$F; tau <- best$tau }
  }

  fit <- build_ms_count_fit(model, em, arm_idx, disp = NULL)
  fit$method <- if (has_field) "nested_laplace" else "laplace"
  if (has_field) {
    fit$spatial <- spatial
    # spatial_field is the per-node intercept field; count_field_offset is the
    # per-site eta contribution (node field mapped through the site -> cell map).
    fit$spatial_field <- as.numeric(Ffield[, 1L])
    fit$model$count_field_offset <- site_field_off(Ffield)
    sigma_field <- 1 / sqrt(tau)
    fit$spatial_hyper <- list(type = ptype, tau = tau, sigma = sigma_field,
                              rho = if (is_car) rho else NA_real_,
                              field_labels = field_labels)
    fit$means <- c(fit$means, stats::setNames(sigma_field,
                              paste0("sigma_field_", field_labels)))
    if (K > 1L) {
      fit$trend_fields <- lapply(2:K, function(k) as.numeric(Ffield[, k]))
      names(fit$trend_fields) <- field_labels[-1L]
      fit$trend_field <- fit$trend_fields[[1L]]
    }
  }
  if (has_factor) {
    fit$latent <- latent
    Sigma_res <- tcrossprod(lambda); dimnames(Sigma_res) <-
      list(model$species_names, model$species_names)
    rownames(lambda) <- model$species_names
    colnames(lambda) <- paste0("factor", seq_len(Qk))
    fit$ms_factor <- list(
      n_factors = Qk, loadings = lambda, factors = eta,
      residual_cov = Sigma_res,
      residual_cor = stats::cov2cor(Sigma_res + diag(1e-10, S)))
    fit$model$count_factor_offset <- tcrossprod(eta, lambda)
  }
  fit
}

# Back-compat thin wrappers -> the unified latent fitter.
.tobs_fit_ms_count_spatial <- function(model, spatial, ...) {
  .tobs_fit_ms_count_latent(model, spatial = spatial, latent = NULL, ...)
}
