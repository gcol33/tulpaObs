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
                                  max_iter = 50L, tol = 1e-8) {
  Ns <- nrow(offset_mat); K <- ncol(W)
  ywt <- vapply(seq_len(K), function(k) as.numeric(W[, k] * y_rowsum), numeric(Ns))
  build_H <- function(w) {
    blocks <- vector("list", K * K)
    for (k in seq_len(K)) for (l in seq_len(K)) {
      D <- Matrix::Diagonal(x = W[, k] * W[, l] * w)
      blocks[[(k - 1L) * K + l]] <- if (k == l) D + tau[k] * Q else D
    }
    do.call(rbind, lapply(seq_len(K), function(k)
      do.call(cbind, blocks[((k - 1L) * K + 1L):(k * K)])))
  }
  for (it in seq_len(max_iter)) {
    field_off <- rowSums(W * F)                          # per-site field eta
    w   <- rowSums(exp(pmin(offset_mat + field_off, 700)))
    g   <- unlist(lapply(seq_len(K), function(k)
      ywt[, k] - W[, k] * w - tau[k] * as.numeric(Q %*% F[, k])))
    H   <- build_H(w)
    step <- as.numeric(Matrix::solve(H, g))
    F   <- F + matrix(step, Ns, K)
    for (k in seq_len(K)) F[, k] <- F[, k] - mean(F[, k])
    if (max(abs(step)) < tol) break
  }
  field_off <- rowSums(W * F)
  w    <- rowSums(exp(pmin(offset_mat + field_off, 700)))
  Cov  <- Matrix::solve(build_H(w))
  tau_new <- numeric(K)
  for (k in seq_len(K)) {
    idx  <- (k - 1L) * Ns + seq_len(Ns)
    quad <- as.numeric(t(F[, k]) %*% (Q %*% F[, k])) +
            sum(Matrix::diag(Q %*% Cov[idx, idx, drop = FALSE]))
    tau_new[k] <- (Ns - 1) / max(quad, 1e-8)
  }
  list(F = F, tau = tau_new)
}


# Fit the community-spatial count model. `model` is the (natural-scale) ms_count
# model; `spatial` the resolved icar term. Block coordinate ascent between the
# community EM (field as offset) and the Poisson-ICAR field update.
.tobs_fit_ms_count_spatial <- function(model, spatial,
                                       max.iter = 200L, tol = 1e-4,
                                       sigma.beta = 5, priors = NULL,
                                       max.outer = 20L, verbose = FALSE, ...) {
  if (!identical(model$response %||% "poisson", "poisson")) {
    stop("Community-spatial count (ms_count + areal field) is Poisson-only: ",
         "with one field node per site an overdispersed community count is not ",
         "identified against the shared field (gcol33/tulpaObs#117).",
         call. = FALSE)
  }
  X   <- model$X
  P   <- ncol(X)
  S   <- model$n_species
  Ns  <- model$n_sites
  su  <- model$summaries
  # Every species must be observed at every site for the shared-field row
  # aggregation (one field node per site). Missing species-site cells are not
  # yet wired.
  if (any(!model$valid)) {
    stop("ms_count() areal field needs a complete y (no NA species-site cells); ",
         "missing cells with a shared field are not yet wired ",
         "(gcol33/tulpaObs#117).", call. = FALSE)
  }
  y_mat    <- matrix(as.numeric(model$y), Ns, S)
  y_rowsum <- rowSums(y_mat)

  # Resolve the field(s): a plain icar() is one intercept field; the bar
  # spatial(~ 1 + w || cell, graph) is an intercept field + one varying-
  # coefficient (SVC) field per covariate (the svcMsAbund case). All share one
  # graph; each contributes a column to the per-site weight matrix W (all-ones
  # for the intercept, the covariate values for a trend field). One field node
  # per site (identity map).
  fields <- .tobs_resolve_occu_spatial_fields(spatial, model)
  A <- fields[[1L]]$graph
  if (is.null(A)) {
    stop("ms_count() areal field needs the adjacency graph on the icar() term.",
         call. = FALSE)
  }
  if (nrow(A) != Ns) {
    stop(sprintf("icar graph has %d nodes but the model has %d sites; one field ",
                 "node per site is required for the community field.",
                 nrow(A), Ns), call. = FALSE)
  }
  if (!all(vapply(fields, function(fd) identical(fd$type %||% "icar", "icar"),
                  logical(1)))) {
    stop("ms_count() community field supports icar() only (bym2()/car_proper() ",
         "are follow-ups, gcol33/tulpaObs#117).", call. = FALSE)
  }
  K <- length(fields)
  W <- matrix(1, Ns, K)
  field_labels <- character(K)
  for (k in seq_len(K)) {
    wk <- fields[[k]]$weight
    if (!is.null(wk)) {
      if (length(wk) != Ns) {
        stop("ms_count() varying-coefficient field weight must be one value ",
             "per site.", call. = FALSE)
      }
      W[, k] <- as.numeric(wk)
      field_labels[k] <- fields[[k]]$weight_label %||% paste0("trend", k - 1L)
    } else {
      field_labels[k] <- "intercept"
    }
  }
  Q <- Matrix::Diagonal(x = rowSums(A)) - methods::as(A, "CsparseMatrix")

  arm_idx <- list(mu = seq_len(P))
  Ffield  <- matrix(0, Ns, K)
  tau     <- rep(1, K)
  mu0 <- numeric(P); mu0[1L] <- log(max(mean(y_rowsum / S), 0.1))
  em  <- NULL

  for (outer in seq_len(max.outer)) {
    field_off <- rowSums(W * Ffield)                   # per-site eta offset
    # (a) community EM given the current field offset.
    sp_ll <- function(s, theta, global) {
      eta <- as.numeric(su[[s]]$X %*% theta) + field_off
      sum(stats::dpois(su[[s]]$y, exp(pmin(eta, 700)), log = TRUE))
    }
    sp_grad <- function(s, theta, global) {
      mu_s <- exp(pmin(as.numeric(su[[s]]$X %*% theta) + field_off, 700))
      as.numeric(crossprod(su[[s]]$X, su[[s]]$y - mu_s))
    }
    em <- .tobs_community_em(
      S = S, P = P, arm_idx = arm_idx, sp_ll = sp_ll, sp_grad = sp_grad,
      init_mu = if (is.null(em)) mu0 else em$mu, init_global = numeric(0),
      penalize_global = FALSE, sigma_beta = sigma.beta, priors = priors,
      sigma_init = 0.3, max_iter = min(as.integer(max.iter), 60L),
      tol = as.numeric(tol), newton_max = 30L, verbose = FALSE)

    # (b) field update given the coefficients.
    offset_mat <- vapply(seq_len(S),
                         function(s) as.numeric(X %*% (em$mu + em$b_list[[s]])),
                         numeric(Ns))
    fu <- .ms_count_field_solve(offset_mat, Q, Ffield, tau, W, y_rowsum)
    delta <- max(abs(fu$F - Ffield))
    Ffield <- fu$F; tau <- fu$tau
    if (isTRUE(verbose)) {
      message(sprintf("[ms_count spatial %d] field delta=%.2e tau=%s",
                      outer, delta, paste(round(tau, 3), collapse = ",")))
    }
    if (outer > 2L && delta < tol) break
  }

  fit <- build_ms_count_fit(model, em, arm_idx, disp = NULL)
  fit$method  <- "nested_laplace"
  fit$spatial <- spatial
  # The intercept field is the headline spatial_field; any varying-coefficient
  # (SVC) fields are the trend field(s). fitted() / WAIC use the combined
  # per-site offset sum_k W[,k] F[,k].
  fit$spatial_field <- as.numeric(Ffield[, 1L])
  fit$model$count_field_offset <- rowSums(W * Ffield)
  sigma_field <- 1 / sqrt(tau)
  fit$spatial_hyper <- list(type = "icar", tau = tau, sigma = sigma_field,
                            field_labels = field_labels)
  fit$means <- c(fit$means,
                 stats::setNames(sigma_field,
                                 paste0("sigma_field_", field_labels)))
  if (K > 1L) {
    fit$trend_fields <- lapply(2:K, function(k) as.numeric(Ffield[, k]))
    names(fit$trend_fields) <- field_labels[-1L]
    fit$trend_field <- fit$trend_fields[[1L]]
  }
  fit
}
