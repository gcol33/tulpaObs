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


# Poisson-ICAR field update. Given the per-(species, site) coefficient offsets
# `offset_mat` [n_sites x n_species] (X_i.(mu + b_s)), refine the shared field f
# and its precision tau. The field enters every species at its site; per site i
# the score aggregates over species. Newton on f with the ICAR prior tau*Q, then
# a closed-form tau M-step (rank / (f'Q f + tr(Q Cov_f)); ICAR rank = n - 1). f
# is demeaned each step (the ICAR null space is the constant; the community
# intercept absorbs the level).
.ms_count_field_solve <- function(offset_mat, Q, f, tau, y_rowsum,
                                  max_iter = 50L, tol = 1e-8) {
  Ns <- nrow(offset_mat)
  for (it in seq_len(max_iter)) {
    mu  <- exp(pmin(offset_mat + f, 700))       # n_sites x n_species (recycled f)
    w   <- rowSums(mu)
    g   <- y_rowsum - w - tau * as.numeric(Q %*% f)
    H   <- Matrix::Diagonal(x = w) + tau * Q
    step <- as.numeric(Matrix::solve(H, g))
    f   <- f + step
    f   <- f - mean(f)
    if (max(abs(step)) < tol) break
  }
  w    <- rowSums(exp(pmin(offset_mat + f, 700)))
  H    <- Matrix::Diagonal(x = w) + tau * Q
  Covf <- Matrix::solve(H)
  quad <- as.numeric(t(f) %*% (Q %*% f)) + sum(Matrix::diag(Q %*% Covf))
  tau_new <- (Ns - 1) / max(quad, 1e-8)
  list(f = f, tau = tau_new)
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
  ptype <- spatial$type %||% "icar"
  if (!identical(ptype, "icar")) {
    stop(sprintf(paste0(
      "ms_count() areal field supports icar(); got '%s'. bym2()/car_proper() ",
      "are not yet wired for the community count field (gcol33/tulpaObs#117)."),
      ptype), call. = FALSE)
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

  # ICAR precision Q = D - W from the adjacency graph.
  A <- spatial$graph
  if (is.null(A)) {
    stop("ms_count() areal field needs the adjacency graph on the icar() term.",
         call. = FALSE)
  }
  if (nrow(A) != Ns) {
    stop(sprintf("icar graph has %d nodes but the model has %d sites; one field ",
                 "node per site is required.", nrow(A), Ns), call. = FALSE)
  }
  Q <- Matrix::Diagonal(x = rowSums(A)) - methods::as(A, "CsparseMatrix")

  arm_idx <- list(mu = seq_len(P))
  f   <- numeric(Ns)
  tau <- 1
  mu0 <- numeric(P); mu0[1L] <- log(max(mean(y_rowsum / S), 0.1))
  em  <- NULL

  for (outer in seq_len(max.outer)) {
    # (a) community EM given the current field offset f.
    sp_ll <- function(s, theta, global) {
      eta <- as.numeric(su[[s]]$X %*% theta) + f
      sum(stats::dpois(su[[s]]$y, exp(pmin(eta, 700)), log = TRUE))
    }
    sp_grad <- function(s, theta, global) {
      mu_s <- exp(pmin(as.numeric(su[[s]]$X %*% theta) + f, 700))
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
    fu <- .ms_count_field_solve(offset_mat, Q, f, tau, y_rowsum)
    delta <- max(abs(fu$f - f))
    f <- fu$f; tau <- fu$tau
    if (isTRUE(verbose)) {
      message(sprintf("[ms_count spatial %d] field delta=%.2e tau=%.3f",
                      outer, delta, tau))
    }
    if (outer > 2L && delta < tol) break
  }

  fit <- build_ms_count_fit(model, em, arm_idx, disp = NULL)
  # Attach the shared field (demeaned) + its SD, and make fitted() / WAIC
  # field-aware via a per-site offset on the model (the field is a per-site eta
  # offset, one node per site).
  fit$method        <- "nested_laplace"
  fit$spatial       <- spatial
  fit$spatial_field <- as.numeric(f)
  fit$spatial_hyper <- list(type = "icar", tau = tau, sigma = 1 / sqrt(tau))
  fit$means         <- c(fit$means, sigma_field = 1 / sqrt(tau))
  fit$model$count_field_offset <- as.numeric(f)
  fit
}
