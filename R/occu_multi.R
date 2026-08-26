# occu_multi.R - Multi-species co-occurrence occupancy (Rota et al. 2016;
# unmarked occuMulti). The joint occupancy state of S species at a site is a
# vector z in {0,1}^S drawn from a log-linear (Ising-type) model whose natural
# parameters are first-order (per species) and second-order (per species pair)
# terms; conditional on presence each species is detected with its own
# probability. Species interactions are the second-order natural parameters
# (positive = co-occur more than independent, negative = avoid).
#
#   P(z) = exp(f(z)) / sum_{z'} exp(f(z')),
#   f(z) = sum_k a_k z_k + sum_{k<l} b_kl z_k z_l,
#   a_k = X_state . beta_a_k,   b_kl = X_state . beta_b_kl,   logit p_k = X_det . beta_p_k
#
# and given z, species k contributes prod_j p_kj^{y_kj}(1-p_kj)^{1-y_kj} when
# z_k = 1, or 1{no detection of k} when z_k = 0. The latent z integrates out by
# enumerating the 2^S states, so the exact marginal is maximised (optim BFGS)
# with an observed-information vcov -- the royle_nichols() recipe generalised to
# a joint multi-species state. Detection is site-level; first + second order;
# a shared covariate design across the natural parameters (each still carries its
# own coefficients). Per-natural-parameter formulas, higher-order terms, and
# visit-level detection are documented follow-ups.
#
#   .tobs_build_occu_multi()   data binder -> model_type = "occu_multi"
#   .tobs_fit_occu_multi()     optim over the 2^S-state marginal
#   .dispatch_occu_multi()     tobs() entry (bind + fit + assemble)

# All 2^S binary state vectors as rows of a [2^S x S] matrix (species in columns).
.occu_multi_states <- function(S) {
  as.matrix(expand.grid(rep(list(c(0L, 1L)), S))[, rev(seq_len(S)), drop = FALSE])
}

# Species pairs (k < l) for the second-order natural parameters.
.occu_multi_pairs <- function(S) {
  if (S < 2L) return(list())
  utils::combn(S, 2L, simplify = FALSE)
}

# ---------------------------------------------------------------------------
# Marginal log-likelihood (per site), reused by fit / ploglik / simulate.
# ---------------------------------------------------------------------------

# natpar1 [N x S], natpar2 [N x n_pairs], p_mat [N x S] detection probs, k_mat /
# n_mat [N x S] per-species detections / valid visits, Z [2^S x S] states, pairs
# the k<l index pairs. Returns length-N marginal log-likelihood.
.occu_multi_site_loglik <- function(natpar1, natpar2, p_mat, k_mat, n_mat,
                                    Z, pairs) {
  N        <- nrow(p_mat); S <- ncol(Z); n_states <- nrow(Z)
  p_mat    <- pmin(pmax(p_mat, 1e-12), 1 - 1e-12)
  logem1   <- k_mat * log(p_mat) + (n_mat - k_mat) * log(1 - p_mat)   # z_k = 1
  anydet   <- k_mat > 0                                               # [N x S]

  # Second-order state indicators z_k z_l per state.
  Zpair <- matrix(0, n_states, length(pairs))
  for (m in seq_along(pairs))
    Zpair[, m] <- Z[, pairs[[m]][1L]] * Z[, pairs[[m]][2L]]

  # f(z) per site x state.
  Fmat <- natpar1 %*% t(Z)
  if (length(pairs) > 0L) Fmat <- Fmat + natpar2 %*% t(Zpair)         # [N x n_states]

  # Emission log-prob per site x state: sum_k (z_k=1: logem1) ; a species with
  # z_k = 0 but a detection makes the state impossible (-Inf).
  logG <- matrix(0, N, n_states)
  for (s in seq_len(n_states)) {
    z  <- Z[s, ]
    em <- as.numeric(logem1 %*% z)                       # sum_k z_k logem1_ik
    invalid <- as.numeric(anydet %*% (1 - z)) > 0        # z_k=0 but detected
    em[invalid] <- -Inf
    logG[, s] <- em
  }

  lse <- function(r) { m <- max(r); if (!is.finite(m)) return(-Inf)
                       m + log(sum(exp(r - m))) }
  logZ   <- apply(Fmat, 1L, lse)                         # partition function
  lognum <- apply(Fmat + logG, 1L, lse)                  # weighted emission
  lognum - logZ
}

# ---------------------------------------------------------------------------
# Data binder
# ---------------------------------------------------------------------------

# `y` is a length-S list of n_sites x max_visits 0/1/NA detection matrices (or a
# 3D array [n_sites x max_visits x S]). `species` names the arms. `state_formula`
# is the shared occupancy covariate design (each natural parameter carries its
# own coefficients); `det_formula` the shared site-level detection design.
.tobs_build_occu_multi <- function(state_formula, det_formula, data, y, species) {
  if (is.array(y) && length(dim(y)) == 3L) {
    S <- dim(y)[3L]
    y <- lapply(seq_len(S), function(s) y[, , s])
  }
  if (!is.list(y)) stop("occu_multi() y must be a list of S detection matrices ",
                        "or a 3D [sites x visits x species] array.", call. = FALSE)
  S <- length(y)
  if (S < 2L) stop("occu_multi() needs at least 2 species.", call. = FALSE)
  if (is.null(species)) species <- names(y) %||% paste0("sp", seq_len(S))

  y <- lapply(y, function(m) matrix(as.integer(round(m)), nrow(m), ncol(m)))
  n_sites    <- nrow(y[[1L]]); max_visits <- ncol(y[[1L]])
  .tobs_check_site_count(n_sites, nrow(data), "sites")

  k_mat <- n_mat <- matrix(0, n_sites, S)
  for (s in seq_len(S)) {
    vs <- !is.na(y[[s]])
    if (any(y[[s]][vs] != 0L & y[[s]][vs] != 1L))
      stop(sprintf("occu_multi() species '%s' y must be 0/1/NA.", species[s]),
           call. = FALSE)
    k_mat[, s] <- rowSums(y[[s]] == 1L & vs)
    n_mat[, s] <- rowSums(vs)
  }

  bind    <- .tobs_bind_formulas(list(state = state_formula, det = det_formula),
                                 data)
  X_state <- stats::model.matrix(bind$fe$state, data)
  X_det   <- stats::model.matrix(bind$fe$det,   data)

  pairs <- .occu_multi_pairs(S)
  Z     <- .occu_multi_states(S)

  # process_info: one entry per natural parameter (first + second order) then one
  # per species detection arm, so coef() / vcov() name every coordinate.
  pi_list <- list()
  for (s in seq_len(S))
    pi_list[[length(pi_list) + 1L]] <- list(
      name = paste0("f_", species[s]), p = ncol(X_state),
      coef_names = colnames(X_state), link = "identity")
  for (pr in pairs)
    pi_list[[length(pi_list) + 1L]] <- list(
      name = paste0("f_", species[pr[1L]], "_", species[pr[2L]]),
      p = ncol(X_state), coef_names = colnames(X_state), link = "identity")
  for (s in seq_len(S))
    pi_list[[length(pi_list) + 1L]] <- list(
      name = paste0("p_", species[s]), p = ncol(X_det),
      coef_names = colnames(X_det), link = "logit")

  structure(list(
    model_type  = "occu_multi",
    y           = y,
    species     = species,
    S           = S,
    pairs       = pairs,
    Z           = Z,
    k_mat       = k_mat,
    n_mat       = n_mat,
    X_state     = X_state,
    X_det       = X_det,
    formulas    = list(state = bind$fe$state, det = bind$fe$det),
    structured_terms = bind$terms,
    data        = data,
    n_sites     = n_sites,
    max_visits  = max_visits,
    process_info = pi_list
  ), class = "tobs_model")
}

# Split a packed theta into per-natural-parameter state betas, per-pair betas,
# and per-species detection betas, and form the per-site natural parameters +
# detection probabilities.
.occu_multi_unpack <- function(theta, model) {
  S <- model$S; pairs <- model$pairs
  p_state <- ncol(model$X_state); p_det <- ncol(model$X_det)
  n_first <- S; n_second <- length(pairs)
  off <- 0L
  natpar1 <- matrix(0, model$n_sites, S)
  for (s in seq_len(S)) {
    natpar1[, s] <- as.vector(model$X_state %*% theta[off + seq_len(p_state)])
    off <- off + p_state
  }
  natpar2 <- matrix(0, model$n_sites, max(n_second, 1L))
  if (n_second > 0L) {
    natpar2 <- matrix(0, model$n_sites, n_second)
    for (m in seq_len(n_second)) {
      natpar2[, m] <- as.vector(model$X_state %*% theta[off + seq_len(p_state)])
      off <- off + p_state
    }
  }
  p_mat <- matrix(0, model$n_sites, S)
  for (s in seq_len(S)) {
    p_mat[, s] <- stats::plogis(as.vector(model$X_det %*% theta[off + seq_len(p_det)]))
    off <- off + p_det
  }
  list(natpar1 = natpar1, natpar2 = natpar2, p_mat = p_mat)
}

# ---------------------------------------------------------------------------
# Fitter
# ---------------------------------------------------------------------------

.tobs_fit_occu_multi <- function(model, verbose = TRUE,
                                 max.iter = NULL, tol = NULL, ...) {
  S <- model$S; pairs <- model$pairs
  p_state <- ncol(model$X_state); p_det <- ncol(model$X_det)
  n_theta <- (S + length(pairs)) * p_state + S * p_det
  k_mat <- model$k_mat; n_mat <- model$n_mat; Z <- model$Z

  nll <- function(theta) {
    up  <- .occu_multi_unpack(theta, model)
    ll  <- .occu_multi_site_loglik(up$natpar1, up$natpar2, up$p_mat,
                                   k_mat, n_mat, Z, pairs)
    val <- -sum(ll)
    if (is.finite(val)) val else 1e10
  }

  # Init: first-order natural params from each species' occupancy logit, second
  # order at 0 (independence), detection from the per-species detection rate.
  init <- numeric(n_theta)
  off  <- 0L
  for (s in seq_len(S)) {
    occ_hat <- mean(model$k_mat[, s] > 0)
    init[off + 1L] <- stats::qlogis(min(max(occ_hat, 0.05), 0.95))
    off <- off + p_state
  }
  off <- off + length(pairs) * p_state                 # second order at 0
  for (s in seq_len(S)) {
    p_hat <- sum(model$k_mat[, s]) / max(sum(model$n_mat[, s]), 1)
    init[off + 1L] <- stats::qlogis(min(max(p_hat, 0.05), 0.95))
    off <- off + p_det
  }

  par_names <- unlist(lapply(model$process_info, function(pp)
    paste0(pp$name, "_", pp$coef_names)))

  .tobs_bfgs_marginal_fit(nll, init, par_names, model, N = model$n_sites,
                          max.iter = max.iter %||% 800L, tol = tol)
}

# ---------------------------------------------------------------------------
# tobs() dispatcher
# ---------------------------------------------------------------------------

.dispatch_occu_multi <- function(formula, data, family, detection, y, visits,
                                 engine, priors, control,
                                 approx = "gaussian_laplace",
                                 correction = "none", species = NULL, ...) {
  if (is.null(detection))
    stop("occu_multi() requires a `detection` formula (shared site-level ",
         "per-species detection).", call. = FALSE)
  if (is.null(y))
    stop("occu_multi() requires `y` (a list of S detection matrices or a 3D ",
         "[sites x visits x species] array).", call. = FALSE)
  if (!is.null(visits))
    stop("occu_multi() detection is site-level; visit-level detection ",
         "covariates (`visits`) are not yet supported.", call. = FALSE)
  model <- .tobs_build_occu_multi(
    state_formula = formula, det_formula = detection, data = data, y = y,
    species = species)
  .tobs_fit_occu_multi(model, verbose = isTRUE(control$verbose),
                       max.iter = control$max.iter, tol = control$tol)
}

# ---------------------------------------------------------------------------
# fitted / predict / residuals / pointwise log-likelihood / simulate
# ---------------------------------------------------------------------------

# Per-site marginal occupancy probability psi_k = P(z_k = 1) per species, and the
# per-species detection probability.
.tobs_fitted_occu_multi <- function(object) {
  model <- object$model
  up    <- .occu_multi_unpack(object$means, model)
  Z <- model$Z; pairs <- model$pairs
  Zpair <- matrix(0, nrow(Z), length(pairs))
  for (m in seq_along(pairs)) Zpair[, m] <- Z[, pairs[[m]][1L]] * Z[, pairs[[m]][2L]]
  Fmat <- up$natpar1 %*% t(Z)
  if (length(pairs) > 0L) Fmat <- Fmat + up$natpar2 %*% t(Zpair)
  W <- exp(Fmat - apply(Fmat, 1L, max))
  W <- W / rowSums(W)                                  # P(z) per site x state
  psi <- W %*% Z                                       # marginal P(z_k = 1)
  colnames(psi) <- model$species
  colnames(up$p_mat) <- model$species
  list(psi = psi, p = up$p_mat)
}

# predict(): marginal occupancy (default) or detection per species at the fitted
# or new design.
.tobs_predict_occu_multi <- function(object, newdata = NULL,
                                     type = c("state", "detection")) {
  type  <- match.arg(type)
  model <- object$model
  if (is.null(newdata)) {
    fv <- .tobs_fitted_occu_multi(object)
    return(if (identical(type, "detection")) fv$p else fv$psi)
  }
  # New design: rebuild a shim model on newdata (same coefficients), then reduce.
  shim <- model
  shim$X_state <- stats::model.matrix(model$formulas$state, newdata)
  shim$X_det   <- stats::model.matrix(model$formulas$det,   newdata)
  shim$n_sites <- nrow(shim$X_state)
  obj2 <- object; obj2$model <- shim
  fv <- .tobs_fitted_occu_multi(obj2)
  if (identical(type, "detection")) fv$p else fv$psi
}

# residuals(): per-(site, species) response residual on the observed
# any-detection indicator against the marginal detection probability
# psi_k * (1 - (1 - p_k)^{n_ik}).
.tobs_residuals_occu_multi <- function(object, type) {
  model <- object$model
  fv    <- .tobs_fitted_occu_multi(object)
  psi <- fv$psi; p <- fv$p
  pdet <- psi * (1 - (1 - p)^model$n_mat)
  obs  <- (model$k_mat > 0) * 1
  eps  <- 1e-10; pc <- pmin(pmax(pdet, eps), 1 - eps)
  res <- switch(type,
    response = obs - pdet,
    pearson  = (obs - pdet) / sqrt(pc * (1 - pc) + eps),
    deviance = sign(obs - pdet) * sqrt(2 * abs(
      ifelse(obs > 0, obs * log(obs / pc), 0) +
      ifelse(obs < 1, (1 - obs) * log((1 - obs) / (1 - pc)), 0))))
  list(occ = res, det = NULL)
}

# Pointwise log-likelihood [n_draws x n_sites] over the posterior draws.
.tobs_ploglik_occu_multi <- function(object, n.draws = 1000L, n.threads = 1L) {
  model <- object$model
  draws <- object$draws
  if (!is.null(n.draws) && as.integer(n.draws) < nrow(draws)) {
    draws <- draws[seq_len(as.integer(n.draws)), , drop = FALSE]
  }
  k_mat <- model$k_mat; n_mat <- model$n_mat; Z <- model$Z; pairs <- model$pairs
  t(vapply(seq_len(nrow(draws)), function(d) {
    up <- .occu_multi_unpack(draws[d, ], model)
    .occu_multi_site_loglik(up$natpar1, up$natpar2, up$p_mat, k_mat, n_mat, Z, pairs)
  }, numeric(model$n_sites)))
}

# Posterior replicate detection histories: draw a coefficient vector, a joint
# occupancy state per site from P(z), then per-species detections at the observed
# visit pattern. Returns a list of S matrices (or a list of those for nsim > 1).
.tobs_simulate_occu_multi <- function(object, nsim = 1) {
  model <- object$model; S <- model$S; Z <- model$Z; pairs <- model$pairs
  valid <- lapply(model$y, function(m) !is.na(m))
  draw_one <- function() {
    idx <- sample.int(nrow(object$draws), 1L)
    up  <- .occu_multi_unpack(object$draws[idx, ], model)
    Zpair <- matrix(0, nrow(Z), length(pairs))
    for (m in seq_along(pairs)) Zpair[, m] <- Z[, pairs[[m]][1L]] * Z[, pairs[[m]][2L]]
    Fmat <- up$natpar1 %*% t(Z)
    if (length(pairs) > 0L) Fmat <- Fmat + up$natpar2 %*% t(Zpair)
    W <- exp(Fmat - apply(Fmat, 1L, max)); W <- W / rowSums(W)
    out <- lapply(seq_len(S), function(s)
      matrix(NA_integer_, model$n_sites, model$max_visits))
    for (i in seq_len(model$n_sites)) {
      st <- sample.int(nrow(Z), 1L, prob = W[i, ])
      z  <- Z[st, ]
      for (s in seq_len(S)) {
        vj <- which(valid[[s]][i, ])
        if (!length(vj)) next
        out[[s]][i, vj] <- if (z[s] == 1L)
          stats::rbinom(length(vj), 1L, up$p_mat[i, s]) else 0L
      }
    }
    names(out) <- model$species
    out
  }
  if (nsim == 1L) return(draw_one())
  lapply(seq_len(nsim), function(s) draw_one())
}

# ---------------------------------------------------------------------------
# Simulator for recovery tests
# ---------------------------------------------------------------------------

#' Simulate a multi-species co-occurrence occupancy data set
#'
#' Draws from the [occu_multi()] model: a joint occupancy state for `S` species
#' from the log-linear (first + second order) occupancy model, then per-species
#' site-level detections at the observed visit pattern.
#'
#' @param S Number of species (default 2).
#' @param N Number of sites (default 300).
#' @param J Number of replicate visits (default 4).
#' @param n_state_covs Number of shared occupancy covariates (default 1).
#' @param beta_first Length-`S` list of first-order natural-parameter
#'   coefficients `c(intercept, slopes...)`. Default moderate occupancy.
#' @param beta_second Length-`choose(S, 2)` list of second-order (interaction)
#'   coefficients (in `combn(S, 2)` order). Default a single positive
#'   interaction for `S = 2`, else 0.
#' @param beta_p Length-`S` list of detection coefficients (logit). Default
#'   `c(qlogis(0.5), ...)`.
#' @param seed Optional random seed.
#' @return A list with `y` (a length-`S` list of `N x J` matrices), `data`,
#'   `species`, and `truth`.
#' @export
simulate_occu_multi <- function(S = 2, N = 300, J = 4, n_state_covs = 1,
                                beta_first = NULL, beta_second = NULL,
                                beta_p = NULL, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  pairs <- .occu_multi_pairs(S)
  if (is.null(beta_first))
    beta_first <- lapply(seq_len(S), function(s)
      c(stats::qlogis(0.45), stats::runif(n_state_covs, -0.4, 0.4)))
  if (is.null(beta_second))
    beta_second <- lapply(seq_along(pairs), function(m)
      c(if (S == 2L) 1.0 else 0.0, rep(0, n_state_covs)))
  if (is.null(beta_p))
    beta_p <- lapply(seq_len(S), function(s) c(stats::qlogis(0.5)))

  if (n_state_covs > 0L) {
    state_covs <- data.frame(matrix(stats::rnorm(N * n_state_covs), N, n_state_covs))
    names(state_covs) <- paste0("scov", seq_len(n_state_covs))
    data    <- state_covs
    X_state <- stats::model.matrix(~ ., state_covs)
  } else {
    data    <- data.frame(row.names = seq_len(N))
    X_state <- stats::model.matrix(~ 1, data.frame(row.names = seq_len(N)))
  }
  X_det   <- stats::model.matrix(~ 1, data.frame(row.names = seq_len(N)))

  natpar1 <- vapply(seq_len(S), function(s)
    as.vector(X_state %*% beta_first[[s]]), numeric(N))            # [N x S]
  natpar2 <- if (length(pairs))
    vapply(seq_along(pairs), function(m)
      as.vector(X_state %*% beta_second[[m]]), numeric(N)) else matrix(0, N, 0)
  p_mat <- vapply(seq_len(S), function(s)
    plogis(as.vector(X_det %*% beta_p[[s]])), numeric(N))          # [N x S]

  Z <- .occu_multi_states(S)
  Zpair <- matrix(0, nrow(Z), length(pairs))
  for (m in seq_along(pairs)) Zpair[, m] <- Z[, pairs[[m]][1L]] * Z[, pairs[[m]][2L]]

  species <- paste0("sp", seq_len(S))
  y <- lapply(seq_len(S), function(s) matrix(0L, N, J))
  z_all <- matrix(0L, N, S)
  for (i in seq_len(N)) {
    f <- as.numeric(natpar1[i, ] %*% t(Z))
    if (length(pairs)) f <- f + as.numeric(natpar2[i, ] %*% t(Zpair))
    w  <- exp(f - max(f)); w <- w / sum(w)
    st <- sample.int(nrow(Z), 1L, prob = w)
    z  <- Z[st, ]; z_all[i, ] <- z
    for (s in seq_len(S)) {
      y[[s]][i, ] <- if (z[s] == 1L) stats::rbinom(J, 1L, p_mat[i, s]) else 0L
    }
  }
  names(y) <- species
  list(y = y, data = data, species = species,
       truth = list(beta_first = beta_first, beta_second = beta_second,
                    beta_p = beta_p, z = z_all, pairs = pairs))
}
