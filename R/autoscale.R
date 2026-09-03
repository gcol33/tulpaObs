# =============================================================================
# autoscale.R - Internal design-matrix autoscaling for tobs() engines
#
# Numeric design columns on a large additive scale (calendar year ~ 2000,
# elevation in m, area in m^2, ...) form near-singular pairs with the
# intercept in the MAP Hessian. The optimizer can satisfy the data with a
# large negative intercept compensated by a small slope along the
# high-magnitude column. Coefficients are mathematically consistent with the
# prediction but numerically nonsensical, and any user-specified prior on the
# intercept is silently retargeted.
#
# Fix: center+scale numeric columns inside the engine before optimization,
# then transform betas / SEs / draws / Hessians back to the user-facing
# natural scale before returning. Mirrors glmnet / lme4's internal
# standardization and is invisible to the user.
# =============================================================================


# Decide which columns of a design matrix to standardize.
# Skip the intercept, skip binary indicators / factor contrasts (<=2 unique
# values), skip columns with sd == 0. Returns 1-based column indices.
.autoscale_pickcols <- function(X) {
  if (is.null(X) || ncol(X) == 0L) return(integer(0))
  cn <- colnames(X)
  intercept_cols <- if (!is.null(cn)) which(cn == "(Intercept)") else integer(0)
  candidates <- setdiff(seq_len(ncol(X)), intercept_cols)
  keep <- integer(0)
  for (j in candidates) {
    col <- X[, j]
    if (!is.numeric(col)) next
    finite <- col[is.finite(col)]
    if (length(finite) == 0L) next
    u <- unique(finite)
    if (length(u) <= 2L) next
    s <- stats::sd(col)
    if (!is.finite(s) || s <= 0) next
    keep <- c(keep, j)
  }
  keep
}

# Compute mean/sd of the columns flagged by .autoscale_pickcols. Returns a
# `scale` list that .apply_scale / .unscale_* consume. `intercept_col` is the
# 1-based intercept index (NA_integer_ if no intercept). `cols`, `means`,
# `sds` are aligned vectors.
.scale_meta <- function(X) {
  if (is.null(X) || ncol(X) == 0L) {
    return(list(cols = integer(0), means = numeric(0), sds = numeric(0),
                intercept_col = NA_integer_, ncol = if (is.null(X)) 0L else ncol(X)))
  }
  cn <- colnames(X)
  intercept_col <- if (!is.null(cn)) {
    ix <- which(cn == "(Intercept)")
    if (length(ix) == 0L) NA_integer_ else ix[1L]
  } else NA_integer_
  cols <- .autoscale_pickcols(X)
  if (length(cols) == 0L) {
    return(list(cols = integer(0), means = numeric(0), sds = numeric(0),
                intercept_col = intercept_col, ncol = ncol(X)))
  }
  means <- vapply(cols, function(j) mean(X[, j]), numeric(1))
  sds   <- vapply(cols, function(j) stats::sd(X[, j]), numeric(1))
  list(cols = cols, means = means, sds = sds,
       intercept_col = intercept_col, ncol = ncol(X))
}

# Apply a scale spec to a (possibly new) design matrix X. Used both at fit
# time (on the training X) and could be reused for predict() at new X with
# the cached mean/sd. Columns indexed by `scale$cols` are replaced by
# `(X[, j] - mean_j) / sd_j`. Other columns (intercept, factor contrasts,
# binary indicators) are untouched.
.apply_scale_to_X <- function(X, scale) {
  if (is.null(scale) || length(scale$cols) == 0L) return(X)
  for (k in seq_along(scale$cols)) {
    j <- scale$cols[k]
    X[, j] <- (X[, j] - scale$means[k]) / scale$sds[k]
  }
  X
}

# Scale + return the (X, scale) pair in one shot. Convenience wrapper.
.autoscale_design <- function(X) {
  scale <- .scale_meta(X)
  list(X = .apply_scale_to_X(X, scale), scale = scale)
}

# Linear transform matrix T such that beta_natural = T %*% beta_scaled.
# Identity except:
#   T[j, j]            = 1 / sd_j         for each scaled slope j
#   T[intercept, j]    = -mean_j / sd_j   for each scaled slope j
# Non-scaled non-intercept columns get T[k, k] = 1.
.scale_transform <- function(scale) {
  p <- as.integer(scale$ncol)
  T <- diag(1, p)
  if (length(scale$cols) == 0L) return(T)
  i0 <- scale$intercept_col
  for (k in seq_along(scale$cols)) {
    j <- scale$cols[k]
    T[j, j] <- 1 / scale$sds[k]
    if (!is.na(i0)) {
      T[i0, j] <- -scale$means[k] / scale$sds[k]
    }
  }
  T
}

# Transform a vector of betas from scaled space to natural space.
.unscale_beta_vec <- function(beta_sc, scale) {
  if (length(scale$cols) == 0L) return(beta_sc)
  as.numeric(.scale_transform(scale) %*% as.numeric(beta_sc))
}

# Inverse of `.unscale_beta_vec`: map natural-scale betas back into the
# centered+scaled parameterization. Used for round-trip tests and for
# evaluating arm log-likelihoods (`.loglik_cover_*`) at user-facing
# natural-scale coefficients when the encoded design X is scaled.
.scale_beta_vec <- function(beta_nat, scale) {
  if (length(scale$cols) == 0L) return(beta_nat)
  beta_sc <- beta_nat
  i0 <- scale$intercept_col
  for (k in seq_along(scale$cols)) {
    j <- scale$cols[k]
    beta_sc[j] <- scale$sds[k] * beta_nat[j]
  }
  if (!is.na(i0)) {
    shift <- 0
    for (k in seq_along(scale$cols)) {
      j <- scale$cols[k]
      shift <- shift + scale$means[k] * beta_nat[j]
    }
    beta_sc[i0] <- beta_nat[i0] + shift
  }
  beta_sc
}

# Transform a SINGLE scale-block's scaled-space variance-covariance matrix to
# natural space: V_natural = T %*% V_sc %*% t(T). (The multi-process,
# block-diagonal generalisation is `.unscale_vcov()` in R/occu_fit.R; this
# single-block helper is what the cover-hurdle decode applies per arm.)
.unscale_vcov_block <- function(V_sc, scale) {
  if (is.null(V_sc) || length(scale$cols) == 0L) return(V_sc)
  T <- .scale_transform(scale)
  T %*% V_sc %*% t(T)
}

# Transform per-coefficient SEs from scaled space to natural space, under
# the diagonal-Var approximation (Var(beta_sc) ~ diag(sds_sc^2)). This is
# what `build_laplace_fit()` uses today (independent normal pseudo-draws),
# so propagating with the same approximation keeps the natural-scale SEs
# internally consistent with the natural-scale draws.
#
# For scaled slope j: sds_natural[j] = sds_sc[j] / sd_j.
# For the intercept:  Var_natural[0] = sds_sc[0]^2
#                                    + sum_k (mean_k/sd_k)^2 * sds_sc[k]^2.
# All other columns are unchanged.
.unscale_sds_vec <- function(sds_sc, scale) {
  if (length(scale$cols) == 0L) return(sds_sc)
  sds_nat <- sds_sc
  i0 <- scale$intercept_col
  add_var_intercept <- 0
  for (k in seq_along(scale$cols)) {
    j <- scale$cols[k]
    sds_nat[j] <- sds_sc[j] / scale$sds[k]
    add_var_intercept <- add_var_intercept +
      (scale$means[k] / scale$sds[k])^2 * sds_sc[j]^2
  }
  if (!is.na(i0)) {
    sds_nat[i0] <- sqrt(sds_sc[i0]^2 + add_var_intercept)
  }
  sds_nat
}

# Marginal SEs from a covariance block, `NA` when the block is unavailable.
# `p` is the width to report, so a caller that could not invert its Hessian
# still gets a vector of the right length.
.sd_from_vcov <- function(V, p) {
  if (is.null(V)) return(rep(NA_real_, p))
  .tobs_sds_from_vcov(V)
}

# Transform an `n_draws x p` matrix of pseudo-draws (rows = posterior
# samples, columns = coefficients) from scaled space to natural space.
.unscale_draws_mat <- function(draws_sc, scale) {
  if (is.null(draws_sc) || length(scale$cols) == 0L) return(draws_sc)
  T <- .scale_transform(scale)
  out <- draws_sc %*% t(T)
  colnames(out) <- colnames(draws_sc)
  out
}


# ---------------------------------------------------------------------------
# Multi-process helpers for the occu / dyn_occu / ms_occu / int_occu / jsdm
# pipeline. These iterate over the per-process design matrices in
# `model$X_processes` and apply the scaling independently for each.
# ---------------------------------------------------------------------------

# Build a shallow copy of `model` with each `X_processes[[k]]` centered and
# scaled by its own (mean, sd) cache, and return the list of per-process
# scale specs alongside. The original model is unmodified.
.autoscale_model_X <- function(model) {
  if (is.null(model$X_processes) || length(model$X_processes) == 0L) {
    return(list(model = model, scales = NULL))
  }
  X_list <- vector("list", length(model$X_processes))
  scales <- vector("list", length(model$X_processes))
  for (k in seq_along(model$X_processes)) {
    sc <- .autoscale_design(model$X_processes[[k]])
    X_list[[k]] <- sc$X
    scales[[k]] <- sc$scale
  }
  scaled <- model
  scaled$X_processes <- X_list
  list(model = scaled, scales = scales)
}

# Walk a fit's `means` / `sds` / `draws` and unscale the per-process slices
# in place. Process layout follows `process_info` (one block per element,
# `p` columns each). Any tail beyond the process_info coefficients (e.g.
# `p_visit_*` columns from `X_det_visit`) is left untouched -- those
# matrices aren't passed through `.autoscale_model_X()` and so are still
# on natural scale.
.unscale_fit_per_process <- function(fit, scales, process_info) {
  if (is.null(scales) || is.null(process_info)) return(fit)
  if (is.null(fit$means)) return(fit)
  off <- 0L
  for (k in seq_along(process_info)) {
    p_k <- as.integer(process_info[[k]]$p)
    if (p_k == 0L) next
    sc <- scales[[k]]
    if (!is.null(sc) && length(sc$cols) > 0L) {
      idx <- off + seq_len(p_k)
      mu_sc <- as.numeric(fit$means[idx])
      sd_sc <- if (!is.null(fit$sds)) as.numeric(fit$sds[idx]) else NULL
      fit$means[idx] <- .unscale_beta_vec(mu_sc, sc)
      if (!is.null(sd_sc)) {
        fit$sds[idx] <- .unscale_sds_vec(sd_sc, sc)
      }
      if (!is.null(fit$draws) && ncol(fit$draws) >= max(idx)) {
        fit$draws[, idx] <- .unscale_draws_mat(
          fit$draws[, idx, drop = FALSE], sc
        )
      }
    }
    off <- off + p_k
  }
  fit
}
