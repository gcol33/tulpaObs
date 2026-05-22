# =============================================================================
# Local EM + Laplace driver for tulpaObs, with coefficient priors.
#
# tulpa::tulpa_em_laplace() does not currently route priors on fixed effects
# through to its inner tulpa_laplace() fits, so for the prior-aware path we
# run a tulpaObs-local EM driver that calls .fit_block_penalized() per block.
# The structure is intentionally a near-copy of tulpa_em_laplace() so that
# any future upstream prior support drops in cleanly.
#
# All convergence semantics, weight-damping rules, and return shape are
# preserved so build_laplace_fit() and the EM correction post-step
# (MI / Gibbs) continue to work without changes downstream.
# =============================================================================


# Helpers cloned in spirit from tulpa:::.fits_to_param_vec / .max_rel_change.
# Kept local so this driver compiles on a tulpa that doesn't export them.
.fits_to_param_vec_local <- function(fits) {
  parts <- lapply(fits, function(f) {
    if (is.null(f) || is.null(f$mode)) return(numeric(0))
    as.numeric(f$mode)
  })
  unlist(parts, use.names = FALSE)
}

.max_rel_change_local <- function(prev, curr) {
  if (length(prev) == 0 || length(curr) == 0 ||
      length(prev) != length(curr)) return(Inf)
  denom <- pmax(abs(prev), 1e-8)
  max(abs(curr - prev) / denom, na.rm = TRUE)
}


# Attach a per-submodel prior to every block returned by m_step_encode.
# The encoded blocks already know their submodel name via `names(blocks)`
# (set by the callbacks: c("occ", "det"), c("occ", "det", "col", "ext"), ...).
# We map those names to the user-facing prior buckets via
# .prior_for_submodel() in penalized_irls.R.
.attach_priors_to_blocks <- function(blocks, model, prior_spec) {
  if (is.null(prior_spec)) return(blocks)
  pi_list <- model$process_info
  # Build a name -> process_info mapping. The submodel-block keys in the
  # callbacks are ("occ", "det", "col", "ext", "det1", "det2", ...) while
  # process_info names are ("psi", "p", "gamma", "epsilon", ...) or
  # ("psi", source_name_1, source_name_2, ...) for integrated.
  block_to_pi <- function(block_name, k) {
    # Single-season / community: block "occ" -> pi "psi", "det" -> pi "p"
    if (block_name == "occ" && length(pi_list) >= 1L) return(pi_list[[1]])
    if (block_name == "det" && length(pi_list) >= 2L) return(pi_list[[2]])
    # Dynamic: "occ"->psi1, "det"->p, "col"->gamma, "ext"->epsilon (indices 1..4)
    name_to_idx <- c(occ = 1L, det = 2L, col = 3L, ext = 4L)
    if (block_name %in% names(name_to_idx)) {
      idx <- name_to_idx[[block_name]]
      if (idx <= length(pi_list)) return(pi_list[[idx]])
    }
    # Integrated: "det1", "det2", ... map to pi 2, 3, ...
    if (grepl("^det[0-9]+$", block_name)) {
      idx <- as.integer(sub("^det", "", block_name)) + 1L
      if (idx <= length(pi_list)) return(pi_list[[idx]])
    }
    # Fall back to positional.
    if (k <= length(pi_list)) return(pi_list[[k]])
    NULL
  }

  bn <- names(blocks); if (is.null(bn)) bn <- rep("", length(blocks))
  has_visit <- !is.null(model$det_visit_names) &&
               length(model$det_visit_names) > 0L
  for (k in seq_along(blocks)) {
    pi <- block_to_pi(bn[k], k)
    if (is.null(pi)) next
    # Map block name -> submodel-prior bucket via the process_info name.
    pr <- .prior_for_submodel(prior_spec, pi$name, pi$coef_names)
    # Detection block carries visit-level slopes on the tail of its X matrix;
    # extend the prior to match. Visit-level coefs always route through the
    # `p_slope` bucket (they are slopes by construction, not intercepts).
    if (has_visit && identical(bn[k], "det")) {
      pr_visit <- .prior_for_submodel(
        prior_spec, "p",
        coef_names = paste0("visit_", model$det_visit_names)
      )
      pr <- list(
        mean = c(pr$mean, pr_visit$mean),
        sd   = c(pr$sd,   pr_visit$sd)
      )
    }
    if (is.null(pr) || all(!is.finite(pr$sd))) next   # no penalty
    blocks[[k]]$prior <- list(mean = pr$mean, sd = pr$sd)
  }
  blocks
}


# Local EM + Laplace driver (no MI / Gibbs correction in this variant: the
# point estimate plus the penalised Hessian SE is what build_laplace_fit()
# consumes for tobs occupancy).
.tobs_em_laplace_penalized <- function(model, callbacks, prior_spec,
                                       max_iter = 50L, tol = 1e-4,
                                       damping = 0.3, verbose = TRUE) {
  e_step <- callbacks$e_step
  m_step_encode <- callbacks$m_step_encode

  fits <- list()
  weights <- NULL
  prev_params <- numeric(0)
  history <- data.frame(iter = integer(0), delta = double(0))
  converged <- FALSE
  iter <- 0L

  for (iter in seq_len(max_iter)) {
    # ---- E-step ----
    e_result <- e_step(fits)
    if (!is.list(e_result) || is.null(e_result$weights)) {
      stop(".tobs_em_laplace_penalized: e_step must return list with `weights`",
           call. = FALSE)
    }
    weights_new <- e_result$weights
    weights <- if (is.null(weights)) {
      weights_new
    } else {
      (1 - damping) * weights_new + damping * weights
    }

    # ---- M-step (with prior attachment) ----
    blocks <- m_step_encode(weights)
    if (!is.list(blocks) || length(blocks) == 0L) {
      stop(".tobs_em_laplace_penalized: m_step_encode returned no blocks",
           call. = FALSE)
    }
    block_names <- names(blocks)
    if (is.null(block_names)) block_names <- paste0("submodel_", seq_along(blocks))

    blocks <- .attach_priors_to_blocks(blocks, model, prior_spec)

    new_fits <- vector("list", length(blocks))
    names(new_fits) <- block_names
    for (k in seq_along(blocks)) {
      # Warm-start beta_init from previous iteration's mode (if shapes match).
      if (!is.null(fits[[block_names[k]]]) &&
          !is.null(fits[[block_names[k]]]$mode)) {
        prev_mode <- fits[[block_names[k]]]$mode
        if (length(prev_mode) >= ncol(blocks[[k]]$X)) {
          blocks[[k]]$beta_init <- prev_mode[seq_len(ncol(blocks[[k]]$X))]
        }
      }
      new_fits[[k]] <- .fit_block_penalized(blocks[[k]])
    }
    fits <- new_fits

    # ---- Convergence ----
    curr_params <- .fits_to_param_vec_local(fits)
    delta <- .max_rel_change_local(prev_params, curr_params)
    prev_params <- curr_params
    history <- rbind(history, data.frame(iter = iter, delta = delta))
    if (verbose) cat(sprintf("  EM iter %d: delta = %.6g\n", iter, delta))
    if (is.finite(delta) && delta < tol) {
      converged <- TRUE
      if (verbose) cat(sprintf("  EM converged after %d iterations\n", iter))
      break
    }
  }

  list(
    fits       = fits,
    weights    = weights,
    n_iter     = iter,
    converged  = converged,
    convergence = list(converged = converged, n_iter = iter, history = history),
    history    = history,
    correction = "none",
    pooled     = NULL
  )
}
