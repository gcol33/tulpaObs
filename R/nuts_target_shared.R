# nuts_target_shared.R - shared pieces for the single-arm-vector marginal NUTS
# families (dyn_abun, fp_occu, ...).
#
# Each such family has a flat coefficient vector split into k contiguous arm
# blocks and a joint log-posterior that is the marginal log-likelihood (from the
# family's C++ kernel) plus a weak isotropic Gaussian prior. The per-family
# wrappers supply only their arm names, design matrices, and the matching
# `grad_eta_*` field selectors. The C++ FullGradFn for each family is the actual
# sampler target; these R functions are the cross-check oracle.

# Contiguous coefficient blocks for `sizes` (in flat-vector order). Returns the
# per-arm index vectors named by `arms`.
.tobs_nuts_arm_idx <- function(arms, sizes) {
  off <- cumsum(c(0L, sizes[-length(sizes)]))
  stats::setNames(Map(function(o, n) o + seq_len(n), off, sizes), arms)
}

# Joint log-posterior for a uniform k-arm marginal target. `eval_out` is the
# kernel return carrying `$log_lik` plus one `grad_eta_*` field per arm; `arms`
# lists, in flat-vector order, each arm's theta `idx`, design matrix `X`, and the
# name of its `grad` field in `eval_out`. Weak Gaussian prior N(0, sigma.beta^2)
# on every coefficient.
.tobs_nuts_logpost_k <- function(theta, eval_out, arms, total, sigma.beta = 10) {
  grad <- numeric(total)
  for (a in arms) grad[a$idx] <- as.numeric(crossprod(a$X, eval_out[[a$grad]]))
  ib2 <- 1 / sigma.beta^2
  list(lp   = eval_out$log_lik - 0.5 * ib2 * sum(theta^2),
       grad = grad - ib2 * theta)
}
