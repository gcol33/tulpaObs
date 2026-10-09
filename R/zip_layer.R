# zip_layer.R - the structural-zero (ZIP / ZINB) additive layer over a per-site
# marginal, with optional analytic per-arm gradients. A structural-zero site is
# never occupied, so every count at it is zero; the observed per-site marginal is
#   L_i = omega * 1{all y_i = 0} + (1 - omega) * L_i^fam,
# with L_i^fam the family's own exact marginal (the Royle sum for abun(), the
# Dail-Madsen forward HMM for dyn_abun(), HMM + band multinomials for
# distsamp_open()). omega is an intercept-only structural-zero probability on
# the logit scale, parameter `zi_logit` at position `izi` of the family's theta
# vector.
#
# The family supplies one or two closures:
#   eval_site(theta) -> NULL on a failed evaluation, else a list carrying `llr`
#     (the per-site log marginal, length(az)) plus whatever its gradient needs;
#   grad_arms(theta, ev, w) -> numeric(length(theta)) with the family's own arm
#     entries filled from the eta gradients in `ev`, each weighted by the
#     posterior "not-a-structural-zero" weight w_i (the arms enter L_i only
#     through L_i^fam). NULL when the family optimises on the likelihood alone;
#     `neg_grad` is then NULL too.
# The layer owns the mixture, the ZI-logit score (1 - om) - w_i summed over
# sites, the NB dispersion score (the kernel returns that one summed across
# sites, so it cannot be ZIP-weighted analytically and is central-differenced
# on the exact ZIP objective, two extra marginal evaluations), and the seed for
# the structural-zero logit: `zi_seed[1]` times the all-zero share of sites,
# clamped to [`zi_seed[2]`, `zi_seed[3]`].
.tobs_zip_layer <- function(eval_site, grad_arms = NULL, az, izi, ir = NA_integer_,
                            zi_seed = c(0.5, 0.05, 0.7)) {
  is_nb   <- !is.na(ir)
  n_sites <- length(az)

  # ZIP marginal log-lik and the posterior weight w_i on the family component: a
  # detected site is certainly not a structural zero (w = 1); an all-zero site
  # mixes the structural point mass in, w_i = (1 - om) L_fam_i / L_i.
  pieces <- function(theta, llr) {
    om <- stats::plogis(theta[izi]); log1m <- log1p(-om)
    ll <- numeric(n_sites); w <- rep(1, n_sites)
    ll[!az] <- log1m + llr[!az]
    a <- log1m + llr[az]; b <- log(om); mx <- pmax(a, b)
    Li <- mx + log(exp(a - mx) + exp(b - mx))
    ll[az] <- Li
    w[az] <- exp(a - Li)
    list(log_lik = sum(ll), w = w, om = om)
  }

  neg_ll <- function(theta) {
    ev <- eval_site(theta)
    if (is.null(ev) || any(!is.finite(ev$llr))) return(1e10)
    val <- -pieces(theta, ev$llr)$log_lik
    if (is.finite(val)) val else 1e10
  }

  neg_grad <- if (is.null(grad_arms)) NULL else function(theta) {
    ev <- eval_site(theta)
    if (is.null(ev)) return(rep(0, length(theta)))
    zp <- pieces(theta, ev$llr); w <- zp$w; om <- zp$om
    g <- grad_arms(theta, ev, w)
    g[izi] <- sum((1 - om) - w)
    if (is_nb) {
      h <- 1e-4; th <- theta
      th[ir] <- theta[ir] + h; fp <- -neg_ll(th)
      th[ir] <- theta[ir] - h; fm <- -neg_ll(th)
      g[ir] <- (fp - fm) / (2 * h)
    }
    -g
  }

  list(neg_ll = neg_ll, neg_grad = neg_grad,
       zi_logit0 = stats::qlogis(min(max(mean(az) * zi_seed[1L], zi_seed[2L]),
                                     zi_seed[3L])))
}
