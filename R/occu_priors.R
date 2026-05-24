# =============================================================================
# Weakly-informative priors on fixed-effect coefficients for the EM-Laplace
# engine. See R/penalized_irls.R for the math.
#
# Defaults (2026-05-15):
#   * `p_(Intercept)`          ~ Normal(0, 1.5)
#   * `p` slopes               ~ Normal(0, 2.5)
#   * `psi_(Intercept)`        ~ Normal(0, 2)
#   * `psi` slopes             ~ Normal(0, 5)
#
# Rationale. INLAabun D1 (N=600, J=6, 30 seeds) showed `engine='laplace'`
# had a ~1.5-unit logit bias on the psi intercept while NUTS recovered the
# same fixture cleanly. The MAP sits on the psi-p identifiability ridge that
# NUTS escapes through its priors; the unpenalised MAP doesn't. The
# detection-arm prior is the load-bearing one: a Normal(0, 1.5) intercept
# is wide enough to be uninformative across detection probabilities ~5 - 95%
# but informative enough to break the ridge. The occupancy-arm prior is
# wider (sd = 2 / 5) so it barely affects N >= 600 fits but prevents the
# same identifiability disaster in even smaller-N regimes.
#
# Set any `sd` to `Inf` to disable that prior. This is the same single
# objective: the penalty term is `(beta - mu)^2 / (2 * Inf^2) = 0`.
# =============================================================================


#' Weakly-informative priors for occupancy Laplace fits
#'
#' Constructs a prior specification consumed by [tobs()] when a Laplace
#' method (`method = "laplace"` etc.) is used for the occupancy family. The
#' penalty is applied as
#' a quadratic term `+ sum((beta_j - mu_j)^2 / (2 * sd_j^2))` added to the
#' negative-log-posterior in each M-step submodel block. See
#' `R/penalized_irls.R` for the full derivation.
#'
#' The defaults are weakly informative: they pull the detection intercept
#' toward `p = 0.5` strongly enough to break the psi-p identifiability
#' ridge at small `J`, but stay wide enough not to bias estimates at
#' N >= 600 with informative covariates. Setting any `sd` to `Inf` disables
#' that component of the prior (same penalised objective, just `1/Inf^2 = 0`).
#'
#' @param p_intercept Prior on the detection intercept (logit scale).
#'   A list `list(mean = numeric(1), sd = numeric(1))`.
#'   Default `list(mean = 0, sd = 1.5)`.
#' @param p_slope Prior on detection slopes (every non-intercept detection
#'   coefficient). Default `list(mean = 0, sd = 2.5)`.
#' @param beta_occ_intercept Prior on the occupancy (psi / psi1) intercept.
#'   Default `list(mean = 0, sd = 2)`.
#' @param beta_occ_slope Prior on occupancy slopes. Default
#'   `list(mean = 0, sd = 5)`.
#' @return An `occu_priors` object, ready to pass to `tobs(..., priors = ...)`.
#' @examples
#' \dontrun{
#' # default weakly-informative priors
#' fit <- tobs(~ x, data = d, family = occu(), detection = ~ z, y = y,
#'             method = "laplace", priors = occu_priors())
#'
#' # disable detection-slope penalty
#' priors <- occu_priors(p_slope = list(mean = 0, sd = Inf))
#' }
#' @export
occu_priors <- function(p_intercept       = list(mean = 0, sd = 1.5),
                        p_slope           = list(mean = 0, sd = 2.5),
                        beta_occ_intercept = list(mean = 0, sd = 2),
                        beta_occ_slope    = list(mean = 0, sd = 5)) {
  .check_one <- function(p, nm) {
    if (is.null(p)) return(NULL)
    if (!is.list(p) || is.null(p$mean) || is.null(p$sd)) {
      stop(sprintf("`%s` must be a list with `mean` and `sd` (or NULL).", nm),
           call. = FALSE)
    }
    if ((!is.numeric(p$mean) && !is.logical(p$mean)) ||
        (!is.numeric(p$sd)   && !is.logical(p$sd))) {
      stop(sprintf("`%s$mean` and `%s$sd` must be numeric.", nm, nm),
           call. = FALSE)
    }
    if (any(is.na(p$mean)) || any(!is.finite(suppressWarnings(as.numeric(p$mean))))) {
      stop(sprintf("`%s$mean` must be finite (got NA / -Inf / +Inf).", nm),
           call. = FALSE)
    }
    if (any(p$sd <= 0)) {
      stop(sprintf("`%s$sd` must be positive (Inf is allowed = no penalty).", nm),
           call. = FALSE)
    }
    p
  }
  out <- list(
    p_intercept        = .check_one(p_intercept,        "p_intercept"),
    p_slope            = .check_one(p_slope,            "p_slope"),
    beta_occ_intercept = .check_one(beta_occ_intercept, "beta_occ_intercept"),
    beta_occ_slope     = .check_one(beta_occ_slope,     "beta_occ_slope")
  )
  structure(out, class = c("occu_priors", "tobs_priors_spec"))
}


#' @export
print.occu_priors <- function(x, ...) {
  cat("occu_priors (weakly-informative for method='laplace'):\n")
  fmt <- function(p, label) {
    if (is.null(p)) {
      cat(sprintf("  %s: <unset>\n", label))
    } else {
      cat(sprintf("  %s ~ Normal(%.2g, %.2g)\n", label, p$mean, p$sd))
    }
  }
  fmt(x$p_intercept,        "p_(Intercept)")
  fmt(x$p_slope,            "p slopes      ")
  fmt(x$beta_occ_intercept, "psi_(Intercept)")
  fmt(x$beta_occ_slope,     "psi slopes    ")
  invisible(x)
}


# Resolve a user-supplied prior into the internal prior_spec list used by
# .attach_priors_to_blocks(). Accepts:
#   * NULL   -> default occu_priors()
#   * an `occu_priors` object -> used as-is
#   * a plain list with the same field names -> coerced via occu_priors()
#
# Returns NULL to mean "no penalty" only if the user explicitly passes
# `FALSE` or `"none"` — this is the escape hatch for tests that want
# the historical unpenalised MAP behavior.
.resolve_occu_priors <- function(priors) {
  if (identical(priors, FALSE) || identical(priors, "none")) {
    return(NULL)
  }
  if (is.null(priors)) {
    return(occu_priors())
  }
  if (inherits(priors, "occu_priors")) {
    return(priors)
  }
  if (is.list(priors)) {
    args <- priors[intersect(names(priors),
                             c("p_intercept", "p_slope",
                               "beta_occ_intercept", "beta_occ_slope"))]
    return(do.call(occu_priors, args))
  }
  stop("`priors` must be NULL, FALSE, an `occu_priors` object, or a named list.",
       call. = FALSE)
}
