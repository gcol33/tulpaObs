# =============================================================================
# Weakly-informative priors on fixed-effect coefficients for the EM-Laplace
# engine.
#
# The penalty itself is implemented in tulpa: each M-step block carries a
# per-coefficient `beta_prior = list(mean, sd)` (attached by
# `.attach_priors_to_blocks()` below) that tulpa::tulpa_laplace() adds to the
# negative-log-posterior as `sum((beta_j - mu_j)^2 / (2 * sd_j^2))`. The same
# prior flows through tulpa's MI / Gibbs correction refits, so the penalised
# estimate and the stochastic-correction draws regularise identically. `sd_j
# = Inf` sets that coefficient's precision to 0 (no penalty).
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


# Validate one `list(mean, sd)` prior bucket (shared by occu_priors() and
# cover_priors()). Returns the bucket unchanged, or NULL for an unset bucket.
.check_prior_bucket <- function(p, nm) {
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


# Validate each named bucket and wrap the result in a priors-spec object of the
# given class. Shared assembly for occu_priors() and cover_priors(); each keeps
# its own argument surface and only hands the named buckets and class here.
.build_prior_spec <- function(buckets, class) {
  out <- Map(.check_prior_bucket, buckets, names(buckets))
  structure(out, class = c(class, "tobs_priors_spec"))
}


#' Weakly-informative priors for occupancy Laplace fits
#'
#' Constructs a prior specification consumed by [tobs()] when a Laplace
#' method (`method = "laplace"` etc.) is used for the occupancy family. The
#' penalty is applied as
#' a quadratic term `+ sum((beta_j - mu_j)^2 / (2 * sd_j^2))` added to the
#' negative-log-posterior in each M-step submodel block by
#' [tulpa::tulpa_laplace()] (via the per-block `beta_prior`).
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
  .build_prior_spec(
    list(p_intercept        = p_intercept,
         p_slope            = p_slope,
         beta_occ_intercept = beta_occ_intercept,
         beta_occ_slope     = beta_occ_slope),
    "occu_priors")
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


# =============================================================================
# Prior-spec -> per-submodel -> per-block `beta_prior` plumbing.
#
# These map the user-facing occu_priors() buckets onto the M-step submodel
# blocks the callbacks in R/laplace.R emit, so tulpa's penalized Laplace fits
# every block with the right per-coefficient prior. The penalty maths lives in
# tulpa (tulpa_laplace's `beta_prior`); here we only resolve which bucket
# applies to which coefficient.
# =============================================================================


# Normalise a prior block to a (mean vector, sd vector) of length p.
# `prior` is a list with $mean (scalar or length-p) and $sd (scalar or
# length-p), where sd may include +Inf entries (no penalty on that coef).
.expand_prior <- function(prior, p) {
  if (is.null(prior)) {
    return(list(mean = rep(0, p), sd = rep(Inf, p)))
  }
  m <- prior$mean
  s <- prior$sd
  if (is.null(m)) m <- 0
  if (is.null(s)) s <- Inf
  if (length(m) == 1L) m <- rep(m, p)
  if (length(s) == 1L) s <- rep(s, p)
  if (length(m) != p) {
    stop(sprintf(
      ".expand_prior(): prior$mean has length %d, expected %d", length(m), p),
      call. = FALSE)
  }
  if (length(s) != p) {
    stop(sprintf(
      ".expand_prior(): prior$sd has length %d, expected %d", length(s), p),
      call. = FALSE)
  }
  if (any(!is.finite(m))) {
    stop(".expand_prior(): prior$mean must be finite (got NA / +-Inf)",
         call. = FALSE)
  }
  if (any(s <= 0)) {
    stop(".expand_prior(): prior$sd must be positive (Inf is allowed for no penalty)",
         call. = FALSE)
  }
  list(mean = as.numeric(m), sd = as.numeric(s))
}


# Expand an (intercept-bucket, slope-bucket) pair into per-coefficient
# (mean, sd) vectors aligned to `coef_names`. The intercept column is the one
# named "(Intercept)" (model.matrix convention), falling back to the first
# column. A NULL bucket means "no penalty" (mean 0, sd Inf) for that role.
# Shared by the occupancy (.prior_for_submodel) and cover (.cover_arm_prior)
# prior builders so the intercept/slope split lives in one place.
.expand_int_slope_prior <- function(int_b, slp_b, coef_names) {
  p <- length(coef_names)
  if (p == 0L) return(list(mean = numeric(0), sd = numeric(0)))
  is_intercept <- coef_names == "(Intercept)"
  if (!any(is_intercept)) {
    is_intercept <- c(TRUE, rep(FALSE, p - 1L))[seq_len(p)]
  }
  mean_vec <- numeric(p)
  sd_vec   <- numeric(p)
  for (j in seq_len(p)) {
    b <- if (is_intercept[j]) int_b else slp_b
    mean_vec[j] <- if (is.null(b)) 0   else b$mean %||% 0
    sd_vec[j]   <- if (is.null(b)) Inf else b$sd   %||% Inf
  }
  list(mean = mean_vec, sd = sd_vec)
}


# Build a prior block for the occupancy / detection submodel given the
# user-facing tobs prior list and the coefficient names. The prior spec
# accepts four named buckets:
#   $p_intercept       - list(mean, sd) applied to the detection intercept
#   $p_slope           - list(mean, sd) applied to every non-intercept detection coef
#   $beta_occ_intercept - list(mean, sd) applied to the occupancy / psi intercept
#   $beta_occ_slope    - list(mean, sd) applied to every non-intercept occupancy coef
#
# The "intercept" column is detected by name `(Intercept)`, which is what
# model.matrix() uses; falls back to the first column otherwise.
.prior_for_submodel <- function(prior_spec, sub_name, coef_names) {
  p <- length(coef_names)
  if (p == 0L) return(.expand_prior(NULL, 0L))
  if (is.null(prior_spec)) return(.expand_prior(NULL, p))

  if (sub_name %in% c("p")) {
    int_b <- prior_spec$p_intercept
    slp_b <- prior_spec$p_slope
  } else if (sub_name %in% c("psi", "psi1")) {
    int_b <- prior_spec$beta_occ_intercept
    slp_b <- prior_spec$beta_occ_slope
  } else if (sub_name %in% c("gamma")) {
    int_b <- prior_spec$gamma_intercept %||% prior_spec$beta_occ_intercept
    slp_b <- prior_spec$gamma_slope %||% prior_spec$beta_occ_slope
  } else if (sub_name %in% c("epsilon")) {
    int_b <- prior_spec$epsilon_intercept %||% prior_spec$beta_occ_intercept
    slp_b <- prior_spec$epsilon_slope %||% prior_spec$beta_occ_slope
  } else if (grepl("^det[0-9]+$", sub_name) || sub_name == "det") {
    int_b <- prior_spec$p_intercept
    slp_b <- prior_spec$p_slope
  } else {
    int_b <- prior_spec$p_intercept
    slp_b <- prior_spec$p_slope
  }

  .expand_int_slope_prior(int_b, slp_b, coef_names)
}


# Attach a per-submodel prior to every block returned by m_step_encode as a
# per-block `beta_prior = list(mean, sd)`. tulpa's EM block fitter
# (.fit_block_via_laplace) reads `block$beta_prior` and threads it through
# tulpa_laplace() in every phase (EM iterations and the MI/Gibbs correction
# refits), so the penalised estimate and the corrected draws regularise the
# same way.
#
# The encoded blocks know their submodel name via `names(blocks)` (set by the
# callbacks: c("occ", "det"), c("occ", "det", "col", "ext"), ...); we map those
# to the user-facing prior buckets via .prior_for_submodel(). Spatial blocks
# are skipped: tulpa_laplace() rejects beta_prior on the SPDE/NNGP path (which
# carries its own fixed-effect prior), so they stay unpenalised.
.attach_priors_to_blocks <- function(blocks, model, prior_spec) {
  if (is.null(prior_spec)) return(blocks)
  pi_list <- model$process_info
  # Map a block key ("occ", "det", "col", "ext", "det1", ...) to its
  # process_info entry ("psi", "p", "gamma", "epsilon", source names, ...).
  block_to_pi <- function(block_name, k) {
    if (block_name == "occ" && length(pi_list) >= 1L) return(pi_list[[1]])
    if (block_name == "det" && length(pi_list) >= 2L) return(pi_list[[2]])
    name_to_idx <- c(occ = 1L, det = 2L, col = 3L, ext = 4L)
    if (block_name %in% names(name_to_idx)) {
      idx <- name_to_idx[[block_name]]
      if (idx <= length(pi_list)) return(pi_list[[idx]])
    }
    if (grepl("^det[0-9]+$", block_name)) {
      idx <- as.integer(sub("^det", "", block_name)) + 1L
      if (idx <= length(pi_list)) return(pi_list[[idx]])
    }
    if (k <= length(pi_list)) return(pi_list[[k]])
    NULL
  }

  bn <- names(blocks); if (is.null(bn)) bn <- rep("", length(blocks))
  has_visit <- !is.null(model$det_visit_names) &&
               length(model$det_visit_names) > 0L
  for (k in seq_along(blocks)) {
    # Spatial blocks carry their own prior in the SPDE/NNGP solver; leave them
    # for the unpenalised spatial path (tulpa_laplace errors on beta_prior +
    # spatial). The non-spatial detection block of a spatial fit is also left
    # alone because attachment is gated on `is.null(spatial)` upstream.
    if (!is.null(blocks[[k]]$spatial)) next
    pi <- block_to_pi(bn[k], k)
    if (is.null(pi)) next
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
    if (length(pr$sd) == 0L) next
    blocks[[k]]$beta_prior <- list(mean = pr$mean, sd = pr$sd)
  }
  blocks
}


# =============================================================================
# Cover-hurdle fixed-effect priors.
#
# The cover hurdle has two linear predictors with their own fixed-effect
# coefficients: an occurrence arm (Bernoulli presence) and a positive-cover arm
# (beta or lognormal). Unlike occupancy, the two arms are separately identified
# (presence from 0/>0, positive cover from the positive values), so there is no
# psi-p-style ridge that demands a default prior. Cover priors are therefore
# OPT-IN: `tobs(..., family = cover(), priors = NULL)` fits unpenalised (the
# historical behaviour the recovery tests assume). Passing `cover_priors()`
# adds the same quadratic `beta_prior` penalty tulpa_laplace() applies on the
# occupancy path, mainly to tame perfect separation in the occurrence arm at
# small N. The prior is specified on natural-scale coefficients and applied on
# the autoscaled design, matching the occupancy convention.
# =============================================================================


#' Weakly-informative priors for cover-hurdle Laplace fits
#'
#' Builds an opt-in fixed-effect prior for [cover()] models fit with a Laplace
#' method. The cover hurdle has two arms with their own coefficients -- an
#' occurrence arm (Bernoulli presence) and a positive-cover arm (beta or
#' lognormal) -- each with an intercept and a slope bucket. The penalty is the
#' quadratic `+ sum((beta_j - mu_j)^2 / (2 sd_j^2))` added to each arm's
#' negative log-likelihood, the same mechanism [occu_priors()] uses.
#'
#' Cover priors are opt-in: without `priors`, a cover fit is unpenalised. The
#' main use is regularising perfect separation in the occurrence arm at small
#' `N`. Set any `sd = Inf` to drop that component.
#'
#' @section Coverage: Applied on the separate-Laplace path
#'   (`method = "laplace"` / `"laplace_sla"`) when the formula has no spatial
#'   term. Both arms are penalised: the occurrence and lognormal-positive arms
#'   through [tulpa::tulpa_laplace()], and the beta-positive arm through
#'   [tulpa::tulpa_laplace_beta()]. The joint nested-Laplace path
#'   (`method = "nested_laplace"` / `"nested_laplace_sla"`, used for spatial
#'   cover formulas) also threads the priors, as a per-arm `beta_prior_mean` /
#'   `beta_prior_prec` on the joint engine's occurrence and positive responses.
#'   The separate-Laplace path keeps rejecting a prior alongside a spatial term
#'   (that solver carries its own); add the spatial term through the
#'   nested-Laplace path to combine the two.
#'
#' @param occ_intercept Prior on the occurrence (presence) intercept, logit
#'   scale. `list(mean, sd)`. Default `list(mean = 0, sd = 2)`.
#' @param occ_slope Prior on occurrence slopes. Default
#'   `list(mean = 0, sd = 2.5)`.
#' @param pos_intercept Prior on the positive-cover intercept (logit scale for
#'   `"beta"`, log scale for `"lognormal"`). Default `list(mean = 0, sd = 3)`.
#' @param pos_slope Prior on positive-cover slopes. Default
#'   `list(mean = 0, sd = 2.5)`.
#' @return A `cover_priors` object, ready to pass to `tobs(..., priors = ...)`.
#' @seealso [occu_priors()]
#' @examples
#' \dontrun{
#' # regularise the occurrence arm, leave the positive arm unpenalised
#' fit <- tobs(~ x, data = d, family = cover(response = "lognormal"),
#'             y = cov, method = "laplace",
#'             priors = cover_priors(pos_intercept = list(mean = 0, sd = Inf),
#'                                   pos_slope     = list(mean = 0, sd = Inf)))
#' }
#' @export
cover_priors <- function(occ_intercept = list(mean = 0, sd = 2),
                         occ_slope     = list(mean = 0, sd = 2.5),
                         pos_intercept = list(mean = 0, sd = 3),
                         pos_slope     = list(mean = 0, sd = 2.5)) {
  .build_prior_spec(
    list(occ_intercept = occ_intercept,
         occ_slope     = occ_slope,
         pos_intercept = pos_intercept,
         pos_slope     = pos_slope),
    "cover_priors")
}


#' @export
print.cover_priors <- function(x, ...) {
  cat("cover_priors (opt-in, for cover() Laplace fits):\n")
  fmt <- function(p, label) {
    if (is.null(p)) cat(sprintf("  %s: <unset>\n", label))
    else cat(sprintf("  %s ~ Normal(%.2g, %.2g)\n", label, p$mean, p$sd))
  }
  fmt(x$occ_intercept, "occ_(Intercept)")
  fmt(x$occ_slope,     "occ slopes     ")
  fmt(x$pos_intercept, "pos_(Intercept)")
  fmt(x$pos_slope,     "pos slopes     ")
  invisible(x)
}


# Resolve a user `priors` argument for the cover hurdle. Cover priors are
# opt-in, so NULL / FALSE / "none" all mean "no penalty" (returns NULL). A
# `cover_priors` object is used as-is; a plain list with the cover bucket names
# is coerced via cover_priors(). An `occu_priors` object is rejected (the
# buckets differ: occ/pos vs psi/p).
.resolve_cover_priors <- function(priors) {
  if (is.null(priors) || identical(priors, FALSE) || identical(priors, "none")) {
    return(NULL)
  }
  if (inherits(priors, "cover_priors")) return(priors)
  if (inherits(priors, "occu_priors")) {
    stop("cover() takes cover_priors(), not occu_priors() ",
         "(the arms are occ/pos, not psi/p).", call. = FALSE)
  }
  if (is.list(priors)) {
    args <- priors[intersect(names(priors),
                             c("occ_intercept", "occ_slope",
                               "pos_intercept", "pos_slope"))]
    return(do.call(cover_priors, args))
  }
  stop("`priors` for cover() must be NULL, FALSE, a cover_priors object, ",
       "or a named list.", call. = FALSE)
}


# Per-coefficient (mean, sd) beta_prior for one cover arm, or NULL if the arm
# is unpenalised (no bucket active -> all sd Inf). `arm` is "occ" or "pos";
# `coef_names` are the arm's design-matrix column names.
.cover_arm_prior <- function(cprior, arm, coef_names) {
  if (is.null(cprior)) return(NULL)
  pr <- .expand_int_slope_prior(cprior[[paste0(arm, "_intercept")]],
                                cprior[[paste0(arm, "_slope")]],
                                coef_names)
  if (length(pr$sd) == 0L || all(!is.finite(pr$sd))) return(NULL)
  pr
}
