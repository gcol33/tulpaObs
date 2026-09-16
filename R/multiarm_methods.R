# ============================================================================
# S3 surface for multi-arm fits (cover_fit, occu_categorical_fit).
#
# These families store a fit as independent per-arm coefficient blocks (a
# presence arm plus a positive-cover or nominal-class arm) rather than the flat
# `$means` / `$draws` / `$model` layout that the inherited `tulpa_fit` and
# `tobs_fit` methods assume. Tagging them `tobs_multiarm_fit` (between the
# concrete class and `tobs_fit`) routes the generic extractors through this one
# arm-aware implementation instead of the flat-layout methods, which otherwise
# error on the missing structure.
#
# A new multi-arm family joins by (1) adding "tobs_multiarm_fit" to its class
# vector and (2) providing an `.tobs_arm_blocks()` branch. Nothing else.
#
# The arms are fit independently, so the joint coefficient covariance is exactly
# block-diagonal (no cross-arm term to drop). When a fit already carries a full
# posterior (`$draws`, the NUTS cover path), the inference methods defer to the
# draw-based `tulpa_fit` methods via NextMethod().
# ============================================================================

# --- per-family arm accessor -------------------------------------------------
# Returns an ordered named list of arms; each arm is
#   list(estimate = <named numeric>, se = <named numeric>, vcov = <matrix>).
# `estimate`/`se` names are the coefficient labels within the arm.
.tobs_arm_blocks <- function(object) {
  if (inherits(object, "occu_categorical_fit"))
    return(.tobs_arm_blocks_categorical(object))
  if (inherits(object, "cover_fit"))
    return(.tobs_arm_blocks_cover(object))
  stop("no arm-block accessor for class ", class(object)[1L], call. = FALSE)
}

# Exact per-arm covariance from a stored value when its dimension matches the
# coefficient vector, else the SE-diagonal (independent-coefficient) fallback.
.tobs_arm_vcov <- function(V, se) {
  p <- length(se)
  if (!is.null(V)) {
    V <- as.matrix(V)
    if (nrow(V) >= p && ncol(V) >= p) {
      Vb <- V[seq_len(p), seq_len(p), drop = FALSE]
      dimnames(Vb) <- list(names(se), names(se))
      return(Vb)
    }
  }
  Vb <- diag(se^2, nrow = p)
  dimnames(Vb) <- list(names(se), names(se))
  Vb
}

.tobs_arm_blocks_cover <- function(object) {
  b_occ <- object[["beta_occ"]]; s_occ <- object[["se_occ"]]
  b_pos <- object[["beta_pos"]]; s_pos <- object[["se_pos"]]
  V_occ <- tryCatch(solve(as.matrix(object[["occ"]][["H_beta"]])),
                    error = function(e) NULL)
  V_pos <- tryCatch(solve(as.matrix(object[["pos"]][["H_beta"]])),
                    error = function(e) NULL)
  list(
    presence = list(estimate = b_occ, se = s_occ,
                    vcov = .tobs_arm_vcov(V_occ, s_occ)),
    positive = list(estimate = b_pos, se = s_pos,
                    vcov = .tobs_arm_vcov(V_pos, s_pos))
  )
}

.tobs_arm_blocks_categorical <- function(object) {
  b_occ <- object[["beta_occ"]]; s_occ <- object[["se_occ"]]
  # Class arm: flatten the (p x (K-1)) matrix column-major to align with
  # `vcov_class` (its diagonal was reshaped column-major into `se_class`).
  Beta <- object[["beta_class"]]; SE <- object[["se_class"]]
  cls_labels <- colnames(Beta); coef_labels <- rownames(Beta)
  nm <- as.vector(outer(coef_labels, cls_labels,
                        function(co, cl) paste0(cl, ":", co)))
  b_cls <- as.vector(Beta); names(b_cls) <- nm
  s_cls <- as.vector(SE);   names(s_cls) <- nm
  list(
    presence = list(estimate = b_occ, se = s_occ,
                    vcov = .tobs_arm_vcov(object[["vcov_occ"]], s_occ)),
    class    = list(estimate = b_cls, se = s_cls,
                    vcov = .tobs_arm_vcov(object[["vcov_class"]], s_cls))
  )
}

# Joint log-likelihood field (cover: c(occ, pos); categorical: c(occ, class);
# spatial joint: c(joint)); arms are independent so the finite entries sum.
.tobs_multiarm_loglik_val <- function(object) {
  ll <- object[["log_marginal"]]
  if (is.null(ll)) ll <- object[["loglik"]]
  sum(ll[is.finite(ll)])
}

.tobs_multiarm_npar <- function(object) {
  p_occ  <- length(object[["beta_occ"]])
  p_arm2 <- length(object[["beta_pos"]])
  if (p_arm2 == 0L) p_arm2 <- length(as.vector(object[["beta_class"]]))
  n_disp <- sum(is.finite(c(object[["phi_pos"]], object[["sigma_pos"]])))
  as.integer(p_occ + p_arm2 + n_disp)
}

.tobs_multiarm_unsupported <- function(what, object) {
  fam <- object[["family"]][["name"]] %||% "multi-arm"
  stop(sprintf(
    "%s() is not defined for a %s fit. Use predict() for expected values on new data, or simulate_%s() to simulate a fresh data set.",
    what, fam, fam), call. = FALSE)
}

# --- number of observations --------------------------------------------------
#' S3 methods for multi-arm fits
#'
#' Coefficient, covariance, and summary accessors for fits that store
#' independent per-arm coefficient blocks (`cover_fit`, `occu_categorical_fit`)
#' rather than a flat posterior. The arms are fit independently, so the joint
#' covariance is block-diagonal. When a full posterior is present (the NUTS
#' cover path) the inference methods defer to the draw-based `tulpa_fit`
#' methods.
#'
#' @param object,x A `tobs_multiarm_fit` (a `cover_fit` or
#'   `occu_categorical_fit`).
#' @param parm Ignored (present for `confint()` generic compatibility).
#' @param level Confidence level for `confint()` and `summary()` (default 0.95).
#' @param ... Ignored, or forwarded to `NextMethod()` on the posterior path.
#' @param arm Optional arm name (`"presence"`, `"positive"` or `"class"`) for
#'   `coef()`; see [coef.tobs_fit()].
#' @return The layout every `tobs_fit` returns: `nobs()` an integer; `coef()` a
#'   named vector (`presence_(Intercept)`, `positive_x`, `class_2:x`, ...);
#'   `vcov()` the block-diagonal covariance and `confint()` a two-column matrix,
#'   both named as `coef()`; `logLik()` a `logLik` object; `tidy()` the
#'   `arm` / `term` / `estimate` / `std.error` / `conf.low` / `conf.high` table;
#'   `summary()` the estimate / std.error / interval data frame;
#'   `glance()` the one-row column set of [glance.tobs_fit()]; `plot()` one
#'   Gaussian-density panel per coefficient (up to 4), invisibly.
#' @name tobs_multiarm_methods
#' @export
nobs.tobs_multiarm_fit <- function(object, ...) {
  as.integer(object[["n_total"]] %||% NA_integer_)
}

# --- coefficients (per arm) --------------------------------------------------
# On the posterior (NUTS) path a flat named coefficient vector already exists;
# defer to the draw-based tulpa_fit method so the flat surface is preserved.
#' @rdname tobs_multiarm_methods
#' @export
coef.tobs_multiarm_fit <- function(object, arm = NULL, ...) {
  if (!is.null(object[["draws"]])) return(NextMethod())
  flat <- .tobs_multiarm_flat(object, "estimate")
  if (is.null(arm)) flat else .tobs_arm_coef(flat, object, arm)
}

# One per-arm block field ("estimate" or "se") as a flat `<arm>_<term>` vector.
.tobs_multiarm_flat <- function(object, field) {
  blocks <- .tobs_arm_blocks(object)
  unlist(lapply(names(blocks), function(a)
    stats::setNames(blocks[[a]][[field]],
                    paste0(a, "_", names(blocks[[a]][[field]])))))
}

# --- block-diagonal covariance ----------------------------------------------
#' @rdname tobs_multiarm_methods
#' @export
vcov.tobs_multiarm_fit <- function(object, ...) {
  if (!is.null(object[["draws"]])) return(NextMethod())
  blocks <- .tobs_arm_blocks(object)
  arms   <- names(blocks)
  labels <- names(.tobs_multiarm_flat(object, "estimate"))
  p  <- length(labels)
  V  <- matrix(0, p, p, dimnames = list(labels, labels))
  at <- 0L
  for (a in arms) {
    pa <- length(blocks[[a]][["estimate"]])
    idx <- at + seq_len(pa)
    V[idx, idx] <- blocks[[a]][["vcov"]]
    at <- at + pa
  }
  V
}

# --- Wald confidence intervals ----------------------------------------------
#' @rdname tobs_multiarm_methods
#' @export
confint.tobs_multiarm_fit <- function(object, parm, level = 0.95, ...) {
  if (!is.null(object[["draws"]])) return(NextMethod())
  est <- .tobs_multiarm_flat(object, "estimate")
  se  <- .tobs_multiarm_flat(object, "se")
  z   <- stats::qnorm(1 - (1 - level) / 2)
  ci  <- cbind(est - z * se, est + z * se)
  colnames(ci) <- paste0(format(100 * c((1 - level) / 2, 1 - (1 - level) / 2),
                                 trim = TRUE), " %")
  if (!missing(parm) && !is.null(parm)) ci <- ci[parm, , drop = FALSE]
  ci
}

# --- log-likelihood ----------------------------------------------------------
# What .tobs_multiarm_loglik_val() actually computed, so AIC()/BIC()'s
# `quantity == "log_likelihood"` gate reads a stated quantity instead of an
# unrecorded one. occu_categorical's two Laplace arms are an unpenalized MLE
# only when neither carries a caller-supplied prior (`fit_occu_categorical()`
# stamps `penalized`); a prior-carrying fit reports a log posterior mode
# instead, which AIC/BIC do not accept. Other multiarm families (cover) are
# not established either way, so they keep declining as before.
.tobs_multiarm_loglik_quantity <- function(object) {
  if (inherits(object, "occu_categorical_fit")) {
    if (isTRUE(object[["penalized"]])) {
      return(list(value = NA_character_, declined = "prior_attached"))
    }
    return(list(value = "log_likelihood", declined = NA_character_))
  }
  list(value = NA_character_, declined = "unrecorded")
}

#' @rdname tobs_multiarm_methods
#' @export
logLik.tobs_multiarm_fit <- function(object, ...) {
  if (!is.null(object[["draws"]])) return(NextMethod())
  val <- .tobs_multiarm_loglik_val(object)
  q   <- .tobs_multiarm_loglik_quantity(object)
  ll  <- structure(val, df = .tobs_multiarm_npar(object),
                   nobs = as.integer(object[["n_total"]] %||% NA_integer_),
                   class = "logLik")
  attr(ll, "quantity") <- q$value
  if (!is.na(q$declined)) attr(ll, "declined") <- q$declined
  ll
}

# --- one-row model summary ---------------------------------------------------
#' @rdname tobs_multiarm_methods
#' @export
glance.tobs_multiarm_fit <- function(x, ...) {
  if (!is.null(x[["draws"]])) return(NextMethod())
  .tobs_glance_layout(x, list(n_fixed = length(.tobs_multiarm_flat(x, "estimate"))))
}

# --- tidy coefficient table --------------------------------------------------
#' @rdname tobs_multiarm_methods
#' @param conf.level Interval level for `tidy()` (default 0.95).
#' @export
tidy.tobs_multiarm_fit <- function(x, conf.level = 0.95, ...) {
  if (!is.null(x[["draws"]])) return(NextMethod())
  est <- .tobs_multiarm_flat(x, "estimate")
  se  <- .tobs_multiarm_flat(x, "se")
  ci  <- stats::confint(x, level = conf.level)
  sp  <- .tobs_split_terms(names(est), .tobs_fit_arms(x))
  data.frame(sp, estimate = unname(est), std.error = unname(se),
             conf.low = unname(ci[, 1L]), conf.high = unname(ci[, 2L]),
             row.names = NULL, stringsAsFactors = FALSE)
}

# --- summary -------------------------------------------------------------------
#' @rdname tobs_multiarm_methods
#' @export
summary.tobs_multiarm_fit <- function(object, level = 0.95, ...) {
  if (!is.null(object[["draws"]])) return(NextMethod())
  est <- .tobs_multiarm_flat(object, "estimate")
  ci  <- stats::confint(object, level = level)
  out <- data.frame(estimate  = unname(est),
                    std.error = unname(.tobs_multiarm_flat(object, "se")),
                    lower     = unname(ci[, 1L]), upper = unname(ci[, 2L]),
                    row.names = names(est))
  names(out)[3:4] <- colnames(ci)
  out
}

# --- coefficient plot --------------------------------------------------------
# Same no-draws branch as plot.tulpa_fit (tulpa's `.fit_fixed_table()` /
# Gaussian-density panel), read off the flat multiarm estimate/se instead of
# a single fixed-effect table, so the generic doesn't fall through to
# plot.tulpa_fit and hit its "carries no $draws, $mode/$H_beta, $cov, or
# $means" refusal.
#' @rdname tobs_multiarm_methods
#' @export
plot.tobs_multiarm_fit <- function(x, ...) {
  if (!is.null(x[["draws"]])) return(NextMethod())
  est <- .tobs_multiarm_flat(x, "estimate")
  se  <- .tobs_multiarm_flat(x, "se")
  np  <- length(est)
  n_panel <- min(np, 4L)
  old_par <- graphics::par(mfrow = c(n_panel, 1L), mar = c(4, 4, 1, 1))
  on.exit(graphics::par(old_par))
  for (j in seq_len(n_panel)) {
    m <- est[j]; s <- se[j]
    if (!is.finite(s) || s <= 0) {
      plot(m, 0, type = "p", xlab = names(est)[j], ylab = "density")
      next
    }
    xs <- seq(m - 4 * s, m + 4 * s, length.out = 200)
    plot(xs, stats::dnorm(xs, m, s), type = "l", xlab = names(est)[j], ylab = "density")
    graphics::abline(v = m, col = "red", lty = 2)
  }
  invisible(x)
}

# --- response-scale methods that need a design / RNG decision ----------------
#' @rdname tobs_multiarm_methods
#' @export
fitted.tobs_multiarm_fit <- function(object, ...)
  .tobs_multiarm_unsupported("fitted", object)

#' @rdname tobs_multiarm_methods
#' @export
residuals.tobs_multiarm_fit <- function(object, ...)
  .tobs_multiarm_unsupported("residuals", object)
