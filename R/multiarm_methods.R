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
#' @param level Confidence level for `confint()`.
#' @param nsim,seed Present for the `simulate()` generic; `simulate()` is not
#'   supported (use the family's `simulate_*()` function).
#' @param ... Ignored, or forwarded to `NextMethod()` on the posterior path.
#' @return `nobs()` an integer; `coef()` a named list of per-arm estimates;
#'   `vcov()` the block-diagonal covariance; `confint()` a two-column matrix;
#'   `logLik()` a `logLik` object; `glance()`/`tidy()` a data frame. On a joint
#'   nested-Laplace fit `glance()` also carries `outer_grid_placement` and
#'   `outer_grid_recenter_declined`.
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
coef.tobs_multiarm_fit <- function(object, ...) {
  if (!is.null(object[["draws"]])) return(NextMethod())
  blocks <- .tobs_arm_blocks(object)
  lapply(blocks, `[[`, "estimate")
}

# --- block-diagonal covariance ----------------------------------------------
#' @rdname tobs_multiarm_methods
#' @export
vcov.tobs_multiarm_fit <- function(object, ...) {
  if (!is.null(object[["draws"]])) return(NextMethod())
  blocks <- .tobs_arm_blocks(object)
  arms   <- names(blocks)
  labels <- unlist(lapply(arms, function(a)
    paste0(a, ":", names(blocks[[a]][["estimate"]]))), use.names = FALSE)
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
  blocks <- .tobs_arm_blocks(object)
  est <- unlist(lapply(names(blocks), function(a)
    stats::setNames(blocks[[a]][["estimate"]],
                    paste0(a, ":", names(blocks[[a]][["estimate"]])))))
  se  <- unlist(lapply(names(blocks), function(a)
    stats::setNames(blocks[[a]][["se"]],
                    paste0(a, ":", names(blocks[[a]][["se"]])))))
  z   <- stats::qnorm(1 - (1 - level) / 2)
  ci  <- cbind(est - z * se, est + z * se)
  colnames(ci) <- paste0(format(100 * c((1 - level) / 2, 1 - (1 - level) / 2),
                                 trim = TRUE), " %")
  if (!missing(parm) && !is.null(parm)) ci <- ci[parm, , drop = FALSE]
  ci
}

# --- log-likelihood ----------------------------------------------------------
#' @rdname tobs_multiarm_methods
#' @export
logLik.tobs_multiarm_fit <- function(object, ...) {
  if (!is.null(object[["draws"]])) return(NextMethod())
  val <- .tobs_multiarm_loglik_val(object)
  structure(val, df = .tobs_multiarm_npar(object),
            nobs = as.integer(object[["n_total"]] %||% NA_integer_),
            class = "logLik")
}

# --- one-row model summary ---------------------------------------------------
#' @rdname tobs_multiarm_methods
#' @export
glance.tobs_multiarm_fit <- function(x, ...) {
  ll <- logLik(x)
  g <- data.frame(
    n         = as.integer(x[["n_total"]] %||% NA_integer_),
    logLik    = as.numeric(ll),
    df        = attr(ll, "df"),
    converged = isTRUE(x[["convergence"]][["converged"]] %||% x[["converged"]]),
    stringsAsFactors = FALSE
  )
  # This method is terminal for cover_fit, so the joint outer-grid placement
  # has to be added here as well as in glance.tobs_fit().
  .tobs_glance_outer_grid(g, x)
}

# --- tidy coefficient table --------------------------------------------------
#' @rdname tobs_multiarm_methods
#' @export
tidy.tobs_multiarm_fit <- function(x, ...) {
  if (!is.null(x[["draws"]])) return(NextMethod())
  blocks <- .tobs_arm_blocks(x)
  do.call(rbind, lapply(names(blocks), function(a) {
    b <- blocks[[a]]
    z <- stats::qnorm(0.975)
    data.frame(
      arm       = a,
      term      = names(b[["estimate"]]),
      estimate  = as.numeric(b[["estimate"]]),
      std.error = as.numeric(b[["se"]]),
      conf.low  = as.numeric(b[["estimate"]] - z * b[["se"]]),
      conf.high = as.numeric(b[["estimate"]] + z * b[["se"]]),
      row.names = NULL, stringsAsFactors = FALSE
    )
  }))
}

# --- summary (serves occu_categorical_fit; cover_fit has its own) ------------
#' @rdname tobs_multiarm_methods
#' @export
summary.tobs_multiarm_fit <- function(object, ...) {
  if (!is.null(object[["draws"]])) return(NextMethod())
  td  <- tidy(object)
  fam <- object[["family"]][["name"]] %||% "multi-arm"
  cat(sprintf("<%s: %d observations>\n", fam,
              as.integer(object[["n_total"]] %||% NA_integer_)))
  for (a in unique(td$arm)) {
    sub <- td[td$arm == a, , drop = FALSE]
    tab <- data.frame(
      Estimate   = sub$estimate,
      `Std.Error` = sub$std.error,
      z          = sub$estimate / sub$std.error,
      `Pr(>|z|)` = 2 * stats::pnorm(-abs(sub$estimate / sub$std.error)),
      row.names  = sub$term, check.names = FALSE)
    cat(sprintf("\n  %s arm:\n", a))
    print(round(tab, 4))
  }
  invisible(td)
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

#' @rdname tobs_multiarm_methods
#' @export
simulate.tobs_multiarm_fit <- function(object, nsim = 1, seed = NULL, ...)
  .tobs_multiarm_unsupported("simulate", object)
