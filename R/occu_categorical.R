# =============================================================================
# occu_categorical.R -- presence + nominal (unordered) K-class hurdle
# (gcol33/tulpaObs#106).
#
# Observed-state sub-model: each unit is either absent (y = 0) or present in one
# of K nominal classes (y in 1..K). Presence is a Bernoulli arm; the class given
# present is a baseline-category multinomial logit (class K the baseline). The two
# arms factorise the likelihood exactly:
#   P(y = 0) = 1 - psi,   P(y = k) = psi * p_k,   p = softmax(X beta_class).
# This is the categorical analogue of the cover() hurdle (presence + continuous
# magnitude): here the positive part is an unordered class, not a magnitude, so it
# uses a multinomial logit rather than beta / lognormal.
#
# The multinomial math (softmax, score, coupled Hessian) is the FD-validated
# tulpa kernel `cpp_multinomial_logit_terms` (src/multinomial_logit.h); the
# non-spatial Laplace fit below is the vectorised R Newton over the same closed
# forms. The native multi-process LikelihoodSpec path (spatial fields / NUTS) is
# the documented follow-up.
# =============================================================================


# Vectorised multinomial-logit Laplace (Newton) over a coefficient matrix
# Beta [p x (K-1)], baseline class K. `cls` is the observed class (1..K) per row
# of `X`. A weak ridge `prior_prec` (precision of an N(0, 1/prior_prec) prior on
# every coefficient) keeps the Hessian non-singular under a rare / empty class.
# Returns the mode `Beta`, the joint vcov of vec(Beta) (column-major: all p
# coefficients of class 1, then class 2, ...), the log-likelihood and convergence.
.tobs_mlogit_fit <- function(X, cls, K, max_iter = 100L, tol = 1e-8,
                             prior_prec = 1e-4) {
  n <- nrow(X); p <- ncol(X); Km1 <- K - 1L
  if (length(cls) != n) stop("length(cls) must equal nrow(X).", call. = FALSE)
  if (any(cls < 1L | cls > K)) stop("`cls` must be in 1..K.", call. = FALSE)

  # Class indicator over the non-baseline classes (baseline K -> all-zero row).
  Y <- matrix(0, n, Km1)
  for (j in seq_len(Km1)) Y[cls == j, j] <- 1

  Beta <- matrix(0, p, Km1)
  Ip   <- diag(prior_prec, p * Km1)
  converged <- FALSE; ll <- NA_real_
  for (it in seq_len(max_iter)) {
    eta <- X %*% Beta                                   # n x (K-1)
    E   <- exp(pmin(eta, 700))                          # overflow guard
    denom <- 1 + rowSums(E)
    P   <- E / denom                                    # n x (K-1) class probs
    grad <- crossprod(X, Y - P) - prior_prec * Beta     # p x (K-1)

    # Block Hessian: -d2 ll / dBeta_j dBeta_l = X' diag(P_j (1{j=l} - P_l)) X.
    H <- matrix(0, p * Km1, p * Km1)
    for (j in seq_len(Km1)) for (l in seq_len(Km1)) {
      w <- P[, j] * ((j == l) - P[, l])
      H[((j - 1L) * p + 1L):(j * p), ((l - 1L) * p + 1L):(l * p)] <-
        crossprod(X, w * X)
    }
    H <- H + Ip
    step <- solve(H, as.numeric(grad))                  # column-major vec(grad)
    Beta <- Beta + matrix(step, p, Km1)
    if (max(abs(step)) < tol) { converged <- TRUE; break }
  }
  eta <- X %*% Beta; E <- exp(pmin(eta, 700)); denom <- 1 + rowSums(E)
  pmat <- cbind(E / denom, 1 / denom)                   # n x K (baseline last)
  ll <- sum(log(pmat[cbind(seq_len(n), cls)]))
  V  <- tryCatch(solve(H), error = function(e) NULL)
  list(Beta = Beta, vcov = V, loglik = ll, converged = converged, K = K,
       n_iter = it)
}


# Encode an occu_categorical response: split y in {0, 1..K} into a presence
# indicator and the class on the present subset, with the shared fixed-effects
# design for each arm.
.encode_occu_categorical <- function(formula, data, y, K) {
  if (!is.numeric(y)) y <- as.integer(y)
  if (any(is.na(y))) stop("`y` must not contain NA.", call. = FALSE)
  if (any(y != round(y)) || any(y < 0) || any(y > K)) {
    stop("`y` must be integers in 0..K (0 = absent, k = class k).", call. = FALSE)
  }
  fe <- .tobs_bind_formulas(list(state = formula), data)$fe$state
  X  <- stats::model.matrix(fe, data)
  present <- as.integer(y > 0)
  is_pos  <- present == 1L
  list(present = present, X_occ = X, X_class = X[is_pos, , drop = FALSE],
       cls = as.integer(y[is_pos]), idx_pos = which(is_pos), K = K,
       formula = fe, N = length(y))
}


# ---------------------------------------------------------------------------
# Dispatcher (called from tobs())
# ---------------------------------------------------------------------------

.dispatch_occu_categorical <- function(formula, data, family, detection, y,
                                       visits, engine, priors, control,
                                       approx = "gaussian_laplace",
                                       correction = "none", ...) {
  if (!identical(engine, "laplace")) {
    stop("occu_categorical() currently supports method = 'laplace' only (got '",
         engine, "'). Spatial / NUTS paths are not yet wired.", call. = FALSE)
  }
  if (!is.null(detection)) {
    stop("occu_categorical() does not use a `detection` formula ",
         "(one observation per unit). Drop it.", call. = FALSE)
  }
  if (is.null(y)) {
    stop("occu_categorical() requires `y` (a length-N integer vector in 0..K: ",
         "0 = absent, k = class k).", call. = FALSE)
  }
  classes <- family$params$classes
  K <- if (!is.null(classes)) length(classes) else as.integer(max(y))
  if (K < 2L) {
    stop("occu_categorical() needs at least two classes (max(y) = ", K, ").",
         call. = FALSE)
  }
  class_labels <- classes %||% as.character(seq_len(K))

  enc <- .encode_occu_categorical(formula, data, y, K)
  .tobs_check_site_count(length(y), nrow(data), "values")
  if (length(enc$cls) < K) {
    stop("occu_categorical(): only ", length(enc$cls), " present unit(s) for ",
         K, " classes; need at least one present unit per class to identify the ",
         "multinomial.", call. = FALSE)
  }
  fit_occu_categorical(enc, class_labels, priors, control, family)
}


# ---------------------------------------------------------------------------
# Fit: Bernoulli presence + multinomial class, assembled into a tobs_fit
# ---------------------------------------------------------------------------

fit_occu_categorical <- function(enc, class_labels, priors, control, family) {
  max_iter   <- control$max.iter  %||% 100L
  tol        <- control$tol       %||% 1e-8
  prior_prec <- control$prior.prec %||% 1e-4

  # Presence arm: a plain binomial Laplace on (y > 0) ~ X.
  m_occ <- tulpa::tulpa_laplace(
    y = enc$present, n_trials = rep(1L, enc$N), X = enc$X_occ,
    family = "binomial", max_iter = max_iter, tol = tol)
  p_occ    <- ncol(enc$X_occ)
  beta_occ <- m_occ$mode[seq_len(p_occ)]
  V_occ    <- tryCatch(solve(m_occ$H_beta), error = function(e) NULL)
  se_occ   <- if (is.null(V_occ)) rep(NA_real_, p_occ) else
    sqrt(pmax(diag(as.matrix(V_occ)), 0))
  names(beta_occ) <- names(se_occ) <- colnames(enc$X_occ)

  # Class arm: baseline-category multinomial logit on the present subset.
  mc <- .tobs_mlogit_fit(enc$X_class, enc$cls, enc$K,
                         max_iter = max_iter, tol = tol, prior_prec = prior_prec)
  Beta <- mc$Beta
  rownames(Beta) <- colnames(enc$X_class)
  colnames(Beta) <- class_labels[seq_len(enc$K - 1L)]   # non-baseline classes
  se_class <- if (is.null(mc$vcov)) matrix(NA_real_, nrow(Beta), ncol(Beta)) else
    matrix(sqrt(pmax(diag(mc$vcov), 0)), nrow(Beta), ncol(Beta),
           dimnames = dimnames(Beta))

  structure(
    list(
      beta_occ    = beta_occ,
      se_occ      = se_occ,
      vcov_occ    = V_occ,
      beta_class  = Beta,
      se_class    = se_class,
      vcov_class  = mc$vcov,
      K           = enc$K,
      class_labels = class_labels,
      baseline    = class_labels[enc$K],
      n_total     = enc$N,
      n_present   = length(enc$cls),
      loglik      = c(occ = as.numeric(m_occ$log_marginal %||% NA_real_),
                      class = mc$loglik),
      converged   = isTRUE(m_occ$converged) && isTRUE(mc$converged),
      convergence = list(converged = isTRUE(m_occ$converged) && isTRUE(mc$converged),
                         n_iter = mc$n_iter),
      encoding    = enc,
      family      = family
    ),
    class = c("occu_categorical_fit", "tobs_multiarm_fit", "tobs_fit", "tulpa_fit"))
}


# ---------------------------------------------------------------------------
# S3 + predict
# ---------------------------------------------------------------------------

#' @export
print.occu_categorical_fit <- function(x, ...) {
  cat("<occu_categorical: presence + nominal class hurdle>\n")
  cat(sprintf("  classes      : %d (%s; baseline = %s)\n", x$K,
              paste(x$class_labels, collapse = ", "), x$baseline))
  cat(sprintf("  N total      : %d\n", x$n_total))
  cat(sprintf("  N present    : %d (%.1f%%)\n", x$n_present,
              100 * x$n_present / x$n_total))
  cat(sprintf("  converged    : %s\n", x$convergence$converged))
  cat("\n  presence (logit psi):\n")
  print(round(x$beta_occ, 4))
  cat("\n  class (baseline-category multinomial logit):\n")
  print(round(x$beta_class, 4))
  invisible(x)
}

#' @export
coef.occu_categorical_fit <- function(object, ...) {
  list(presence = object$beta_occ, class = object$beta_class)
}

#' Predict from an occu_categorical fit
#'
#' Returns, per row of `newdata`, the presence probability `psi`, the
#' conditional class probabilities `p_1..p_K` (given present), and the
#' unconditional class probabilities `psi * p_k` plus the absent probability
#' `1 - psi`.
#'
#' @param object an `occu_categorical_fit`.
#' @param newdata a data frame of predictors for the shared formula.
#' @param ... ignored.
#' @return A list with `psi` (length n), `cond` (n x K conditional class
#'   probabilities), and `joint` (n x (K+1): absent, then class 1..K).
#' @export
predict.occu_categorical_fit <- function(object, newdata, ...) {
  if (missing(newdata) || is.null(newdata)) {
    stop("`newdata` is required.", call. = FALSE)
  }
  X <- stats::model.matrix(object$encoding$formula, newdata)
  psi  <- as.numeric(stats::plogis(X %*% object$beta_occ))
  eta  <- X %*% object$beta_class                   # n x (K-1)
  E    <- exp(pmin(eta, 700)); denom <- 1 + rowSums(E)
  cond <- cbind(E / denom, 1 / denom)               # n x K (baseline last)
  colnames(cond) <- object$class_labels
  joint <- cbind(absent = 1 - psi, psi * cond)
  list(psi = psi, cond = cond, joint = joint)
}


# ---------------------------------------------------------------------------
# Simulate
# ---------------------------------------------------------------------------

#' Simulate from an occu_categorical model
#'
#' @param N number of units.
#' @param beta_occ length-2 presence coefficients (intercept, slope on `x`).
#' @param beta_class a `2 x (K-1)` matrix of class coefficients (rows
#'   intercept / slope on `x`; columns the non-baseline classes), or `NULL`
#'   for a built-in 3-class default.
#' @param seed optional RNG seed.
#' @return A list with `data` (a data frame with `x`), `y` (length-N integer in
#'   `0..K`), and `truth`.
#' @export
simulate_occu_categorical <- function(N = 500L, beta_occ = c(0.2, 0.8),
                                      beta_class = NULL, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  if (is.null(beta_class)) {
    beta_class <- matrix(c(0.4, 1.0, -0.5, -0.8, 0.2, 0.6), nrow = 2)  # K = 4
  }
  K <- ncol(beta_class) + 1L
  x <- stats::rnorm(N)
  X <- cbind(1, x)
  psi   <- stats::plogis(as.numeric(X %*% beta_occ))
  present <- stats::rbinom(N, 1L, psi)
  eta <- X %*% beta_class
  E   <- exp(eta); denom <- 1 + rowSums(E)
  P   <- cbind(E / denom, 1 / denom)                       # N x K
  cls <- apply(P, 1L, function(pr) sample.int(K, 1L, prob = pr))
  y   <- ifelse(present == 1L, cls, 0L)
  list(data = data.frame(x = x), y = as.integer(y),
       truth = list(beta_occ = beta_occ, beta_class = beta_class, K = K))
}

