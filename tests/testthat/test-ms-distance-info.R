# The assembled per-species marginal observed information for ms_distance()
# (gcol33/tulpaObs#161). The community EM's fallback central-differences this
# matrix at 2U sweeps per species per Newton step, and every distance sweep sums
# over the latent N, which made this family the most expensive file in the
# recovery tier. The kernel already returns every piece of the Louis (1982) block
#
#   B_i = diag(info_lam_i, info_sig_i) - Var(N_i|y) v_i v_i',
#         v_i = (-swl_i, -vN_sig_i)
#
# so .tobs_ms_distance_info_block() assembles it instead. What is asserted here
# is that the assembled matrix IS the finite difference it replaces, and that the
# fit is unchanged by the substitution.

.msdi_cut <- c(0, 25, 50, 75, 100)

# Model + closures for a small community distance fixture: the engine, the two
# designs, and the gradient the finite difference is taken of. Single source for
# every block below.
.msdi_setup <- function(n_species = 3L, N = 25L, seed = 11L) {
  d <- simulate_ms_distance(n_species = n_species, N = N,
                            cutpoints = .msdi_cut, seed = seed)
  fit <- tobs(~ abund_cov1, detection = ~ 1, data = d$data,
              family = ms_distance(cutpoints = .msdi_cut), y = d$y,
              species = d$species, method = "laplace",
              control = list(max.iter = 3L, progress = FALSE))
  model <- fit$model
  eng   <- tulpaObs:::.tobs_ms_distance_engine(model, K_max = NULL)
  X_lam <- model$X_processes[[1L]]; X_sig <- model$X_processes[[2L]]
  p_lam <- ncol(X_lam); p_sig <- ncol(X_sig)
  P <- p_lam + p_sig
  lam_idx <- seq_len(p_lam); sig_idx <- p_lam + seq_len(p_sig)
  eta_of <- function(theta) list(
    lam = as.numeric(X_lam %*% theta[lam_idx]),
    sig = as.numeric(X_sig %*% theta[sig_idx]))
  list(
    d = d, fit = fit, model = model, eng = eng,
    X_lam = X_lam, X_sig = X_sig, P = P,
    lam_idx = lam_idx, sig_idx = sig_idx, S = model$n_species,
    theta0 = as.numeric(fit$means[seq_len(P)]),
    sweep = function(s, theta) {
      e <- eta_of(theta); eng$sweep(s, e$lam, e$sig, 0)
    },
    grad = function(s, theta) {
      e  <- eta_of(theta); sw <- eng$sweep(s, e$lam, e$sig, 0)
      c(as.numeric(crossprod(X_lam, sw$grad_lam)),
        as.numeric(crossprod(X_sig, sw$grad_sig)))
    })
}

# The community EM's own fallback, verbatim: central difference of the gradient,
# symmetrized (R/community_em.R).
.msdi_fd <- function(st, s, theta, h = 1e-4) {
  J <- matrix(0, st$P, st$P)
  for (k in seq_len(st$P)) {
    up <- theta; up[k] <- up[k] + h
    dn <- theta; dn[k] <- dn[k] - h
    J[, k] <- (st$grad(s, up) - st$grad(s, dn)) / (2 * h)
  }
  -0.5 * (J + t(J))
}

.msdi_block <- function(st, s, theta) {
  tulpaObs:::.tobs_ms_distance_info_block(
    st$sweep(s, theta), st$X_lam, st$X_sig, st$lam_idx, st$sig_idx, st$P)
}

test_that("the assembled Louis block reproduces the finite-difference Hessian", {
  st <- .msdi_setup()
  # Evaluated at the mode and away from it: at the mode the score is ~0, which is
  # where a wrong information matrix is easiest to hide.
  offsets <- list(numeric(st$P),
                  c(0.3, -0.2, 0.25)[seq_len(st$P)],
                  c(-0.4, 0.3, -0.35)[seq_len(st$P)])
  for (off in offsets) {
    theta <- st$theta0 + off
    for (s in seq_len(st$S)) {
      A  <- .msdi_block(st, s, theta)
      FD <- .msdi_fd(st, s, theta)
      # Relative to the matrix scale: the entries run to ~1e3 here, so an
      # absolute tolerance would be reporting the scale rather than the fit.
      expect_lt(max(abs(A - FD)) / max(abs(FD)), 1e-5)
    }
  }
})

test_that("the cross block is large enough for that agreement to mean something", {
  # The lambda-sigma block is the only part the sign convention inside v moves,
  # and the whole of it comes from the rank-1 score-variance subtraction. If it
  # were ~0 on this fixture the test above would pass on a wrong assembly.
  st <- .msdi_setup()
  theta <- st$theta0 + c(0.3, -0.2, 0.25)[seq_len(st$P)]
  A <- .msdi_block(st, 1L, theta)
  cross <- max(abs(A[st$lam_idx, st$sig_idx]))
  diagb <- max(abs(c(A[st$lam_idx, st$lam_idx], A[st$sig_idx, st$sig_idx])))
  expect_gt(cross / diagb, 0.1)
})

test_that("the sign inside v is pinned, not inherited from the N-mixture", {
  # nmix_site_marginal()'s block uses v = (-swl, +p); the distance kernel stores
  # vN_d already negated, so the second component takes a further flip. Getting
  # that wrong leaves the diagonal untouched and the cross block wrong by a
  # factor of ~2, so assert the wrong convention is REJECTED by the same
  # comparison the block above passes.
  st <- .msdi_setup()
  theta <- st$theta0 + c(0.3, -0.2, 0.25)[seq_len(st$P)]
  sw <- st$sweep(1L, theta)
  wrong <- {
    b11 <- sw$info_lam     - sw$var_N * sw$swl^2
    b22 <- sw$info_sig_obs - sw$var_N * sw$vN_sig^2
    b12 <-                   sw$var_N * sw$swl * sw$vN_sig   # sign flipped
    I <- matrix(0, st$P, st$P)
    I[st$lam_idx, st$lam_idx] <- crossprod(st$X_lam, st$X_lam * b11)
    I[st$sig_idx, st$sig_idx] <- crossprod(st$X_sig, st$X_sig * b22)
    cr <- crossprod(st$X_lam, st$X_sig * b12)
    I[st$lam_idx, st$sig_idx] <- cr
    I[st$sig_idx, st$lam_idx] <- t(cr)
    (I + t(I)) / 2
  }
  FD <- .msdi_fd(st, 1L, theta)
  expect_gt(max(abs(wrong - FD)) / max(abs(FD)), 1e-2)
})

test_that("the hazard key keeps the finite-difference fallback and still fits", {
  # grad_b / info_b come back already summed over sites and the per-site
  # detection cross terms are not exported, so the shared log-shape global cannot
  # be sandwiched. The fitter passes sp_info = NULL there.
  d <- simulate_ms_distance(n_species = 3, N = 25, cutpoints = .msdi_cut,
                            seed = 4)
  fit <- tobs(~ abund_cov1, detection = ~ 1, data = d$data,
              family = ms_distance(cutpoints = .msdi_cut, key = "hazard"),
              y = d$y, species = d$species, method = "laplace",
              control = list(max.iter = 3L, progress = FALSE))
  expect_s3_class(fit, "tobs_fit")
  expect_true(all(is.finite(coef(fit)$lambda)))
})

test_that("supplying the block leaves the community fit unchanged", {
  # Same closures, same start, same iteration budget: the only difference is
  # which Hessian the Newton step uses. Both must land on the same mode, the same
  # community covariances and the same standard errors -- the FD Hessian was
  # accurate, just expensive.
  st <- .msdi_setup(n_species = 5L, N = 40L, seed = 7L)
  sp_ll <- function(s, theta, global) sum(st$sweep(s, theta)$log_lik)
  sp_grad <- function(s, theta, global) st$grad(s, theta)
  sp_info <- function(s, theta, global) .msdi_block(st, s, theta)
  run <- function(info) tulpaObs:::.tobs_community_em(
    S = st$S, P = st$P,
    arm_idx = list(lambda = st$lam_idx, sigma = st$sig_idx),
    sp_ll = sp_ll, sp_grad = sp_grad, sp_info = info,
    init_mu = st$theta0, init_global = numeric(0), penalize_global = FALSE,
    sigma_beta = 5, priors = NULL, sigma_init = 0.3,
    max_iter = 30L, tol = 1e-4, newton_max = 20L, verbose = FALSE)
  a <- run(sp_info)
  b <- run(NULL)
  expect_lt(max(abs(a$mu - b$mu)), 1e-8)
  for (arm in names(a$Sigma))
    expect_lt(max(abs(a$Sigma[[arm]] - b$Sigma[[arm]])), 1e-7)
  expect_lt(max(abs(unlist(a$b_list) - unlist(b$b_list))), 1e-8)
  se_a <- sqrt(diag(a$Vf)); se_b <- sqrt(diag(b$Vf))
  expect_lt(max(abs(se_a - se_b) / se_b), 1e-6)
})
