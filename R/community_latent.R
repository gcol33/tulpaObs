# =============================================================================
# community_latent.R - shared latent-structure engine for the community families
#
# A community model carrying a shared areal field and/or latent factors on one
# arm has the linear predictor
#
#   eta_{s,i} = X_i . (mu + b_s) + sum_k W[i,k] F[u(i),k]
#                                + sum_q lambda_{s,q} zeta_{q,i}
#
# with b_s ~ N(0, Sigma) the per-species coefficient deviations, F the areal
# field(s) (an all-ones weight column for the intercept field, covariate values
# for a varying-coefficient field), and (zeta, lambda) the per-site latent
# factors and their per-species loadings.
#
# The fit is a block coordinate ascent. The latent structure enters each
# species' log-likelihood as a fixed per-(site, species) offset, so a
# coefficient update is an ordinary community Laplace-EM (.tobs_community_em,
# used unchanged); the field and factor updates given the coefficients are
# self-contained Newton solves. Alternating the blocks converges to the joint
# mode. This sidesteps the community EM's finite-difference Hessian, which does
# not scale to an O(n_sites) field, and needs no C++.
#
# Every family enters through ONE interface. A family supplies an `oracle`, and
# the only thing the field solve and the factor update read from it is
#
#   working(eta_mat) -> list(score = [n_sites x n_species],
#                            curv  = [n_sites x n_species])
#
# the score and curvature of that family's marginal log-likelihood with respect
# to an additive offset on the structured arm's linear predictor. Poisson gives
# (y - mu, mu); the occupancy two-state marginal gives its own pair; a Bernoulli
# JSDM gives (y - psi, psi(1-psi)). Adding a family to every latent route is
# therefore one callback, not a new fitter.
#
# Oracle fields:
#   working(eta_mat)   (required) score / curvature, as above.
#   data_ll(eta_mat)   (required) the summed data log-likelihood, used to compare
#                      field hyperparameters across a grid.
#   cov_curv(eta_mat)  (optional) the curvature used for the Laplace covariance
#                      in the tau M-step, when a family stabilises it away from
#                      the Newton curvature; defaults to working()$curv.
#   n_sites, n_species (required) dimensions.
# =============================================================================


# ---------------------------------------------------------------------------
# Loadings
# ---------------------------------------------------------------------------

# The factors are a LATENT, not a parameter: zeta_i ~ N(0, I_Q). Everything in
# this section follows from taking that seriously, and the estimator arrived at
# it in two steps.
#
# The factor update holds the factors at their joint mode and treats
# zeta t(lambda) as a known offset. That penalised objective is UNBOUNDED in the
# magnitude direction -- lambda is free to rescale, so sending zeta -> 0 and
# lambda -> Inf at fixed product drives the N(0, I) prior cost to zero (measured
# on lfJSDM with the unit-variance anchor removed: loadings to 48x truth). The
# anchor masks that rather than fixing it, leaving the magnitude at pure maximum
# likelihood, which over-fits and inflates the community coefficients to
# compensate (gcol33/tulpaObs#153).
#
# Setting the magnitude from the joint site marginal fixed the worst of it, but
# not the rest: the DIRECTION handed to that search is itself a joint-mode
# estimate over Ns * Q incidental parameters, so it carries the factors'
# estimation error and no single scalar can take that back out
# (gcol33/tulpaObs#156). Integrating the factors out of the loadings entirely is
# what actually closes it.
#
# So the pieces below are, in the order the fit uses them: the shared quadrature
# grid, the 1-D magnitude search (now an INITIALISER only), the marginal-
# likelihood loading EM, and the offset the coefficient block conditions on.

# Gauss-Hermite nodes and weights for the standard normal, by Golub-Welsch on the
# probabilists' Hermite Jacobi matrix. Weights sum to one, so
# integral f(u) phi(u) du = sum_k w_k f(x_k).
.tobs_gh_nodes <- function(n = 5L) {
  n <- as.integer(n)
  k <- seq_len(n - 1L)
  J <- matrix(0, n, n)
  J[cbind(k, k + 1L)] <- sqrt(k)
  J[cbind(k + 1L, k)] <- sqrt(k)
  e <- eigen(J, symmetric = TRUE)
  o <- order(e$values)
  list(x = e$values[o], w = (e$vectors[1L, o])^2)
}

# The adaptive Gauss-Hermite grid for the JOINT site marginal, integrating zeta_i
# over all Q dimensions with the species at a site kept together:
#
#   L_i = integral prod_s f(y_is | eta_is + lambda_s' z) N(z; 0, I_Q) dz
#
# on a tensor Gauss-Hermite grid. A node is a fixed Q-vector, so it shifts cell
# (i, s) by lambda_s' z_k -- the same shift at every site -- and one node costs a
# single ll_cell call.
#
# Returns the per-site node positions and the unnormalised per-node log terms.
# Log-sum-exp over the nodes gives the marginal (.tobs_latent_joint_marginal);
# normalising them over the nodes gives the posterior p(z_i | y_i) that the
# loading M-step averages the complete-data score over
# (.tobs_latent_factor_mmle). Both readers want the same grid, so it is built
# once here rather than twice.
#
# The rule is ADAPTIVE: the nodes are placed at the mode and curvature of each
# site's own integrand rather than at the prior's scale. Fixed prior-scale nodes
# are accurate when the per-site likelihood is flat and the N(0, I) prior
# dominates -- which is the Bernoulli case, and is the only case the node-count
# stability in test-community-latent-quad.R ever exercised. A Poisson site with
# counts of a few carries a sharply peaked, shifted integrand that five
# prior-scale nodes do not resolve, and the resulting quadrature error does NOT
# cancel in the argmax: measured on the clean fixture (true eta, known loadings
# handed in inflated), the implied magnitude wanders 0.69-0.93 across node counts
# 5-35 on a Poisson oracle while sitting stable to under 0.4% on a Bernoulli one.
# That is the ~25% low magnitude on the count routes in gcol33/tulpaObs#154.
#
# Non-adaptive is the zhat = 0, A = I special case of what follows, so there is
# one path here rather than two: when the mode-find cannot proceed (a non-finite
# score far from the mode, which the count marginals can return) it falls back to
# exactly that and reproduces the previous estimator.
.tobs_latent_joint_grid <- function(oracle, eta_base, lambda, gh,
                                    adapt = TRUE) {
  Ns <- nrow(eta_base); Qk <- ncol(lambda)
  nd <- as.matrix(expand.grid(rep(list(gh$x), Qk)))
  lw <- rowSums(as.matrix(expand.grid(rep(list(log(gh$w)), Qk))))

  zhat <- matrix(0, Ns, Qk)
  A <- array(0, dim = c(Ns, Qk, Qk))
  for (a in seq_len(Qk)) A[, a, a] <- 1
  logdetA <- numeric(Ns)

  if (isTRUE(adapt)) {
    # Mode of h_i(z) = sum_s ll_cell(eta_is + lambda_s' z) - z'z / 2 by Newton.
    # The gradient and curvature are the same quantities the factor update's
    # zeta half-pass uses, assembled across sites at once: the Hessian entry
    # (a, b) is sum_s curv_is lambda_sa lambda_sb + delta_ab.
    # h_i(z) = sum_s ll_cell(eta_is + lambda_s' z) - z'z / 2, the objective whose
    # mode and curvature the nodes adapt to.
    hval <- function(z) {
      v <- rowSums(oracle$ll_cell(eta_base + tcrossprod(z, lambda))) -
        0.5 * rowSums(z^2)
      ifelse(is.finite(v), v, -Inf)
    }
    # A Newton step is not guaranteed to ascend, and the count marginals are
    # exp-linked, so an undamped overshoot drives exp(eta) past overflow and the
    # mode lands somewhere arbitrary -- which then places every node there and
    # makes the quadrature worse than not adapting at all. Backtrack per site:
    # halve only the sites whose step did not improve, keep the ones that did.
    ok <- TRUE
    h0 <- hval(zhat)
    if (any(!is.finite(h0))) ok <- FALSE
    for (it in seq_len(12L)) {
      if (!ok) break
      wk <- oracle$working(eta_base + tcrossprod(zhat, lambda))
      if (any(!is.finite(wk$score)) || any(!is.finite(wk$curv))) { ok <- FALSE; break }
      G <- wk$score %*% lambda - zhat
      step <- matrix(0, Ns, Qk)
      for (i in seq_len(Ns)) {
        Hi <- matrix(0, Qk, Qk)
        for (a in seq_len(Qk)) for (b in seq_len(Qk)) {
          Hi[a, b] <- sum(wk$curv[i, ] * lambda[, a] * lambda[, b]) +
            (if (a == b) 1 else 0)
        }
        step[i, ] <- tryCatch(solve(Hi, G[i, ]),
                              error = function(e) rep(NA_real_, Qk))
      }
      step[!is.finite(step)] <- 0
      tt <- rep(1, Ns)
      moved <- rep(FALSE, Ns)
      for (bt in seq_len(25L)) {
        pend <- !moved & (rowSums(abs(step)) > 0)
        if (!any(pend)) break
        cand <- zhat + step * tt
        hc <- hval(cand)
        take <- pend & is.finite(hc) & (hc > h0)
        if (any(take)) {
          zhat[take, ] <- cand[take, , drop = FALSE]
          h0[take] <- hc[take]
          moved[take] <- TRUE
        }
        tt[pend & !take] <- tt[pend & !take] / 2
      }
      if (max(abs(step * tt)) < 1e-8) break
    }
    if (ok) {
      wk <- oracle$working(eta_base + tcrossprod(zhat, lambda))
      ok <- all(is.finite(wk$score)) && all(is.finite(wk$curv))
    }
    if (ok) {
      # A_i = H_i^{-1/2} via the Cholesky H = L L', A = L^{-T}: then A A' =
      # H^{-1} and log|det A| = -sum(log(diag(L))).
      for (i in seq_len(Ns)) {
        Hi <- matrix(0, Qk, Qk)
        for (a in seq_len(Qk)) for (b in seq_len(Qk)) {
          Hi[a, b] <- sum(wk$curv[i, ] * lambda[, a] * lambda[, b]) +
            (if (a == b) 1 else 0)
        }
        Li <- tryCatch(chol(Hi), error = function(e) NULL)
        if (is.null(Li)) { ok <- FALSE; break }
        A[i, , ] <- backsolve(Li, diag(Qk))       # L^{-T}, since chol() is upper
        logdetA[i] <- -sum(log(diag(Li)))
      }
    }
    if (!ok) {
      zhat[] <- 0
      A[] <- 0
      for (a in seq_len(Qk)) A[, a, a] <- 1
      logdetA[] <- 0
    }
  }

  cphi <- -0.5 * Qk * log(2 * pi)
  Z <- vector("list", nrow(nd))
  L <- vector("list", nrow(nd))
  for (k in seq_len(nrow(nd))) {
    Zk <- zhat
    for (b in seq_len(Qk)) Zk <- Zk + A[, , b] * nd[k, b]
    Z[[k]] <- Zk
    L[[k]] <- lw[k] + rowSums(oracle$ll_cell(eta_base + tcrossprod(Zk, lambda))) +
      (cphi - 0.5 * rowSums(Zk^2)) - (cphi - 0.5 * sum(nd[k, ]^2))
  }
  list(Z = Z, logterm = L, logdetA = logdetA)
}

# The joint site marginal itself: log-sum-exp of the grid's per-node log terms.
.tobs_latent_joint_marginal <- function(oracle, eta_base, lambda, gh,
                                        adapt = TRUE) {
  g <- .tobs_latent_joint_grid(oracle, eta_base, lambda, gh, adapt = adapt)
  m <- Reduce(pmax, g$logterm)
  sum(g$logdetA + m +
        log(Reduce(`+`, lapply(g$logterm, function(z) exp(z - m)))))
}

# The factor SCALE: the scalar c maximising the joint marginal at loadings
# c * lambda. This runs ONCE, on the first outer pass, to put the joint-mode
# direction at a sane magnitude before the loading EM refines all of it. Its
# value there is that it searches GLOBALLY over an expanding bracket, which a
# local ascent cannot do, and the joint-mode iterate it is handed can be an
# order of magnitude out.
#
# It has to be the JOINT marginal. Species s' own marginal integrates zeta in
# closed form to one dimension and depends on lambda only through ||lambda_s||,
# which looks like it identifies the magnitude -- but it does not. Species s
# contributes a single Bernoulli per site, so its own marginal carries no
# replication with which to see overdispersion: a normal-mixed logit is very
# nearly a rescaled logit, plogis(eta / sqrt(1 + 0.346 sigma^2)), leaving sigma
# and the coefficient scale confounded along a ridge. Held at the true
# coefficients that profile peaks correctly; fitted jointly it slides to whichever
# end of the grid it starts toward (measured: to the floor, 0.22x truth). What
# identifies the magnitude is that zeta_i is SHARED across species at a site, so
# it is the cross-species co-occurrence -- visible only in the joint integral --
# that pins it.
#
# One shared scalar rather than a magnitude per species, because a per-species
# search at this stage is weakly identified (estimates scattering 0.00-3.20
# against known 0.16-1.91, and a species pushed to zero loses its row of the
# residual correlation). Per-species magnitudes ARE recoverable, but from the
# marginal likelihood with the factors integrated out, which is the EM below --
# not from a profile over a joint-mode direction.
.tobs_latent_factor_scale <- function(oracle, eta_base, lambda, gh) {
  # The magnitude step reads per-cell densities, which `working()` does not
  # carry. Name the missing callback here: without this the first family wired
  # in without one dies inside rowSums() as "attempt to apply non-function".
  if (!is.function(oracle$ll_cell)) {
    stop("the latent factor magnitude needs the oracle's `ll_cell(eta)` ",
         "returning an [n_sites x n_species] matrix of per-cell marginal ",
         "log-likelihoods; this oracle supplies none.", call. = FALSE)
  }
  if (all(sqrt(rowSums(lambda^2)) < 1e-8)) return(1)

  prof_at <- function(g)
    .tobs_latent_joint_marginal(oracle, eta_base, g * lambda, gh)

  # The bracket EXPANDS rather than clamping. A fixed [0.2, 1.5] window returned
  # its own boundary as though it were an optimum, so a run in which the update's
  # lambda had drifted far from the marginal's answer reported a plausible number
  # instead of a saturated one -- and because the reported loadings feed the next
  # EM offset, a saturated scale and a regrowing update ratchet against each
  # other (gcol33/tulpaObs#154, seed 305: c pinned at the 0.2 floor while the
  # loadings ran to 10.7x truth). Expanding makes the scale an actual argmax; the
  # attribute below reports the cases where even the widened bracket did not
  # close, so saturation is visible rather than silent.
  lo <- 0.2; hi <- 1.5
  grid <- exp(seq(log(lo), log(hi), length.out = 14L))
  prof <- vapply(grid, prof_at, numeric(1))
  saturated <- FALSE
  for (expand in seq_len(6L)) {
    j <- which.max(prof)
    if (j != 1L && j != length(grid)) break
    if (j == 1L) {
      lo <- lo / 4
      add <- exp(seq(log(lo), log(grid[1L]), length.out = 5L))[-5L]
      grid <- c(add, grid); prof <- c(vapply(add, prof_at, numeric(1)), prof)
    } else {
      hi <- hi * 4
      add <- exp(seq(log(grid[length(grid)]), log(hi), length.out = 5L))[-1L]
      grid <- c(grid, add); prof <- c(prof, vapply(add, prof_at, numeric(1)))
    }
    # A magnitude this far from the update's own iterate means the two blocks
    # have come apart; widen a bounded number of times, then say so.
    if (expand == 6L) saturated <- TRUE
  }
  j <- which.max(prof)
  if (j == 1L || j == length(grid)) saturated <- TRUE
  out <- if (j == 1L || j == length(grid)) {
    grid[j]
  } else {
    x <- grid[(j - 1L):(j + 1L)]; y <- prof[(j - 1L):(j + 1L)]
    d <- (y[1] - y[2]) * (x[2] - x[3]) - (y[2] - y[3]) * (x[1] - x[2])
    if (!is.finite(d) || abs(d) < 1e-12) {
      grid[j]
    } else {
      v <- x[2] - 0.5 * ((x[2] - x[1])^2 * (y[2] - y[3]) -
                         (x[2] - x[3])^2 * (y[2] - y[1])) / d
      if (is.finite(v) && v >= x[1] && v <= x[3]) v else grid[j]
    }
  }
  attr(out, "saturated") <- saturated
  out
}

# The loadings by MARGINAL maximum likelihood: the same joint site marginal the
# scale step searches in one dimension, ascended over all S * Q loadings.
#
# The scale step alone is not enough, because the DIRECTION it is handed is
# already over-fit. The factor update holds zeta at its joint mode, so the pair
# (zeta, lambda) is a joint-likelihood estimate with Ns * Q incidental parameters
# growing with the sample -- the Neyman-Scott situation, in which the joint mode
# is inconsistent. The site factors are estimated, their estimation error lands
# in the fitted co-occurrence, and lambda absorbs it. A single scalar cannot
# remove noise variance already baked into the direction: it finds the magnitude
# that is right FOR THE FITTED DIRECTION and so too high for the true one.
#
# Measured on the ms_count fixture (N = 160, Poisson, 6 seeds per cell), mean
# ||lambda_hat||_F / ||lambda_true||_F under the joint-mode estimator, against
# Q / S -- the loadings per species-worth of data, which is what that argument
# says the over-fit should grow with:
#
#   S =  7, Q = 2  (Q/S = 0.286)   1.435   (worst seed 2.229)
#   S = 14, Q = 3  (Q/S = 0.214)   1.076
#   S = 14, Q = 2  (Q/S = 0.143)   1.064
#   S = 28, Q = 2  (Q/S = 0.071)   1.057
#   S = 14, Q = 1  (Q/S = 0.071)   1.014
#
# Monotone in Q / S, so the direction is where the over-fit lives, but NOT
# proportional to it: it is flat around 1.06 across the middle of that range and
# most of the S = 7 excess is two seeds of six. Growing, not a clean law.
#
# Integrating zeta out instead removes the incidental parameters, which is what
# this does. The E-step is the posterior p(z_i | y_i) on the grid the marginal is
# already evaluated on; the M-step maximises the expected complete-data
# log-likelihood, which separates across species into a Qk-dimensional weighted
# Newton per species with the nodes as its design rows. Ascending that ascends
# the marginal (Dempster, Laird & Rubin 1977).
#
# Warm-started at the scale step's answer rather than run from the update's raw
# iterate: the marginal is not concave in lambda, and the 1-D bracket search is a
# cheap GLOBAL pass over the magnitude, which is the direction the update's
# iterate is grossly wrong in. This EM is then a local refinement of both
# magnitude and direction. One objective, two search phases -- not two estimators.
#
# `center` carries the factor update's cross-species loading constraint (used
# when a shared field is also present, so the field owns the shared spatial mean
# and the factors the between-species differences). It is applied after each
# M-step, so the E-step re-adapts to the constrained loadings and the sequence
# stays a projected ascent rather than a projection tacked onto the answer.
.tobs_latent_factor_mmle <- function(oracle, eta_base, lambda, gh,
                                     em.iter = 10L, newton = 2L, tol = 1e-5,
                                     center = FALSE) {
  Ns <- nrow(eta_base); S <- ncol(eta_base); Qk <- ncol(lambda)
  zeta0 <- matrix(0, Ns, Qk)
  # Collapsed loadings carry no latent structure, so the offset is zero rather
  # than absent -- the driver differences it against the previous pass and would
  # fail on a NULL.
  if (all(sqrt(rowSums(lambda^2)) < 1e-8)) {
    return(list(lambda = lambda, zeta = zeta0, converged = TRUE,
                offset = matrix(0, Ns, S), loglik = NA_real_))
  }

  # The expected complete-data log-likelihood at `lam`, over the E-step's fixed
  # nodes and weights. This, not the marginal, is what the M-step line search
  # compares -- the nodes are adapted at the E-step's lambda and held through the
  # M-step, so it is the quantity the Newton step is a step on.
  qobj <- function(lam, Z, W) {
    v <- 0
    for (k in seq_along(Z)) {
      v <- v + sum(W[[k]] *
                     rowSums(oracle$ll_cell(eta_base + tcrossprod(Z[[k]], lam))))
    }
    if (is.finite(v)) v else -Inf
  }

  # The grid is re-adapted every E-step, so successive marginals are not
  # evaluated on identical quadrature and the sequence is only near-monotone.
  # Carry the best rather than the last, so a late quadrature wobble cannot
  # return a worse lambda than one already seen.
  best <- list(ll = -Inf, lambda = lambda, zeta = zeta0)
  ll_prev <- -Inf
  converged <- FALSE

  for (it in seq_len(em.iter)) {
    # ---- E-step: posterior weights over the nodes, and E[z_i | y_i] ----
    g   <- .tobs_latent_joint_grid(oracle, eta_base, lambda, gh)
    m   <- Reduce(pmax, g$logterm)
    tot <- Reduce(`+`, lapply(g$logterm, function(z) exp(z - m)))
    ll  <- sum(g$logdetA + m + log(tot))
    if (!is.finite(ll)) break
    W <- lapply(g$logterm, function(z) exp(z - m) / tot)
    if (!all(vapply(W, function(w) all(is.finite(w)), logical(1)))) break

    zeta <- matrix(0, Ns, Qk)
    for (k in seq_along(W)) zeta <- zeta + W[[k]] * g$Z[[k]]
    if (ll > best$ll) best <- list(ll = ll, lambda = lambda, zeta = zeta)
    if (is.finite(ll_prev) && abs(ll - ll_prev) < tol * (abs(ll) + 1)) {
      converged <- TRUE
      break
    }
    ll_prev <- ll

    # ---- M-step: weighted Newton on lambda_s, species by species ----
    q0 <- qobj(lambda, g$Z, W)
    for (nw in seq_len(newton)) {
      if (!is.finite(q0)) break
      grad <- matrix(0, S, Qk)
      H    <- array(0, dim = c(Qk, Qk, S))
      ok   <- TRUE
      for (k in seq_along(W)) {
        Zk <- g$Z[[k]]; wk <- W[[k]]
        wr <- oracle$working(eta_base + tcrossprod(Zk, lambda))
        if (any(!is.finite(wr$score)) || any(!is.finite(wr$curv))) {
          ok <- FALSE
          break
        }
        # grad[s, a]    = sum_ik w_ik score_isk Z_ika
        # H[a, b, s]    = sum_ik w_ik curv_isk Z_ika Z_ikb
        grad <- grad + crossprod(wr$score, Zk * wk)
        for (a in seq_len(Qk)) for (b in a:Qk) {
          h <- as.numeric(crossprod(wr$curv, wk * Zk[, a] * Zk[, b]))
          H[a, b, ] <- H[a, b, ] + h
          if (b != a) H[b, a, ] <- H[b, a, ] + h
        }
      }
      if (!ok) break

      D <- matrix(0, S, Qk)
      for (s in seq_len(S)) {
        Hs <- matrix(H[, , s], Qk, Qk)
        D[s, ] <- tryCatch(as.numeric(solve(Hs, grad[s, ])), error = function(e)
          tryCatch(as.numeric(solve(Hs + diag(1e-6, Qk), grad[s, ])),
                   error = function(e2) rep(NA_real_, Qk)))
      }
      if (!all(is.finite(D))) break

      # Exp-linked families overshoot, so backtrack on the M-step objective.
      t <- 1; taken <- FALSE
      while (t > 1e-4) {
        cand <- lambda + t * D
        q <- qobj(cand, g$Z, W)
        if (q > q0) { lambda <- cand; q0 <- q; taken <- TRUE; break }
        t <- t / 2
      }
      if (!taken) break
    }
    if (isTRUE(center)) lambda <- sweep(lambda, 2L, colMeans(lambda), "-")
  }

  # The achieved joint marginal at the returned loadings (gcol33/tulpaObs#157):
  # `best$ll` is not necessarily the last E-step's `ll` (the carried-best rule
  # above), so it is surfaced here rather than recomputed, and lets a caller
  # compare two fits (e.g. two directions, or a fit against one reference-
  # started at known truth) on the objective the estimator actually ascends.
  list(lambda = best$lambda, zeta = best$zeta, converged = converged,
       offset = .tobs_latent_factor_offset(oracle, eta_base, best$lambda, gh),
       loglik = best$ll)
}

# The per-cell offset the community EM conditions on.
#
# The driver hands the coefficient update ONE [n_sites x n_species] offset and
# lets it maximise sum_is ll_cell(eta_is(beta) + off_is). The quantity that
# should be maximised is the integrated sum_is E_z[ll_cell(eta_is(beta) +
# lambda_s' z_i)], and NO plug-in point offset reproduces it through a nonlinear
# link: zeta t(lambda) at the posterior mean carries too little latent variance,
# and rescaling the scores to the prior's spread carries too much. Both were
# measured, and each trades one bias for the other -- posterior means put +0.165
# on the community intercept (Jensen: a log link absorbs the missing variance),
# unit-variance scores put the loading magnitude back up to 1.70x with a 4.92x
# tail.
#
# Match the SCORE instead. The two objectives have stationary conditions
#
#   sum_is  score(eta_is + off_is)        d eta_is / d beta = 0     (plug-in)
#   sum_is  E_z[score(eta_is + lam_s'z)]  d eta_is / d beta = 0     (integrated)
#
# so choosing off_is to solve score(eta_is + off_is) = E_z[score(...)] makes them
# the SAME condition, for any family, with no knowledge of the link. The
# expected score comes off the E-step grid already built; the solve is a scalar
# Newton per cell (d score / d eta = -curv, which the oracle also supplies), and
# the block-coordinate loop re-solves it at the current coefficients each pass,
# so it is exact at convergence rather than only near the start.
#
# For a Poisson log link this reduces to off = lambda' zhat + v / 2, the Jensen
# correction, which is the check that the general construction is right.
.tobs_latent_factor_offset <- function(oracle, eta_base, lambda, gh) {
  g <- .tobs_latent_joint_grid(oracle, eta_base, lambda, gh)
  m   <- Reduce(pmax, g$logterm)
  tot <- Reduce(`+`, lapply(g$logterm, function(z) exp(z - m)))
  W   <- lapply(g$logterm, function(z) exp(z - m) / tot)

  # The posterior-mean plug-in offset, E[lambda_s' z_i | y_i]. It needs only the
  # grid, so it is available as a fallback even where the family's working()
  # cannot be evaluated -- and it is exactly the plug-in this step improves on,
  # not a degenerate zero.
  off <- matrix(0, nrow(eta_base), ncol(eta_base))
  for (k in seq_along(W)) off <- off + W[[k]] * tcrossprod(g$Z[[k]], lambda)

  sbar <- matrix(0, nrow(eta_base), ncol(eta_base))
  for (k in seq_along(W)) {
    wr <- oracle$working(eta_base + tcrossprod(g$Z[[k]], lambda))
    if (any(!is.finite(wr$score))) return(off)
    sbar <- sbar + W[[k]] * wr$score
  }

  # Newton on f(off) = score(eta_base + off) - sbar, started at the plug-in mean
  # offset. Steps are capped so a near-flat curvature cannot throw a cell into
  # the region where an exp link overflows.
  for (it in seq_len(30L)) {
    wr <- oracle$working(eta_base + off)
    if (any(!is.finite(wr$score)) || any(!is.finite(wr$curv))) break
    step <- (wr$score - sbar) / pmax(wr$curv, 1e-8)
    step[!is.finite(step)] <- 0
    step <- pmax(pmin(step, 1), -1)
    off <- off + step
    if (max(abs(step)) < 1e-8) break
  }
  off
}


# ---------------------------------------------------------------------------
# Field geometry
# ---------------------------------------------------------------------------

# SPDE (continuous Matern) precision on the mesh at unit scale, for the
# stationary Matern with alpha = nu + d/2 = 2 (nu = 1, d = 2):
#
#   Q(kappa) = kappa^4 C0 + 2 kappa^2 G1 + G1 C0^-1 G1
#
# with C0 the lumped FEM mass diagonal and G1 the stiffness matrix (Lindgren,
# Rue & Lindstrom 2011). The driver's per-field `tau` carries the overall scale
# (the SPDE tau^2), so the field prior is tau Q(kappa) exactly as an areal field's
# is tau Q; kappa is the extra hyperparameter, gridded like proper-CAR's rho.
.tobs_latent_spde_Q <- function(spec, kappa) {
  C0    <- Matrix::Diagonal(x = spec$C0_diag)
  C0inv <- Matrix::Diagonal(x = 1 / spec$C0_diag)
  G1    <- spec$G
  kappa^4 * C0 + 2 * kappa^2 * G1 + G1 %*% C0inv %*% G1
}

# Candidate kappa grid from the observed-coordinate extent. The Matern range is
# rho = sqrt(8 nu) / kappa, so a grid of ranges spanning a tenth to most of the
# domain becomes a kappa grid.
.tobs_latent_spde_kappa_grid <- function(spec) {
  co  <- spec$obs_coords
  ext <- max(apply(co, 2L, function(z) diff(range(z))))
  nu  <- spec$nu %||% 1
  sqrt(8 * nu) / (ext * c(0.1, 0.2, 0.35, 0.5, 0.8))
}

# log|Q| for a proper (full-rank) precision -- the SPDE case. The areal
# counterpart drops the constrained constant direction instead.
.tobs_latent_spde_logdet <- function(Q) {
  as.numeric(Matrix::determinant(Q, logarithm = TRUE)$modulus)
}

# Continuous Matern (spde) field geometry. The mesh nodes carry the field and a
# barycentric projector A [n_sites x n_mesh] maps them onto sites, so the driver
# reads Mt = A (node -> site) and M = t(A) (the adjoint that scatters site scores
# onto nodes) -- the same linear map slot the areal group_var incidence uses.
.tobs_latent_spde_setup <- function(fields, model, what) {
  if (length(fields) > 1L || !is.null(fields[[1L]]$weight)) {
    stop(what, " spde() is the single shared intercept field only (no ",
         "varying-coefficient bar); use icar()/car_proper() for those.",
         call. = FALSE)
  }
  spec <- fields[[1L]]$tulpa_spec
  if (is.null(spec) || is.null(spec$A) || is.null(spec$G) ||
      is.null(spec$C0_diag)) {
    stop(what, " spde() needs the tulpa mesh spec (the projector + FEM ",
         "matrices).", call. = FALSE)
  }
  Ns <- model$n_sites
  if (nrow(spec$A) != Ns) {
    stop(sprintf(paste0("%s spde() projector has %d rows but the model has %d ",
                        "sites; the mesh must be built on the model's ",
                        "coordinates."), what, nrow(spec$A), Ns), call. = FALSE)
  }
  list(graph = NULL, n_nodes = as.integer(spec$n_mesh),
       M = Matrix::t(spec$A), type = "spde", K = 1L,
       W = matrix(1, Ns, 1L), field_labels = "intercept", spec = spec)
}


# Resolve the spatial term(s) on a community formula into the field geometry:
# the per-field per-site weights W [n_sites x K], the field labels, the field
# kind, the adjacency graph, and the node <- site linear map M -- a group_var
# incidence for an areal field with sites > cells, or the barycentric projector
# transpose for a continuous spde mesh. M is NULL when an areal field has one
# node per site. Shared by every community latent fitter.
.tobs_latent_field_setup <- function(spatial, model, what,
                                     allow = c("icar", "car_proper", "bym2")) {
  Ns     <- model$n_sites
  fields <- .tobs_resolve_occu_spatial_fields(spatial, model)
  ptype  <- fields[[1L]]$type %||% "icar"
  if (!all(vapply(fields, function(fd) identical(fd$type %||% "icar", ptype),
                  logical(1))) || !ptype %in% allow) {
    stop(what, " community field supports ",
         paste(sprintf("%s()", allow), collapse = " / "),
         " (one field kind per formula).", call. = FALSE)
  }
  # Continuous Matern field: mesh nodes + projector, no adjacency graph.
  if (identical(ptype, "spde")) {
    return(.tobs_latent_spde_setup(fields, model, what))
  }
  A <- fields[[1L]]$graph
  if (is.null(A)) {
    stop(what, " areal field needs the adjacency graph on the icar() term.",
         call. = FALSE)
  }
  n_nodes <- nrow(A)

  # group_var maps several sites to one field cell; the field has one node per
  # graph cell and the incidence M aggregates the sites onto it.
  Mmap <- NULL
  gv   <- fields[[1L]]$group_var
  if (!is.null(gv)) {
    if (is.null(model$data) || !gv %in% names(model$data)) {
      stop(sprintf("spatial group_var '%s' is not a column of the data.", gv),
           call. = FALSE)
    }
    node_of_site <- as.integer(model$data[[gv]])
    if (length(node_of_site) != Ns || anyNA(node_of_site) ||
        min(node_of_site) < 1L || max(node_of_site) > n_nodes) {
      stop(sprintf(paste0("spatial group_var '%s' must be an integer cell index in ",
                          "1..%d, one per site (%d sites)."), gv, n_nodes, Ns),
           call. = FALSE)
    }
    if (n_nodes != Ns) {
      Mmap <- Matrix::sparseMatrix(i = node_of_site, j = seq_len(Ns), x = 1,
                                   dims = c(n_nodes, Ns))
    }
  } else if (n_nodes != Ns) {
    stop(sprintf(paste0("icar graph has %d nodes but the model has %d sites; add ",
                        "group_var = \"<cell>\" to map sites to cells, or use one ",
                        "node per site."), n_nodes, Ns), call. = FALSE)
  }

  if (identical(ptype, "bym2") &&
      (length(fields) > 1L || !is.null(fields[[1L]]$weight) || !is.null(Mmap))) {
    stop(what, " bym2() field is the single shared intercept field only (no ",
         "varying-coefficient bar / group_var); use icar()/car_proper() for ",
         "those.", call. = FALSE)
  }

  K <- length(fields)
  W <- matrix(1, Ns, K)
  field_labels <- character(K)
  for (k in seq_len(K)) {
    wk <- fields[[k]]$weight
    if (!is.null(wk)) {
      if (length(wk) != Ns) {
        stop(what, " varying-coefficient field weight must be one value per ",
             "site.", call. = FALSE)
      }
      W[, k] <- as.numeric(wk)
      field_labels[k] <- fields[[k]]$weight_label %||% paste0("trend", k - 1L)
    } else {
      field_labels[k] <- "intercept"
    }
  }

  list(graph = A, n_nodes = n_nodes, M = Mmap, type = ptype, K = K, W = W,
       field_labels = field_labels)
}


# ---------------------------------------------------------------------------
# Areal field solve (generic over the family working weights)
# ---------------------------------------------------------------------------

# Newton update for K covariate-weighted areal fields on the structured arm.
# `eta_base` [n_sites x n_species] is the arm's linear predictor from everything
# except these fields (the coefficients, plus the factor part when both are
# present). The per-site field contribution to eta is sum_k W[i,k] F[u(i),k]; per
# field the score aggregates the working residual over species, weighted by
# W[,k]. Joint Newton over the K*n_nodes field vector (a K x K block system, each
# (k,l) block M diag(W[,k] W[,l] curv) M' + (k==l) tau_k Q), then a per-field
# closed-form tau M-step. Each field is demeaned: the intrinsic null space is the
# constant and the fixed effects own the level.
.tobs_latent_field_solve <- function(oracle, eta_base, Q, F, tau, W, M = NULL,
                                     max_iter = 50L, tol = 1e-8,
                                     constrain_mean = TRUE, rankdef = 1L) {
  Ns <- nrow(eta_base); K <- ncol(W); Nn <- nrow(F)
  if (is.null(M)) M <- Matrix::Diagonal(Ns)
  Mt <- Matrix::t(M)
  build_H <- function(cw) {
    blocks <- vector("list", K * K)
    for (k in seq_len(K)) for (l in seq_len(K)) {
      Dw <- M %*% Matrix::Diagonal(x = W[, k] * W[, l] * cw) %*% Mt
      blocks[[(k - 1L) * K + l]] <- if (k == l) Dw + tau[k] * Q else Dw
    }
    do.call(rbind, lapply(seq_len(K), function(k)
      do.call(cbind, blocks[((k - 1L) * K + 1L):(k * K)])))
  }
  field_off_of <- function(Ff) rowSums(W * as.matrix(Mt %*% Ff))

  # The data term is M diag(w) M', whose rank is at most the number of SITES. A
  # continuous mesh can carry MORE nodes than sites, so it is the prior tau Q
  # that makes H full rank -- a small tau then leaves H near-singular and the
  # factorization fails. Retry with a tiny relative ridge (the same guard the
  # BYM2 solve uses for its intrinsic null direction).
  safe_solve <- function(H, g) {
    tryCatch(as.numeric(Matrix::solve(H, g)), error = function(e) {
      d <- mean(Matrix::diag(H))
      r <- if (is.finite(d) && d > 0) 1e-8 * d else 1e-8
      tryCatch(as.numeric(Matrix::solve(H + Matrix::Diagonal(nrow(H), r), g)),
               error = function(e2)
                 as.numeric(Matrix::solve(
                   H + Matrix::Diagonal(nrow(H), max(1e-4 * d, 1e-6)), g)))
    })
  }
  for (it in seq_len(max_iter)) {
    wk <- oracle$working(eta_base + field_off_of(F))
    sc <- rowSums(wk$score); cw <- rowSums(wk$curv)
    g  <- unlist(lapply(seq_len(K), function(k)
      as.numeric(M %*% (W[, k] * sc)) - tau[k] * as.numeric(Q %*% F[, k])))
    step <- safe_solve(build_H(cw), g)
    F <- F + matrix(step, Nn, K)
    if (isTRUE(constrain_mean)) {
      for (k in seq_len(K)) F[, k] <- F[, k] - mean(F[, k])
    }
    if (max(abs(step)) < tol) break
  }

  # tau M-step. The Laplace covariance of the field uses cov_curv when the family
  # stabilises it away from the Newton curvature (an occupancy marginal can drive
  # the curvature to ~0 at undetected sites, which would blow the covariance up).
  eta_full <- eta_base + field_off_of(F)
  cw  <- rowSums(if (is.null(oracle$cov_curv)) oracle$working(eta_full)$curv
                 else oracle$cov_curv(eta_full))
  Hc  <- build_H(cw)
  Cov <- tryCatch(Matrix::solve(Hc), error = function(e) {
    d <- mean(Matrix::diag(Hc))
    Matrix::solve(Hc + Matrix::Diagonal(nrow(Hc),
                                        if (is.finite(d) && d > 0) 1e-8 * d
                                        else 1e-8))
  })
  tau_new <- numeric(K)
  df <- Nn - as.integer(rankdef)
  for (k in seq_len(K)) {
    idx  <- (k - 1L) * Nn + seq_len(Nn)
    quad <- as.numeric(t(F[, k]) %*% (Q %*% F[, k])) +
            sum(Matrix::diag(Q %*% Cov[idx, idx, drop = FALSE]))
    # Floor tau away from zero: with more field nodes than sites the data term
    # alone is rank-deficient, so a collapsing tau would leave the next Newton's
    # H singular.
    tau_new[k] <- max(df / max(quad, 1e-8), 1e-6)
  }
  list(F = F, tau = tau_new)
}


# Marginal (Laplace) objective of the field solve at its converged mode, used to
# pick the proper-CAR rho over a small grid: the family data log-likelihood plus
# the field GMRF prior (0.5 (df log tau + log|Q|) - 0.5 tau f'Q f) minus the
# Laplace normaliser 0.5 log|H|. `logdetQ` is log|Q(rho)| (per unit tau).
.tobs_latent_field_marginal <- function(oracle, eta_base, Q, F, tau, W, logdetQ,
                                        rankdef = 0L, M = NULL) {
  Ns <- nrow(eta_base); K <- ncol(W); Nn <- nrow(F)
  if (is.null(M)) M <- Matrix::Diagonal(Ns)
  Mt <- Matrix::t(M)
  field_off <- rowSums(W * as.matrix(Mt %*% F))
  eta_full  <- eta_base + field_off
  ll  <- oracle$data_ll(eta_full)
  cw  <- rowSums(oracle$working(eta_full)$curv)
  blocks <- vector("list", K * K)
  for (k in seq_len(K)) for (l in seq_len(K)) {
    Dw <- M %*% Matrix::Diagonal(x = W[, k] * W[, l] * cw) %*% Mt
    blocks[[(k - 1L) * K + l]] <- if (k == l) Dw + tau[k] * Q else Dw
  }
  H <- do.call(rbind, lapply(seq_len(K), function(k)
    do.call(cbind, blocks[((k - 1L) * K + 1L):(k * K)])))
  df <- Nn - as.integer(rankdef)
  prior <- 0
  for (k in seq_len(K)) {
    prior <- prior + 0.5 * (df * log(tau[k]) + logdetQ) -
             0.5 * tau[k] * as.numeric(t(F[, k]) %*% (Q %*% F[, k]))
  }
  ll + prior - 0.5 * as.numeric(Matrix::determinant(H, logarithm = TRUE)$modulus)
}


# ---------------------------------------------------------------------------
# BYM2 field solve (generic over the family working weights)
# ---------------------------------------------------------------------------

# Riebler BYM2 scale factor: the geometric mean of the marginal variances of the
# intrinsic ICAR field (its generalised-inverse diagonal), so a unit-scale
# structured component v has geometric-mean marginal variance 1 and rho is
# interpretable as the structured share.
.tobs_latent_bym2_scale <- function(Qic) {
  e   <- eigen(as.matrix(Qic), symmetric = TRUE)
  pos <- e$values > max(e$values) * 1e-10
  vdiag <- rowSums((e$vectors[, pos, drop = FALSE]^2) %*%
                   diag(1 / e$values[pos], sum(pos)))
  exp(mean(log(pmax(vdiag, 1e-12))))
}

# BYM2 field update (one shared field). The combined field phi = a v + b u enters
# eta, with v ~ ICAR (unit precision, sum-to-zero) the structured part and
# u ~ N(0, I) the unstructured part; a = sigma sqrt(rho/scale),
# b = sigma sqrt(1 - rho). Joint Newton over (v, u) at fixed (sigma, rho); the
# hyperparameters are chosen over a grid by the returned field marginal `m`.
.tobs_latent_bym2_solve <- function(oracle, eta_base, Qic, v, u, a, b,
                                    max_iter = 50L, tol = 1e-8) {
  Ns <- nrow(eta_base); Imat <- Matrix::Diagonal(Ns)
  # A tiny ridge regularises the ICAR null direction (a rank-n-1 Qic makes the v
  # block singular when the structured weight a is small); v is demeaned each
  # step, so the ridged constant direction is removed anyway.
  ridge <- Matrix::Diagonal(2L * Ns, 1e-6)
  build_H <- function(cw, a, b) {
    Dw <- Matrix::Diagonal(x = cw)
    rbind(cbind(a * a * Dw + Qic, a * b * Dw),
          cbind(a * b * Dw,       b * b * Dw + Imat)) + ridge
  }
  for (it in seq_len(max_iter)) {
    phi <- a * v + b * u
    wk  <- oracle$working(eta_base + phi)
    sc  <- rowSums(wk$score); cw <- rowSums(wk$curv)
    gv  <- a * sc - as.numeric(Qic %*% v)
    gu  <- b * sc - u
    H   <- build_H(cw, a, b)
    step <- tryCatch(as.numeric(Matrix::solve(H, c(gv, gu))),
                     error = function(e)
                       as.numeric(Matrix::solve(
                         H + Matrix::Diagonal(2L * Ns, 1e-3), c(gv, gu))))
    v <- v + step[seq_len(Ns)]; u <- u + step[Ns + seq_len(Ns)]
    v <- v - mean(v)
    if (max(abs(step)) < tol) break
  }
  # Laplace marginal at the mode (drop constant normalisers): the data
  # log-likelihood in phi, the (v, u) prior quadratics, and -0.5 log|H|.
  phi <- a * v + b * u
  eta_full <- eta_base + phi
  cw  <- rowSums(oracle$working(eta_full)$curv)
  H   <- build_H(cw, a, b)
  m   <- oracle$data_ll(eta_full) -
         0.5 * as.numeric(t(v) %*% (Qic %*% v)) - 0.5 * sum(u * u) -
         0.5 * as.numeric(Matrix::determinant(H, logarithm = TRUE)$modulus)
  list(v = v, u = u, m = m)
}


# log pseudo-determinant of the CAR precision on the sum-to-zero subspace: the
# sum of log of the n-1 largest eigenvalues (the smallest eigenvalue's
# eigenvector is ~ the constant, which the sum-to-zero constraint removes). Used
# to compare the field marginal across a proper-CAR rho grid.
.tobs_latent_car_logdet <- function(Qr) {
  ev <- sort(eigen(as.matrix(Qr), symmetric = TRUE, only.values = TRUE)$values,
             decreasing = TRUE)
  sum(log(pmax(ev[-length(ev)], 1e-10)))
}


# ---------------------------------------------------------------------------
# Latent factor update (generic over the family working weights)
# ---------------------------------------------------------------------------

# Factor update: refine the per-site factors zeta [n_sites x Q] and per-species
# loadings lambda [n_species x Q] against `eta_base` [n_sites x n_species], the
# arm's linear predictor from everything except the factors. Alternating Newton
# (zeta_i | lambda; lambda_s | zeta) with an N(0, I) factor prior and a weak
# loading ridge; the factor columns are centred and scaled to unit variance each
# pass (the scale folds into lambda), a standard factor-analysis anchor.
#
# Site i's working weights depend on zeta[i, ] alone and species s's on
# lambda[s, ] alone, so each inner loop's updates are mutually independent: one
# working() call per half-pass is exact, not a Jacobi approximation.
.tobs_latent_factor_update <- function(oracle, eta_base, zeta, lambda,
                                       inner = 15L, center = FALSE) {
  Ns <- nrow(eta_base); S <- ncol(eta_base); Qk <- ncol(zeta)
  ridge <- 1e-3

  # The penalised objective both half-passes ascend: the data log-likelihood, the
  # N(0, I) factor prior, and the loading ridge. An eta far enough out to overflow
  # an exp link reads as -Inf, so the line search below rejects it.
  obj <- function(z, l) {
    v <- oracle$data_ll(eta_base + tcrossprod(z, l)) -
      0.5 * sum(z * z) - ridge / 2 * sum(l * l)
    if (is.finite(v)) v else -Inf
  }

  # Newton step against the curvature the family supplies, bumped by a small
  # ridge when that curvature is singular.
  nstep <- function(H, g) {
    tryCatch(as.numeric(solve(H, g)), error = function(e)
      tryCatch(as.numeric(solve(H + diag(1e-6, nrow(H)), g)),
               error = function(e2) rep(NA_real_, length(g))))
  }

  # A Newton step is not guaranteed to ascend: the count marginals (N-mixture,
  # distance) are exp-linked, so an overshoot drives exp(eta) past overflow and
  # the family returns a non-finite score, which then propagates into the next
  # solve. Halve the step until it ascends; hold the previous iterate if it never
  # does. Within a half-pass the updates are mutually independent, so one shared
  # scale is exact for whichever step it accepts. A full-length step ascends in
  # the well-behaved case and is applied unchanged.
  ascend <- function(f0, build) {
    t <- 1
    while (t > 1e-4) {
      cand <- build(t)
      f <- obj(cand$zeta, cand$lambda)
      if (f > f0) return(c(cand, list(f = f)))
      t <- t / 2
    }
    NULL
  }

  for (it in seq_len(inner)) {
    # Recomputed per pass: the rescale and the loading centring below both move
    # the objective, so a value carried over from the previous pass is stale.
    f0 <- obj(zeta, lambda)
    if (!is.finite(f0)) break

    wk <- oracle$working(eta_base + tcrossprod(zeta, lambda))
    if (any(!is.finite(wk$score)) || any(!is.finite(wk$curv))) break
    Dz <- matrix(0, Ns, Qk)
    for (i in seq_len(Ns)) {
      g <- as.numeric(crossprod(lambda, wk$score[i, ])) - zeta[i, ]
      H <- crossprod(lambda * sqrt(wk$curv[i, ])) + diag(Qk)
      Dz[i, ] <- nstep(H, g)
    }
    if (all(is.finite(Dz))) {
      cand <- ascend(f0, function(t) list(zeta = zeta + t * Dz, lambda = lambda))
      if (!is.null(cand)) { zeta <- cand$zeta; f0 <- cand$f }
    }

    wk <- oracle$working(eta_base + tcrossprod(zeta, lambda))
    if (any(!is.finite(wk$score)) || any(!is.finite(wk$curv))) break
    Dl <- matrix(0, S, Qk)
    for (s in seq_len(S)) {
      g <- as.numeric(crossprod(zeta, wk$score[, s])) - ridge * lambda[s, ]
      H <- crossprod(zeta * sqrt(wk$curv[, s])) + ridge * diag(Qk)
      Dl[s, ] <- nstep(H, g)
    }
    if (all(is.finite(Dl))) {
      cand <- ascend(f0, function(t) list(zeta = zeta, lambda = lambda + t * Dl))
      if (!is.null(cand)) { lambda <- cand$lambda; f0 <- cand$f }
    }

    for (q in seq_len(Qk)) {
      zeta[, q] <- zeta[, q] - mean(zeta[, q])
      sdq <- stats::sd(zeta[, q])
      if (is.finite(sdq) && sdq > 1e-6) {
        zeta[, q] <- zeta[, q] / sdq; lambda[, q] <- lambda[, q] * sdq
      }
    }
    # When composed with a shared areal field, centre the loadings across species
    # so the field owns the shared spatial mean and the factors own the between-
    # species differences (otherwise a near-constant factor trades off with the
    # field). No-op for the factor-only model.
    if (isTRUE(center)) lambda <- sweep(lambda, 2L, colMeans(lambda), "-")
  }
  list(zeta = zeta, lambda = lambda)
}


# Initial factors / loadings. The factors start at deterministic, mutually
# distinct smooth columns (a fixed seed-free start keeps a fit reproducible); the
# loadings at a small constant, which the first Newton pass spreads.
.tobs_latent_factor_init <- function(latent, Ns, S) {
  Qk <- as.integer(latent$n_factors %||% 1L)
  if (Qk < 1L) stop("latent(): n_factors must be >= 1.", call. = FALSE)
  if (Qk > S - 1L) {
    stop(sprintf("latent(): n_factors (%d) must be < n_species (%d).", Qk, S),
         call. = FALSE)
  }
  zeta <- matrix(0, Ns, Qk)
  for (q in seq_len(Qk)) zeta[, q] <- scale(cos(seq_len(Ns) * q))[, 1]
  list(n_factors = Qk, zeta = zeta, lambda = matrix(0.1, S, Qk))
}

# A second candidate DIRECTION for the outer==1 search (gcol33/tulpaObs#157):
# the top-Q eigenvectors of the coefficient-only working-residual covariance
# across species. This is the classical "principal factor" starting value for
# factor analysis (Lawley & Maxwell 1971) and the residual-correlation-eigen
# start used to initialise latent-factor JSDMs in practice (e.g. the
# per-species-GLM-residual start behind Warton et al. 2015's factor-analytic
# approach, and HMSC's PCA start) -- not a new estimator, a better place to
# hand the existing joint-mode ascent off from.
#
# `working(eta_base)` at zero factor offset is the family-generic quantity
# every oracle already exposes (the Newton step direction score/curv), so this
# needs no family-specific code. Its species x species covariance
# `Sigma_res = V diag(ev) V'` truncated to the top Q components and split as
# `lambda = V[,1:Q] diag(sqrt(ev))`, `zeta = resid V[,1:Q] diag(1/sqrt(ev))`
# reproduces the rank-Q SVD approximation of the residual at this starting
# point -- an ordinary PCA start, deterministic given the data (no random
# rotation to seed). Returns NULL when the decomposition is not well posed
# (fewer than Qk positive eigenvalues, non-finite residuals), so the caller
# falls back to the existing cosine start alone.
.tobs_latent_factor_eigen_init <- function(oracle, eta_base, Qk) {
  wk <- oracle$working(eta_base)
  if (any(!is.finite(wk$score)) || any(!is.finite(wk$curv))) return(NULL)
  # PEARSON residuals (score / sqrt(curv)), not the Newton working-response
  # residual (score / curv). The Newton residual's denominator is the family's
  # curvature itself (mu, for a Poisson log link), so a near-zero fitted mean
  # sends it to +-1/mu and a handful of low-count cells dominate the whole
  # covariance estimate -- measured: it drove the profile magnitude search on
  # the resulting direction to 119x truth (non-convergent), a worse start than
  # the plain cosine one it was meant to improve on. Pearson residuals divide
  # by sqrt(curv) instead (the family's own standard-deviation scale, matching
  # the classical Pearson-residual co-occurrence start used to initialise
  # latent-factor JSDMs), and stay bounded the way an ordinary residual should.
  resid <- wk$score / sqrt(pmax(wk$curv, 1e-8))
  if (!all(is.finite(resid)) || nrow(resid) < 2L) return(NULL)
  resid <- scale(resid, center = TRUE, scale = FALSE)
  Sigma_res <- crossprod(resid) / (nrow(resid) - 1)
  ee <- tryCatch(eigen(Sigma_res, symmetric = TRUE), error = function(e) NULL)
  if (is.null(ee) || length(ee$values) < Qk) return(NULL)
  ev <- ee$values[seq_len(Qk)]
  if (any(!is.finite(ev)) || any(ev <= 1e-8)) return(NULL)
  V <- ee$vectors[, seq_len(Qk), drop = FALSE]
  # `lambda = V`, `zeta = resid %*% V`: the rank-Q PCA scores/loadings pair for
  # the Pearson-residual matrix (V is already orthonormal, so no extra
  # eigenvalue scaling is needed here). This is a DIRECTION only -- the
  # subsequent global magnitude search (.tobs_latent_factor_scale) sets the
  # overall scale on the actual marginal, so the candidate is handed to it
  # directly rather than through the joint-mode ascent (which is the same
  # inconsistent estimator #156 replaced, and re-optimizing a good direction
  # against it can throw the direction away).
  lambda <- V
  zeta <- resid %*% V
  for (q in seq_len(Qk)) {
    zeta[, q] <- zeta[, q] - mean(zeta[, q])
    sdq <- stats::sd(zeta[, q])
    if (is.finite(sdq) && sdq > 1e-6) {
      zeta[, q] <- zeta[, q] / sdq
      lambda[, q] <- lambda[, q] * sdq
    }
  }
  if (!all(is.finite(zeta)) || !all(is.finite(lambda))) return(NULL)
  list(zeta = zeta, lambda = lambda)
}

# K pseudo-random restart DIRECTIONS for the outer==1 candidate search
# (gcol33/tulpaObs#157). A single "smarter" deterministic direction is not
# reliable -- the principal-factor start above helps on some data and is
# measured WORSE than the plain cosine start on others -- so several i.i.d.
# candidates are tried and the honest selector (the converged loading-EM
# marginal, at the call site) picks whichever direction the ascent actually
# reaches a better mode from. Verified against an EM started at the literal
# simulated truth on the fixture that motivated this (ms_count seed 215): one
# of these restarts reached the truth-quality basin (1.28x vs the truth
# start's 1.03x) where the cosine start alone settled at 1.65x.
#
# The seed is FIXED, not random or data-derived, so the set of candidates --
# and so the fit -- stays exactly reproducible run to run. The caller's own
# RNG stream is saved and restored around the draw (the same pattern
# `withr::with_seed` uses), so a `set.seed()` before `tobs()` still reproduces
# everything else about the fit.
.tobs_latent_factor_random_starts <- function(Ns, S, Qk, k = 6L,
                                              seed = 20250157L) {
  has_seed <- exists(".Random.seed", envir = .GlobalEnv)
  if (has_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit({
    if (has_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  })
  set.seed(seed)
  lapply(seq_len(k), function(i) {
    list(zeta = matrix(stats::rnorm(Ns * Qk), Ns, Qk),
         lambda = matrix(stats::rnorm(S * Qk, 0, 0.3), S, Qk),
         ascend = TRUE)
  })
}


# ---------------------------------------------------------------------------
# Block coordinate driver
# ---------------------------------------------------------------------------

# Fit a community model with a shared latent structure: a shared areal field
# (`spatial`), latent factors (`latent`), or BOTH (the spatial-factor case). One
# block coordinate ascent over (a) the community EM with the combined latent as a
# per-species offset, (b) the multi-field areal Newton, (c) the factor update
# (with centred loadings when a field is also present, so the field owns the
# shared spatial mean). Single source of truth for every community latent route
# across every family; the field-only / factor-only models are the special cases.
#
# `oracle_for(eta_base_fn)` is a family callback returning the working oracle;
# `em_fit(site_off, fac_off, em_prev)` runs the family's community EM at the
# given offsets; `offset_of(em)` returns the arm's [n_sites x n_species]
# coefficient predictor.
.tobs_community_latent_ascent <- function(spatial, latent, model, what,
                                          make_oracle, em_fit, offset_of,
                                          allow = c("icar", "car_proper", "bym2"),
                                          tol = 1e-4, max.outer = NULL,
                                          factor.outer = 25L,
                                          n.quad = 5L, verbose = FALSE) {
  S  <- model$n_species
  Ns <- model$n_sites
  has_field  <- !is.null(spatial)
  has_factor <- !is.null(latent)

  # The two blocks converge at very different rates, so the iteration budget is
  # resolved per model rather than shared. A field block reaches `tol` and breaks
  # out early, so its budget is only a cap. The FACTOR block does not: it
  # alternates with the coefficient block along a slow mode (the per-species
  # intercepts and the offset's per-species level both absorb the same latent
  # level), and the per-cell offset change decays about 2% per pass, so `tol`
  # 1e-4 would need roughly 300 of them. Stopping at 25 leaves a real bias in
  # the reported community mean, measured on the ms_count fixture (6 seeds,
  # deviation from the seed's realized mean):
  #
  #   max.outer     25       60      150      400
  #   intercept  +0.0613  +0.0237  +0.0001  -0.0038
  #   slope      +0.0074  +0.0087  +0.0033  +0.0016
  #
  # so a factor model gets `factor.outer`, which each family sets from its OWN
  # measurement rather than inheriting a number measured elsewhere. The curve
  # above is ms_count's; ms_occu recovers at the same budget. The remaining
  # factor families keep 25 until someone measures them, because the cost is not
  # transferable either: ms_count also gained a 2.5x warm start that pays back
  # part of the longer loop, whereas ms_abun_latent already warm-started and so
  # takes the full 6x -- enough to push test-ms-abun-factor.R past 85 minutes on
  # one file. An explicit `control$max.outer` still wins over both.
  if (is.null(max.outer)) max.outer <- if (has_factor) factor.outer else 25L
  max.outer <- as.integer(max.outer)

  # ---- field setup ----
  geom <- NULL; Ffield <- NULL; tau <- NULL; Q <- NULL
  is_bym2 <- FALSE; is_car <- FALSE; is_spde <- FALSE
  rho <- NA_real_; kappa <- NA_real_
  constrain_mean <- TRUE; rankdef <- 1L
  if (has_field) {
    geom <- .tobs_latent_field_setup(spatial, model, what = what, allow = allow)
    is_bym2 <- identical(geom$type, "bym2")
    is_car  <- identical(geom$type, "car_proper")
    is_spde <- identical(geom$type, "spde")
    if (is_spde) {
      # Continuous Matern field on a mesh. Its precision tau Q(kappa) is PROPER
      # (full rank), so unlike an intrinsic areal field it has no constant null
      # space: the field is not demeaned and the M-step keeps the full rank.
      # kappa (the Matern range) is gridded by the field marginal, exactly as
      # proper-CAR's rho is.
      kappa_grid <- .tobs_latent_spde_kappa_grid(geom$spec)
      kappa <- stats::median(kappa_grid)
      Q <- .tobs_latent_spde_Q(geom$spec, kappa)
      constrain_mean <- FALSE; rankdef <- 0L
      Ffield <- matrix(0, geom$n_nodes, geom$K); tau <- rep(1, geom$K)
    } else {
      Dg   <- Matrix::Diagonal(x = rowSums(geom$graph))
      Wadj <- methods::as(geom$graph, "CsparseMatrix")
      rho_grid <- if (is_car) c(0.5, 0.8, 0.95, 0.99) else NA_real_
      rho <- if (is_car) 0.95 else NA_real_
      # The areal field is sum-to-zero (it captures spatial DEVIATIONS; the
      # intercept in X owns the level). icar fixes the dependence at the
      # intrinsic limit (Q = D - W); car_proper estimates the dependence
      # strength rho over a small grid (Q(rho) = D - rho W) by the field
      # marginal likelihood.
      Q <- if (is_car) Dg - rho * Wadj else Dg - Wadj
      Ffield <- matrix(0, geom$n_nodes, geom$K); tau <- rep(1, geom$K)
      if (is_bym2) {
        Qic       <- Dg - Wadj
        bym_scale <- .tobs_latent_bym2_scale(Qic)
        bym_v <- numeric(Ns); bym_u <- numeric(Ns)
        bym_sigma_grid <- c(0.3, 0.6, 1.0, 1.6)
        bym_rho_grid   <- c(0.25, 0.5, 0.75, 0.9)
        bym_sigma <- 0.6; bym_rho <- 0.5; bym_phi <- numeric(Ns)
      }
    }
  }
  Mt_field <- if (has_field && !is.null(geom$M)) Matrix::t(geom$M) else NULL
  site_field_off <- function(Ff) {
    Fs <- if (is.null(Mt_field)) Ff else as.matrix(Mt_field %*% Ff)
    rowSums(geom$W * Fs)
  }
  cur_site_off <- function() {
    if (!has_field) numeric(Ns) else if (is_bym2) bym_phi else site_field_off(Ffield)
  }

  # ---- factor setup ----
  fac <- if (has_factor) .tobs_latent_factor_init(latent, Ns, S) else NULL
  zeta <- fac$zeta; lambda <- fac$lambda
  # `(zeta, lambda)` seeds the first outer pass only. `(zeta_hat, lambda_hat)` is
  # the marginal-likelihood pair, and is the loop's ONLY carried factor state:
  # what the community EM conditions on, what the next pass warm-starts from, and
  # what the fit reports.
  #
  # Running the joint-mode update every pass ALONGSIDE it does not work, in
  # either direction. Writing the refined answer back over the update's state
  # made its next Newton regrow the magnitude to refit the residual while the
  # refinement shrank it again (loadings to 5.3e3 x truth, residual correlation
  # to 0.01). Keeping the two separate instead ratchets: the EM conditions on
  # E[z | y], whose spread is attenuated by the shrinkage a posterior mean
  # carries, the update is not told that and grows lambda to explain the variance
  # the offset did not deliver, and each pass compounds the last (measured on the
  # S = 7 fixture: |lambda| 1.98 -> 2.94 -> 4.29 -> ... -> 925, and because the
  # community EM absorbs the inflated offset into the coefficients, the runaway
  # pair is locally self-consistent and nothing downstream rejects it). One
  # estimator, one state.
  lambda_hat <- lambda
  zeta_hat   <- zeta
  # The score-matched per-cell offset the community EM conditions on, and the
  # loop's convergence target. Zero until the first factor pass fills it.
  fac_offset <- matrix(0, Ns, S)
  gh_joint <- if (has_factor) .tobs_gh_nodes(n.quad) else NULL
  fac_saturated <- FALSE
  fac_em_converged <- FALSE
  fac_loglik <- NA_real_

  em <- NULL
  for (outer in seq_len(max.outer)) {
    site_off <- cur_site_off()
    fac_off  <- fac_offset
    # (a) community EM given the combined latent offset.
    em <- em_fit(site_off, fac_off, em)
    eta_coef <- offset_of(em)
    oracle   <- make_oracle(em)
    delta <- 0
    # (b) field update (its base holds the coefficients + the factor part).
    if (has_field && is_bym2) {
      a  <- bym_sigma * sqrt(bym_rho / bym_scale); b <- bym_sigma * sqrt(1 - bym_rho)
      br <- .tobs_latent_bym2_solve(oracle, eta_coef + fac_off, Qic,
                                    bym_v, bym_u, a, b)
      phi_new <- a * br$v + b * br$u
      delta <- max(delta, max(abs(phi_new - bym_phi)))
      bym_v <- br$v; bym_u <- br$u; bym_phi <- phi_new
    } else if (has_field) {
      fu <- .tobs_latent_field_solve(oracle, eta_coef + fac_off, Q, Ffield, tau,
                                     geom$W, M = geom$M,
                                     constrain_mean = constrain_mean,
                                     rankdef = rankdef)
      delta <- max(delta, max(abs(fu$F - Ffield))); Ffield <- fu$F; tau <- fu$tau
    }
    # (c) factor update (its base holds the coefficients + the field part).
    if (has_factor) {
      base_fac <- eta_coef + matrix(cur_site_off(), Ns, S)
      if (outer == 1L) {
        # Initialise only. The marginal's gradient in lambda vanishes at
        # lambda = 0, so the loading EM cannot start from the zero init; the
        # joint-mode update supplies a non-degenerate direction cheaply, and the
        # 1-D scale search puts its magnitude in the right basin with a global
        # bracket the local EM has no way to perform.
        #
        # gcol33/tulpaObs#157: that ascent is a LOCAL search and inherits
        # whatever basin its starting DIRECTION falls into -- measured, 1 seed
        # in 16 on the ms_count fixture (seed 215) settled 31 nats below what
        # the same ascent reaches from a better direction, with the residual
        # correlation reading 0.90 throughout (row-normalised, blind to it).
        # Rescaling the magnitude of a bad direction does not escape it
        # (measured: 0.001 movement at 3x cost). Escaping needs a different
        # DIRECTION, and a single "smarter" one is not reliable enough on its
        # own: a principal-factor start off the coefficient-only
        # working-residual covariance (.tobs_latent_factor_eigen_init) beat
        # the cosine start's basin on some seeds but was measured WORSE on
        # seed 215 itself (final 7.50x truth against the cosine start's own
        # 1.65x) -- a single alternative direction can just as easily be a
        # worse one. What actually escaped seed 215's basin, verified against
        # an EM started at the literal simulated truth (which reaches 1.03x,
        # the reachable optimum the issue measured), was trying several fixed
        # pseudo-random restarts alongside it and keeping whichever the SAME
        # ascent + scale search + loading EM lands highest on the marginal --
        # the "multi-start over directions, selected on the marginal" the
        # issue names as the honest alternative. On seed 215 the winning
        # random restart reached 701212 against the cosine start's 701195 and
        # the truth start's 701211 -- AT the truth-quality basin, at 1.28x
        # rather than 1.65x.
        #
        # The comparison has to be the CONVERGED loading-EM marginal, not the
        # raw ascent + scale-search value: on seed 215 the candidates' raw
        # scale-search marginals differed by under 0.04% of their scale (all
        # far from any mode) and did not rank the same as their converged
        # values. Running each candidate to its OWN loading-EM convergence and
        # comparing THAT marginal is the same comparison the issue itself used
        # to diagnose the basin (an EM started at truth reaching -4676.7
        # against the seed's -4720.9), so it is what "selected on the
        # marginal" has to mean here. This costs one extra loading-EM run per
        # candidate, on the first outer pass only -- negligible against the
        # `factor.outer` outer passes of community-EM refitting that follow
        # (measured: +30s against an ~90s fit, for 8 extra candidates).
        #
        # The restarts use a FIXED internal seed, not a random or data-derived
        # one, so a fit stays exactly reproducible run to run; the caller's own
        # RNG stream is saved and restored around the draw so a `set.seed()`
        # before `tobs()` still reproduces everything else about the fit.
        starts <- list(list(zeta = zeta, lambda = lambda, ascend = TRUE))
        eig <- .tobs_latent_factor_eigen_init(oracle, base_fac, ncol(lambda))
        if (!is.null(eig)) starts <- c(starts, list(c(eig, list(ascend = FALSE))))
        starts <- c(starts,
                   .tobs_latent_factor_random_starts(Ns, S, ncol(lambda)))

        cand <- NULL
        for (st in starts) {
          gu  <- if (isTRUE(st$ascend)) {
            .tobs_latent_factor_update(oracle, base_fac, st$zeta, st$lambda,
                                       center = has_field)
          } else {
            st
          }
          sc   <- .tobs_latent_factor_scale(oracle, base_fac, gu$lambda, gh_joint)
          lam0 <- gu$lambda * as.numeric(sc)
          sat  <- isTRUE(attr(sc, "saturated"))
          mm_st <- .tobs_latent_factor_mmle(oracle, base_fac, lam0, gh_joint,
                                            center = has_field)
          m <- mm_st$loglik
          if (is.null(cand) ||
              (is.finite(m) && (!is.finite(cand$m) || m > cand$m))) {
            cand <- list(m = m, mm = mm_st, saturated = sat)
          }
        }
        if (isTRUE(cand$saturated)) fac_saturated <- TRUE
        mm <- cand$mm
      } else {
        mm <- .tobs_latent_factor_mmle(oracle, base_fac, lambda_hat, gh_joint,
                                       center = has_field)
      }
      lambda_hat <- mm$lambda
      zeta_hat   <- mm$zeta
      fac_offset <- mm$offset
      fac_em_converged <- isTRUE(mm$converged)
      fac_loglik <- mm$loglik
      delta <- max(delta, max(abs(fac_offset - fac_off)))
    }
    if (isTRUE(verbose)) {
      # The loading magnitude and the spread of the posterior factor scores. The
      # pair is what shows the two blocks ratcheting against each other, which is
      # the failure mode this loop is arranged to avoid.
      message(sprintf("[%s latent %d] delta=%.2e%s", what, outer, delta,
                      if (has_factor)
                        sprintf("  |lambda_hat|=%.3f  sd(zeta_hat)=%.3f",
                                sqrt(sum(lambda_hat^2)),
                                stats::sd(as.numeric(zeta_hat)))
                      else ""))
    }
    if (outer > 2L && delta < tol) break
  }

  # ---- field hyperparameter selection at the converged coefficients ----
  if (has_field) {
    eta_coef <- offset_of(em)
    oracle   <- make_oracle(em)
    base_f   <- eta_coef + (if (has_factor) fac_offset else 0)
    if (is_car) {
      best <- list(m = -Inf)
      for (rr in rho_grid) {
        Qr <- Dg - rr * Wadj
        ld <- .tobs_latent_car_logdet(Qr)
        fr <- .tobs_latent_field_solve(oracle, base_f, Qr, Ffield, tau, geom$W,
                                       M = geom$M, constrain_mean = TRUE,
                                       rankdef = 1L)
        mm <- .tobs_latent_field_marginal(oracle, base_f, Qr, fr$F, fr$tau,
                                          geom$W, ld, rankdef = 1L, M = geom$M)
        if (is.finite(mm) && mm > best$m) {
          best <- list(m = mm, rho = rr, F = fr$F, tau = fr$tau)
        }
      }
      if (is.finite(best$m)) { rho <- best$rho; Ffield <- best$F; tau <- best$tau }
    }
    if (is_spde) {
      # Pick the Matern range (kappa) on the grid by the field marginal. Q(kappa)
      # is proper, so its log-determinant is the full one and the M-step keeps
      # rank n_mesh (no constrained constant direction to drop).
      best <- list(m = -Inf)
      for (kk in kappa_grid) {
        Qk <- .tobs_latent_spde_Q(geom$spec, kk)
        ld <- .tobs_latent_spde_logdet(Qk)
        fr <- tryCatch(
          .tobs_latent_field_solve(oracle, base_f, Qk, Ffield, tau, geom$W,
                                   M = geom$M, constrain_mean = FALSE,
                                   rankdef = 0L),
          error = function(e) NULL)
        if (is.null(fr)) next
        mm <- .tobs_latent_field_marginal(oracle, base_f, Qk, fr$F, fr$tau,
                                          geom$W, ld, rankdef = 0L, M = geom$M)
        if (is.finite(mm) && mm > best$m) {
          best <- list(m = mm, kappa = kk, F = fr$F, tau = fr$tau)
        }
      }
      if (is.finite(best$m)) {
        kappa <- best$kappa; Ffield <- best$F; tau <- best$tau; Q <- Qk
      }
    }
    if (is_bym2) {
      best <- list(m = -Inf)
      for (sg in bym_sigma_grid) for (rr in bym_rho_grid) {
        a  <- sg * sqrt(rr / bym_scale); b <- sg * sqrt(1 - rr)
        br <- tryCatch(.tobs_latent_bym2_solve(oracle, base_f, Qic, bym_v, bym_u,
                                               a, b),
                       error = function(e) NULL)
        if (!is.null(br) && is.finite(br$m) && br$m > best$m) {
          best <- c(br, list(sigma = sg, rho = rr, a = a, b = b))
        }
      }
      if (is.finite(best$m)) {
        bym_v <- best$v; bym_u <- best$u
        bym_sigma <- best$sigma; bym_rho <- best$rho
        bym_phi <- best$a * best$v + best$b * best$u
      }
    }
  }

  list(em = em, geom = geom,
       field = if (!has_field) NULL else list(
         F = Ffield, tau = tau, type = geom$type, rho = rho, kappa = kappa,
         site_off = cur_site_off(),
         bym2 = if (is_bym2) list(sigma = bym_sigma, rho = bym_rho,
                                  scale = bym_scale, phi = bym_phi) else NULL),
       factor = if (!has_factor) NULL else list(
         n_factors = fac$n_factors, zeta = zeta_hat, lambda = lambda_hat,
         offset = fac_offset,
         # TRUE when the initialising magnitude search hit the end of even its
         # widened bracket, i.e. the loadings the EM was started from were a
         # boundary and not an argmax. Every quantity the family reports off the
         # factors -- the loadings and so the residual COVARIANCE -- is then
         # suspect in scale, though the residual CORRELATION is row-normalised
         # and survives.
         magnitude_saturated = fac_saturated,
         # TRUE when the last pass's loading EM met its own tolerance rather than
         # running out of iterations.
         loading_em_converged = fac_em_converged,
         # The joint marginal at (lambda_hat, zeta_hat) from the last factor
         # pass (gcol33/tulpaObs#157): lets a caller compare two fits (e.g. two
         # starting directions) on the objective the estimator ascends, which
         # is otherwise invisible -- the residual correlation is row-normalised
         # and a magnitude regression does not move it.
         marginal_loglik = fac_loglik))
}


# ---------------------------------------------------------------------------
# Result assembly
# ---------------------------------------------------------------------------

# Attach the shared-field summary of a block-coordinate latent fit onto `fit`.
# `offset_slot` names the model slot holding the per-site eta contribution that
# fitted() / WAIC add back (e.g. "count_field_offset").
.tobs_latent_attach_field <- function(fit, res, spatial, offset_slot) {
  fd <- res$field
  if (is.null(fd)) return(fit)
  fit$spatial <- spatial
  fit$method  <- "nested_laplace"
  if (!is.null(fd$bym2)) {
    # The reported field is the combined BYM2 field phi = a v + b u.
    fit$spatial_field <- as.numeric(fd$bym2$phi)
    fit$model[[offset_slot]] <- as.numeric(fd$bym2$phi)
    fit$spatial_hyper <- list(type = "bym2", sigma = fd$bym2$sigma,
                              rho = fd$bym2$rho, scale = fd$bym2$scale,
                              field_labels = "intercept")
    fit$means <- c(fit$means, sigma_field_intercept = fd$bym2$sigma)
    return(fit)
  }
  labels <- res$geom$field_labels
  # spatial_field is the per-node intercept field; the offset slot is the
  # per-site eta contribution (the node field mapped through the site -> cell map).
  fit$spatial_field <- as.numeric(fd$F[, 1L])
  fit$model[[offset_slot]] <- fd$site_off
  sigma_field <- 1 / sqrt(fd$tau)
  fit$spatial_hyper <- list(type = fd$type, tau = fd$tau, sigma = sigma_field,
                            rho = fd$rho, field_labels = labels)
  # A continuous Matern (spde) field additionally reports its range: the mesh
  # field lives on n_mesh nodes, and rho = sqrt(8 nu) / kappa is the practical
  # correlation range in the coordinate units.
  if (identical(fd$type, "spde") && is.finite(fd$kappa %||% NA_real_)) {
    nu <- res$geom$spec$nu %||% 1
    fit$spatial_hyper$kappa <- fd$kappa
    fit$spatial_hyper$range <- sqrt(8 * nu) / fd$kappa
    fit$means <- c(fit$means, spde_range = sqrt(8 * nu) / fd$kappa)
  }
  fit$means <- c(fit$means,
                 stats::setNames(sigma_field, paste0("sigma_field_", labels)))
  if (res$geom$K > 1L) {
    fit$trend_fields <- lapply(2:res$geom$K, function(k) as.numeric(fd$F[, k]))
    names(fit$trend_fields) <- labels[-1L]
    fit$trend_field <- fit$trend_fields[[1L]]
  }
  fit
}

# Attach the latent-factor summary. The loadings and factors are identified only
# up to rotation, so the recoverable target is the residual species covariance
# Sigma_res = lambda lambda' (reported alongside the implied correlation).
.tobs_latent_attach_factor <- function(fit, res, latent, model, offset_slot) {
  fc <- res$factor
  if (is.null(fc)) return(fit)
  fit$latent <- latent
  lambda <- fc$lambda
  Sigma_res <- tcrossprod(lambda)
  dimnames(Sigma_res) <- list(model$species_names, model$species_names)
  rownames(lambda) <- model$species_names
  colnames(lambda) <- paste0("factor", seq_len(fc$n_factors))
  S <- nrow(lambda)
  fit$ms_factor <- list(
    n_factors = fc$n_factors, loadings = lambda, factors = fc$zeta,
    residual_cov = Sigma_res,
    residual_cor = stats::cov2cor(Sigma_res + diag(1e-10, S)),
    # Scale diagnostics. residual_cor is row-normalised and so cannot show a
    # magnitude problem; these can.
    magnitude_saturated = isTRUE(fc$magnitude_saturated),
    loading_em_converged = isTRUE(fc$loading_em_converged),
    marginal_loglik = fc$marginal_loglik)
  # The score-matched offset the fit conditioned on, not zeta t(lambda): fitted()
  # and WAIC read this slot, and they should see the predictor the coefficients
  # were estimated against.
  fit$model[[offset_slot]] <- fc$offset
  fit
}
