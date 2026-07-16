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
      stop(sprintf("spatial group_var '%s' must be an integer cell index in ",
                   "1..%d, one per site (%d sites).", gv, n_nodes, Ns),
           call. = FALSE)
    }
    if (n_nodes != Ns) {
      Mmap <- Matrix::sparseMatrix(i = node_of_site, j = seq_len(Ns), x = 1,
                                   dims = c(n_nodes, Ns))
    }
  } else if (n_nodes != Ns) {
    stop(sprintf("icar graph has %d nodes but the model has %d sites; add ",
                 "group_var = \"<cell>\" to map sites to cells, or use one ",
                 "node per site.", n_nodes, Ns), call. = FALSE)
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
                                          tol = 1e-4, max.outer = 25L,
                                          verbose = FALSE) {
  S  <- model$n_species
  Ns <- model$n_sites
  has_field  <- !is.null(spatial)
  has_factor <- !is.null(latent)

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

  em <- NULL
  for (outer in seq_len(max.outer)) {
    site_off <- cur_site_off()
    fac_off  <- if (has_factor) tcrossprod(zeta, lambda) else matrix(0, Ns, S)
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
      gu <- .tobs_latent_factor_update(
        oracle, eta_coef + matrix(cur_site_off(), Ns, S), zeta, lambda,
        center = has_field)
      delta <- max(delta, max(abs(tcrossprod(gu$zeta, gu$lambda) - fac_off)))
      zeta <- gu$zeta; lambda <- gu$lambda
    }
    if (isTRUE(verbose)) {
      message(sprintf("[%s latent %d] delta=%.2e", what, outer, delta))
    }
    if (outer > 2L && delta < tol) break
  }

  # ---- field hyperparameter selection at the converged coefficients ----
  if (has_field) {
    eta_coef <- offset_of(em)
    oracle   <- make_oracle(em)
    base_f   <- eta_coef + (if (has_factor) tcrossprod(zeta, lambda) else 0)
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
         n_factors = fac$n_factors, zeta = zeta, lambda = lambda))
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
    residual_cor = stats::cov2cor(Sigma_res + diag(1e-10, S)))
  fit$model[[offset_slot]] <- tcrossprod(fc$zeta, lambda)
  fit
}
