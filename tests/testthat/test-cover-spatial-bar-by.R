# Replicated CAR via spatial(by=) on the cover hurdle (gcol33/tulpaObs#82).
#
# `by = "grp"` on a spatial bar replicates the whole field across the levels of a
# factor: the graph becomes the block-diagonal Kronecker I_L (x) Q (L disjoint
# copies) and each observation's node is offset into its level's copy, with the
# field hyperparameters SHARED across levels (one sigma[, rho_car]; one Sigma for
# a correlated `|`). It is orthogonal to the bar character and to `to`, so it
# composes with all three cover bar forms -- shared (both-arm `||`), correlated
# (`|`), and arm-specific (single-arm `||`). The engine remap is
# tulpa::tulpa_bar_field_replicate(); these tests pin the consumer wiring that
# threads it through each cover path.

# Rook-adjacency on a g x g grid (self-contained).
.by_grid_adj <- function(g) {
  n <- g * g
  co <- expand.grid(r = seq_len(g), c = seq_len(g))
  adj <- matrix(0L, n, n)
  for (i in seq_len(n)) for (j in seq_len(n)) {
    if (i < j && abs(co$r[i] - co$r[j]) + abs(co$c[i] - co$c[j]) == 1L) {
      adj[i, j] <- 1L; adj[j, i] <- 1L
    }
  }
  adj
}

# A captured bar spec, exactly as the tobs formula machinery builds it: the bar
# defaults to both arms, and the arm tag (spec$to) is set afterwards the way
# placement does.
.by_spec <- function(bar, adj, to = c("presence", "positive"), by = NULL,
                     model = "icar") {
  rest <- list(graph = adj)
  if (!is.null(by)) rest$by <- by
  spec <- .tobs_spatial_bar_spec(bar, rest, model = model)
  spec$to <- to
  spec
}


# ---- capture: `by` is taken as a string column name ------------------------

test_that("spatial(<bar>, by=) captures the factor column name", {
  adj <- .by_grid_adj(3L)
  spec <- .by_spec(~ 1 + x || cell, adj, by = "grp")
  expect_identical(spec$by_var, "grp")
  # Absent by is NULL (the no-replication identity).
  expect_null(.by_spec(~ 1 + x || cell, adj)$by_var)
})

test_that("spatial(<bar>, by=) rejects a non-string `by`", {
  adj <- .by_grid_adj(3L)
  expect_error(.by_spec(~ 1 || cell, adj, by = 1L),
               "must be a single replication-factor column name")
  expect_error(.by_spec(~ 1 || cell, adj, by = c("a", "b")),
               "must be a single replication-factor column name")
})


# ---- wiring: each cover path replicates the graph + offsets the index -------

test_that("arm-specific `||` + by builds the I_L (x) Q graph and offset index", {
  g <- 3L; n <- g * g; L <- 2L
  adj <- .by_grid_adj(g)
  N <- 40L
  d <- data.frame(cell = rep_len(seq_len(n), N),
                  grp  = rep_len(c("a", "b"), N),
                  x    = rnorm(N))
  spec <- .by_spec(~ 1 + x || cell, adj, to = "positive", by = "grp")
  af <- .tobs_armspecific_bar_fields(spec, d)

  expect_equal(dim(af$graph), c(L * n, L * n))      # block-diagonal I_2 (x) Q
  expect_identical(af$by$n_levels, L)
  expect_identical(af$by$n_nodes, n)
  # Level-b observations sit in the second graph copy (node + n); level-a in the
  # first. The offset index never exceeds L*n and is >= 1.
  lvl <- as.integer(factor(d$grp))
  expect_equal(af$idx_obs, as.integer(d$cell) + (lvl - 1L) * n)
  expect_true(min(af$idx_obs) >= 1L && max(af$idx_obs) <= L * n)
  # A single-level by is the identity: base graph, raw index.
  d1 <- d; d1$grp <- "only"
  af1 <- .tobs_armspecific_bar_fields(.by_spec(~ 1 || cell, adj, to = "positive",
                                               by = "grp"), d1)
  expect_equal(dim(af1$graph), c(n, n))
  expect_equal(af1$idx_obs, as.integer(d1$cell))
})

test_that("correlated `|` + by builds the replicated MCAR graph + offset index", {
  g <- 3L; n <- g * g; L <- 2L
  adj <- .by_grid_adj(g)
  N <- 40L
  d <- data.frame(cell = rep_len(seq_len(n), N),
                  grp  = rep_len(c("a", "b"), N),
                  x    = rnorm(N))
  spec <- .by_spec(~ 1 + x | cell, adj, by = "grp", model = "icar")
  mc <- .cover_build_mcar_spec(spec, d)

  expect_identical(mc$n_spatial_units, as.integer(L * n))
  expect_equal(dim(mc$graph), c(L * n, L * n))
  expect_identical(mc$by$n_levels, L)
  lvl <- as.integer(factor(d$grp))
  expect_equal(mc$idx_occ, as.integer(d$cell) + (lvl - 1L) * n)
  # Sigma is the cross-FIELD covariance (p x p), unchanged by replication.
  expect_identical(mc$n_fields, 2L)
})

test_that("shared `||` + by carries by_var onto the desugared + tulpa specs", {
  adj <- .by_grid_adj(3L)
  d <- data.frame(cell = rep_len(seq_len(9L), 30L),
                  grp  = rep_len(c("a", "b"), 30L), x = rnorm(30L))
  spec <- .by_spec(~ 1 + x || cell, adj, by = "grp")
  terms <- .tobs_expand_spatial_bar(spec, d)
  # Every desugared term carries by_var (the fit replicates with the data in
  # scope); the intercept term's tulpa_spatial spec propagates it.
  expect_true(all(vapply(terms, function(t) identical(t$by_var, "grp"), logical(1))))
  tsp <- .tobs_term_to_tulpa_spatial(terms[[1L]])
  expect_identical(tsp$by_var, "grp")
})


# ---- recovery: per-level fields recover, sharing one sigma ------------------

test_that("shared `||` + by recovers independent per-level fields", {
  skip_on_cran()
  skip_if_fast()
  set.seed(11)
  g <- 6L; n <- g * g
  adj <- .by_grid_adj(g)
  # Two genuinely different per-level fields over the same graph; one shared
  # amplitude. A by fit must recover BOTH (a single non-replicated field cannot).
  smooth <- function(seed) {
    set.seed(seed); z <- rnorm(n)
    for (it in 1:60) {
      zn <- z
      for (i in seq_len(n)) {
        nb <- which(adj[i, ] == 1L); zn[i] <- 0.35 * z[i] + 0.65 * mean(z[nb])
      }
      z <- zn
    }
    z <- z - mean(z); z / stats::sd(z)
  }
  z_a <- smooth(1); z_b <- smooth(2)
  N <- 6000L
  cell <- sample.int(n, N, replace = TRUE)
  grp  <- sample(c("a", "b"), N, replace = TRUE)
  x    <- as.numeric(scale(rnorm(N)))
  sigma_field <- 0.9
  z_obs <- ifelse(grp == "a", z_a[cell], z_b[cell])
  beta_occ <- c(0.2, 0.3); beta_pos <- c(-0.8, 0.2)
  eta_occ <- beta_occ[1] + beta_occ[2] * x + sigma_field * z_obs
  occur   <- rbinom(N, 1, plogis(eta_occ))
  eta_pos <- beta_pos[1] + beta_pos[2] * x + sigma_field * z_obs
  cover   <- ifelse(occur == 1L,
                    pmin(exp(eta_pos + rnorm(N, 0, 0.35)), 1 - 1e-6), 0)
  dat <- data.frame(cell = cell, grp = grp, x = x, cover = cover)

  ctrl <- list(verbose = FALSE, progress = FALSE, integration = "grid",
               sigma.grid = exp(seq(log(0.3), log(2.5), length.out = 6)))
  fit <- tobs(formula = ~ x + spatial(~ 1 || cell, graph = adj, by = "grp"),
              data = dat, family = cover(response = "lognormal"),
              y = dat$cover, method = "nested_laplace", control = ctrl)

  expect_s3_class(fit, "cover_fit")
  bundle <- .tobs_joint_draws(fit, n = 400L)
  z_hat  <- colMeans(bundle$blocks[[1L]]$z)            # over L*n = 2n nodes
  expect_length(z_hat, 2L * n)
  za_hat <- z_hat[seq_len(n)]; zb_hat <- z_hat[n + seq_len(n)]
  za_hat <- za_hat - mean(za_hat); zb_hat <- zb_hat - mean(zb_hat)
  expect_gt(cor(za_hat, z_a - mean(z_a)), 0.8)
  expect_gt(cor(zb_hat, z_b - mean(z_b)), 0.8)
  # Independent realizations stay separate (no cross-level leakage).
  expect_lt(abs(cor(za_hat, z_b - mean(z_b))), 0.4)
})
