# SLA on the joint nested-Laplace cover-hurdle path.
#
# These tests are written ahead of the implementation so the API surface
# is pinned. They mirror the existing standalone-Laplace SLA tests in
# test-sla-cover-hurdle.R but exercise the method = "nested_laplace" path
# (BYM2-spatial joint fit). See dev_notes/sla_joint_proposal.md for the
# spec these tests pin.

suppressPackageStartupMessages({
    library(testthat)
    library(tulpaObs)
    library(tulpa)
})


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

# Chain adjacency matrix on n_s spatial units (1-2-3-...-n_s).
chain_adj <- function(n_s) {
    adj <- matrix(0L, n_s, n_s)
    if (n_s < 2L) return(adj)
    for (s in seq_len(n_s - 1L)) {
        adj[s, s + 1L] <- 1L
        adj[s + 1L, s] <- 1L
    }
    adj
}


# Simulate a small cover-hurdle dataset with a BYM2 spatial field on a
# chain graph. The phi / theta fields are demeaned so they live in the
# sum-to-zero subspace assumed by the model's prior — this keeps the
# truth aligned with the constrained inner posterior, so the validation
# is apples-to-apples (mirrors the demean fix used in dev_notes d3 sims).
.make_cover_data <- function(seed,
                             alpha_true = 1.0,
                             N          = 200,
                             n_s        = 25,
                             sigma      = 0.5,
                             rho        = 0.7,
                             beta_occ   = c(-0.3, 0.7),
                             beta_pos   = c(0.4, -0.5),
                             phi_disp   = 30) {
    set.seed(seed)
    spatial_idx <- sample.int(n_s, N, replace = TRUE)

    phi_f   <- rnorm(n_s, 0, 1); phi_f   <- phi_f   - mean(phi_f)
    theta_f <- rnorm(n_s, 0, 1); theta_f <- theta_f - mean(theta_f)
    w_s     <- sigma * (sqrt(rho) * phi_f + sqrt(1 - rho) * theta_f)

    x <- rnorm(N)
    eta_occ <- beta_occ[1] + beta_occ[2] * x + w_s[spatial_idx]
    occur   <- rbinom(N, 1, plogis(eta_occ))

    eta_pos <- beta_pos[1] + beta_pos[2] * x +
                   alpha_true * w_s[spatial_idx]
    mu_pos  <- plogis(eta_pos)
    y       <- numeric(N)
    is_pos  <- occur == 1L
    if (any(is_pos)) {
        y[is_pos] <- rbeta(
            sum(is_pos),
            mu_pos[is_pos] * phi_disp,
            (1 - mu_pos[is_pos]) * phi_disp
        )
    }
    y[!is_pos] <- 0
    y <- pmin(pmax(y, 0), 1 - 1e-6)

    list(
        data = data.frame(x = x, region = factor(spatial_idx,
                                                  levels = seq_len(n_s))),
        y    = y,
        adj  = chain_adj(n_s),
        truth = list(beta_occ = beta_occ, beta_pos = beta_pos,
                     sigma = sigma, rho = rho, alpha = alpha_true,
                     phi_disp = phi_disp)
    )
}


# ---------------------------------------------------------------------------
# 1. No fallback message under simplified_laplace
# ---------------------------------------------------------------------------

test_that("joint SLA path no longer falls back via message", {
  skip_if_fast()
    set.seed(101)
    sim <- .make_cover_data(seed = 101, N = 200, n_s = 25)

    msgs <- character(0)
    fit <- withCallingHandlers(
        tobs(
            formula = ~ x + bym2(graph = sim$adj, group_var = "region"),
            data    = sim$data,
            family  = cover("beta"),
            y       = sim$y,
            method  = "nested_laplace_sla",
            control = list(
                sigma.grid     = c(0.4, 0.8),
                rho.grid       = c(0.5, 0.9)
            )
        ),
        message = function(m) {
            msgs <<- c(msgs, conditionMessage(m))
            invokeRestart("muffleMessage")
        }
    )

    expect_s3_class(fit, "cover_fit")
    expect_false(any(grepl("falling back", msgs, fixed = TRUE)),
                 info = paste(msgs, collapse = " | "))
    expect_false(any(grepl("currently wired only", msgs, fixed = TRUE)),
                 info = paste(msgs, collapse = " | "))
})


# ---------------------------------------------------------------------------
# 2. API surface returns expected fields
# ---------------------------------------------------------------------------

test_that("SLA joint fit exposes skew + draws fields", {
  skip_if_fast()
    set.seed(102)
    sim <- .make_cover_data(seed = 102, N = 200, n_s = 25)

    fit <- suppressMessages(tobs(
        formula = ~ x + bym2(graph = sim$adj, group_var = "region"),
        data    = sim$data,
        family  = cover("beta"),
        y       = sim$y,
        method  = "nested_laplace_sla",
        control = list(
            sigma.grid     = c(0.4, 0.8),
            rho.grid       = c(0.5, 0.9)
        )
    ))

    # Status field present and well-typed.
    expect_true(is.character(fit$sla_status))
    expect_equal(length(fit$sla_status), 1L)
    expect_true(grepl("^simplified_laplace|^fallback_gaussian",
                       fit$sla_status),
                info = paste0("sla_status = '", fit$sla_status, "'"))

    # Per-coefficient skewness vectors, length = ncol(X) of each arm.
    p_occ <- ncol(fit$encoding$occ_data$X)
    p_pos <- ncol(fit$encoding$pos_data$X)

    expect_true(is.numeric(fit$skew_occ))
    expect_equal(length(fit$skew_occ), p_occ)
    expect_true(all(is.finite(fit$skew_occ)))

    expect_true(is.numeric(fit$skew_pos))
    expect_equal(length(fit$skew_pos), p_pos)
    expect_true(all(is.finite(fit$skew_pos)))

    # Draws matrices: 1000 rows x p_arm columns.
    expect_true(is.matrix(fit$draws_occ))
    expect_equal(nrow(fit$draws_occ), 1000L)
    expect_equal(ncol(fit$draws_occ), p_occ)

    expect_true(is.matrix(fit$draws_pos))
    expect_equal(nrow(fit$draws_pos), 1000L)
    expect_equal(ncol(fit$draws_pos), p_pos)
})


# ---------------------------------------------------------------------------
# 3. Gaussian-mode unchanged: default approx leaves SLA fields off
# ---------------------------------------------------------------------------

test_that("approx='gaussian_laplace' leaves SLA fields off", {
  skip_if_fast()
    set.seed(103)
    sim <- .make_cover_data(seed = 103, N = 200, n_s = 25)

    fit <- suppressMessages(tobs(
        formula = ~ x + bym2(graph = sim$adj, group_var = "region"),
        data    = sim$data,
        family  = cover("beta"),
        y       = sim$y,
        method  = "nested_laplace",
        # Default approx (or explicit gaussian_laplace) — both should leave
        # the SLA fields un-populated.
        control = list(
            sigma.grid     = c(0.4, 0.8),
            rho.grid       = c(0.5, 0.9)
        )
    ))

    # Match the standalone-Laplace convention (decode_cover_hurdle): the
    # off-mode sets sla_status = "off" and leaves skew / draws as NULL.
    expect_identical(fit$sla_status %||% "off", "off")
    expect_null(fit$skew_occ)
    expect_null(fit$skew_pos)
    expect_null(fit$draws_occ)
    expect_null(fit$draws_pos)
})


# ---------------------------------------------------------------------------
# 4. Gaussian-limit sanity: gamma should be near zero at large N
# ---------------------------------------------------------------------------

test_that("SLA gamma near zero at large N", {
    skip_on_cran()
    skip_if_fast()

    set.seed(104)
    sim <- .make_cover_data(seed = 104, N = 1000, n_s = 25,
                            alpha_true = 1.0)

    fit <- suppressMessages(tobs(
        formula = ~ x + bym2(graph = sim$adj, group_var = "region"),
        data    = sim$data,
        family  = cover("beta"),
        y       = sim$y,
        method  = "nested_laplace_sla",
        control = list(
            sigma.grid     = c(0.4, 0.8),
            rho.grid       = c(0.5, 0.9)
        )
    ))

    # Inner Laplace is near-Gaussian at large N for Bernoulli/Beta with
    # moderately disperse data — gamma should be small in magnitude.
    expect_true(all(is.finite(fit$skew_occ)))
    expect_true(all(is.finite(fit$skew_pos)))
    expect_lt(max(abs(fit$skew_occ)), 0.2)
    expect_lt(max(abs(fit$skew_pos)), 0.2)
})


# ---------------------------------------------------------------------------
# 5. Cross-check vs separate-Laplace SLA at near-zero spatial amplitude
# ---------------------------------------------------------------------------

test_that("joint SLA matches separate SLA at vanishing sigma", {
    skip_on_cran()
    skip_if_fast()

    set.seed(105)
    # Use a near-zero alpha and small sigma so the donor field has
    # vanishing amplitude; the joint kernel should then essentially
    # collapse to two independent Laplace arms.
    sim <- .make_cover_data(seed = 105, N = 200, n_s = 25,
                            alpha_true = 0.01, sigma = 0.05)

    fit_joint <- suppressMessages(tobs(
        formula = ~ x + bym2(graph = sim$adj, group_var = "region"),
        data    = sim$data,
        family  = cover("beta"),
        y       = sim$y,
        method  = "nested_laplace_sla",
        control = list(
            sigma.grid     = c(0.01, 0.02, 0.03),
            rho.grid       = c(0.5)
        )
    ))

    # Separate single-Laplace SLA: no spatial spec, no engine override.
    fit_sep <- suppressMessages(tobs(
        formula = ~ x,
        data    = sim$data,
        family  = cover("beta"),
        y       = sim$y,
        method  = "laplace_sla"
    ))

    expect_true(is.numeric(fit_joint$skew_pos))
    expect_true(is.numeric(fit_sep$skew_pos))
    expect_equal(length(fit_joint$skew_pos), length(fit_sep$skew_pos))

    # The joint path still integrates a 3-point (sigma, sigma_pos) grid
    # at scale ~0.02, so per-grid skewness contributions don't fully
    # cancel even at near-zero amplitude. What this test defends is the
    # qualitative collapse: both arms' skewness stay close to zero (well
    # under 0.1 absolute) and the joint stays within the SLA grid-mixing
    # noise floor of the separate path. Empirically the residual is
    # ~0.03-0.07 per coefficient; 1e-1 is the honest tolerance.
    expect_lt(max(abs(fit_joint$skew_pos)), 0.1)
    expect_lt(max(abs(fit_sep$skew_pos)),   0.1)
    expect_equal(fit_joint$skew_pos, fit_sep$skew_pos, tolerance = 1e-1,
                 ignore_attr = TRUE)
})
