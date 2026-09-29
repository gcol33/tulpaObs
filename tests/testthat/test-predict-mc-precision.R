# Monte Carlo precision of draw-based trajectory summaries (#386): the
# order-statistic MC error estimate, and the nsim = "auto" stopping loop.

test_that("the MC error of a quantile matches its sampling SD across replicates", {
    set.seed(11)
    n <- 2000L; R <- 1000L
    # Each row is one replicate of n draws with unit posterior SD, so the
    # spread of a row quantile across rows is its MC standard error in SDs --
    # the quantity the estimate claims to measure.
    gen <- list(gaussian = function(k) stats::rnorm(k),
                lognormal = function(k) {
                    x <- stats::rlnorm(k); x / sqrt((exp(1) - 1) * exp(1)) })
    for (g in names(gen)) {
        m <- matrix(gen[[g]](R * n), R, n)
        for (p in c(0.025, 0.5, 0.975)) {
            err  <- tulpaObs:::.tobs_mc_quantile_error(m, p)
            q    <- apply(m, 1L, stats::quantile, probs = p, names = FALSE)
            spread <- stats::sd(q) / mean(apply(m, 1L, stats::sd))
            expect_equal(mean(err), spread, tolerance = 0.12,
                         label = sprintf("%s p = %s", g, p))
        }
    }
    # Gaussian 0.975 quantile against its asymptotic SE.
    m <- matrix(stats::rnorm(R * n), R, n)
    se_asym <- sqrt(0.975 * 0.025 / n) / stats::dnorm(stats::qnorm(0.975))
    expect_equal(mean(tulpaObs:::.tobs_mc_quantile_error(m, 0.975)), se_asym,
                 tolerance = 0.12)
    # A bracket clipped at 0 keeps a finite, positive error.
    e <- tulpaObs:::.tobs_mc_quantile_error(matrix(stats::rnorm(3 * 60), 3L), 0.005)
    expect_true(all(is.finite(e) & e > 0))
    # A row with no spread carries no MC error.
    expect_identical(
        tulpaObs:::.tobs_mc_quantile_error(matrix(1, 2L, 50L), c(0.5, 0.975)),
        matrix(0, 2L, 2L))
})

.mc_bundle <- function(n) list(n = n, m = matrix(stats::rnorm(3L * n), 3L, n))
.mc_state  <- function(b) function(nd) list(x = b$m)
.mc_run <- function(nsim, redraw = .mc_bundle, mc.tol = 0.05,
                    nsim.max = 10000L, first = 250L) {
    tulpaObs:::.tobs_mc_trajectory(
        .mc_state, .mc_bundle(first), redraw, nds = list(NULL, NULL),
        times = c(0, 1), cell = 1:3, level = 0.95, aggregate = FALSE,
        draws = FALSE, nsim = nsim, mc.tol = mc.tol, nsim.max = nsim.max)
}

test_that("nsim = \"auto\" draws until every bound meets mc.tol", {
    set.seed(12)
    tr <- .mc_run("auto")
    expect_lte(attr(tr, "mc_se_max"), 0.05)
    expect_false(attr(tr, "nsim_capped"))
    expect_identical(attr(tr, "mc_tol"), 0.05)
    # About 7.2 / 0.05^2 = 2,900 draws for a Gaussian quantity at level 0.95.
    expect_gt(attr(tr, "nsim_used"), 2000L)
    expect_lte(attr(tr, "nsim_used"), 6500L)
    expect_null(attr(tr, ".mc_se"))
})

test_that("nsim.max stops the draws and says so", {
    set.seed(13)
    tr <- .mc_run("auto", nsim.max = 600L)
    expect_identical(attr(tr, "nsim_used"), 600L)
    expect_true(attr(tr, "nsim_capped"))
    expect_gt(attr(tr, "mc_se_max"), 0.05)
})

test_that("an integer nsim builds the table once and reports its precision", {
    set.seed(14)
    tr <- .mc_run(300L, redraw = function(n) stop("redrawn"), first = 300L)
    expect_identical(attr(tr, "nsim_used"), 300L)
    expect_true(is.na(attr(tr, "mc_tol")))
    expect_false(attr(tr, "nsim_capped"))
    expect_gt(attr(tr, "mc_se_max"), 0.05)
})

test_that("nsim, mc.tol and nsim.max are validated", {
    chk <- tulpaObs:::.tobs_check_nsim
    expect_identical(chk(500, "change", 0.05, 1e4), 500L)
    expect_identical(chk("auto", "trajectory", 0.05, 1e4), "auto")
    expect_error(chk("auto", "change", 0.05, 1e4), "trajectory\" only")
    expect_error(chk("auto", "trajectory", 0, 1e4), "mc.tol")
    expect_error(chk("auto", "trajectory", 0.05, 1.5), "nsim.max")
    expect_error(chk(10.5, "trajectory", 0.05, 1e4), "whole number")
    expect_error(chk("many", "trajectory", 0.05, 1e4), "whole number")
})
