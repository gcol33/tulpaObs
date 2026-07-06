# Generate the bundled example datasets for tulpaObs.
#
# These are SYNTHETIC datasets drawn from documented generative models with
# fixed seeds. They exist so `?occu` / `?abun` / `?cover` examples and the
# quickstart run without a `simulate_*()` call, and so every dataset carries a
# known `truth` for the parameter-recovery story the package tells. They are
# not field data.
#
# Run from the package root:
#   Rscript data-raw/make_datasets.R
# Regenerates data/peatland_occu.rda, data/foray_counts.rda, data/meadow_cover.rda.

set.seed(20260706)

`%||%` <- function(a, b) if (is.null(a)) b else a

# ---------------------------------------------------------------------------
# 1. peatland_occu -- single-season occupancy (MacKenzie et al. 2002)
#
#   z_i    ~ Bernoulli(psi_i),  logit psi_i = 0.20 + 0.90 wetness - 0.55 elevation
#   y_ij|z ~ Bernoulli(z_i p_ij), logit p_ij = 0.10 + 0.65 effort_ij
#
# 120 sites, 4 visits. A wetland amphibian survey: occupancy tracks site
# wetness (positively) and elevation (negatively); detection tracks per-visit
# survey effort.
# ---------------------------------------------------------------------------
make_peatland_occu <- function() {
  N <- 120L; J <- 4L
  elevation <- as.numeric(scale(runif(N, 200, 1400)))
  wetness   <- as.numeric(scale(rbeta(N, 2, 2)))
  b_psi <- c(`(Intercept)` = 0.20, wetness = 0.90, elevation = -0.55)
  psi   <- plogis(b_psi[1] + b_psi["wetness"] * wetness +
                    b_psi["elevation"] * elevation)
  z     <- rbinom(N, 1L, psi)

  effort <- matrix(as.numeric(scale(rnorm(N * J))), N, J)
  b_p <- c(`(Intercept)` = 0.10, effort = 0.65)
  p   <- plogis(b_p[1] + b_p["effort"] * effort)
  y   <- matrix(rbinom(N * J, 1L, as.vector(z) * p), N, J)
  # a scatter of NA visits, as real detection histories carry
  y[cbind(sample(N, 12), sample(J, 12, replace = TRUE))] <- NA

  coords <- cbind(x = runif(N), y = runif(N))
  list(
    y        = y,
    occ.covs = data.frame(elevation = elevation, wetness = wetness),
    det.covs = list(effort = effort),
    coords   = coords,
    truth    = list(beta_psi = b_psi, beta_p = b_p, z = z, psi = psi)
  )
}

# ---------------------------------------------------------------------------
# 2. foray_counts -- binomial N-mixture abundance (Royle 2004)
#
#   N_i    ~ Poisson(lambda_i), log lambda_i = 0.80 + 0.70 shrub - 0.30 elevation
#   y_ij|N ~ Binomial(N_i, p_ij), logit p_ij = 0.00 + 0.55 effort_ij
#
# 100 sites, 3 visits. Repeated point counts of a shrub-associated bird.
# ---------------------------------------------------------------------------
make_foray_counts <- function() {
  N <- 100L; J <- 3L
  elevation <- as.numeric(scale(runif(N, 100, 900)))
  shrub     <- as.numeric(scale(rgamma(N, 2, 1)))
  b_lambda <- c(`(Intercept)` = 0.80, shrub = 0.70, elevation = -0.30)
  lambda   <- exp(b_lambda[1] + b_lambda["shrub"] * shrub +
                    b_lambda["elevation"] * elevation)
  Ni       <- rpois(N, lambda)

  effort <- matrix(as.numeric(scale(rnorm(N * J))), N, J)
  b_p <- c(`(Intercept)` = 0.00, effort = 0.55)
  p   <- plogis(b_p[1] + b_p["effort"] * effort)
  y   <- matrix(rbinom(N * J, as.vector(Ni), p), N, J)

  list(
    y        = y,
    occ.covs = data.frame(elevation = elevation, shrub = shrub),
    det.covs = list(effort = effort),
    truth    = list(beta_lambda = b_lambda, beta_p = b_p, N = Ni, lambda = lambda)
  )
}

# ---------------------------------------------------------------------------
# 3. meadow_cover -- vegetation cover hurdle (occurrence + conditional cover)
#
#   occur_i ~ Bernoulli(pi_i),  logit pi_i = -0.30 + 1.05 moisture
#   cover_i | occur = 1 ~ Beta(mu_i phi, (1 - mu_i) phi),
#                        logit mu_i = -0.55 + 0.80 moisture - 0.45 grazing,  phi = 6
#   cover_i | occur = 0 = 0
#
# 150 grassland plots across 3 survey years (a resurvey panel; `year` is
# centred so `within_between()` can split cross-plot baseline from within-plot
# trend). Returned as a long data frame -- one row per plot -- ready for
# `cover(response = "beta")`.
# ---------------------------------------------------------------------------
make_meadow_cover <- function() {
  P <- 150L
  plot     <- factor(sprintf("P%03d", seq_len(P)))
  year     <- sample(2018:2020, P, replace = TRUE)
  moisture <- as.numeric(scale(rbeta(P, 2, 2)))
  grazing  <- as.numeric(scale(runif(P)))

  b_occ <- c(`(Intercept)` = -0.30, moisture = 1.05)
  pi    <- plogis(b_occ[1] + b_occ["moisture"] * moisture)
  occur <- rbinom(P, 1L, pi)

  b_pos <- c(`(Intercept)` = -0.55, moisture = 0.80, grazing = -0.45)
  phi   <- 6
  mu    <- plogis(b_pos[1] + b_pos["moisture"] * moisture +
                    b_pos["grazing"] * grazing)
  cover <- ifelse(occur == 1L,
                  rbeta(P, mu * phi, (1 - mu) * phi),
                  0)

  data.frame(
    plot     = plot,
    year     = year,
    year_c   = year - mean(year),
    moisture = moisture,
    grazing  = grazing,
    cover    = cover,
    row.names = NULL,
    stringsAsFactors = FALSE
  )
}

peatland_occu <- make_peatland_occu()
foray_counts  <- make_foray_counts()
meadow_cover  <- make_meadow_cover()

# attach the generative truth to the cover data frame as an attribute so the
# data-frame return stays tidy but the recovery target is still reachable.
attr(meadow_cover, "truth") <- list(
  beta_occ = c(`(Intercept)` = -0.30, moisture = 1.05),
  beta_pos = c(`(Intercept)` = -0.55, moisture = 0.80, grazing = -0.45),
  phi = 6
)

dir.create("data", showWarnings = FALSE)
save(peatland_occu, file = "data/peatland_occu.rda", compress = "xz")
save(foray_counts,  file = "data/foray_counts.rda",  compress = "xz")
save(meadow_cover,  file = "data/meadow_cover.rda",  compress = "xz")

cat(sprintf("peatland_occu: %d sites x %d visits, naive occ %.2f\n",
            nrow(peatland_occu$y), ncol(peatland_occu$y),
            mean(rowSums(peatland_occu$y, na.rm = TRUE) > 0)))
cat(sprintf("foray_counts:  %d sites x %d visits, max count %d\n",
            nrow(foray_counts$y), ncol(foray_counts$y), max(foray_counts$y)))
cat(sprintf("meadow_cover:  %d plots, %.0f%% with cover > 0\n",
            nrow(meadow_cover), 100 * mean(meadow_cover$cover > 0)))
