## Minimal repro: tobs(family = cover(positive = "beta")) loses intercept
## sanity when a design column lives on a large additive scale (here,
## calendar year ~ 2000 from within_between()).
##
## Mechanism: the intercept column (1) and a column with mean ~ 2000 form
## a near-singular pair in the design Hessian. The MAP optimizer satisfies
## the data with a large negative intercept compensated by a small slope
## along the high-magnitude column. Coefficients are mathematically
## consistent with the prediction but numerically nonsensical and the
## prior on the intercept (which the user specifies on the natural scale)
## is no longer doing what the user thinks.
##
## Standard fix used by glmnet / lme4 / brms: center+scale numeric design
## columns internally, optimize in scaled space, transform betas + SEs
## back to the user-facing scale before returning. Priors specified by
## the user on the original scale must be transformed to the scaled space
## (Normal(0, sd) on beta_j becomes Normal(0, sd * sd_j) on beta_sc_j).
##
## This script shows
##   (a) raw within_between() output -> intercept blows up
##   (b) manual scale() of *_btw and *_wtn columns -> sensible recovery
## across two truths (no time effect, small time effect).

suppressPackageStartupMessages({
  library(tulpaObs)
  library(dplyr)
})

set.seed(2026)

n_groups        <- 60L
visits_per_grp  <- 6L

rs_groups <- data.frame(
  group_id  = paste0("G_", seq_len(n_groups)),
  base_year = sample(1985:2005, n_groups, replace = TRUE)
)

dat <- rs_groups[rep(seq_len(n_groups), each = visits_per_grp), ] %>%
  dplyr::mutate(
    visit_no = rep(seq_len(visits_per_grp), n_groups),
    year     = base_year + 3 * (visit_no - 1)
  )

## Truth: psi (per group, constant across visits), p constant, cover_pos
## has a small time slope so the within-group year_wtn covariate has
## something to estimate.
psi_true <- 0.4
p_true   <- 0.6
mu_int   <- -1.2
mu_slope <-  0.02  # cover_pos logit slope per year (small but identifiable)

z_per_grp  <- stats::rbinom(n_groups, 1, psi_true)
dat$occur_true <- z_per_grp[match(dat$group_id, rs_groups$group_id)]
dat$detect     <- ifelse(dat$occur_true == 1,
                         stats::rbinom(nrow(dat), 1, p_true), 0)

mu_cover <- plogis(mu_int + mu_slope * (dat$year - 2000))
phi_beta <- 8
cov_pos  <- stats::rbeta(nrow(dat), mu_cover * phi_beta,
                         (1 - mu_cover) * phi_beta)
dat$y    <- ifelse(dat$detect == 1, cov_pos, 0)

dat <- tulpaObs::within_between(dat, group = "group_id", vars = "year")

## (A) Raw within_between output -- year_btw lives at ~ 2000
fit_raw <- tulpaObs::tobs(
  formula = ~ year_btw + year_wtn,
  data    = dat,
  family  = tulpaObs::cover(positive = "beta"),
  y       = dat$y,
  engine  = "laplace"
)

cat("(A) raw within_between (year_btw ~ 2000):\n")
cat(sprintf("  occ_(Intercept)            %+.3f\n",  fit_raw$beta_occ[["(Intercept)"]]))
cat(sprintf("  occ_year_btw               %+.5f\n", fit_raw$beta_occ[["year_btw"]]))
cat(sprintf("  occ_year_wtn               %+.5f\n", fit_raw$beta_occ[["year_wtn"]]))
cat(sprintf("  pos_(Intercept)            %+.3f\n",  fit_raw$beta_pos[["(Intercept)"]]))
cat(sprintf("  pos_year_btw               %+.5f\n", fit_raw$beta_pos[["year_btw"]]))
cat(sprintf("  pos_year_wtn               %+.5f\n", fit_raw$beta_pos[["year_wtn"]]))

## (B) User-side workaround: scale both decomposed columns
dat_sc <- dat %>% dplyr::mutate(
  year_btw_sc = as.numeric(scale(year_btw)),
  year_wtn_sc = as.numeric(scale(year_wtn))
)

fit_sc <- tulpaObs::tobs(
  formula = ~ year_btw_sc + year_wtn_sc,
  data    = dat_sc,
  family  = tulpaObs::cover(positive = "beta"),
  y       = dat_sc$y,
  engine  = "laplace"
)

cat("\n(B) user-scaled inputs:\n")
cat(sprintf("  occ_(Intercept)            %+.3f\n",  fit_sc$beta_occ[["(Intercept)"]]))
cat(sprintf("  occ_year_btw_sc            %+.3f\n",  fit_sc$beta_occ[["year_btw_sc"]]))
cat(sprintf("  occ_year_wtn_sc            %+.3f\n",  fit_sc$beta_occ[["year_wtn_sc"]]))
cat(sprintf("  pos_(Intercept)            %+.3f\n",  fit_sc$beta_pos[["(Intercept)"]]))
cat(sprintf("  pos_year_btw_sc            %+.3f\n",  fit_sc$beta_pos[["year_btw_sc"]]))
cat(sprintf("  pos_year_wtn_sc            %+.3f\n",  fit_sc$beta_pos[["year_wtn_sc"]]))

cat("\ntruth:\n")
cat(sprintf("  occ_(Intercept)  ~ logit(psi*p) = %+.3f\n",
            qlogis(psi_true * p_true)))
cat(sprintf("  pos_(Intercept)  = %+.3f at year=2000\n", mu_int))
cat(sprintf("  pos slope on raw year = %+.5f /year\n", mu_slope))

cat("\nsessionInfo() tulpaObs:\n")
print(packageVersion("tulpaObs"))
