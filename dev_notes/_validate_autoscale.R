# Compare:
#   - (A) raw formula  ~ year_btw + year_wtn        (autoscale ON internally)
#   - (B) user-scaled  ~ scale(year_btw) + scale(year_wtn)
# These two MUST produce equivalent predictions and the autoscale (A) must
# converge to the same MAP as (B) (just expressed in natural-scale coefs).
devtools::load_all()
suppressPackageStartupMessages({library(dplyr)})

set.seed(2026)
n_groups <- 60L
visits_per_grp <- 6L
rs_groups <- data.frame(
  group_id = paste0("G_", seq_len(n_groups)),
  base_year = sample(1985:2005, n_groups, replace = TRUE)
)
dat <- rs_groups[rep(seq_len(n_groups), each = visits_per_grp), ] %>%
  dplyr::mutate(visit_no = rep(seq_len(visits_per_grp), n_groups),
                year = base_year + 3 * (visit_no - 1))
psi_true <- 0.4; p_true <- 0.6; mu_int <- -1.2; mu_slope <- 0.02
z_per_grp <- stats::rbinom(n_groups, 1, psi_true)
dat$occur_true <- z_per_grp[match(dat$group_id, rs_groups$group_id)]
dat$detect <- ifelse(dat$occur_true == 1,
                     stats::rbinom(nrow(dat), 1, p_true), 0)
mu_cover <- plogis(mu_int + mu_slope * (dat$year - 2000))
phi_beta <- 8
cov_pos <- stats::rbeta(nrow(dat), mu_cover * phi_beta,
                       (1 - mu_cover) * phi_beta)
dat$y <- ifelse(dat$detect == 1, cov_pos, 0)
dat <- within_between(dat, group = "group_id", vars = "year")

# (A) Raw formula — autoscale should kick in internally.
fit_a <- tobs(~ year_btw + year_wtn, data = dat,
              family = cover(positive = "beta"), y = dat$y,
              engine = "laplace")

# (B) User-side scaled formula, no autoscale needed (sd already ~1).
dat_b <- dat %>% mutate(year_btw_sc = as.numeric(scale(year_btw)),
                        year_wtn_sc = as.numeric(scale(year_wtn)))
fit_b <- tobs(~ year_btw_sc + year_wtn_sc, data = dat_b,
              family = cover(positive = "beta"), y = dat_b$y,
              engine = "laplace")

cat("--- (A) raw formula, autoscale internal ---\n")
print(fit_a$beta_occ)
print(fit_a$beta_pos)

cat("\n--- (B) user-scaled formula ---\n")
print(fit_b$beta_occ)
print(fit_b$beta_pos)

# Predictions on the same training data should match closely between A and B.
p_a <- predict(fit_a, newdata = dat, type = "expected")
p_b <- predict(fit_b, newdata = dat_b, type = "expected")
cat(sprintf("\nMax |pred_a - pred_b| = %.3e\n", max(abs(p_a - p_b))))

# Cross-check: (A)'s natural-scale slope * sd(year_btw) should equal (B)'s scaled slope.
sd_ybtw <- sd(dat$year_btw)
sd_ywtn <- sd(dat$year_wtn)
cat(sprintf("\nSlope cross-check:\n"))
cat(sprintf("  occ: a$year_btw * sd = %+.4f  vs  b$year_btw_sc = %+.4f\n",
            fit_a$beta_occ[["year_btw"]] * sd_ybtw,
            fit_b$beta_occ[["year_btw_sc"]]))
cat(sprintf("  occ: a$year_wtn * sd = %+.4f  vs  b$year_wtn_sc = %+.4f\n",
            fit_a$beta_occ[["year_wtn"]] * sd_ywtn,
            fit_b$beta_occ[["year_wtn_sc"]]))
cat(sprintf("  pos: a$year_btw * sd = %+.4f  vs  b$year_btw_sc = %+.4f\n",
            fit_a$beta_pos[["year_btw"]] * sd_ybtw,
            fit_b$beta_pos[["year_btw_sc"]]))

cat("\nTruth:\n")
cat(sprintf("  psi*p   ~ %.3f  -> occ intercept (at mean) ~ %+.3f\n",
            psi_true * p_true, qlogis(psi_true * p_true)))
cat(sprintf("  pos intercept at year=2000: %+.3f\n", mu_int))
cat(sprintf("  pos slope (raw year): %+.5f\n", mu_slope))
cat(sprintf("  pos slope (sd-scaled year, sd=%.2f): %+.5f\n",
            sd_ybtw, mu_slope * sd_ybtw))
