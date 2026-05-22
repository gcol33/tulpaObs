suppressPackageStartupMessages({
  devtools::load_all(quiet = TRUE)
})

set.seed(42)
n_sites <- 60; max_visits <- 4
elev <- rnorm(n_sites)
psi <- plogis(-0.3 + 0.8 * elev)
z <- rbinom(n_sites, 1, psi)
effort_mat <- matrix(runif(n_sites * max_visits, 0.2, 1.5),
                     n_sites, max_visits)
p_mat <- plogis(qlogis(0.4) + 1.5 * (effort_mat - mean(effort_mat)))
y <- matrix(0L, n_sites, max_visits)
for (i in seq_len(n_sites)) {
  if (z[i] == 1) {
    for (j in seq_len(max_visits)) y[i, j] <- rbinom(1, 1, p_mat[i, j])
  }
}

long <- data.frame(
  site_id = rep(paste0("s", seq_len(n_sites)), each = max_visits),
  year    = rep(seq_len(max_visits), times = n_sites),
  occur   = as.vector(t(y)),
  effort  = as.vector(t(effort_mat))
)
od <- tobs_data(long, y = "occur", site = "site_id", visit = "year",
                det.covs = "effort")
site_df <- data.frame(site_id = paste0("s", seq_len(n_sites)))

fit <- tobs(
  formula    = ~ 1,
  data       = site_df,
  y          = od$y,
  detection  = ~ effort,
  visit_data = od$det.covs,
  family     = occu(),
  engine     = "laplace",
  control    = list(verbose = FALSE)
)

cat("class(fit):\n"); print(class(fit))
cat("names(fit):\n"); print(names(fit))
cat("draws cols:\n"); print(colnames(fit$draws))
if (!is.null(fit$means)) {
  cat("means:\n"); print(fit$means)
}
cat("model$X_det_visit dim:\n"); print(dim(fit$model$X_det_visit))
cat("model$det_visit_names:\n"); print(fit$model$det_visit_names)
cat("model$process_info:\n"); print(fit$model$process_info)
