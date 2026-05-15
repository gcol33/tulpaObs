devtools::load_all(".", quiet = TRUE)
chain_adj <- function(n_s) {adj <- matrix(0L, n_s, n_s); for (s in seq_len(n_s)) for (j in setdiff(c(s - 1L, s + 1L), c(0L, n_s + 1L))) adj[s, j] <- 1L; adj}
adj <- chain_adj(25L)
set.seed(50101L)
N <- 300L; n_s <- 25L; sigma <- 0.6; rho <- 0.7
region <- sample.int(n_s, N, replace = TRUE)
pf <- rnorm(n_s); tf <- rnorm(n_s); w  <- sigma * (sqrt(rho)*pf + sqrt(1-rho)*tf)
x  <- rnorm(N); alpha_t <- 1.0; s_t <- 0.3
occur <- rbinom(N, 1L, plogis(-0.3 + 0.7*x + w[region]))
eta_pos <- -1.5 + 0.3*x + alpha_t*w[region]
log_y <- rnorm(N, eta_pos, s_t)
y <- numeric(N); is_pos <- occur == 1L; y[is_pos] <- exp(log_y[is_pos])
y <- pmin(y, 1-1e-6)
data <- data.frame(x = x, region = factor(region))
spec <- tulpa::spatial_bym2(adj, level = "group", group_var = "region")

enc <- tulpaObs:::encode_cover_hurdle(~x, data, y, spec, positive="lognormal")
data_obs <- data[enc$obs_keep,,drop=FALSE]
prior <- tulpa::prior_from_spec(spec, data_obs)
spi_full <- prior$spatial_idx
spi_pos <- spi_full[enc$idx_pos]
N_pos <- length(enc$pos_data$y)
arm_occ <- list(y=as.numeric(enc$occ_data$y), n_trials=enc$occ_data$n_trials, X=enc$occ_data$X, spatial_idx=as.integer(spi_full), re_idx=rep(0,length(spi_full)), n_re_groups=0L, sigma_re=1.0, family="binomial", phi=1.0)

prior_j <- prior; prior_j$spatial_idx <- NULL; prior_j$rho_bounds <- NULL
prior_j$sigma_grid <- c(0.3,0.6,0.9); prior_j$rho_grid <- c(0.5,0.7,0.9)

for (test_phi in c(1.0, 0.62, 0.3, 0.1)) {
  arm_pos <- list(y=as.numeric(enc$pos_data$y), n_trials=rep(1L, N_pos), X=enc$pos_data$X, spatial_idx=as.integer(spi_pos), re_idx=rep(0,N_pos), n_re_groups=0L, sigma_re=1.0, family="gaussian", phi=test_phi)
  fit <- tulpa::tulpa_nested_laplace_joint(responses=list(occ=arm_occ, pos=arm_pos), prior=prior_j, copy=list(arm="pos", alpha_grid=c(0,0.5,1,1.5)), max_iter=50L, tol=1e-6, n_threads=1L)
  df <- as.data.frame(cbind(fit$theta_grid, log_marginal=fit$log_marginal, n_iter=fit$n_iter %||% NA))
  cat("\nphi =", test_phi, "\nall grid cells:\n"); print(df)
}
