# Tiny NB+RE problem to see one optim step at a time.
suppressMessages({ devtools::load_all(".", quiet = TRUE) })

set.seed(7)
N <- 40; J <- 4; ngrp <- 8
b <- stats::rnorm(ngrp, sd = 0.5)
grp <- rep(seq_len(ngrp), length.out = N)
data <- data.frame(g = factor(grp))
eta_l <- 0.5 + b[grp]
Nlat  <- stats::rnbinom(N, size = 4, mu = exp(eta_l))
p <- 0.5
y <- matrix(NA_integer_, N, J)
for (i in seq_len(N)) y[i, ] <- stats::rbinom(J, Nlat[i], p)

cat("max(y) =", max(y), "  K_max default =", max(y) + 100, "\n")

# Just exercise .tobs_nmix_re_aghq directly, see how many optim steps.
model <- tulpaObs:::.tobs_build_abun(
  abund_formula = ~ 1 + (1 | g),
  det_formula   = ~ 1,
  data = data, y = y)
# Strip the formula into structured terms (mimics .tobs_structures_from_model)
re_list <- list()
for (t in model$structured_terms) {
  spec  <- t$spec
  procs <- t$processes
  if (inherits(spec, "tobs_re")) {
    spec$shared <- c(1L %in% procs, 2L %in% procs)
    re_list[[length(re_list) + 1L]] <- spec
  }
}
cat("re_list length =", length(re_list), "\n")
cat("shared =", re_list[[1]]$shared, "\n")

# Skip the dispatcher, call the engine path directly with verbose options.
arms <- tulpaObs:::.tobs_nmix_re_split_arms(re_list, model)
design <- arms$lambda
cat("design length =", length(design), "  arm = lambda\n")

t0 <- Sys.time()
ref <- tulpaObs:::.tobs_nmix_re_aghq(
  model, design,
  beta_lambda = c(0.5),
  beta_p      = c(0),
  Sigma_list  = list(matrix(0.25, 1, 1)),
  b           = numeric(ngrp),
  mixture     = "NB", r_init = 4,
  K_max       = NULL, n_quad = 1L, lkj_eta = 1.5,
  theta_prior_sd = 100, max_iter = 30L, verbose = FALSE)
el <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
cat("\n.tobs_nmix_re_aghq took", round(el, 2), "s\n")
str(ref, max.level = 1)
