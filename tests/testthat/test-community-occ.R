test_that("community occupancy model runs", {
  set.seed(42)
  n_sites <- 20
  n_species <- 3
  max_visits <- 3
  sp_effects <- rnorm(n_species, 0, 0.3)

  y_list <- list()
  for (s in seq_len(n_species)) {
    z <- rbinom(n_sites, 1, plogis(sp_effects[s]))
    y_s <- matrix(0L, n_sites, max_visits)
    for (i in seq_len(n_sites)) if (z[i] == 1) y_s[i, ] <- rbinom(max_visits, 1, 0.4)
    y_list[[paste0("sp", s)]] <- y_s
  }

  fit <- tobs(
    formula   = ~ 1,
    data      = data.frame(x = rnorm(n_sites)),
    family    = ms_occu(),
    detection = ~ 1,
    y         = y_list,
    species   = TRUE,
    engine    = "laplace",
    control   = list(verbose = FALSE)
  )
  expect_s3_class(fit, "tobs_fit")
  expect_true(fit$n_params >= 2)
})
