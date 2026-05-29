devtools::load_all(".", quiet = TRUE)
set.seed(20260529)

si <- simulate_int_occu(N_total = 400, n_data = 2, J = c(5, 3), n_shared = 120,
                        beta_occ = c(0.4, 1.0),
                        beta_det = list(c(0.6, 0), c(-0.2, 0)), seed = 42)
cat("--- y structure ---\n")
str(si$y, max.level = 1)
cat("dim s1:", dim(si$y[[1]]), " dim s2:", dim(si$y[[2]]), "\n")
cat("data cols:", names(si$data), "\n")
cat("nrow data:", nrow(si$data), "\n")

fit <- tobs(~ x, data = si$data, family = int_occu(), detection = ~ 1,
            y = list(structured = si$y[[1]], casual = si$y[[2]]),
            method = "laplace", control = list(verbose = FALSE))
cat("\n--- coef ---\n"); print(coef(fit))
cat("\n--- summary ---\n"); print(summary(fit))
cat("\n--- confint ---\n"); print(confint(fit))
cat("\n--- truth beta_occ ---\n"); print(si$truth$beta_occ)

# occu on source 1 alone
src1 <- si$y[[1]]
rows1 <- si$site_maps[[1]]
d1 <- si$data[rows1, , drop = FALSE]
fit1 <- tobs(~ x, data = d1, family = occu(), detection = ~ 1,
             y = src1, method = "laplace", control = list(verbose = FALSE))
cat("\n--- source1-only coef ---\n"); print(coef(fit1))
cat("\n--- source1-only confint ---\n"); print(confint(fit1))

# source-specific detection covariate: add per-site covariate to data
si$data$wind <- rnorm(nrow(si$data))
fit2 <- tobs(~ x, data = si$data, family = int_occu(),
             detection = list(structured = ~ wind, casual = ~ 1),
             y = list(structured = si$y[[1]], casual = si$y[[2]]),
             method = "laplace", control = list(verbose = FALSE))
cat("\n--- per-source detection formula coef ---\n"); print(coef(fit2))

# NUTS smoke
fitn <- tobs(~ x, data = si$data, family = int_occu(), detection = ~ 1,
             y = list(structured = si$y[[1]], casual = si$y[[2]]),
             method = "nuts",
             control = list(n.iter = 300, n.warmup = 150, seed = 1, verbose = FALSE))
cat("\n--- NUTS coef ---\n"); print(coef(fitn))
