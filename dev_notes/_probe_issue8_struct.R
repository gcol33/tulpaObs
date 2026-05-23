suppressMessages(devtools::load_all("C:/Users/Gilles Colling/Documents/dev/tulpaObs"))
set.seed(1)
df <- data.frame(
  site_id = rep(paste0("s", 1:6), each = 3),
  year    = rep(1:3, times = 6),
  occur   = sample(0:1, 18, replace = TRUE),
  effort  = runif(18)
)
od <- tobs_data(df, y = "occur", site = "site_id", visit = "year",
                det.covs = "effort")
fit <- tobs(formula = ~ 1,
            data = data.frame(site_id = unique(df$site_id)),
            y = od$y, detection = ~ effort,
            visit_data = od$det.covs, family = occu(), engine = "laplace")

cat("=== names(fit) ===\n"); print(names(fit))
cat("\n=== fit$process_info ===\n"); print(fit$process_info)
cat("\n=== fit$coefficients ===\n"); print(fit$coefficients)
cat("\n=== fit$model$process_info (names/p) ===\n")
mi <- fit$model$process_info
if (!is.null(mi)) for (pi in mi) cat(sprintf("  %s: p=%s, cols=%s\n",
  pi$name, pi$p, paste(pi$colnames, collapse=",")))
cat("\n=== summary ===\n"); print(summary(fit))
