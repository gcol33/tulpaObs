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
cat(">> class(coef(fit)):", class(coef(fit)), "\n")
cat(">> str(coef(fit)):\n"); str(coef(fit))
cat(">> coef(fit):\n"); print(coef(fit))
