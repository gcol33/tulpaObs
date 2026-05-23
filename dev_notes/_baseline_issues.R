# Baseline check before fixing issues #8 and #10.
# Run: load_all (compiles current src), then exercise both issue paths.
suppressMessages(devtools::load_all("."))
cat("=== load_all OK ===\n")

## ---- Issue #8: tobs_data() output composes with tobs(visit_data=) ----
set.seed(1)
df <- data.frame(
  site_id = rep(paste0("s", 1:5), each = 3),
  year    = rep(1:3, times = 5),
  occur   = sample(0:1, 15, replace = TRUE),
  effort  = runif(15)
)
od <- tobs_data(df, y = "occur", site = "site_id", visit = "year",
                det.covs = "effort")
cat("--- str(od$det.covs) ---\n"); str(od$det.covs)
res8 <- tryCatch({
  fit <- tobs(
    formula    = ~ 1,
    data       = data.frame(site_id = unique(df$site_id)),
    y          = od$y,
    detection  = ~ effort,
    visit_data = od$det.covs,
    family     = occu(),
    engine     = "laplace"
  )
  cat("ISSUE8: OK — fit class:", paste(class(fit), collapse="/"), "\n")
  print(coef(fit))
  "OK"
}, error = function(e) { cat("ISSUE8: ERROR —", conditionMessage(e), "\n"); "ERR" })

## ---- Issue #10: bar-syntax forms ----
bar_case <- function(tag, f) {
  set.seed(2)
  N <- 60; J <- 4
  d <- data.frame(
    g = factor(sample(letters[1:6], N, replace = TRUE)),
    x = rnorm(N), z = rnorm(N)
  )
  y <- matrix(rbinom(N*J, 1, 0.4), N, J)
  out <- tryCatch({
    fit <- tobs(formula = f, data = d, y = y, detection = ~ 1,
                family = occu(), engine = "laplace")
    paste0(tag, ": OK")
  }, error = function(e) paste0(tag, ": ERROR — ", conditionMessage(e)))
  cat(out, "\n")
}
bar_case("(1|g)        ", ~ (1 | g))
bar_case("(x|g)        ", ~ (x | g))
bar_case("(x||g)       ", ~ (x || g))
bar_case("(0+x|g)slope ", ~ (0 + x | g))
bar_case("(1+x+z|g)mult", ~ (1 + x + z | g))
bar_case("(1|g/h)nested", ~ (1 | g/x))   # x as 2nd grouping just to trip parser
cat("=== baseline done ===\n")
