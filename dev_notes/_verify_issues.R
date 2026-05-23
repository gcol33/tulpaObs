# Verify issue #8 (visit_data composition) and issue #10 (bar-syntax forms)
# after the multi-slope / nested / slope-only work + tulpa ABI 22.
suppressMessages(devtools::load_all("C:/Users/Gilles Colling/Documents/dev/tulpaObs"))
cat("=== load_all OK ===\n\n")

## ---- Issue #10 part A: AST desugaring maps each bar to the right re() ----
show_desugar <- function(f) {
  out <- tulpaObs:::.tobs_desugar_bars(f)
  cat(sprintf("  %-16s -> %s\n", deparse(f[[length(f)]]),
              deparse(out[[length(out)]])))
}
cat("--- bar desugaring ---\n")
show_desugar(~ (1 | g))
show_desugar(~ (x | g))
show_desugar(~ (x || g))
show_desugar(~ (1 + x + z | g))
show_desugar(~ (0 + x | g))
show_desugar(~ (1 | g:h))
show_desugar(~ (1 | g/h))
show_desugar(~ (1 + x | g/h))
cat("\n")

## ---- Issue #10 part B: each form fits without error (Laplace) ----
bar_case <- function(tag, f) {
  set.seed(2)
  N <- 80; J <- 5
  d <- data.frame(
    g = factor(sample(letters[1:6], N, replace = TRUE)),
    h = factor(sample(LETTERS[1:3], N, replace = TRUE)),
    x = rnorm(N), z = rnorm(N)
  )
  y <- matrix(rbinom(N * J, 1, 0.4), N, J)
  out <- tryCatch({
    fit <- tobs(formula = f, data = d, y = y, detection = ~ 1,
                family = occu(), engine = "laplace")
    sprintf("%-16s OK", tag)
  }, error = function(e) sprintf("%-16s ERROR -- %s", tag, conditionMessage(e)))
  cat(" ", out, "\n")
}
cat("--- fits ---\n")
bar_case("(1|g)",        ~ (1 | g))
bar_case("(x|g)",        ~ (x | g))
bar_case("(x||g)",       ~ (x || g))
bar_case("(1+x+z|g)",    ~ (1 + x + z | g))
bar_case("(0+x|g)",      ~ (0 + x | g))
bar_case("(1|g:h)",      ~ (1 | g:h))
bar_case("(1|g/h)",      ~ (1 | g/h))
bar_case("(1+x|g/h)",    ~ (1 + x | g/h))
cat("\n")

## ---- Issue #8: tobs_data() output composes with tobs(visit_data=) ----
set.seed(1)
df <- data.frame(
  site_id = rep(paste0("s", 1:6), each = 3),
  year    = rep(1:3, times = 6),
  occur   = sample(0:1, 18, replace = TRUE),
  effort  = runif(18)
)
od <- tobs_data(df, y = "occur", site = "site_id", visit = "year",
                det.covs = "effort")
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
  cat("ISSUE8: OK -- coef:\n"); print(round(coef(fit), 3)); "OK"
}, error = function(e) { cat("ISSUE8: ERROR --", conditionMessage(e), "\n"); "ERR" })

cat("\n=== verify done ===\n")
