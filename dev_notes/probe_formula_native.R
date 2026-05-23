# Probe: exercise the formula-native API end-to-end (no spatial=/re= args).
# Verifies builders -> structured_terms -> .tobs_structures_from_model ->
# engine, for the core occupancy paths.
devtools::load_all(".", quiet = TRUE)

set.seed(1)
n <- 60L
adj <- matrix(0L, n, n)
for (i in seq_len(n - 1L)) { adj[i, i + 1L] <- 1L; adj[i + 1L, i] <- 1L }

dat <- data.frame(
  elev     = rnorm(n),
  lon      = runif(n), lat = runif(n),
  observer = factor(sample(letters[1:3], n, replace = TRUE))
)
psi <- plogis(0.3 + 0.9 * dat$elev)
z   <- rbinom(n, 1, psi)
J   <- 4L
Y   <- matrix(rbinom(n * J, 1, z * 0.5), n, J)

ok <- function(tag, expr) {
  out <- tryCatch({ force(expr); "OK" },
                  error = function(e) paste("FAIL:", conditionMessage(e)))
  cat(sprintf("[%s] %s\n", out, tag))
}

# 1) plain occupancy, laplace (no structured terms)
ok("plain occu (laplace)", {
  f <- tobs(~ elev, dat, occu(), detection = ~ 1, y = Y)
  stopifnot(is.finite(f$intercepts$psi))
})

# 2) icar() spatial field in the occupancy formula, NUTS
ok("icar(graph=adj) on psi (nuts)", {
  f <- tobs(~ elev + icar(graph = adj), dat, occu(), detection = ~ 1, y = Y,
            engine = "nuts", control = list(iter = 60, warmup = 30))
  stopifnot(inherits(f, "tobs_fit"), !is.null(f$spatial))
})

# 3) re() random intercept on detection, NUTS
ok("re(observer) on detection (nuts)", {
  f <- tobs(~ elev, dat, occu(), detection = ~ re(observer), y = Y,
            engine = "nuts", control = list(iter = 60, warmup = 30))
  stopifnot(inherits(f, "tobs_fit"), !is.null(f$re))
})

# 4) copy(): one icar field shared across psi and detection
ok("icar(id) on psi + copy() on det (nuts)", {
  f <- tobs(~ elev + icar(graph = adj, id = "u"), dat, occu(),
            detection = ~ copy("u"), y = Y,
            engine = "nuts", control = list(iter = 60, warmup = 30))
  stopifnot(isTRUE(f$spatial$shared[1]), isTRUE(f$spatial$shared[2]))
})

# 5) structured term on an unsupported process errors clearly
ok("icar on colonization rejected", {
  res <- tryCatch(
    tobs(~ elev, dat, dyn_occu(), detection = ~ 1, y = array(Y, c(n, J, 1L)),
         col_formula = ~ icar(graph = adj), ext_formula = ~ 1, engine = "nuts",
         control = list(iter = 20, warmup = 10)),
    error = function(e) conditionMessage(e))
  stopifnot(is.character(res), grepl("occupancy/state and detection", res))
})

cat("\nprobe complete\n")
