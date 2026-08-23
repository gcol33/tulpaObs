# Gates fire from error paths, so a malformed message in one is invisible to
# every fitting test. Two checks here: the message text a user actually sees on
# the repaired gates, and a package-wide scan for the sprintf() arity defect
# that produced them.

test_that("sprintf() calls in R/ have a matching format and argument count", {
  skip_on_cran()
  r_dir <- normalizePath(file.path(testthat::test_path(), "..", "..", "R"),
                         mustWork = FALSE)
  skip_if_not(dir.exists(r_dir), "package sources not available")

  # A multi-line message written as comma-separated literals inside sprintf()
  # without a paste0() wrapper leaves only the first literal as the format; the
  # rest silently become `...`. Where the format carries no directive the extra
  # args are dropped (truncated message); where a string literal lands on %d,
  # sprintf() itself errors.
  count_directives <- function(fmt) {
    stripped <- gsub("%%", "", fmt, fixed = TRUE)
    m <- gregexpr("%[-+ #0]*[0-9*]*(\\.[0-9*]+)?[disefgGxXoObaA]", stripped)[[1]]
    if (identical(as.integer(m), -1L)) 0L else length(m)
  }
  bad_directives <- function(fmt) {
    stripped <- gsub("%%", "", fmt, fixed = TRUE)
    m <- regmatches(stripped, gregexpr("%[-+ #0]*[0-9*]*(\\.[0-9*]+)?[a-zA-Z]",
                                       stripped))[[1]]
    m[!grepl("[disefgGxXoObaA]$", m)]
  }

  bad <- character()
  scan_call <- function(e, file) {
    if (!is.call(e)) return(invisible(NULL))
    fn <- e[[1]]
    if (is.name(fn) && identical(as.character(fn), "sprintf") && length(e) >= 2) {
      fmt <- e[[2]]
      if (is.character(fmt) && length(fmt) == 1) {
        args <- as.list(e)[-c(1, 2)]
        nms <- names(args)
        if (!is.null(nms)) args <- args[!nzchar(nms)]
        invalid <- bad_directives(fmt)
        if (count_directives(fmt) != length(args) || length(invalid)) {
          bad <<- c(bad, sprintf("%s: %s", basename(file), substr(fmt, 1, 60)))
        }
      }
    }
    for (i in seq_along(e)) {
      # A default-less formal parses to the empty symbol; forcing it errors.
      is_call_slot <- tryCatch(is.call(e[[i]]), error = function(...) FALSE)
      if (is_call_slot) scan_call(e[[i]], file)
    }
    invisible(NULL)
  }

  for (f in list.files(r_dir, pattern = "\\.R$", full.names = TRUE)) {
    for (e in parse(f, keep.source = FALSE)) scan_call(e, f)
  }
  expect_equal(bad, character())
})

test_that("temporal() errors on the Laplace engine instead of being dropped", {
  skip_on_cran()
  s <- simulate_occu(N = 30, J = 3, seed = 1)
  s$data$year <- rep(1:4, length.out = nrow(s$data))
  # .tobs_laplace() has no temporal channel; fitting anyway would omit the field
  # the user asked for. The term is carried by nested_laplace / nuts.
  expect_error(
    tobs(~ occ_cov1 + temporal(year), detection = ~ 1, data = s$data,
         family = occu(), y = s$y, method = "laplace", verbose = FALSE),
    "temporal\\(\\) is not consumed by the Laplace engine")
})

test_that("abun() spatial gates interpolate the term they rejected", {
  skip_on_cran()
  a <- simulate_abun(N = 20, J = 3, seed = 2)
  adj <- matrix(0, 20, 20)
  adj[cbind(1:19, 2:20)] <- 1; adj[cbind(2:20, 1:19)] <- 1

  # The rejected type must reach the message: it is what tells the user which
  # term to swap out.
  expect_error(
    tobs(~ abund_cov1 + car(graph = adj), detection = ~ 1, data = a$data,
         family = abun(), y = a$y, method = "nested_laplace", verbose = FALSE),
    "got 'car'")

  a$data$lon <- runif(nrow(a$data)); a$data$lat <- runif(nrow(a$data))
  expect_error(
    tobs(~ abund_cov1 + gp(lon, lat, prior_range = c(0.3, 0.5)),
         detection = ~ 1, data = a$data, family = abun(), y = a$y,
         method = "nested_laplace", verbose = FALSE),
    "the dense GP term 'gp' is not wired")

  # A node/site count mismatch reports both counts (this one used to die inside
  # sprintf() with a '%d' type error before reaching the user).
  adj5 <- matrix(0, 5, 5)
  adj5[cbind(1:4, 2:5)] <- 1; adj5[cbind(2:5, 1:4)] <- 1
  expect_error(
    tobs(~ abund_cov1 + icar(graph = adj5), detection = ~ 1, data = a$data,
         family = abun(), y = a$y, method = "nested_laplace", verbose = FALSE),
    "has 5 units but the model has 20 sites")
})

test_that("distance() cutpoints gate reports the expected length", {
  skip_on_cran()
  expect_error(
    tobs(~ 1, detection = ~ 1, data = data.frame(z = rep(0, 10)),
         family = distance(key = "halfnorm", transect = "line",
                           cutpoints = c(0, 50)),
         y = matrix(0L, 10, 3), method = "laplace", verbose = FALSE),
    "cutpoints must have length ncol\\(y\\) \\+ 1 = 4")
})

test_that("a family that does not consume a structured term rejects it", {
  skip_on_cran()
  # The binder strips every registry term and lme4 bar out of the fixed-effect
  # formula, so a dispatcher that never reads `structured_terms` returns a
  # fixed-effects fit that looks like the structured one the caller asked for.
  # Each case below is a term the family cannot honour; every one must error.
  chain <- function(n) {
    a <- matrix(0, n, n)
    a[cbind(1:(n - 1), 2:n)] <- 1; a[cbind(2:n, 1:(n - 1))] <- 1
    a
  }
  cp <- c(0, 25, 50, 75, 100)

  rn <- simulate_royle_nichols(N = 20, J = 3, seed = 1)
  rn$data$g <- factor(rep(1:4, length.out = nrow(rn$data)))
  expect_error(
    tobs(~ x + icar(graph = chain(nrow(rn$data))), detection = ~ 1,
         data = rn$data, family = royle_nichols(), y = rn$y,
         method = "laplace", verbose = FALSE),
    "royle_nichols[(][)]: the spatial field icar[(][)] is not wired")
  expect_error(
    tobs(~ x + (1 | g), detection = ~ 1, data = rn$data,
         family = royle_nichols(), y = rn$y, method = "laplace",
         verbose = FALSE),
    "royle_nichols[(][)]: a random effect .* is not wired")

  do <- simulate_double_observer(N = 20, seed = 1)
  do$data$g  <- factor(rep(1:4, length.out = nrow(do$data)))
  do$data$yr <- rep(1:4, length.out = nrow(do$data))
  expect_error(
    tobs(~ abund_cov1 + temporal(yr), detection = ~ 1, data = do$data,
         family = double_observer(), y = do$y, method = "laplace",
         verbose = FALSE),
    "double_observer[(][)]: temporal[(][)] is not wired")
  # A term written in `detection =` is arm-tagged, not a separate channel: it
  # has to be rejected on the detection arm too.
  expect_error(
    tobs(~ abund_cov1, detection = ~ (1 | g), data = do$data,
         family = double_observer(), y = do$y, method = "laplace",
         verbose = FALSE),
    "double_observer[(][)]: a random effect .* is not wired")

  gd <- simulate_gdistremoval(N = 20, seed = 1)
  expect_error(
    tobs(~ abund_cov1 + icar(graph = chain(nrow(gd$data))), detection = ~ 1,
         removal = ~ 1, data = gd$data, family = gdistremoval(cutpoints = cp),
         y = gd$y, y_rem = gd$y_rem, method = "laplace", verbose = FALSE),
    "gdistremoval[(][)]: the spatial field icar[(][)] is not wired")

  ds <- simulate_distsamp_open(N = 20, seed = 1)
  ds$data$g <- factor(rep(1:4, length.out = nrow(ds$data)))
  expect_error(
    tobs(~ abund_cov1 + (1 | g), detection = ~ 1, data = ds$data,
         family = distsamp_open(cutpoints = cp), y = ds$y, method = "laplace",
         verbose = FALSE),
    "distsamp_open[(][)]: a random effect .* is not wired")
})

test_that("ms_int_occu() parses its formulas instead of passing them to model.matrix()", {
  skip_on_cran()
  s  <- simulate_ms_int_occu(seed = 1)
  n  <- nrow(s$data)
  sp <- paste0("sp", seq_len(dim(s$y[[1L]])[3L]))

  # `1 | site` with a numeric `site` is valid base-R that evaluates to a
  # constant TRUE column, so an unparsed formula fits a silently DIFFERENT
  # model rather than dropping the term.
  s$data$site <- seq_len(n)
  expect_error(
    tobs(~ x + (1 | site), detection = ~ 1, data = s$data,
         family = ms_int_occu(), y = s$y, species = sp, method = "laplace",
         verbose = FALSE),
    "ms_int_occu[(][)]: a random effect .* is not wired")

  # A registry term used to die as "could not find function" instead of the
  # family's own pointer.
  adj <- matrix(0, n, n)
  adj[cbind(1:(n - 1), 2:n)] <- 1; adj[cbind(2:n, 1:(n - 1))] <- 1
  expect_error(
    tobs(~ x + icar(graph = adj), detection = ~ 1, data = s$data,
         family = ms_int_occu(), y = s$y, species = sp, method = "laplace",
         verbose = FALSE),
    "ms_int_occu[(][)]: the spatial field icar[(][)] is not wired")

  # The fixed-effects fit is unaffected by the binder: same designs, same fit.
  f <- tobs(~ x, detection = ~ 1, data = s$data, family = ms_int_occu(),
            y = s$y, species = sp, method = "laplace",
            control = list(verbose = FALSE, progress = FALSE))
  expect_s3_class(f, "tobs_fit")
  expect_identical(colnames(f$model$X_psi), c("(Intercept)", "x"))
  expect_true(all(vapply(f$model$X_p, function(X)
    identical(colnames(X), "(Intercept)"), logical(1))))
})
