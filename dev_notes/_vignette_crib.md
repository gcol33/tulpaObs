# tulpaObs vignette crib (verified against current code, 2026-05-29)

All calls below were run under `devtools::load_all(".")` with R 4.6.0 and
recover their simulated truth. Use these as the spine; read the cited R/
source before adding anything not listed here. **Never invent a function.**

## Verification (mandatory)

After writing your `.Rmd`, run:

```
"/c/Program Files/R/R-4.6.0/bin/Rscript.exe" dev_notes/_render_check.R vignettes/<your-file>.Rmd
```

It purls the vignette and executes every chunk under `load_all()`. Iterate
until it prints `==== OK: <file> ====` with no errors. The chunks must run
clean — a single failing chunk fails the vignette.

## Global rules

- Default fit is `method = "laplace"`. The EM prints iteration progress
  unless you pass `control = list(verbose = FALSE)`. Put `verbose = FALSE`
  in the control of **every** laplace/nested_laplace fit so output is clean.
- NUTS control names are dotted: `n.iter`, `n.warmup`, `n.chains`, `seed`,
  `adapt.delta`, `max.treedepth`, `n.thin`, `n.threads`. There is NO `iter`
  or `warmup` — those error.
- `logLik()` returns `NaN` for occupancy Laplace fits — do NOT showcase it.
- `tidy()` / `glance()` are broom generics and are NOT available unless
  `library(broom)` is attached — avoid them.
- `confint()`, `vcov()`, `summary()`, `coef()` all work (inherited from
  `tulpa::tulpa_fit`). `summary()` adds Rhat/ESS columns under NUTS.
- Math in Rmd uses LaTeX `$...$` / `$$...$$` (Unicode is fine in Rmd, unlike
  Rd). Keep figures with `fig.alt = "..."`.
- YAML header + setup chunk: copy exactly from `vignettes/quickstart.Rmd`
  (output `rmarkdown::html_vignette`, `set.seed(20260529)`, collapse/comment).

## Voice (read these first)

Read `vignettes/quickstart.Rmd` (just written) and
`../tulpa/vignettes/spatial-models.Rmd` for the target voice. Collegial,
concrete-first (show the R behaviour, then name the concept), flat
statements, specific numbers not adjectives, 2-4 sentences of prose between
chunks. No AI tells (no "delve/robust/comprehensive", no em-dash overuse,
no "It's important to note", no aphoristic contrast pairs, no "not X but Y").
Max 2 em-dashes per file.

## Verified family calls

### Single-season occupancy — occu()  (R/occu.R, R/data.R)
```r
sim <- simulate_occu(N = 300, J = 6, n_occ_covs = 1, n_det_covs = 1,
                     beta_occ = c(0.3, 1.0), beta_det = c(0.7, 0.6), seed = 1)
# sim$y is N x J matrix; sim$data has occ_cov1, det_cov1; sim$truth has
#   beta_occ, beta_det, psi, p, z
fit <- tobs(~ occ_cov1, data = sim$data, family = occu(),
            detection = ~ det_cov1, y = sim$y, method = "laplace",
            control = list(verbose = FALSE))
coef(fit)              # list: $psi (named vec), $p (named vec)
summary(fit); confint(fit)
fitted(fit)            # list: $psi, $p, $z (per-site)
tobs_waic(fit)         # list: $waic, $elpd, $p_waic, $lppd
tobs_marginal_effect(fit, "occ_cov1", process = "occupancy")  # plot()-able
predict(fit, X.0 = cbind(`(Intercept)`=1, occ_cov1=c(-1,0,1)))
# NUTS:
tobs(~ occ_cov1, data=sim$data, family=occu(), detection=~det_cov1, y=sim$y,
     method="nuts", control=list(n.iter=800, n.warmup=400, seed=1, verbose=FALSE))
```

### N-mixture abundance — abun()  (R/abun.R)
```r
sa <- simulate_abun(N=200, J=5, n_abund_covs=1, n_det_covs=1,
                    beta_lambda=c(log(3),0.5), beta_p=c(0.4,0.3), seed=2)
# sa$truth: beta_lambda, beta_p, lambda, p, N, mixture, size
fa <- tobs(~ abund_cov1, data=sa$data, family=abun(), detection=~det_cov1,
           y=sa$y, method="laplace", control=list(verbose=FALSE))
coef(fa)               # $lambda (log link! intercept back-transforms exp()), $p
fitted(fa)             # nmix fitted (lambda/p/N)
predict(fa)            # nmix prediction; type "abundance"/"detection"
# Negative binomial:
sb <- simulate_abun(N=200, J=5, n_abund_covs=1, n_det_covs=1, mixture="negbin",
                    size=2, beta_lambda=c(log(4),0.4), beta_p=c(0.5,0.2), seed=5)
fb <- tobs(~ abund_cov1, data=sb$data, family=abun(mixture="negbin"),
           detection=~det_cov1, y=sb$y, method="laplace", control=list(verbose=FALSE))
fb$nmix_dispersion     # list: $r, $log_r, $r_sd
# Areal spatial (Poisson or NB): bym2()/icar()/car_proper() on abundance arm,
# method="nested_laplace". adj = dense 0/1 adjacency matrix, one unit per site.
# Verify this call yourself with the render-check before relying on it.
```

### Community occupancy — ms_occu()  (R/occu.R)
```r
sm <- simulate_ms_occu(N=120, J=5, n_species=8, seed=3)   # sm$y is N x J x n_species
fm <- tobs(~ x, data=sm$data, family=ms_occu(), detection=~1, y=sm$y,
           species=paste0("sp",1:8), method="laplace", control=list(verbose=FALSE))
coef(fm)               # community means $psi, $p
tobs_richness(fm)      # data.frame: site, mean, sd, q2.5, q97.5
# NOTE: ranef(fm) is EMPTY for ms_occu (community-RE handled internally) — do
# not showcase ranef() here. sm$truth has beta_comm_mean/_sd, alpha_comm_mean/_sd.
```

### Community N-mixture — ms_abun()  (R/ms_abun.R)
```r
sab <- simulate_ms_abun(n_species=8, N=80, J=4, n_abund_covs=1, n_det_covs=1, seed=7)
fab <- tobs(~ abund_cov1, data=sab$data, family=ms_abun(), detection=~det_cov1,
            y=sab$y, species=paste0("sp",1:8), method="laplace",
            control=list(verbose=FALSE))
coef(fab)              # community means $lambda, $p
ranef(fab)             # per-species coefficient deviations (WORKS for ms_abun)
# sab$truth: mu_lambda, mu_p, sd_lambda, sd_p
```

### Dynamic (multi-season) occupancy — dyn_occu()  (R/occu.R, R/data.R)
```r
sd <- simulate_dyn_occu(N=100, J=4, n_seasons=5, beta_occ=c(0.5), beta_det=c(0),
                        gamma=0.2, epsilon=0.1, seed=1)   # sd$y is N x J x n_seasons array
# dyn_occu REQUIRES col_formula (colonisation) and ext_formula (extinction):
fit <- tobs(~ 1, data=sd$data, family=dyn_occu(), detection=~1, y=sd$y,
            col_formula=~1, ext_formula=~1, method="laplace",
            control=list(verbose=FALSE))
fit$intercepts         # $psi1, and colonisation/extinction
# y MUST be a 3D array [sites x visits x seasons]; a 2D matrix errors with "3D array".
```

### Integrated occupancy — int_occu()  (R/occu.R, R/data.R)
```r
si <- simulate_int_occu(N_total=150, n_data=2, J=c(4,3), n_shared=20,
                        beta_occ=c(0.5,0.3),
                        beta_det=list(c(0.2,-0.4), c(-0.1,0.3)), seed=42)
# si$y is a list of matrices (one per source)
fit <- tobs(~ x, data=si$data, family=int_occu(), detection=~1,
            y=list(s1=si$y[[1]], s2=si$y[[2]]), method="laplace",
            control=list(verbose=FALSE))
# shared psi across sources, source-specific detection.
```

### Random effects — lme4 bar syntax  (R/em_laplace_re.R, R/re_aghq.R)
```r
# (1 | g) random intercept on occupancy, Laplace via variance-component EM:
fit <- tobs(~ x + (1 | g), data=d, y=y, detection=~1, family=occu(),
            method="laplace", control=list(verbose=FALSE))
# fit$means has a sigma_ hyperparameter; ranef(fit) returns per-group BLUPs
#   (data.frame: group, level, term, estimate, std.error).
# RE on detection: detection = ~ (1 | g). Correlated slope: (1 + x | g).
# Uncorrelated: (x || g). Build d with factor(g); see test-re-laplace-recovery.R
# for a working simulator. NUTS fits every RE form.
```

### Areal spatial occupancy — nested_laplace  (R/em_nested_laplace.R)
```r
# Build a dense 0/1 adjacency matrix adj (one row/col per site). Chain or grid.
fit <- tobs(~ x + bym2(graph = adj), data=d, family=occu(), detection=~1, y=y,
            method="nested_laplace", control=list(max.iter=10, verbose=FALSE))
# Terms: icar(graph=adj), bym2(graph=adj), car(graph=adj), car_proper(graph=adj).
# fit$nested_laplace$multi_prior[[1]]$type == "bym2".
# predict(fit, type="state") returns marginalised per-site psi with calibrated
#   psi_lower/psi_upper (and held-out all-NA sites interpolated by the field).
# See test-nested-laplace-occu.R for a panel simulator with a chain graph.
```

### Cover hurdle — cover()  (R/family_cover_hurdle.R, R/sim_cover_hurdle.R)
```r
sc <- simulate_cover(N=200, beta_occ=c(-0.3,0.8), beta_pos=c(-1,0.3),
                     sigma_pos=0.4, seed=4)
# sc$data: cover (response, 0 or positive), x, lon, lat
# There are existing cover vignettes (cover-hurdle.Rmd, cover-hurdle-vs-inla.Rmd,
# cover-hurdle-motivate.Rmd). Do NOT duplicate them.
```

### Diagnostics  (R/diagnostics.R)
```r
tobs_waic(fit)                                   # $waic, $elpd, $p_waic, $lppd
tobs_ppc(fit, fit.stat = "freeman-tukey")        # posterior predictive check
tobs_pit_residuals(fit)
tobs_test_dispersion(fit); tobs_test_zero_inflation(fit)
tobs_test_outliers(fit); tobs_test_uniformity(fit)
# residuals(fit, type = "deviance"|"pearson"|"response") -> list($occ, $det)
# Read R/diagnostics.R for exact return shapes and any Moran's I / variogram
# helpers BEFORE using them. Verify each with render-check.
```

## Data formatting  (R/data.R)
- `tobs_format()`, `tobs_format_ms()` — read R/data.R for signatures and the
  expected `y` structure (matrix for single, list for integrated, array for
  community/dynamic). `tobs_data()` builds/validates the data object;
  `summary()`/`plot()` methods exist on `tobs_data`.
