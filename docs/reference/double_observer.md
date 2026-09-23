# Double-observer abundance family

Abundance from a double-observer protocol (unmarked `multinomPois` with
a double-observer pi-function). Site abundance `N ~ Poisson(lambda)` is
surveyed by two observers with detection `p1` / `p2`. The state
`formula` models `log lambda`; `detection` the shared site-level
per-observer detection design (observers 1 and 2 carry separate
coefficients). By Poisson-multinomial thinning the observable cell
counts are independent Poissons, so the marginal is closed form with no
latent-abundance summation.

## Usage

``` r
double_observer(type = c("independent", "dependent"))
```

## Arguments

- type:

  `"independent"` (default; `N x 3` cell counts) or `"dependent"`
  (removal-style, `N x 2` cell counts, needs `primary =` for
  identifiability).

## Value

A `tobs_family` object for
[`tobs()`](https://gillescolling.com/tulpaObs/reference/tobs.md).

## Independent vs dependent protocol

`type = "independent"` (the default): the two observers detect
independently and each individual is recorded by observer 1 only,
observer 2 only, or both, so `y` is an `N x 3` matrix of cell counts in
that column order and the three cells are
`pi = (p1(1 - p2), (1 - p1)p2, p1 p2)`.

`type = "dependent"`: a removal-style protocol where a *primary*
observer records what it detects and a *secondary* observer records only
what the primary missed, so `y` is an `N x 2` matrix
`(primary-detected, secondary-only)` with cells
`pi = (p_pri, (1 - p_pri)p_sec)`. A single fixed primary observer gives
only two cells for three parameters (`lambda`, `p1`, `p2`), a ridge that
does not separate the two observer detections; observer
**role-swapping** identifies them – pass `primary =` (a length-`N`
vector in `{1, 2}` naming each site's primary observer), and with
observer 1 primary at some sites and observer 2 at others
`(p_pri, p_sec)` alternates between `(p1, p2)` and `(p2, p1)`, giving
the four cell means that recover `(lambda, p1, p2)`. With `p1 = p2` (an
interchangeable pair) the dependent protocol reduces to a two-pass
[`removal()`](https://gillescolling.com/tulpaObs/reference/removal.md)
model.

## Examples

``` r
# \donttest{
sim <- simulate_double_observer(N = 200, seed = 1)
fit <- tobs(~ abund_cov1, data = sim$data, family = double_observer(),
            detection = ~ det_cov1, y = sim$y, control = list(verbose = FALSE))
coef(fit)

# Dependent (role-swapping) protocol:
sd <- simulate_double_observer(N = 300, type = "dependent", seed = 1)
fd <- tobs(~ abund_cov1, data = sd$data, family = double_observer("dependent"),
           detection = ~ det_cov1, y = sd$y, primary = sd$primary,
           control = list(verbose = FALSE))
coef(fd)
# }
```
