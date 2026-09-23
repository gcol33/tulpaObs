# Simulate from an occu_categorical model

Simulate from an occu_categorical model

## Usage

``` r
simulate_occu_categorical(
  N = 500L,
  beta_occ = c(0.2, 0.8),
  beta_class = NULL,
  seed = NULL
)
```

## Arguments

- N:

  number of units.

- beta_occ:

  length-2 presence coefficients (intercept, slope on `x`).

- beta_class:

  a `2 x (K-1)` matrix of class coefficients (rows intercept / slope on
  `x`; columns the non-baseline classes), or `NULL` for a built-in
  3-class default.

- seed:

  optional RNG seed.

## Value

A list with `data` (a data frame with `x`), `y` (length-N integer in
`0..K`), and `truth`.

## Examples

``` r
sim <- simulate_occu_categorical(N = 100, seed = 1)
table(sim$y)
```
