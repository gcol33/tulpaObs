# Third derivative of binomial log-likelihood wrt logit(p)

For y_i ~ Binomial(n_i, plogis(eta_i)), d^3 log p / d eta^3 = -n_i \*
p_i \* (1 - p_i) \* (1 - 2 \* p_i)

## Usage

``` r
.l3_binomial_logit(eta, n_trials = 1L)
```

## Arguments

- eta:

  Linear predictor at sites (numeric vector).

- n_trials:

  Per-site trial counts (integer or numeric vector, recycled to
  length(eta) if scalar).

## Value

Vector of third derivatives (length = length(eta)).
