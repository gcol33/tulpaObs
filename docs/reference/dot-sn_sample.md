# Sample from a skew-normal distribution

Two-component method (Azzalini 1985): Z1, Z2 ~ N(0,1) i.i.d., then Z =
delta \* \|Z1\| + sqrt(1 - delta^2) \* Z2 X = xi + omega \* Z

## Usage

``` r
.sn_sample(n, sn)
```

## Arguments

- n:

  Number of samples.

- sn:

  Skew-normal parameter list from
  [`tulpa::sn_match()`](https://gillescolling.com/tulpa/reference/sn_match.html)
  (elements `xi`, `omega`, `alpha`).

## Value

Numeric vector of length n.
