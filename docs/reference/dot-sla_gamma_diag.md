# Assemble per-coefficient skewness from per-site third derivatives

Implements equation (2.2) of
`dev_notes/simplified_laplace_derivation.md`: gamma_j = sigma_j^(-3) \*
sum_i l3_i \* v\_(i,j)^3, v\_(i,j) = (X Sigma)\_(i,j)

## Usage

``` r
.sla_gamma_diag(l3, X, Sigma)
```

## Arguments

- l3:

  Per-site third derivative (length n).

- X:

  Design matrix (n x p).

- Sigma:

  Posterior covariance of the fixed-effect coefficients (p x p, must be
  symmetric PD).

## Value

Vector of length p with one skewness per coefficient.

## Details

For *diagonal-in-eta* likelihoods only (binomial, Poisson, Gaussian).
