# Generic finite-difference SLA gamma for non-diagonal families

Computes per-coefficient SLA gamma by 5-point central finite difference
of the original observation log-likelihood along the direction
`Sigma[, j]` in beta-space. Use this when the family's third derivative
in eta-space does not decompose into a per-site sum (e.g. HMM forward
likelihood, integrated shared-process, hurdle joints).

## Usage

``` r
.sla_gamma_fd(beta_hat, Sigma, log_lik_fn, h = NULL)
```

## Arguments

- beta_hat:

  Mode of the joint posterior (length p).

- Sigma:

  Posterior covariance at the mode (p x p, symmetric PD).

- log_lik_fn:

  Callable: takes a length-p beta vector, returns scalar log-likelihood
  at that beta (no prior contribution – SLA gamma uses the observation
  likelihood only, per RMC 2009 sec.3.2).

- h:

  Optional step size override (scalar or length-p vector).

## Value

Numeric vector of length p with per-coefficient gamma.

## Details

Step size defaults to `eps^{1/5} * sigma_j / ||v_j||`, which keeps the
displacement in beta-space on the natural scale of the j-th marginal
posterior. Override `h` (scalar or length-p vector) to inspect
truncation/cancellation behaviour.
