# Predict from an occu_categorical fit

Returns, per row of `newdata`, the presence probability `psi`, the
conditional class probabilities `p_1..p_K` (given present), and the
unconditional class probabilities `psi * p_k` plus the absent
probability `1 - psi`.

## Usage

``` r
# S3 method for class 'occu_categorical_fit'
predict(object, newdata, ...)
```

## Arguments

- object:

  an `occu_categorical_fit`.

- newdata:

  a data frame of predictors for the shared formula.

- ...:

  ignored.

## Value

A list with `psi` (length n), `cond` (n x K conditional class
probabilities), and `joint` (n x (K+1): absent, then class 1..K).

## Examples

``` r
# \donttest{
sim <- simulate_occu_categorical(N = 300, seed = 1)
fit <- tobs(~ x, data = sim$data, family = occu_categorical(), y = sim$y)
pr <- predict(fit, newdata = data.frame(x = c(-1, 0, 1)))
round(pr$joint, 3)
# }
```
