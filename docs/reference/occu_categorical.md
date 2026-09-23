# Presence + nominal class hurdle family

A hurdle for a response that is either absent or, when present, one of
`K` **nominal (unordered)** classes – a colour morph, microhabitat /
substrate use, or a cryptic-species / classifier label observed given
the organism is present. The presence/absence part is a Bernoulli arm;
the class given present is a baseline-category multinomial logit (the
last class is the baseline). The two arms factorise the likelihood
exactly: \$\$P(y = 0) = 1 - \psi, \qquad P(y = k) = \psi \\ p_k,\$\$
with \\p = \mathrm{softmax}(X \beta\_{class})\\.

## Usage

``` r
occu_categorical(classes = NULL)
```

## Arguments

- classes:

  optional character vector of class labels (length `K`), used only to
  name the coefficient blocks; when `NULL`, `K` is taken from `max(y)`
  and the classes are labelled `1..K`.

## Value

A `tobs_family` object.

## Details

This is the categorical counterpart of
[`cover()`](https://gillescolling.com/tulpaObs/reference/cover.md)
(presence + a magnitude): here the positive part is an unordered class,
so it uses a multinomial logit rather than beta / lognormal. For an
*ordered* class response (Braun-Blanquet cover bands) use
`cover(response = "ordinal")`, which exploits the ordering; this family
is for classes with no ordering.

## Response

`y` is a length-N integer vector in `0..K`: `0` is absent, `k` is class
`k`. It may sit on the formula left-hand side (`y ~ ...`), so `y =` can
be dropped. The `formula` predictors are shared by both arms.

## Scope

Non-spatial Laplace (`method = "laplace"`). The multinomial math is the
FD-validated tulpa kernel (`multinomial_logit.h`); the non-spatial fit
is the vectorised R Newton over the same closed forms. Spatial fields /
NUTS (the native multi-process likelihood) and the latent-class
*misclassification* variant (the K-class generalisation of
[`fp_occu()`](https://gillescolling.com/tulpaObs/reference/fp_occu.md),
a confusion matrix on the observed label) are documented follow-ups.

## See also

[`cover()`](https://gillescolling.com/tulpaObs/reference/cover.md)
(presence + magnitude),
[`fp_occu()`](https://gillescolling.com/tulpaObs/reference/fp_occu.md)
(two-state false-positive detection).

## Examples

``` r
f <- occu_categorical(classes = c("red", "green", "blue"))
f
```
