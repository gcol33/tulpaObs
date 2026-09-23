# S3 methods for multi-arm fits

Coefficient, covariance, and summary accessors for fits that store
independent per-arm coefficient blocks (`cover_fit`,
`occu_categorical_fit`) rather than a flat posterior. The arms are fit
independently, so the joint covariance is block-diagonal. When a full
posterior is present (the NUTS cover path) the inference methods defer
to the draw-based `tulpa_fit` methods.

## Usage

``` r
# S3 method for class 'tobs_multiarm_fit'
nobs(object, ...)

# S3 method for class 'tobs_multiarm_fit'
coef(object, arm = NULL, ...)

# S3 method for class 'tobs_multiarm_fit'
vcov(object, ...)

# S3 method for class 'tobs_multiarm_fit'
confint(object, parm, level = 0.95, ...)

# S3 method for class 'tobs_multiarm_fit'
logLik(object, ...)

# S3 method for class 'tobs_multiarm_fit'
glance(x, ...)

# S3 method for class 'tobs_multiarm_fit'
tidy(x, conf.level = 0.95, ...)

# S3 method for class 'tobs_multiarm_fit'
summary(object, level = 0.95, ...)

# S3 method for class 'tobs_multiarm_fit'
plot(x, ...)

# S3 method for class 'tobs_multiarm_fit'
fitted(object, ...)

# S3 method for class 'tobs_multiarm_fit'
residuals(object, ...)
```

## Arguments

- object, x:

  A `tobs_multiarm_fit` (a `cover_fit` or `occu_categorical_fit`).

- ...:

  Ignored, or forwarded to
  [`NextMethod()`](https://rdrr.io/r/base/UseMethod.html) on the
  posterior path.

- arm:

  Optional arm name (`"presence"`, `"positive"` or `"class"`) for
  [`coef()`](https://rdrr.io/r/stats/coef.html); see
  [`coef.tobs_fit()`](https://gillescolling.com/tulpaObs/reference/coef.tobs_fit.md).

- parm:

  Ignored (present for
  [`confint()`](https://rdrr.io/r/stats/confint.html) generic
  compatibility).

- level:

  Confidence level for
  [`confint()`](https://rdrr.io/r/stats/confint.html) and
  [`summary()`](https://rdrr.io/r/base/summary.html) (default 0.95).

- conf.level:

  Interval level for
  [`tidy()`](https://generics.r-lib.org/reference/tidy.html) (default
  0.95).

## Value

The layout every `tobs_fit` returns:
[`nobs()`](https://rdrr.io/r/stats/nobs.html) an integer;
[`coef()`](https://rdrr.io/r/stats/coef.html) a named vector
(`presence_(Intercept)`, `positive_x`, `class_2:x`, ...);
[`vcov()`](https://rdrr.io/r/stats/vcov.html) the block-diagonal
covariance and [`confint()`](https://rdrr.io/r/stats/confint.html) a
two-column matrix, both named as
[`coef()`](https://rdrr.io/r/stats/coef.html);
[`logLik()`](https://rdrr.io/r/stats/logLik.html) a `logLik` object;
[`tidy()`](https://generics.r-lib.org/reference/tidy.html) the `arm` /
`term` / `estimate` / `std.error` / `conf.low` / `conf.high` table;
[`summary()`](https://rdrr.io/r/base/summary.html) the estimate /
std.error / interval data frame;
[`glance()`](https://generics.r-lib.org/reference/glance.html) the
one-row column set of
[`glance.tobs_fit()`](https://gillescolling.com/tulpaObs/reference/glance.tobs_fit.md);
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) one
Gaussian-density panel per coefficient (up to 4), invisibly.
