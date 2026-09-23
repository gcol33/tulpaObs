# Grassland vegetation-cover panel (synthetic)

A vegetation cover-hurdle dataset for 150 grassland plots surveyed
across three years, drawn from a known generative model: the species
occurs with a probability driven by `moisture`, and where it occurs, its
conditional cover (a Beta-distributed proportion in `(0, 1)`) increases
with `moisture` and decreases with `grazing`. Roughly half the plots
record zero cover. `year_c` is the centred survey year, so
[`within_between()`](https://gillescolling.com/tulpaObs/reference/within_between.md)
can split the cross-plot baseline from the within-plot temporal
deviation. Synthetic, not field data; the generative coefficients are
attached as `attr(meadow_cover, "truth")`.

## Usage

``` r
meadow_cover
```

## Format

A data frame with 150 rows (one per plot) and columns:

- plot:

  Plot identifier (factor).

- year:

  Survey year (2018-2020).

- year_c:

  Survey year centred on its mean.

- moisture:

  Standardised soil-moisture index.

- grazing:

  Standardised grazing-pressure index.

- cover:

  Recorded cover as a proportion in `[0, 1]` (0 where absent).

## See also

[`cover()`](https://gillescolling.com/tulpaObs/reference/cover.md),
[`tobs()`](https://gillescolling.com/tulpaObs/reference/tobs.md),
[`within_between()`](https://gillescolling.com/tulpaObs/reference/within_between.md),
[`simulate_cover()`](https://gillescolling.com/tulpaObs/reference/simulate_cover.md)

## Examples

``` r
data(meadow_cover)
mean(meadow_cover$cover > 0)
# \donttest{
fit <- tobs(cover ~ moisture, data = meadow_cover,
            family = cover(response = "beta"))
coef(fit)
# }
```
