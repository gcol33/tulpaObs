# Within / between (Mundlak) decomposition for longitudinal data

For each variable in `vars`, adds two new columns to `data`:

## Usage

``` r
within_between(data, group, vars, suffix = c("_btw", "_wtn"), na.rm = TRUE)
```

## Arguments

- data:

  A data frame.

- group:

  Character name (or character vector of names) identifying the grouping
  column(s). With multiple names, grouping is by the interaction.

- vars:

  Character vector of numeric columns to decompose. Must all exist in
  `data`.

- suffix:

  Length-2 character vector giving the suffixes for the between- and
  within-group columns. Default `c("_btw", "_wtn")`.

- na.rm:

  Logical; passed to [`mean()`](https://rdrr.io/r/base/mean.html) when
  computing per-group means. Default `TRUE`.

## Value

The input `data` with `2 * length(vars)` new columns appended.

## Details

- `<var><suffix[1]>` - the per-group mean of `var` (the
  **between**-group component, constant within each group),

- `<var><suffix[2]>` - `var - <var><suffix[1]>` (the **within**-group
  deviation, summing to zero within each group).

The decomposition is exact: `data[[var]] == data[[btw]] + data[[wtn]]`.
Including both columns in a regression (`~ x_btw + x_wtn`) separates the
cross-group association from the within-group association without
imposing that they be the same coefficient (the random-effects
assumption).

Typical use in the MOTIVATE cover-resurvey workflow:

    plots <- within_between(plots, group = "plot_id", vars = "year")
    fit   <- tobs(~ year_btw + year_wtn + bym2(graph = adj), data = plots,
                  family = cover("beta"), y = plots$cover,
                  method = "nested_laplace")

`year_wtn` then carries the within-plot temporal trend; `year_btw`
absorbs across-plot timing heterogeneity (older versus newer baseline
plots).

## References

Mundlak, Y. (1978). On the pooling of time series and cross section
data. *Econometrica*, 46(1), 69-85.

Bell, A. & Jones, K. (2015). Explaining fixed effects: random effects
modeling of time-series cross-sectional and panel data. *Political
Science Research and Methods*, 3(1), 133-153.

## Examples

``` r
set.seed(1)
plots <- data.frame(
  plot_id = rep(1:5, each = 4),
  year    = rep(2000:2003, times = 5) + rep(c(0, 5, 10, 15, 20), each = 4),
  cover   = stats::runif(20)
)
decomposed <- within_between(plots, group = "plot_id", vars = "year")
head(decomposed)
```
