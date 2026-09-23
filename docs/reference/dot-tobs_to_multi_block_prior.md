# Build a multi-block latent prior list for tulpa::tulpa_nested_laplace()

Converts the tobs-level spatial / temporal / RE specs into the
list-of-blocks shape that
[`tulpa::tulpa_nested_laplace()`](https://gillescolling.com/tulpa/reference/tulpa_nested_laplace.html)
expects under its multi-block dispatch (`.is_multi_block_prior`).

## Usage

``` r
.tobs_to_multi_block_prior(
  spatial = NULL,
  temporal = NULL,
  re = NULL,
  model,
  grids = list()
)
```

## Arguments

- spatial:

  Optional `tobs_spatial` term (from an
  [`icar()`](https://gillescolling.com/tulpaObs/reference/tobs_terms.md)
  /
  [`bym2()`](https://gillescolling.com/tulpaObs/reference/tobs_terms.md)
  formula term). Only BYM2 / ICAR are wired through to the multi-block
  engine at present; GP / multiscale_gp / SVC are not yet supported and
  raise.

- temporal:

  Optional `tobs_temporal` term (from a
  [`temporal()`](https://gillescolling.com/tulpaObs/reference/tobs_terms.md)
  formula term). Types `"ar1"`, `"rw1"`, `"rw2"`, `"iid"` are supported.

- re:

  Optional list of `tobs_re` terms (from
  [`re()`](https://gillescolling.com/tulpaObs/reference/tobs_terms.md)
  formula terms). Only `model = "iid"` terms are converted to IID latent
  blocks here; correlated structures (ar1 / rw1 / rw2 on RE groups) are
  passed through as temporal-like blocks.

- model:

  A `tobs_model` from
  [`.tobs_build_model()`](https://gillescolling.com/tulpaObs/reference/dot-tobs_build_model.md).
  The structured terms carry pre-resolved index codes; `model` pins
  `N = model$n_sites` for single-season fits.

## Value

`NULL` when no latent block is supplied; a single-block list when
exactly one is supplied; a list-of-blocks otherwise. Each block is the
minimal field set that the tulpa multi-block dispatch fills defaults
around (`.NL_REGISTRY[[type]]$defaults`).
