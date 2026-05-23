## Summary

`tobs()` now accepts `lme4` bar syntax as sugar over `re()` — `(1 | g)`,
`(x | g)`, `(x || g)` are rewritten into `re()` calls on the formula AST by
`.tobs_desugar_bars()` (`R/formula_parse.R`) before `terms()` runs. Three bar
forms are deliberately rejected with an informative error because the current
`re()`/encoding layer cannot hold them in one block:

1. **Slope without intercept** — `(0 + x | g)`. The engine's random-slope block
   always carries the group intercept (`build_re_spec` hardcodes
   `n_coefs <- 2L`), so a pure slope is not expressible.
2. **Multiple slopes in one bar** — `(1 + x + z | g)`. `re()` takes a single
   `covariate`, and `build_re_spec` again hardcodes `n_coefs <- 2L`.
3. **Nested / interaction grouping** — `(1 | g/h)`, `(1 | g:h)`. The bar
   grouping must currently be a single factor; build the interaction explicitly
   (`re(interaction(a, b))`).

Errors are raised in `.tobs_bar_to_re()` (`R/formula_parse.R`) and point users at
an explicit `re()`.

## Note on scope (not an upstream `tulpa` change)

The C++ populate path already supports arbitrary block sizes: `populate_helpers.h`
loops over `n_coefs[t]` slope columns and sizes the Cholesky factor as
`k*(k+1)/2`. So **multi-slope correlated RE blocks are an R-side enhancement**
to `build_re_spec()` (`R/occu_fit.R`) plus a multi-covariate `re()` signature —
not an engine/`tulpa` change.

## Proposed work

- [ ] `re()` accepts multiple slope covariates (e.g. `covariate = ` taking a
      vector/one-sided formula), and an explicit "no intercept" option for the
      slope-only case.
- [ ] `build_re_spec()` computes `n_coefs` from the covariate count and stacks
      all slope columns into `slope_matrices[[t]]`.
- [ ] `.tobs_bar_to_re()` maps `(1 + x + z | g)` / `(0 + x | g)` once `re()`
      supports them; optionally expand `g/h` to `g + g:h`.
- [ ] recovery test for a 3x3 correlated RE block.

Until then, users write the equivalent `re()` calls directly.
