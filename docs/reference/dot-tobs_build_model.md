# Build a tobs model object

Inferred model type from arguments:

- **Single-season**: no `col_formula`, no `species`

- **Dynamic**: `col_formula` and/or `ext_formula` provided

- **Integrated**: `integrated = TRUE`, `y` a list of matrices

## Usage

``` r
.tobs_build_model(
  occ_formula,
  det_formula = NULL,
  data,
  y,
  col_formula = NULL,
  ext_formula = NULL,
  species = NULL,
  integrated = FALSE,
  abundance = FALSE,
  count = FALSE,
  count_response = "poisson",
  count_trials = NULL,
  det_visit_formula = NULL,
  det_visit_data = NULL,
  site.map = NULL
)
```

## Details

Community families (ms_occu / ms_dyn_occu / ms_int_occu / ms_occu_cover
/ ms_count / jsdm) have their own in-tree binders + Laplace-EM fitter
and do not pass through here.

A structured term enters the process whose formula it is written in. For
an integrated model this makes `detection = ~ x + spde(lon, lat, ...)` a
**detection-arm** field: the occupancy state stays whatever the
occupancy formula says, and the continuous Matern field sits on the
detection logit. The detection arm there is one arm observed through S
sources, so the field is one structure fit once per source block –
source `s` fits its own realization at its own sites, and two sources
covering the same location each carry their own value there, the same
way each source carries its own detection coefficients.
`fit$spatial_field_det` is the per-source named list of mesh fields
(source names as given on `y`); project a source's field to its sites
with `fit$spatial$tulpa_spec$A`.

The arm accepts a continuous
[`spde()`](https://gillescolling.com/tulpaObs/reference/tobs_terms.md)
field under `method = "laplace"`. The areal kinds and the temporal / re
/ svc / latent classes are grid-integrated as latent blocks on the state
arm, so they are written on the occupancy formula; on a detection
formula they error with a pointer rather than being fit against the
wrong arm.

Each source of an integrated response surveys its own subset of the
sites in `data`, in its own order, so a source's rows are joined to the
site-level frame on a site key rather than read in row order. The key is
`site.map` when given – a list of one entry per source, each naming or
indexing one site of `data` per source row – otherwise
`rownames(y[[s]])` matched against `rownames(data)`. A source carrying
neither key is keyed by row position, which requires one row per site; a
shorter one errors with a pointer to `site.map`.
[`simulate_int_occu()`](https://gillescolling.com/tulpaObs/reference/simulate_int_occu.md)
names its response rows, so its output joins on that key.
