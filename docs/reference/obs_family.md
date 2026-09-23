# Construct a tobs family object

Low-level constructor for `tobs_family` objects. End users should call
the specific family functions
([`occu()`](https://gillescolling.com/tulpaObs/reference/occu.md),
[`abun()`](https://gillescolling.com/tulpaObs/reference/abun.md),
[`cover()`](https://gillescolling.com/tulpaObs/reference/cover.md), ...)
rather than this directly.

## Usage

``` r
obs_family(
  name,
  class_long,
  latent,
  observation,
  replicates = c("required", "optional", "single"),
  default_engine = c("laplace", "nested_laplace", "nuts", "pg_gibbs"),
  status = c("working", "planned", "experimental"),
  params = list(),
  control_keys = character(0),
  control_groups = character(0),
  response = c("matrix", "vector")
)
```

## Arguments

- name:

  short slug, e.g. `"occu"`.

- class_long:

  human-readable name, e.g. `"single-season occupancy"`.

- latent:

  latent-state distribution, e.g. `"bernoulli"`, `"poisson"`,
  `"lognormal"`, `"hurdle"`.

- observation:

  observation likelihood, e.g. `"binomial_detection"`, `"binomial_N"`,
  `"beta"`, `"distance_binned"`.

- replicates:

  one of `"required"`, `"optional"`, `"single"`.

- default_engine:

  `"laplace"`, `"nested_laplace"`, `"nuts"`, or `"pg_gibbs"`. This is
  the route `method = "auto"` resolves to.

- status:

  `"working"`, `"planned"`, or `"experimental"`.

- params:

  named list of family-specific parameters carried with the object
  (K_max, positive-part link, etc.).

- control_keys:

  character vector of extra `control` names this family's dispatcher
  accepts beyond the engine/route controls. These are added to the
  allowlist
  [`tobs()`](https://gillescolling.com/tulpaObs/reference/tobs.md)
  validates `control` against, so a family with a bespoke dispatcher
  (e.g. the cover hurdle's grid controls) is not rejected. Admitted on
  every route the family supports.

- control_groups:

  character vector of capability groups (names of
  `.tobs_control_groups`) this family participates in beyond the ones
  its route admits unconditionally. Unlike `control_keys` these stay
  route-gated, so a group tied to the Laplace engines is still rejected
  under `"nuts"`. `"block_coordinate"` is the case this exists for:
  `max.outer` / `factor.starts` mean something only to a family whose
  latent structure is fit by the block-coordinate driver.

- response:

  shape of the family's response: `"vector"` for a plain length-N
  response vector (the cover hurdle), or `"matrix"` (the default) for a
  detection-history / count matrix, 3D array, or list. A
  single-vector-response family accepts the response on the top formula
  left-hand side (`response ~ predictors`) so `y =` may be omitted; a
  matrix-response family takes the response via `y =` only (a matrix
  does not sit on a formula LHS). Consulted by
  [`tobs()`](https://gillescolling.com/tulpaObs/reference/tobs.md).

## Value

A `tobs_family` object.

## Examples

``` r
f <- obs_family("my_occu", "custom occupancy", latent = "bernoulli",
                observation = "binomial_detection")
f
```
