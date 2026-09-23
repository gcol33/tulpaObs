# tulpaObs: Bayesian Occupancy, Abundance and Detection Models via 'tulpa'

Fits Bayesian models for species observed with imperfect detection:
single-season and dynamic site occupancy (MacKenzie et al. (2002)
[doi:10.1890/0012-9658(2002)083\[2248:ESORWD\]2.0.CO;2](https://doi.org/10.1890/0012-9658%282002%29083%5B2248%3AESORWD%5D2.0.CO%3B2)
), N-mixture abundance from replicated counts (Royle (2004)
[doi:10.1111/j.0006-341X.2004.00142.x](https://doi.org/10.1111/j.0006-341X.2004.00142.x)
), distance sampling, removal and double-observer counts, false-positive
and time-to-detection occupancy, integrated models combining several
data sources, multi-species and joint species distribution models, and a
hurdle model for vegetation cover. Linear predictors take fixed, random,
spatial and temporal effects. Each model supplies its observation
likelihood to the 'tulpa' engine, which fits it by Laplace
approximation, by nested Laplace integration over the hyperparameters,
or by Hamiltonian Monte Carlo with the No-U-Turn Sampler (NUTS) using
automatic differentiation. Simulation, prediction and diagnostic
functions (posterior predictive checks, leave-one-out cross-validation,
simulation-based calibration) accompany each model.

## See also

Useful links:

- <https://github.com/gcol33/tulpaObs>

- Report bugs at <https://github.com/gcol33/tulpaObs/issues>

## Author

**Maintainer**: Gilles Colling <gilles.colling051@gmail.com>
([ORCID](https://orcid.org/0000-0003-3070-6066)) \[copyright holder\]

Authors:

- Gilles Colling <gilles.colling051@gmail.com>
  ([ORCID](https://orcid.org/0000-0003-3070-6066)) \[copyright holder\]
