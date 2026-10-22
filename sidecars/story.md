# Story: tulpaObs

## Claim

Bayesian occupancy, abundance and detection models that separate the ecological state from imperfect detection, fitted through the tulpa engine.

## Audience

Ecologists with repeated-visit detection histories, replicated counts, distance or removal data, or vegetation cover. The README assumes familiarity with data frames and lme4-style formulas; it introduces the state/detection split.

## Relation

- to: tulpa
- phrase: each model supplies its observation likelihood to the tulpa engine
- same: inference (Laplace, nested Laplace, NUTS), structured effects (spatial, temporal, random effects), diagnostics generics
- differs: observation likelihoods for imperfect detection, survey data helpers, simulators per family
- home: vignettes/quickstart.Rmd

## Arcs

### A1 Separate occupancy from detection

- question: Is a zero an absence or a missed detection?
- answer: every family writes the latent state and the detection process as separate levels and integrates the state out.
- evidence: vignettes/occupancy.Rmd; vignettes/occupancy-vs-inla.Rmd
- home: vignettes/occupancy.Rmd
- also: README

### A2 Choose the model for the survey design

- question: Which family matches my data?
- answer: one fitter, `tobs()`, and one family constructor per model.
- evidence: _pkgdown.yml reference groups; family help pages
- home: reference index
- also: README

### A3 Structured effects

- question: How do I add spatial, temporal and grouped effects?
- answer: terms inside the formula (`spatial()`, `temporal()`, lme4 bars).
- evidence: vignettes/spatial-occupancy.Rmd; vignettes/occupancy-spatial-spde.Rmd; vignettes/random-effects.Rmd
- home: vignettes/spatial-occupancy.Rmd
- also: README

### A4 Choose inference and check the fit

- question: Laplace or NUTS, and is the fit trustworthy?
- answer: `method` picks the route; simulators with known truth and diagnostics check it.
- evidence: vignettes/diagnostics.Rmd; vignettes/random-effects.Rmd (AGHQ bias numbers)
- home: vignettes/diagnostics.Rmd
- also: README

## Documents

- DESCRIPTION: claim, relation, the model list.
- README.Rmd / README.md: claim, relation, installation, first fit, one section per arc, routing, citation.
- vignettes/quickstart.Rmd: Get started.
- Topic vignettes: one family or structure each; two comparisons with INLA.
- Help pages and API.md: argument and return-value contracts.
