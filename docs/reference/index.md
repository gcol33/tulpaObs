# Package index

## Front door

One fitter; the model is chosen by the family object.

- [`tobs()`](https://gillescolling.com/tulpaObs/reference/tobs.md) : Fit
  a hierarchical latent-state observation model
- [`tobs_data()`](https://gillescolling.com/tulpaObs/reference/tobs_data.md)
  : Convert long-format data to a site x visit observation object
- [`tobs_format()`](https://gillescolling.com/tulpaObs/reference/tobs_format.md)
  : Format occupancy data from matrices/lists
- [`tobs_format_ms()`](https://gillescolling.com/tulpaObs/reference/tobs_format_ms.md)
  : Format multi-species occupancy data
- [`tobs_check_id()`](https://gillescolling.com/tulpaObs/reference/tobs_check_id.md)
  : Check model identifiability
- [`tobs_terms`](https://gillescolling.com/tulpaObs/reference/tobs_terms.md)
  [`icar`](https://gillescolling.com/tulpaObs/reference/tobs_terms.md)
  [`bym2`](https://gillescolling.com/tulpaObs/reference/tobs_terms.md)
  [`car_proper`](https://gillescolling.com/tulpaObs/reference/tobs_terms.md)
  [`gp`](https://gillescolling.com/tulpaObs/reference/tobs_terms.md)
  [`multiscale_gp`](https://gillescolling.com/tulpaObs/reference/tobs_terms.md)
  [`spde`](https://gillescolling.com/tulpaObs/reference/tobs_terms.md)
  [`svc`](https://gillescolling.com/tulpaObs/reference/tobs_terms.md)
  [`re`](https://gillescolling.com/tulpaObs/reference/tobs_terms.md)
  [`temporal`](https://gillescolling.com/tulpaObs/reference/tobs_terms.md)
  [`latent`](https://gillescolling.com/tulpaObs/reference/tobs_terms.md)
  [`spatial`](https://gillescolling.com/tulpaObs/reference/tobs_terms.md)
  : Structured terms inside a tobs() formula
- [`tobs_get()`](https://gillescolling.com/tulpaObs/reference/tobs_get.md)
  : Extract one species' fit from a batched occu_cover fit.
- [`print(`*`<tobs_family>`*`)`](https://gillescolling.com/tulpaObs/reference/print.tobs_family.md)
  : Print method for tobs_family

## Occupancy families

- [`occu()`](https://gillescolling.com/tulpaObs/reference/occu.md) :
  Single-season occupancy family
- [`dyn_occu()`](https://gillescolling.com/tulpaObs/reference/dyn_occu.md)
  : Dynamic (multi-season) occupancy family
- [`int_occu()`](https://gillescolling.com/tulpaObs/reference/int_occu.md)
  : Integrated occupancy family
- [`dyn_int_occu()`](https://gillescolling.com/tulpaObs/reference/dyn_int_occu.md)
  : Multi-season integrated occupancy family
- [`ms_occu()`](https://gillescolling.com/tulpaObs/reference/ms_occu.md)
  : Multispecies (community) single-season occupancy family
- [`ms_dyn_occu()`](https://gillescolling.com/tulpaObs/reference/ms_dyn_occu.md)
  : Community (multispecies) dynamic occupancy family
- [`ms_int_occu()`](https://gillescolling.com/tulpaObs/reference/ms_int_occu.md)
  : Community (multispecies) integrated occupancy family
- [`jsdm()`](https://gillescolling.com/tulpaObs/reference/jsdm.md) :
  Joint species distribution family (no detection process)
- [`t_occu()`](https://gillescolling.com/tulpaObs/reference/t_occu.md) :
  Multi-season occupancy with a temporal (AR1) year effect (tPGOcc)
- [`fp_occu()`](https://gillescolling.com/tulpaObs/reference/fp_occu.md)
  : Multistate false-positive occupancy family
- [`fp_occu_laplace()`](https://gillescolling.com/tulpaObs/reference/fp_occu_laplace.md)
  : Maximum-likelihood fit of the multistate false-positive occupancy
  model
- [`royle_nichols()`](https://gillescolling.com/tulpaObs/reference/royle_nichols.md)
  : Royle-Nichols occupancy family
- [`occu_ttd()`](https://gillescolling.com/tulpaObs/reference/occu_ttd.md)
  : Time-to-detection occupancy family
- [`occu_multi()`](https://gillescolling.com/tulpaObs/reference/occu_multi.md)
  : Multi-species co-occurrence occupancy family
- [`occu_categorical()`](https://gillescolling.com/tulpaObs/reference/occu_categorical.md)
  : Presence + nominal class hurdle family
- [`occu_priors()`](https://gillescolling.com/tulpaObs/reference/occu_priors.md)
  : Weakly-informative priors for occupancy Laplace fits

## Abundance, distance and removal families

- [`abun()`](https://gillescolling.com/tulpaObs/reference/abun.md) :
  N-mixture abundance family
- [`ms_abun()`](https://gillescolling.com/tulpaObs/reference/ms_abun.md)
  : Multispecies N-mixture family
- [`dyn_abun()`](https://gillescolling.com/tulpaObs/reference/dyn_abun.md)
  : Open-population (Dail-Madsen) N-mixture family
- [`dyn_abun_laplace()`](https://gillescolling.com/tulpaObs/reference/dyn_abun_laplace.md)
  : Maximum-likelihood fit of the Dail-Madsen open N-mixture
- [`count()`](https://gillescolling.com/tulpaObs/reference/count.md) :
  Count / relative-abundance GLMM family (no detection process)
- [`ms_count()`](https://gillescolling.com/tulpaObs/reference/ms_count.md)
  : Multispecies (community) relative-abundance GLMM family
- [`distance()`](https://gillescolling.com/tulpaObs/reference/distance.md)
  : Binned distance-sampling family
- [`distance_laplace()`](https://gillescolling.com/tulpaObs/reference/distance_laplace.md)
  : Laplace fit of the binned distance-sampling abundance model
- [`ms_distance()`](https://gillescolling.com/tulpaObs/reference/ms_distance.md)
  : Community (multispecies) binned distance-sampling family
- [`distsamp_open()`](https://gillescolling.com/tulpaObs/reference/distsamp_open.md)
  : Open-population distance-sampling family
- [`gdistremoval()`](https://gillescolling.com/tulpaObs/reference/gdistremoval.md)
  : Joint distance + removal sampling family
- [`removal()`](https://gillescolling.com/tulpaObs/reference/removal.md)
  : Removal-sampling family
- [`removal_laplace()`](https://gillescolling.com/tulpaObs/reference/removal_laplace.md)
  : Laplace fit of the removal-sampling abundance model
- [`double_observer()`](https://gillescolling.com/tulpaObs/reference/double_observer.md)
  : Double-observer abundance family
- [`nmix_laplace()`](https://gillescolling.com/tulpaObs/reference/nmix_laplace.md)
  : Laplace fit of the Royle (2004) N-mixture model
- [`nmix_laplace_bym2()`](https://gillescolling.com/tulpaObs/reference/nmix_laplace_bym2.md)
  : BYM2 Royle (2004) N-mixture model via nested Laplace
- [`nmix_laplace_car_proper()`](https://gillescolling.com/tulpaObs/reference/nmix_laplace_car_proper.md)
  : Proper CAR Royle (2004) N-mixture model via nested Laplace
- [`nmix_laplace_icar()`](https://gillescolling.com/tulpaObs/reference/nmix_laplace_icar.md)
  : Spatial Royle (2004) N-mixture model via nested Laplace
- [`nmix_laplace_re()`](https://gillescolling.com/tulpaObs/reference/nmix_laplace_re.md)
  : Community / multispecies N-mixture by Laplace
- [`nmix_site_marginal()`](https://gillescolling.com/tulpaObs/reference/nmix_site_marginal.md)
  : Per-site N-mixture marginal as a composable random-effect callback

## Cover families

- [`cover()`](https://gillescolling.com/tulpaObs/reference/cover.md) :
  Cover hurdle family (vegetation cover, MOTIVATE pattern)
- [`occu_cover()`](https://gillescolling.com/tulpaObs/reference/occu_cover.md)
  : Joint occupancy-detection + cover hurdle family
- [`occu_cover_inputs()`](https://gillescolling.com/tulpaObs/reference/occu_cover_inputs.md)
  : Build occu_cover() model inputs from a long, plot-level frame
- [`ms_occu_cover()`](https://gillescolling.com/tulpaObs/reference/ms_occu_cover.md)
  : Community (multispecies) joint occupancy-detection + cover family
- [`occu_multiscale_cover()`](https://gillescolling.com/tulpaObs/reference/occu_multiscale_cover.md)
  : Three-level (multiscale) occupancy + cover hurdle family
- [`cover_priors()`](https://gillescolling.com/tulpaObs/reference/cover_priors.md)
  : Weakly-informative priors for cover-hurdle Laplace fits
- [`occu_aggregation_scan()`](https://gillescolling.com/tulpaObs/reference/occu_aggregation_scan.md)
  : Suggest spatial and temporal aggregation for an identifiable
  occupancy model
- [`plot(`*`<tobs_aggregation_scan>`*`)`](https://gillescolling.com/tulpaObs/reference/plot.tobs_aggregation_scan.md)
  : Plot an occupancy aggregation scan

## Fit accessors

S3 methods on `tobs_fit`.

- [`print(`*`<tobs_fit>`*`)`](https://gillescolling.com/tulpaObs/reference/print.tobs_fit.md)
  : Print method for tobs_fit
- [`summary(`*`<tobs_fit>`*`)`](https://gillescolling.com/tulpaObs/reference/summary.tobs_fit.md)
  : Summary for tobs_fit, with skewness column when simplified Laplace
  is used
- [`coef(`*`<tobs_fit>`*`)`](https://gillescolling.com/tulpaObs/reference/coef.tobs_fit.md)
  : Coefficients for a tobs_fit
- [`coef(`*`<tobs_batch>`*`)`](https://gillescolling.com/tulpaObs/reference/coef.tobs_batch.md)
  : Per-species coefficients from a batched occu_cover fit.
- [`print(`*`<tobs_batch>`*`)`](https://gillescolling.com/tulpaObs/reference/print.tobs_batch.md)
  : Print a batched multi-response occu_cover fit.
- [`tidy(`*`<tobs_fit>`*`)`](https://gillescolling.com/tulpaObs/reference/tidy.tobs_fit.md)
  : Coefficient table for a tobs_fit
- [`glance(`*`<tobs_fit>`*`)`](https://gillescolling.com/tulpaObs/reference/glance.tobs_fit.md)
  : One-row model summary for a tobs_fit
- [`fitted(`*`<tobs_fit>`*`)`](https://gillescolling.com/tulpaObs/reference/fitted.tobs_fit.md)
  : Fitted values (occupancy and detection probabilities)
- [`residuals(`*`<tobs_fit>`*`)`](https://gillescolling.com/tulpaObs/reference/residuals.tobs_fit.md)
  : Residuals from occupancy model
- [`ranef(`*`<tobs_fit>`*`)`](https://gillescolling.com/tulpaObs/reference/ranef.tobs_fit.md)
  : Random-effect estimates (BLUPs) for a tobs_fit
- [`logLik(`*`<tobs_fit>`*`)`](https://gillescolling.com/tulpaObs/reference/logLik.tobs_fit.md)
  : Log-likelihood of a tobs fit
- [`nobs(`*`<tobs_fit>`*`)`](https://gillescolling.com/tulpaObs/reference/nobs.tobs_fit.md)
  : Number of observations
- [`update(`*`<tobs_fit>`*`)`](https://gillescolling.com/tulpaObs/reference/update.tobs_fit.md)
  : Update and refit a tobs model
- [`predict(`*`<tobs_fit>`*`)`](https://gillescolling.com/tulpaObs/reference/predict.tobs_fit.md)
  : Predict from occupancy model
- [`predict(`*`<cover_fit>`*`)`](https://gillescolling.com/tulpaObs/reference/predict.cover_fit.md)
  : Predict cover from a cover_fit
- [`predict(`*`<occu_categorical_fit>`*`)`](https://gillescolling.com/tulpaObs/reference/predict.occu_categorical_fit.md)
  : Predict from an occu_categorical fit
- [`simulate(`*`<tobs_fit>`*`)`](https://gillescolling.com/tulpaObs/reference/simulate.tobs_fit.md)
  : Simulate replicate datasets from the posterior
- [`` `$`( ``*`<tobs_fit>`*`)`](https://gillescolling.com/tulpaObs/reference/cash-.tobs_fit.md)
  : Access spOccupancy-compatible fields from tobs fits
- [`tobs_predict_spatial()`](https://gillescolling.com/tulpaObs/reference/tobs_predict_spatial.md)
  : Predict the state process at new spatial locations
- [`tobs_marginal_effect()`](https://gillescolling.com/tulpaObs/reference/tobs_marginal_effect.md)
  : Compute marginal effect of a covariate
- [`tobs_richness()`](https://gillescolling.com/tulpaObs/reference/tobs_richness.md)
  : Estimate species richness from community model
- [`tobs_associations()`](https://gillescolling.com/tulpaObs/reference/tobs_associations.md)
  : Residual species-association matrices from a spatial-factor
  community fit
- [`nobs(`*`<tobs_multiarm_fit>`*`)`](https://gillescolling.com/tulpaObs/reference/tobs_multiarm_methods.md)
  [`coef(`*`<tobs_multiarm_fit>`*`)`](https://gillescolling.com/tulpaObs/reference/tobs_multiarm_methods.md)
  [`vcov(`*`<tobs_multiarm_fit>`*`)`](https://gillescolling.com/tulpaObs/reference/tobs_multiarm_methods.md)
  [`confint(`*`<tobs_multiarm_fit>`*`)`](https://gillescolling.com/tulpaObs/reference/tobs_multiarm_methods.md)
  [`logLik(`*`<tobs_multiarm_fit>`*`)`](https://gillescolling.com/tulpaObs/reference/tobs_multiarm_methods.md)
  [`glance(`*`<tobs_multiarm_fit>`*`)`](https://gillescolling.com/tulpaObs/reference/tobs_multiarm_methods.md)
  [`tidy(`*`<tobs_multiarm_fit>`*`)`](https://gillescolling.com/tulpaObs/reference/tobs_multiarm_methods.md)
  [`summary(`*`<tobs_multiarm_fit>`*`)`](https://gillescolling.com/tulpaObs/reference/tobs_multiarm_methods.md)
  [`plot(`*`<tobs_multiarm_fit>`*`)`](https://gillescolling.com/tulpaObs/reference/tobs_multiarm_methods.md)
  [`fitted(`*`<tobs_multiarm_fit>`*`)`](https://gillescolling.com/tulpaObs/reference/tobs_multiarm_methods.md)
  [`residuals(`*`<tobs_multiarm_fit>`*`)`](https://gillescolling.com/tulpaObs/reference/tobs_multiarm_methods.md)
  : S3 methods for multi-arm fits

## Diagnostics and model comparison

- [`check_model(`*`<tobs_fit>`*`)`](https://gillescolling.com/tulpaObs/reference/check_model.tobs_fit.md)
  : Comprehensive model checking
- [`pit_residuals(`*`<tobs_fit>`*`)`](https://gillescolling.com/tulpaObs/reference/pit_residuals.tobs_fit.md)
  : PIT residuals for a tobs fit
- [`ppc()`](https://gillescolling.com/tulpaObs/reference/ppc.md) :
  Posterior predictive check
- [`sbc(`*`<tobs_fit>`*`)`](https://gillescolling.com/tulpaObs/reference/sbc.tobs_fit.md)
  : Simulation-based calibration for a fitted tobs model
- [`test_dispersion(`*`<tobs_fit>`*`)`](https://gillescolling.com/tulpaObs/reference/tobs_gof_tests.md)
  [`test_zero_inflation(`*`<tobs_fit>`*`)`](https://gillescolling.com/tulpaObs/reference/tobs_gof_tests.md)
  [`test_outliers(`*`<tobs_fit>`*`)`](https://gillescolling.com/tulpaObs/reference/tobs_gof_tests.md)
  : Goodness-of-fit tests for a tobs fit
- [`waic(`*`<tobs_fit>`*`)`](https://gillescolling.com/tulpaObs/reference/tobs_criteria.md)
  [`loo(`*`<tobs_fit>`*`)`](https://gillescolling.com/tulpaObs/reference/tobs_criteria.md)
  [`dic(`*`<tobs_fit>`*`)`](https://gillescolling.com/tulpaObs/reference/tobs_criteria.md)
  [`cpo(`*`<tobs_fit>`*`)`](https://gillescolling.com/tulpaObs/reference/tobs_criteria.md)
  [`pointwise_loglik(`*`<tobs_fit>`*`)`](https://gillescolling.com/tulpaObs/reference/tobs_criteria.md)
  : Model criteria for occupancy / cover models
- [`tobs_stack()`](https://gillescolling.com/tulpaObs/reference/tobs_stack.md)
  [`print(`*`<tobs_stack>`*`)`](https://gillescolling.com/tulpaObs/reference/tobs_stack.md)
  [`fitted(`*`<tobs_stack>`*`)`](https://gillescolling.com/tulpaObs/reference/tobs_stack.md)
  [`predict(`*`<tobs_stack>`*`)`](https://gillescolling.com/tulpaObs/reference/tobs_stack.md)
  : Stack fitted models into a LOO-weighted ensemble
- [`convergence()`](https://gillescolling.com/tulpaObs/reference/convergence.md)
  [`converged()`](https://gillescolling.com/tulpaObs/reference/convergence.md)
  : Convergence record for a fitted model

## Simulation

- [`simulate_abun()`](https://gillescolling.com/tulpaObs/reference/simulate_abun.md)
  : Simulate Royle (2004) N-mixture abundance data
- [`simulate_count()`](https://gillescolling.com/tulpaObs/reference/simulate_count.md)
  : Simulate count / relative-abundance GLMM data
- [`simulate_cover()`](https://gillescolling.com/tulpaObs/reference/simulate_cover.md)
  : Simulate cover-hurdle data (lognormal positive part)
- [`simulate_cover_joint()`](https://gillescolling.com/tulpaObs/reference/simulate_cover_joint.md)
  : Simulate cover-hurdle data with a shared BYM2 spatial field
- [`simulate_distance()`](https://gillescolling.com/tulpaObs/reference/simulate_distance.md)
  : Simulate binned distance-sampling abundance data
- [`simulate_distsamp_open()`](https://gillescolling.com/tulpaObs/reference/simulate_distsamp_open.md)
  : Simulate an open-population distance-sampling data set
- [`simulate_double_observer()`](https://gillescolling.com/tulpaObs/reference/simulate_double_observer.md)
  : Simulate a double-observer abundance data set
- [`simulate_dyn_abun()`](https://gillescolling.com/tulpaObs/reference/simulate_dyn_abun.md)
  : Simulate Dail-Madsen open N-mixture data
- [`simulate_dyn_int_occu()`](https://gillescolling.com/tulpaObs/reference/simulate_dyn_int_occu.md)
  : Simulate a multi-season integrated occupancy data set
- [`simulate_dyn_occu()`](https://gillescolling.com/tulpaObs/reference/simulate_dyn_occu.md)
  : Simulate temporal (multi-season) occupancy data
- [`simulate_fp_occu()`](https://gillescolling.com/tulpaObs/reference/simulate_fp_occu.md)
  : Simulate multistate false-positive occupancy data
- [`simulate_gdistremoval()`](https://gillescolling.com/tulpaObs/reference/simulate_gdistremoval.md)
  : Simulate a joint distance + removal data set
- [`simulate_int_occu()`](https://gillescolling.com/tulpaObs/reference/simulate_int_occu.md)
  : Simulate integrated (multi-source) occupancy data
- [`simulate_jsdm()`](https://gillescolling.com/tulpaObs/reference/simulate_jsdm.md)
  : Simulate joint species distribution (presence/absence) data
- [`simulate_ms_abun()`](https://gillescolling.com/tulpaObs/reference/simulate_ms_abun.md)
  : Simulate community (multispecies) N-mixture abundance data
- [`simulate_ms_count()`](https://gillescolling.com/tulpaObs/reference/simulate_ms_count.md)
  : Simulate community relative-abundance (count / continuous) data
- [`simulate_ms_distance()`](https://gillescolling.com/tulpaObs/reference/simulate_ms_distance.md)
  : Simulate community (multispecies) binned distance-sampling data
- [`simulate_ms_dyn_occu()`](https://gillescolling.com/tulpaObs/reference/simulate_ms_dyn_occu.md)
  : Simulate temporal multi-species occupancy data
- [`simulate_ms_int_occu()`](https://gillescolling.com/tulpaObs/reference/simulate_ms_int_occu.md)
  : Simulate integrated multi-species occupancy data
- [`simulate_ms_occu()`](https://gillescolling.com/tulpaObs/reference/simulate_ms_occu.md)
  : Simulate multi-species occupancy data
- [`simulate_ms_occu_cover()`](https://gillescolling.com/tulpaObs/reference/simulate_ms_occu_cover.md)
  : Simulate community (multispecies) joint occupancy-cover data
- [`simulate_ms_occu_cover_spatial()`](https://gillescolling.com/tulpaObs/reference/simulate_ms_occu_cover_spatial.md)
  : Simulate a reduced-rank spatial-factor community occu_cover data set
  (K = 1)
- [`simulate_occu()`](https://gillescolling.com/tulpaObs/reference/simulate_occu.md)
  : Simulate single-species occupancy data
- [`simulate_occu_categorical()`](https://gillescolling.com/tulpaObs/reference/simulate_occu_categorical.md)
  : Simulate from an occu_categorical model
- [`simulate_occu_cover()`](https://gillescolling.com/tulpaObs/reference/simulate_occu_cover.md)
  : Simulate joint occupancy-detection + cover data
- [`simulate_occu_multi()`](https://gillescolling.com/tulpaObs/reference/simulate_occu_multi.md)
  : Simulate a multi-species co-occurrence occupancy data set
- [`simulate_occu_multiscale_cover()`](https://gillescolling.com/tulpaObs/reference/simulate_occu_multiscale_cover.md)
  : Simulate a three-level occupancy + cover hurdle data set
- [`simulate_occu_ttd()`](https://gillescolling.com/tulpaObs/reference/simulate_occu_ttd.md)
  : Simulate a time-to-detection occupancy data set
- [`simulate_removal()`](https://gillescolling.com/tulpaObs/reference/simulate_removal.md)
  : Simulate removal-sampling abundance data
- [`simulate_royle_nichols()`](https://gillescolling.com/tulpaObs/reference/simulate_royle_nichols.md)
  : Simulate Royle-Nichols occupancy data
- [`simulate_t_occu()`](https://gillescolling.com/tulpaObs/reference/simulate_t_occu.md)
  : Simulate multi-season occupancy with an AR1 year effect (tPGOcc)

## Data and utilities

- [`foray_counts`](https://gillescolling.com/tulpaObs/reference/foray_counts.md)
  : Repeated point-count abundance survey (synthetic)
- [`meadow_cover`](https://gillescolling.com/tulpaObs/reference/meadow_cover.md)
  : Grassland vegetation-cover panel (synthetic)
- [`peatland_occu`](https://gillescolling.com/tulpaObs/reference/peatland_occu.md)
  : Peatland occupancy survey (synthetic)
- [`summary(`*`<tobs_data>`*`)`](https://gillescolling.com/tulpaObs/reference/summary.tobs_data.md)
  : Summarise occupancy data
- [`plot(`*`<tobs_data>`*`)`](https://gillescolling.com/tulpaObs/reference/plot.tobs_data.md)
  : Plot detection history patterns
- [`within_between()`](https://gillescolling.com/tulpaObs/reference/within_between.md)
  : Within / between (Mundlak) decomposition for longitudinal data
