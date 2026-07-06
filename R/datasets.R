# Documentation for the bundled example datasets. The data objects are built by
# data-raw/make_datasets.R (synthetic, fixed-seed generative models with a known
# `truth`); this file only documents them.

#' Peatland occupancy survey (synthetic)
#'
#' A single-season occupancy detection-history dataset for 120 wetland sites
#' surveyed on 4 visits, drawn from a known generative model (MacKenzie et al.
#' 2002). Occupancy increases with site `wetness` and decreases with
#' `elevation`; per-visit detection increases with survey `effort`. Twelve
#' detection cells are set to `NA` to mimic incomplete histories. Synthetic, not
#' field data; the generative coefficients are in `truth`.
#'
#' @format A list with components:
#' \describe{
#'   \item{y}{120 x 4 integer matrix of detections (0/1/`NA`).}
#'   \item{occ.covs}{Data frame of site covariates: `elevation`, `wetness`
#'     (both standardised).}
#'   \item{det.covs}{Named list with one 120 x 4 visit covariate matrix,
#'     `effort` (standardised).}
#'   \item{coords}{120 x 2 matrix of `x` / `y` coordinates in the unit square.}
#'   \item{truth}{List of generative quantities: `beta_psi`, `beta_p`, the
#'     realised latent occupancy `z`, and per-site `psi`.}
#' }
#' @seealso [occu()], [tobs()], [simulate_occu()]
#' @examples
#' data(peatland_occu)
#' str(peatland_occu$y)
#' \donttest{
#' fit <- tobs(~ elevation + wetness, data = peatland_occu$occ.covs,
#'             family = occu(), detection = ~ effort,
#'             y = peatland_occu$y, visits = peatland_occu$det.covs)
#' coef(fit)
#' }
"peatland_occu"

#' Repeated point-count abundance survey (synthetic)
#'
#' A binomial N-mixture dataset (Royle 2004) of repeated counts at 100 sites
#' over 3 visits, drawn from a known generative model. Latent abundance
#' increases with `shrub` cover and decreases with `elevation`; per-visit
#' detection increases with survey `effort`. Synthetic, not field data; the
#' generative coefficients and the realised latent abundances are in `truth`.
#'
#' @format A list with components:
#' \describe{
#'   \item{y}{100 x 3 integer matrix of counts.}
#'   \item{occ.covs}{Data frame of site covariates: `elevation`, `shrub` (both
#'     standardised).}
#'   \item{det.covs}{Named list with one 100 x 3 visit covariate matrix,
#'     `effort` (standardised).}
#'   \item{truth}{List of generative quantities: `beta_lambda`, `beta_p`, the
#'     realised latent abundance `N`, and per-site `lambda`.}
#' }
#' @seealso [abun()], [tobs()], [simulate_abun()]
#' @examples
#' data(foray_counts)
#' summary(as.vector(foray_counts$y))
#' \donttest{
#' fit <- tobs(~ elevation + shrub, data = foray_counts$occ.covs,
#'             family = abun(), detection = ~ effort,
#'             y = foray_counts$y, visits = foray_counts$det.covs)
#' coef(fit)
#' }
"foray_counts"

#' Grassland vegetation-cover panel (synthetic)
#'
#' A vegetation cover-hurdle dataset for 150 grassland plots surveyed across
#' three years, drawn from a known generative model: the species occurs with a
#' probability driven by `moisture`, and where it occurs, its conditional cover
#' (a Beta-distributed proportion in `(0, 1)`) increases with `moisture` and
#' decreases with `grazing`. Roughly half the plots record zero cover. `year_c`
#' is the centred survey year, so [within_between()] can split the cross-plot
#' baseline from the within-plot temporal deviation. Synthetic, not field data;
#' the generative coefficients are attached as `attr(meadow_cover, "truth")`.
#'
#' @format A data frame with 150 rows (one per plot) and columns:
#' \describe{
#'   \item{plot}{Plot identifier (factor).}
#'   \item{year}{Survey year (2018-2020).}
#'   \item{year_c}{Survey year centred on its mean.}
#'   \item{moisture}{Standardised soil-moisture index.}
#'   \item{grazing}{Standardised grazing-pressure index.}
#'   \item{cover}{Recorded cover as a proportion in `[0, 1]` (0 where absent).}
#' }
#' @seealso [cover()], [tobs()], [within_between()], [simulate_cover()]
#' @examples
#' data(meadow_cover)
#' mean(meadow_cover$cover > 0)
#' \donttest{
#' fit <- tobs(cover ~ moisture, data = meadow_cover,
#'             family = cover(response = "beta"))
#' coef(fit)
#' }
"meadow_cover"
