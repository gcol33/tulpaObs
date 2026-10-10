# =============================================================================
# count_routes.R - the route table of the count / observation families that
# `.tobs_fit_model()` serves (abun, removal, distance, dyn_abun, fp_occu).
#
# One row per `model_type`. A family is one row: its fitter on each route it
# supports, the Laplace iteration budget of each Laplace route, and the
# parameters of the gates the router applies before it picks a route.
# Everything not named in a row -- route selection, argument forwarding, the
# sampler knobs, the natural-scale tail -- is the router's.
#
# Routes:
#
#   plain       no structured term, no random effect; Laplace
#   re          grouped random effect; Laplace + AGHQ. A family without one
#               fits its random effect on `plain`
#   field       areal term, optionally composed with temporal() / svc();
#               nested Laplace
#   field_bfgs  temporal() / svc() block with no areal term, on the shared
#               areal-BFGS driver. A family without one runs it on `field`
#   nuts        no structured term; NUTS (a single-intercept RE rides along)
#   nuts_field  structured term; NUTS with a fixed-hyper non-centred field
#
# Route entry:
#
#   fit       the fitter's name, resolved in the namespace at call time
#   max.iter  Laplace iteration cap of this route (Laplace routes only). These
#   tol       are per-ROUTE values, not a sampler profile: the Newton / EM
#             routes converge at 300 / 1e-4, the BFGS routes over an exact
#             marginal at their own tolerance. A caller's control$max.iter /
#             control$tol replaces them on every route
#   args      fitter formal -> argument-pool name, for a fitter that spells an
#             argument differently
#   fixed     fitter formal -> value pinned on this route
#
# The router hands each fitter the arguments of one pool (the structure specs,
# `mixture`, `K_max`, `n_quad`, `lkj_eta`, `integration`, `verbose`,
# `sigma.beta`, `sigma.logr`, the sampler knobs) restricted to the formals the
# fitter declares, so a fitter receives exactly what it reads.
#
# Gate parameters:
#
#   family        family name in user-facing messages
#   structured    structure slots that send a fit to a field route
#   field_terms   slots counted as a field by `.tobs_check_areal_engine()`;
#                 NULL skips that gate
#   temporal_arm  arm named by `.tobs_check_count_temporal()`; NULL skips it
#   select        name of a family hook `(model, structs, method)` run first.
#                 It raises a gate only that family has, and returns NULL to
#                 leave the route to the router or a route name to take it
# =============================================================================

.TOBS_COUNT_ROUTES <- list(

  nmix = list(
    family = "abun", structured = "spatial", field_terms = "spatial",
    temporal_arm = NULL, select = ".tobs_abun_route",
    routes = list(
      plain = list(fit = ".tobs_fit_nmix", max.iter = 300L, tol = 1e-4,
                   args = c(sigma_beta = "sigma.beta"),
                   fixed = list(method = "laplace")),
      field = list(fit = ".tobs_fit_nmix", max.iter = 300L, tol = 1e-4,
                   args = c(sigma_beta = "sigma.beta"),
                   fixed = list(method = "nested_laplace")),
      nuts = list(fit = ".tobs_fit_abun_nuts"),
      nuts_field = list(fit = ".tobs_fit_abun_nuts_spatial"))),

  removal = list(
    family = "removal", structured = c("spatial", "temporal", "svc"),
    field_terms = c("spatial", "temporal"), temporal_arm = "abundance",
    routes = list(
      plain = list(fit = ".tobs_fit_removal", max.iter = 300L, tol = 1e-4),
      re = list(fit = ".tobs_fit_removal_re", max.iter = 300L, tol = 1e-4,
                args = c(theta_prior_sd = "sigma.beta")),
      field = list(fit = ".tobs_fit_removal_spatial", max.iter = 300L,
                   tol = 1e-4),
      field_bfgs = list(fit = ".tobs_fit_removal_spatial_bfgs",
                        max.iter = 300L, tol = 1e-4),
      nuts = list(fit = ".tobs_fit_removal_nuts"),
      nuts_field = list(fit = ".tobs_fit_removal_nuts_spatial"))),

  distance = list(
    family = "distance", structured = c("spatial", "temporal", "svc"),
    field_terms = c("spatial", "temporal"), temporal_arm = "abundance",
    routes = list(
      plain = list(fit = ".tobs_fit_distance", max.iter = 300L, tol = 1e-4),
      re = list(fit = ".tobs_fit_distance_re", max.iter = 300L, tol = 1e-4,
                args = c(theta_prior_sd = "sigma.beta")),
      field = list(fit = ".tobs_fit_distance_spatial", max.iter = 300L,
                   tol = 1e-4),
      nuts = list(fit = ".tobs_fit_distance_nuts"),
      nuts_field = list(fit = ".tobs_fit_distance_nuts_spatial"))),

  dyn_abun = list(
    family = "dyn_abun", structured = c("spatial", "temporal", "svc"),
    field_terms = c("spatial", "temporal"), temporal_arm = "initial-abundance",
    select = ".tobs_dyn_abun_route",
    routes = list(
      plain = list(fit = ".tobs_fit_dyn_abun", max.iter = 300L, tol = 1e-8),
      zip = list(fit = ".tobs_fit_dyn_abun_zip", max.iter = 300L, tol = 1e-8),
      re = list(fit = ".tobs_fit_dyn_abun_re", max.iter = 300L, tol = 1e-8,
                args = c(theta_prior_sd = "sigma.beta")),
      field = list(fit = ".tobs_fit_dyn_abun_spatial", max.iter = 300L,
                   tol = 1e-8),
      nuts = list(fit = ".tobs_fit_dyn_abun_nuts"),
      nuts_field = list(fit = ".tobs_fit_dyn_abun_nuts_spatial"))),

  fp_occu = list(
    family = "fp_occu", structured = c("spatial", "temporal", "svc"),
    field_terms = c("spatial", "temporal"), temporal_arm = "occupancy",
    routes = list(
      plain = list(fit = ".tobs_fit_fp_occu", max.iter = 500L, tol = 1e-8,
                   fixed = list(sigma.beta = NULL)),
      re = list(fit = ".tobs_fit_fp_occu_re", max.iter = 300L, tol = 1e-4),
      field = list(fit = ".tobs_fit_fp_occu_spatial", max.iter = 300L,
                   tol = 1e-8),
      nuts = list(fit = ".tobs_fit_fp_occu_nuts"),
      nuts_field = list(fit = ".tobs_fit_fp_occu_nuts_spatial")))
)


# Pick the route of one count-family fit. The family hook runs first; then the
# shared gates, each parameterised by the row, which hold on a route the hook
# chose as well; then the hook's route, else the route that follows from
# whether a structured term is present, the engine, and a random effect.
.tobs_count_route <- function(row, model, structs, method) {
  hook_route <- if (!is.null(row$select)) {
    get(row$select, envir = topenv(), mode = "function")(model, structs, method)
  }
  present <- function(slots) {
    any(!vapply(structs[slots], is.null, logical(1)))
  }
  if (!is.null(structs$temporal) && !is.null(row$temporal_arm)) {
    .tobs_check_count_temporal(structs$temporal, structs$spatial, method,
                               row$family, row$temporal_arm,
                               allow_temporal_only = TRUE,
                               allow_nuts_temporal = TRUE)
  }
  if (!is.null(row$field_terms)) {
    .tobs_check_areal_engine(method, has_field = present(row$field_terms),
                             family = row$family,
                             has_svc = !is.null(structs$svc))
  }
  if (!is.null(hook_route)) return(hook_route)
  routes <- names(row$routes)
  nuts <- identical(method, "nuts")
  if (present(row$structured)) {
    if (nuts) "nuts_field"
    else if (is.null(structs$spatial) && "field_bfgs" %in% routes) "field_bfgs"
    else "field"
  } else if (nuts) {
    "nuts"
  } else if (!is.null(structs$re) && "re" %in% routes) {
    "re"
  } else {
    "plain"
  }
}


# Call one route's fitter with the argument pool restricted to its formals.
# `max_iter` / `tol` take the caller's request, else the route's own budget.
.tobs_count_route_call <- function(row, route, pool, max_iter_req, tol_req) {
  spec <- row$routes[[route]]
  if (is.null(spec)) {
    stop(sprintf("No \"%s\" route is registered for %s() in .TOBS_COUNT_ROUTES.",
                 route, row$family), call. = FALSE)
  }
  if (!is.null(spec$max.iter)) pool$max_iter <- max_iter_req %||% spec$max.iter
  if (!is.null(spec$tol))      pool$tol      <- tol_req %||% spec$tol
  for (formal in names(spec$args)) pool[formal] <- pool[spec$args[[formal]]]
  pool[names(spec$fixed)] <- spec$fixed
  f <- get(spec$fit, envir = topenv(), mode = "function")
  do.call(f, pool[intersect(names(pool), names(formals(f)))])
}
