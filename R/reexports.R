# Re-exports from tulpaMesh. The SPDE spatial term (`spde()`) builds its Matern
# field on a triangular mesh; `fem_matrices()` is the mesh-assembly entry point,
# re-exported here so the SPDE workflow is reachable without attaching tulpaMesh.

#' @importFrom tulpaMesh fem_matrices
#' @export
tulpaMesh::fem_matrices


# Diagnostic generics. Each name belongs to whichever package owns the concept:
# the engine's own checks come from tulpa, and WAIC / PSIS-LOO from loo, whose
# generics the wider ecosystem (loo_compare, model weights) dispatches on.
# Re-exporting them here means attaching tulpaObs is enough to reach every
# diagnostic door, and defining none of them means no masking for a session
# that also attaches tulpa or loo.

#' @importFrom tulpa sbc
#' @export
tulpa::sbc

#' @importFrom tulpa pit_residuals
#' @export
tulpa::pit_residuals

#' @importFrom tulpa test_uniformity
#' @export
tulpa::test_uniformity

#' @importFrom tulpa test_dispersion
#' @export
tulpa::test_dispersion

#' @importFrom tulpa test_outliers
#' @export
tulpa::test_outliers

#' @importFrom tulpa test_zero_inflation
#' @export
tulpa::test_zero_inflation

#' @importFrom tulpa check_model
#' @export
tulpa::check_model

#' @importFrom tulpa dic
#' @export
tulpa::dic

#' @importFrom tulpa cpo
#' @export
tulpa::cpo

#' @importFrom loo waic
#' @export
loo::waic

#' @importFrom loo loo
#' @export
loo::loo
