# Re-exports from tulpaMesh. The SPDE spatial term (`spde()`) builds its Matern
# field on a triangular mesh; `fem_matrices()` is the mesh-assembly entry point,
# re-exported here so the SPDE workflow is reachable without attaching tulpaMesh.

#' @importFrom tulpaMesh fem_matrices
#' @export
tulpaMesh::fem_matrices
