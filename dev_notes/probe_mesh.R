## Probe the tulpaMesh::tulpa_mesh structure
suppressMessages({
  library(tulpaMesh)
})

set.seed(42)
coords <- cbind(runif(400), runif(400))

cat("--- variant A: raw coords, no max_edge ---\n")
m1 <- tulpa_mesh(coords = coords)
cat("  vertices:", m1$n_vertices, "  triangles:", m1$n_triangles, "\n")

cat("\n--- variant B: max_edge = 0.3 (scalar) ---\n")
m2 <- tulpa_mesh(coords = coords, max_edge = 0.3)
cat("  vertices:", m2$n_vertices, "  triangles:", m2$n_triangles, "\n")

cat("\n--- variant C: max_edge = c(0.3, 0.6) (vector) ---\n")
m3 <- tulpa_mesh(coords = coords, max_edge = c(0.3, 0.6))
cat("  vertices:", m3$n_vertices, "  triangles:", m3$n_triangles, "\n")

cat("\n--- variant D: max_edge = c(0.3, 0.6) + min_angle = 20 (Ruppert) ---\n")
m4 <- tulpa_mesh(coords = coords, max_edge = c(0.3, 0.6), min_angle = 20)
cat("  vertices:", m4$n_vertices, "  triangles:", m4$n_triangles, "\n")

mesh <- m3  # the failing one

cat("class:", class(mesh), "\n")
cat("names:", paste(names(mesh), collapse = ", "), "\n\n")

for (nm in names(mesh)) {
  obj <- mesh[[nm]]
  cat(sprintf("[%s]: class = %s, ", nm, paste(class(obj), collapse="/")))
  if (is.matrix(obj) || is.data.frame(obj)) {
    cat(sprintf("dim = %dx%d\n", nrow(obj), ncol(obj)))
  } else if (is.numeric(obj) || is.integer(obj)) {
    cat(sprintf("len = %d, first = %s\n", length(obj),
                paste(head(obj, 3), collapse=",")))
  } else {
    cat("\n")
  }
}

cat("\n--- fem_matrices probe ---\n")
fem <- fem_matrices(mesh, obs_coords = coords)
cat("fem names:", paste(names(fem), collapse = ", "), "\n")
cat("fem$C class:", class(fem$C), "\n")
cat("fem$C dim:", dim(fem$C), "\n")
cat("fem$C nnz:", length(fem$C@x), "\n")
cat("fem$C@x summary:\n"); print(summary(fem$C@x))
cat("fem$G@x summary:\n"); print(summary(fem$G@x))
cat("fem$n_mesh:", fem$n_mesh, "\n")
