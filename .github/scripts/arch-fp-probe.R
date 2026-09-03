#!/usr/bin/env Rscript

# Measures how far apart two independent implementations of the same Euclidean
# distance land on the running architecture.
#
# compute_nngp_neighbors() fills the pairwise neighbour block with
# stats::dist(), a compiled C loop (stats/src/distance.c) accumulating
# `dist += dev * dev`. test-nngp-neighbors.R rebuilds the same block from the
# explicit two-column formula, evaluated by the R interpreter. The two are
# separate implementations in separate languages and nothing makes them agree
# to the bit -- on x86_64 they happen to, because the baseline instruction set
# has no fused multiply-add for the compiler to contract that accumulation
# into, and on aarch64 FMA is baseline and gcc contracts it by default.
#
# Run on both architectures. The gap is reported in ULPs of the value itself,
# which is the scale a tolerance for it has to be written in.

fixture <- function() {
  set.seed(11L)
  n <- 30L
  co <- cbind(stats::runif(n), stats::runif(n))
  ord <- order(co[, 1], co[, 2])
  co[ord, , drop = FALSE]
}

explicit_pairs <- function(nb) {
  m <- nrow(nb)
  outer(seq_len(m), seq_len(m), function(a, b) {
    sqrt((nb[a, 1] - nb[b, 1])^2 + (nb[a, 2] - nb[b, 2])^2)
  })
}

cs <- fixture()
k <- 5L

worst_abs <- 0
worst_ulp <- 0
n_exact <- 0L
n_seen <- 0L

for (i in c(2L, 6L, 20L, 30L)) {
  m <- min(i - 1L, k)
  if (m <= 1L) next
  d <- sqrt((cs[seq_len(i - 1L), 1] - cs[i, 1])^2 +
            (cs[seq_len(i - 1L), 2] - cs[i, 2])^2)
  o <- order(d)[seq_len(m)]
  nb <- cs[o, , drop = FALSE]

  from_dist <- unname(as.matrix(stats::dist(nb)))
  from_formula <- explicit_pairs(nb)

  gap <- abs(from_dist - from_formula)
  # ULP of the value, not of 1: these distances are O(0.1), so an absolute
  # tolerance written against 1 would be far looser than it reads.
  scale <- pmax(abs(from_formula), .Machine$double.xmin)
  ulp <- gap / (scale * .Machine$double.eps)

  n_seen <- n_seen + 1L
  if (identical(from_dist, from_formula)) n_exact <- n_exact + 1L
  worst_abs <- max(worst_abs, max(gap))
  worst_ulp <- max(worst_ulp, max(ulp))

  cat(sprintf("i = %2d  m = %d  bit-identical = %-5s  max|gap| = %.3e  max ULP = %.2f\n",
              i, m, identical(from_dist, from_formula), max(gap), max(ulp)))
}

cat("\n")
cat(sprintf("arch                : %s\n", R.version$arch))
cat(sprintf("platform            : %s\n", R.version$platform))
cat(sprintf("blocks compared     : %d\n", n_seen))
cat(sprintf("bit-identical blocks: %d\n", n_exact))
cat(sprintf("worst absolute gap  : %.6e\n", worst_abs))
cat(sprintf("worst gap in ULPs   : %.4f\n", worst_ulp))
