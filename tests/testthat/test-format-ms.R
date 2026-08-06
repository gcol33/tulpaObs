# =============================================================================
# test-format-ms.R - tobs_format_ms(), the multi-species data constructor
# (gcol33/tulpaObs#179).
#
# The contract the community binders rely on is that the two entry points --
# a [sites x visits x species] array and a named list of [sites x visits]
# matrices -- describe the same data, so the same data given either way must
# produce the same object. The species-major fill in the list branch is where
# that can break silently: a transposed or mis-strided write still returns an
# array of the right dimension holding the right multiset of values, so a shape
# check passes while species 1's detections have become species 2's.
#
# Every assertion below therefore names the expected VALUES (the exact array,
# the exact species labels, the exact printed line), not the shape.
# =============================================================================

.fms_fixture <- function(seed = 1L, n_sites = 5L, n_visits = 3L,
                         n_species = 3L) {
  set.seed(seed)
  arr <- array(stats::rbinom(n_sites * n_visits * n_species, 1L, 0.4),
               dim = c(n_sites, n_visits, n_species))
  list(arr = arr,
       lst = stats::setNames(
         lapply(seq_len(n_species), function(s) arr[, , s]),
         paste0("sp_", letters[seq_len(n_species)])))
}


test_that("an array and a named list of matrices give the same object", {
  f <- .fms_fixture()
  nms <- names(f$lst)

  from_array <- tobs_format_ms(f$arr, species_names = nms)
  from_list  <- tobs_format_ms(f$lst)

  expect_identical(from_list, from_array)

  # And the values are the input's, not merely equal in aggregate: a
  # species-major fill that strided wrongly would still produce an array of the
  # right dimension holding the same multiset of 0/1.
  expect_identical(from_list$y, f$arr)
  for (s in seq_along(nms)) {
    expect_identical(from_list$y[, , s], f$arr[, , s])
  }
})


test_that("species labels come from the list names, the argument, or a default", {
  f <- .fms_fixture()
  nms <- names(f$lst)

  expect_identical(tobs_format_ms(f$lst)$species_names, nms)

  # An explicit argument wins over the list names.
  expect_identical(
    tobs_format_ms(f$lst, species_names = c("A", "B", "C"))$species_names,
    c("A", "B", "C"))

  # No names anywhere: positional defaults, in order.
  expect_identical(tobs_format_ms(unname(f$lst))$species_names,
                   c("sp1", "sp2", "sp3"))
  expect_identical(tobs_format_ms(f$arr)$species_names,
                   c("sp1", "sp2", "sp3"))

  # The labels stay aligned with the slices they name.
  d <- tobs_format_ms(f$lst)
  expect_identical(d$y[, , match("sp_b", d$species_names)], f$arr[, , 2])
})


test_that("dimensions and the species count are carried through", {
  f <- .fms_fixture(n_sites = 7L, n_visits = 4L, n_species = 2L)

  for (d in list(tobs_format_ms(f$arr), tobs_format_ms(f$lst))) {
    expect_s3_class(d, "tobs_data")
    expect_identical(dim(d$y), c(7L, 4L, 2L))
    expect_identical(d$n_species, 2L)
    expect_length(d$species_names, 2L)
  }

  # One species is a degenerate case the binders still have to handle: the
  # third dimension must survive rather than dropping to a matrix.
  one <- tobs_format_ms(list(only = f$arr[, , 1]))
  expect_identical(dim(one$y), c(7L, 4L, 1L))
  expect_identical(one$n_species, 1L)
  expect_identical(one$species_names, "only")
  expect_identical(one$y[, , 1], f$arr[, , 1])
})


test_that("missing visits survive both entry points as NA", {
  # An unvisited site-visit is NA, and it must stay NA rather than becoming a
  # fabricated zero -- a zero there is a non-detection the model would score.
  f <- .fms_fixture()
  arr <- f$arr
  arr[1, 2, 1] <- NA
  arr[4, 3, 3] <- NA
  lst <- stats::setNames(lapply(seq_len(3), function(s) arr[, , s]),
                         names(f$lst))

  for (d in list(tobs_format_ms(arr, species_names = names(f$lst)),
                 tobs_format_ms(lst))) {
    expect_true(is.na(d$y[1, 2, 1]))
    expect_true(is.na(d$y[4, 3, 3]))
    expect_equal(sum(is.na(d$y)), 2L)
    expect_identical(d$y, arr)
  }
})


test_that("covariates and coordinates are carried unchanged", {
  f <- .fms_fixture()
  occ <- data.frame(elev = c(1.5, -2, 0, 3, 0.25), grp = letters[1:5])
  det <- list(effort = matrix(seq_len(15) / 10, 5, 3))
  xy  <- cbind(lon = seq(0, 1, length.out = 5), lat = seq(1, 2, length.out = 5))

  d <- tobs_format_ms(f$lst, occ.covs = occ, det.covs = det, coords = xy)
  expect_identical(d$occ.covs, occ)
  expect_identical(d$det.covs, det)
  expect_identical(d$coords, xy)

  # A plain list of site-level covariates is coerced to a data.frame, keeping
  # its names and values.
  dl <- tobs_format_ms(f$lst, occ.covs = list(elev = occ$elev, grp = occ$grp))
  expect_s3_class(dl$occ.covs, "data.frame")
  expect_identical(names(dl$occ.covs), c("elev", "grp"))
  expect_identical(dl$occ.covs$elev, occ$elev)

  # Absent slots stay absent rather than becoming empty structures.
  bare <- tobs_format_ms(f$lst)
  expect_null(bare$occ.covs)
  expect_null(bare$det.covs)
  expect_null(bare$coords)
})


test_that("the object prints its dimensions and its populated slots", {
  f <- .fms_fixture()
  occ <- data.frame(elev = stats::rnorm(5))
  det <- list(effort = matrix(0, 5, 3))
  xy  <- cbind(0, seq_len(5))

  expect_output(print(tobs_format_ms(f$lst)), "tobs_data: 5 sites, 3 visits")

  out <- utils::capture.output(
    print(tobs_format_ms(f$lst, occ.covs = occ, det.covs = det, coords = xy)))
  expect_true(any(grepl("Occupancy covariates: elev", out, fixed = TRUE)))
  expect_true(any(grepl("Detection covariates: effort", out, fixed = TRUE)))
  expect_true(any(grepl("Coordinates: yes", out, fixed = TRUE)))

  # The covariate lines are reported only when the slot is populated.
  bare <- utils::capture.output(print(tobs_format_ms(f$lst)))
  expect_false(any(grepl("covariates", bare)))
  expect_false(any(grepl("Coordinates", bare)))
})


test_that("a list with no species is refused", {
  # Nothing downstream can consume a zero-species object, so this must stop
  # rather than build one.
  expect_error(tobs_format_ms(list()), "empty list")
})


test_that("a species off the shared sites x visits grid is refused", {
  # The array slice is filled by assignment, so a matrix whose length divides
  # the first's is RECYCLED rather than rejected: before this guard, species 2
  # below produced a full 6 x 4 slice built by tiling its 12 values twice, and
  # nothing downstream could tell. Every shape below is a real way to get this
  # wrong, and each must stop.
  ok <- matrix(rbinom(24, 1L, 0.5), 6L, 4L)

  half <- matrix(rbinom(12, 1L, 0.5), 6L, 2L)      # length divides 24 exactly
  expect_error(tobs_format_ms(list(a = ok, b = half)), "must share one")

  transposed <- t(ok)                               # 4 x 6, same length
  expect_error(tobs_format_ms(list(a = ok, b = transposed)), "must share one")

  ragged <- matrix(rbinom(20, 1L, 0.5), 5L, 4L)     # one site short
  expect_error(tobs_format_ms(list(a = ok, b = ragged)), "must share one")

  expect_error(tobs_format_ms(list(a = ok, b = as.vector(ok))), "not a matrix")

  # The species is named in the message when the list is named, so the caller
  # is told which one to look at.
  expect_error(tobs_format_ms(list(setophaga = ok, dendroica = half)),
               "dendroica")

  # A matrix on the shared grid still goes through untouched.
  expect_equal(dim(tobs_format_ms(list(a = ok, b = ok))$y), c(6L, 4L, 2L))
})


test_that("a non-integer response is refused rather than truncated", {
  # The slice is integer storage, so 2.7 would land as 2 with no warning.
  ok   <- matrix(rbinom(24, 1L, 0.5), 6L, 4L)
  frac <- matrix(as.numeric(ok), 6L, 4L); frac[1L, 1L] <- 2.7
  expect_error(tobs_format_ms(list(a = ok, b = frac)), "non-integer")

  # Whole-number doubles are lossless, so they are accepted and stored as
  # integer -- and NA is not mistaken for a fractional value.
  whole <- matrix(as.numeric(ok), 6L, 4L); whole[2L, 2L] <- NA_real_
  res <- tobs_format_ms(list(a = ok, b = whole))
  expect_type(res$y, "integer")
  expect_true(is.na(res$y[2L, 2L, 2L]))
})
