# ---------------------------------------------------------------------------
# share(residual = ): a shared field plus an arm-specific deviation
#
# `share(spatial())` puts the occurrence field on the cover arm as `alpha * w`:
# one amplitude, the same pattern. That is cheap -- the latent problem stays K,
# the node count -- and it is the whole model whenever the two arms really do
# respond to one surface. Where they do not, the copy has no way to bend, and
# the mismatch lands wherever it can: the amplitude drifts toward whichever arm
# carries more information, and a cover-arm pattern the occurrence field does
# not have is simply not fitted.
#
# The alternative already in the package is a free field on the cover arm -- a
# `spatial(~ 1 || cell, graph = adj)` term placed in the positive formula, the
# arm being chosen by placement. It fits any pattern and costs a second full
# field: the latent problem is 2K, and the joint Laplace
# solve degrades worse than 2x in K. It also leaves the amplitude weakly
# identified -- the free field can reproduce `alpha * w` for any alpha, which is
# what #110 measured for the static intercept field.
#
# `residual =` is the middle: one shared field plus a deviation,
#
#     w_occ(s)   = w(s)
#     w_cover(s) = alpha * w(s) + delta(s)
#
# with `delta` orthogonal to the shared component, so the two halves cannot
# trade against each other and the amplitude means something. Written
# `residual = r`, the deviation is carried on r basis functions and the latent
# problem is K + r; written `residual = "full"` it is carried on the node set
# itself, one latent per node, which spans the same surfaces as the free
# arm-specific field and differs from it only in coordinates -- the same fit,
# reported as (amplitude, deviation) rather than as one free field.
#
# Two doors, one meaning: a field PLACED in the positive formula is an
# INDEPENDENT field on that arm and keeps its current, non-orthogonalized
# reading; `residual =` says the arm's field is a deviation FROM the shared one.
# A fit writes one or the other.
# ---------------------------------------------------------------------------


# Pull the residual request off the positive arm's share() specs.
#
# The deviation belongs to the cover arm as a whole, not to one field block, so
# at most one share() in a fit may name it -- two would be two answers to one
# question, and silently taking the last is how a fit ends up carrying a
# structure nobody wrote. Returns NULL when no share() asks for one.
.occu_cover_residual_spec <- function(copies) {
  asked <- Filter(function(cp) !is.null(cp$residual), copies)
  if (length(asked) == 0L) return(NULL)
  if (length(asked) > 1L) {
    stop("occu_cover(): residual = is the cover arm's deviation from the ",
         "shared field, so one share() names it; got ", length(asked),
         ". Write it on a single share().", call. = FALSE)
  }
  cp <- asked[[1L]]
  # A copy pinned at alpha = 0 carries no shared component, so there is nothing
  # for its deviation to deviate FROM: that model is an independent field on the
  # arm, which has its own spelling.
  if (!is.null(cp$alpha_grid) && !isTRUE(cp$alpha_integrate) &&
      length(cp$alpha_grid) == 1L && cp$alpha_grid == 0) {
    stop("share(alpha = 0, residual = ): an amplitude pinned at 0 shares ",
         "nothing, so there is no shared component to deviate from. Write an ",
         "independent field on the arm instead, by placing it in the positive ",
         "formula: positive = ~ ... + spatial(~ 1 || <node>, graph = adj).",
         call. = FALSE)
  }
  c(cp$residual, list(id = cp$id))
}


# Gate the residual against the fit it was written in. Every rejection here is a
# configuration the arm-specific block itself does not carry, so the residual
# inherits the arm-specific field's scope rather than declaring its own.
.occu_cover_residual_check <- function(res, spatial_info, engine) {
  if (is.null(res)) return(invisible(NULL))
  if (!identical(engine, "nested_laplace")) {
    stop("share(residual = ) rides the joint nested-Laplace engine (the ",
         "deviation is a latent block beside the shared field); got method = ",
         "\"", engine, "\". Use method = \"nested_laplace\".", call. = FALSE)
  }
  if (!is.null(spatial_info$armspec[["pos"]])) {
    stop("occu_cover(): the cover arm already carries a field of its own (a ",
         "spatial() term placed in the positive formula), and ",
         "share(residual = ) is the other way of writing one. Keep the ",
         "independent field, or drop it and write the deviation on the ",
         "share().", call. = FALSE)
  }
  if (isTRUE(spatial_info$correlated)) {
    stop("occu_cover(): share(residual = ) does not compose with a correlated ",
         "spatial bar (`|`, free-Sigma MCAR), which already spans the whole ",
         "coupled field structure with its own copy. Use the independent ",
         "spelling `||`.", call. = FALSE)
  }
  invisible(NULL)
}


# Build the deviation as the arm-specific field block the joint fitter already
# takes, on the shared field's own graph and site -> node map. `residual` says
# the arm's field IS a deviation from the shared one, so it is never a second
# geometry: it reuses the occurrence field's nodes exactly, which is also why it
# needs no graph of its own here (the fitter's CSR is built from the shared
# graph and every arm-specific block reads it).
#
# Shape matches `.tobs_armspecific_bar_fields()` -- the same list a placed
# positive-arm bar produces -- so the fitter has ONE arm-specific block builder
# and this door cannot drift from it.
.occu_cover_residual_armspec <- function(res, site_cell, n_sites) {
  if (is.null(res)) return(NULL)
  list(
    arm     = "pos",
    type    = "icar",
    graph   = NULL,
    node    = NULL,
    idx_obs = as.integer(site_cell),
    by      = NULL,
    fields  = list(list(is_intercept = TRUE,
                        column_name  = "(Intercept)",
                        weight       = rep(1.0, n_sites))))
}


# The basis a rank-r deviation is carried on.
#
# Which r directions? Not the r smoothest Laplacian modes on their own: those
# are the directions the ICAR prior puts most of its variance in, so they are
# also where the copied field lives, and a basis built from them lands on the
# subspace the copy already spans -- the maximally confounded one. #110 measured
# that confounding for the full-rank version, where the free cover field and the
# alpha copy trade against each other.
#
# So the basis is the smooth modes PROJECTED off the shared field. Writing
# `Q = D - W` for the ICAR precision with eigenpairs `(lambda_j, v_j)` ordered by
# increasing eigenvalue and dropping the null direction of each connected
# component,
#
#     A = [v_1 / sqrt(lambda_1), ..., v_r / sqrt(lambda_r)]
#     B = (I - w w' / <w, w>) A
#
# `A A'` is the ICAR covariance truncated to its r highest-variance directions,
# so `u ~ N(0, sigma^2 I_r)` makes `A u` a reduced-rank ICAR field, and the
# projection removes the one direction the copy already carries: `B u` is
# orthogonal to `w` for every `u`, exactly and by construction, not on average.
#
# The columns are then scaled so the geometric mean of the implied per-node
# marginal variance is 1 (Sorbye-Rue, the convention the engine's areal blocks
# report `sigma` in), so a rank-r deviation's SD is on the same scale as the
# full-rank field's and the two can be read against each other -- which is what
# a rank sweep needs.
.tobs_residual_basis <- function(graph, w_ref, rank, tol = 1e-8) {
  n <- nrow(graph)
  Q <- diag(rowSums(graph), n, n) - graph
  ev  <- eigen(Q, symmetric = TRUE)
  ord <- rev(seq_len(n))                       # eigen() returns descending
  lam <- ev$values[ord]
  V   <- ev$vectors[, ord, drop = FALSE]
  keep <- which(lam > tol * max(1, max(lam)))  # drop one null per component
  r_max <- length(keep)
  if (r_max < 1L) {
    stop("share(residual = ): the areal graph has no non-null direction to ",
         "carry a deviation.", call. = FALSE)
  }
  if (rank > r_max) {
    stop(sprintf(paste0(
      "share(residual = %d): the graph carries at most %d basis functions ",
      "(%d nodes, %d connected component(s)). Ask for at most %d, or ",
      "residual = \"full\"."), rank, r_max, n, n - r_max, r_max), call. = FALSE)
  }
  idx <- keep[seq_len(rank)]
  A   <- V[, idx, drop = FALSE] %*% diag(1 / sqrt(lam[idx]), rank, rank)

  ww <- sum(w_ref * w_ref)
  B  <- if (is.finite(ww) && ww > 0)
          A - outer(as.numeric(w_ref), as.numeric(crossprod(w_ref, A)) / ww)
        else A
  # Sorbye-Rue scaling on the implied marginal variances, so `sigma` means the
  # same thing here as on a full areal block.
  mvar <- rowSums(B * B)
  gm   <- exp(mean(log(pmax(mvar, .Machine$double.eps))))
  if (is.finite(gm) && gm > 0) B <- B / sqrt(gm)
  B
}


# What a rank-r deviation needs before the fit that carries it can be built: the
# shared field to be orthogonal TO, and the scale to be pinned AT. Both come
# from the same warm fit -- the model without the deviation, which is the fit the
# deviation is a refinement of.
#
# One warm fit per fit, not per rank: a rank sweep over the same data reuses one
# reference, which is also what makes the ranks comparable to each other (they
# differ in how many directions the deviation has, not in where it is anchored
# or how big it is allowed to be).
#
# `residual = "full"` needs neither -- it is the arm-specific block, whose scale
# rides its own outer axis -- so it prepares nothing and pays for no warm fit.
.occu_cover_residual_prepare <- function(res, fit_args, graph) {
  if (is.null(res) || identical(res$rank, "full")) return(NULL)
  warm <- do.call(.tobs_fit_occu_cover_joint, fit_args)
  w    <- warm$spatial_field
  if (is.null(w)) {
    stop("occu_cover(): the warm fit for share(residual = ) reported no ",
         "shared field to orthogonalize against.", call. = FALSE)
  }
  # The pin is the shared field's own marginal SD, stated in the SAME convention
  # the basis is built in: `.tobs_residual_basis()` scales its columns to unit
  # geometric-mean marginal variance, so a coefficient SD of `s` makes the
  # deviation a field of geo-mean marginal SD `s` -- which is what
  # `field_sd_mean` reports for the shared field, and not what its raw amplitude
  # `sigma_mean` reports (those differ by the graph's ICAR scale factor). Pinning
  # on the raw number would silently mis-scale the deviation by that factor.
  sigma <- warm$spatial$field_sd_mean %||% warm$spatial$sigma_mean
  if (is.null(sigma) || !is.finite(sigma) || sigma <= 0) sigma <- 1
  list(basis = .tobs_residual_basis(graph, w, res$rank),
       sigma = sigma,
       warm  = warm)
}


# The deviation's posterior mean surface, read off the fit's latent vector.
#
# The r basis coefficients are the LAST r blocks of the joint latent (they trail
# the field and RE blocks, which is what keeps every other reading of the layout
# valid), each block one scalar. Their grid-weighted mean `u` maps back to the
# node surface as `B %*% u` -- the deviation in the same per-node units the
# shared field is reported in.
.occu_cover_residual_field <- function(fit, basis) {
  jf <- fit$joint_fit
  if (is.null(jf) || is.null(jf$modes) || is.null(jf$arm_layout)) return(NULL)
  bstart <- jf$arm_layout$block_start
  bsize  <- jf$arm_layout$block_size
  r      <- ncol(basis)
  if (is.null(bstart) || length(bstart) < r) return(NULL)
  blocks <- seq.int(length(bstart) - r + 1L, length(bstart))
  if (!all(bsize[blocks] == 1L)) return(NULL)
  oc    <- .tobs_joint_ok_cells(jf, "occu_cover residual")
  modes <- oc$fit$modes[oc$ok_cells, bstart[blocks] + 1L, drop = FALSE]
  u     <- as.numeric(crossprod(oc$w, modes))
  as.numeric(basis %*% u)
}


# The identified decomposition of the cover arm's surface.
#
# The fit reports two pieces on the cover arm: the copied field at its
# amplitude, and the deviation. In the raw parameterization those two trade --
# the deviation can carry any multiple of the shared field, so the amplitude the
# outer grid integrates is only weakly pinned, which is the #110 measurement.
# Their SUM is not weakly pinned: it is the cover arm's surface, which the data
# sees.
#
# So the reported amplitude is read off that sum, as its projection onto the
# occurrence arm's own surface,
#
#     alpha = <s_cover, s_occ> / <s_occ, s_occ>,   delta = s_cover - alpha s_occ
#
# leaving `delta` orthogonal to the shared component by construction. Both
# surfaces are taken on the LINEAR-PREDICTOR scale, which is the only scale on
# which they are commensurable: the engine stores every field block as a
# unit-scale latent and keeps its amplitude on the outer grid, so `spatial_field`
# and the deviation's latent are each a shape, not a contribution, and adding
# them unscaled would compare two different units. `alpha` is then dimensionless
# and directly the simulator's / the copy's amplitude.
#
# This is a change of coordinates on the fitted surface, NOT a refit: the model
# is the one the engine solved, and what changes is which half of it is called
# the copy. The raw pieces stay on the fit (`alpha_mean`, `pos_field`) so the two
# readings can be compared rather than one silently replacing the other.
.occu_cover_residual_attach <- function(fit, res) {
  if (is.null(res)) return(fit)
  w     <- fit$spatial_field
  sigma <- fit$spatial$sigma_mean
  if (is.null(w) || is.null(sigma) || !is.finite(sigma)) return(fit)

  # Full rank rides the arm-specific block, so the deviation IS `pos_field`, at
  # that block's own amplitude; a rank-r deviation lives in the trailing
  # basis-coefficient blocks, at the SD it was pinned at, and is read back
  # through its basis.
  #
  # `sigma_d` is the amplitude that multiplies the latent to give the arm's eta
  # contribution; `sigma_field` is the same field's geometric-mean marginal SD,
  # which is what a reader compares across ranks. For a rank-r deviation the two
  # coincide, because the basis is built in that convention. For the full-rank
  # block they differ by the graph's ICAR scale factor, exactly as the shared
  # field's `sigma` and `field_sd` do -- so both are converted here rather than
  # reported in whichever convention their block happened to use.
  scale_q <- sqrt(.occu_cover_icar_scale(fit$spatial$graph))
  if (identical(res$rank, "full")) {
    d_raw   <- fit$pos_field
    sig_nm  <- grep("^sigma_pos_field", names(fit$means), value = TRUE)
    sigma_d <- if (length(sig_nm)) fit$means[[sig_nm[[1L]]]] else NA_real_
    sigma_field <- sigma_d * scale_q
  } else {
    d_raw   <- .occu_cover_residual_field(fit, res$basis)
    sigma_d <- res$sigma
    sigma_field <- sigma_d
  }
  if (is.null(d_raw) || length(d_raw) != length(w)) return(fit)
  if (is.null(sigma_d) || !is.finite(sigma_d)) sigma_d <- 0

  alpha0 <- fit$spatial$alpha_mean
  if (is.null(alpha0) || !is.finite(alpha0)) alpha0 <- 0

  s_occ   <- sigma * w
  s_cover <- alpha0 * s_occ + sigma_d * d_raw
  ss      <- sum(s_occ * s_occ)
  alpha   <- if (ss > 0) sum(s_cover * s_occ) / ss else NA_real_
  delta   <- if (is.finite(alpha)) s_cover - alpha * s_occ else s_cover

  fit$residual <- list(
    rank        = res$rank,
    n_basis     = if (identical(res$rank, "full")) length(w) else res$rank,
    sigma       = sigma_field,
    alpha       = alpha,
    alpha_raw   = alpha0,
    field       = delta,
    cover_field = s_cover,
    occ_field   = s_occ,
    reference   = "spatial_field")
  fit
}
