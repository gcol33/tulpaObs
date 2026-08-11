# Compact-but-PRINCIPLED outer-grid for a correlated-slope `miid` block's free
# Sigma, in the engine's column-major lower-triangular log-Cholesky coordinates
# (gcol33/tulpa#114). A p = 2 (intercept + one slope) block -- the common
# `(1 + x | g)` -- gets a (sigma_0, sigma_1, rho) tensor sized to compose with the
# shared field + copy amplitude under the engine's outer-grid cap. The grid is a
# coarsened version of the engine's `.mcar_default_logchol_grid`: SYMMETRIC
# correlation nodes that include 0 and reach strong +/- (so the marginal
# correlation is not forced into a lop-sided range), and log-spaced SD nodes
# spanning small to large. The slope covariate is standardized
# (.occu_cover_obs_re_design), so a fixed SD bracket is meaningful for any
# covariate scale. Users widen it via `control$re.logchol.grid.p` /
# `re.logchol.grid.pos`. For p != 2 (multi-slope correlated, rare) return NULL so
# the engine fills its own default; that design is grid-heavy and usually needs
# an explicit coarse grid.
.occu_cover_miid_logchol_grid <- function(p, sig_grid = NULL, rho_grid = NULL) {
  if (!identical(as.integer(p), 2L)) return(NULL)
  sig_grid <- sig_grid %||% c(0.35, 0.8, 1.6)
  rho_grid <- rho_grid %||% c(-0.7, -0.3, 0, 0.3, 0.7)
  g <- expand.grid(s1 = sig_grid, s2 = sig_grid, rho = rho_grid,
                   KEEP.OUT.ATTRS = FALSE)
  out <- cbind(L11 = log(g$s1),
               L21 = g$rho * g$s2,
               L22 = log(g$s2 * sqrt(1 - g$rho^2)))
  colnames(out) <- c("L11", "L21", "L22")
  as.matrix(out)
}

# Public hyperparameter name for each RE block, aligned with the descriptor
# list (gcol33/tulpaObs#103). A lone term on an arm keeps the legacy bare name
# (sigma_re / sigma_re_p / sigma_re_pos for psi / detection / positive cover);
# crossed / nested terms sharing an arm are disambiguated by the grouping var
# (sigma_re_p_<var>), so every block's variance gets a distinct, stable name.
.occu_cover_re_sigma_names <- function(re_descs) {
  if (!length(re_descs)) return(character(0))
  base <- c(psi = "sigma_re", p = "sigma_re_p", pos = "sigma_re_pos")
  arms <- vapply(re_descs, `[[`, character(1), "arm")
  counts <- table(arms)
  vapply(seq_along(re_descs), function(i) {
    arm <- re_descs[[i]]$arm
    nm  <- base[[arm]]
    if (counts[[arm]] > 1L) {
      var <- re_descs[[i]]$var
      nm  <- paste0(nm, "_", make.names(if (is.na(var)) as.character(i) else var))
    }
    nm
  }, character(1))
}

# `fit$re` key for each RE block, on the same rule as the hyperparameter names: a
# lone term on an arm is keyed by the arm ("psi" / "p" / "pos"), crossed / nested
# terms sharing an arm by "<arm>:<var>". ranef() stacks the list and predict()
# sums every term on the arm it predicts, so both read these keys.
.occu_cover_re_keys <- function(arms, vars) {
  counts <- table(arms)
  vapply(seq_along(arms), function(i)
    if (counts[[arms[[i]]]] > 1L) paste0(arms[[i]], ":", vars[[i]])
    else arms[[i]], character(1))
}

# Post-process an occu_cover joint-coupled engine fit into a tobs_fit. `fit` is
# the tulpa_nested_laplace_joint return (single-species) or a per-species slice
# of a batched fused fit assembled to the same shape (gcol33/tulpa#66); `ctx`
# carries the Part-A context the shaping needs. Marginalisation here is a
# weighted sum over outer-grid cells (order-invariant), so a fused fixed-grid
# slice and an adaptive single-species fit shape identically given the same
# cells.
.occu_cover_jc_postprocess <- function(fit, ctx) {
  adj            <- ctx$adj
  is_latent      <- ctx$is_latent
  pi_list        <- ctx$pi_list
  n_cells        <- ctx$n_cells
  disp2_fixed    <- ctx$disp2_fixed
  n_quad_latent  <- ctx$n_quad_latent
  sigma_pos_init <- ctx$sigma_pos_init
  has_trend      <- ctx$has_trend
  n_trend        <- ctx$n_trend
  coupled_trends <- ctx$coupled_trends
  model          <- ctx$model

  # Unpack per-arm posterior means + SDs from the joint modes.
  # arm_layout$beta_start[k] is the 0-based offset of arm k's betas in the
  # latent vector; arm_layout$p[k] is the count.
  layout <- fit$arm_layout
  p_psi  <- layout$p[1L]
  p_p    <- layout$p[2L]
  p_pos  <- layout$p[3L]
  bpsi_idx <- layout$beta_start[1L] + seq_len(p_psi)
  bp_idx   <- layout$beta_start[2L] + seq_len(p_p)
  bpos_idx <- layout$beta_start[3L] + seq_len(p_pos)

  # Converged outer-grid cells, their softmax weights and the reconciled
  # fit$weights (joint_postprocess_shared.R -- shared with .occu_jc_postprocess).
  oc       <- .tobs_joint_ok_cells(fit, "occu_cover joint")
  ok_cells <- oc$ok_cells
  w        <- oc$w
  fit      <- oc$fit

  modes <- fit$modes[ok_cells, , drop = FALSE]
  beta_psi_m  <- as.numeric(crossprod(w, modes[, bpsi_idx, drop = FALSE]))
  beta_p_m    <- as.numeric(crossprod(w, modes[, bp_idx,   drop = FALSE]))
  beta_pos_m  <- as.numeric(crossprod(w, modes[, bpos_idx, drop = FALSE]))

  # Joint (betas + field) posterior covariance by the law of total covariance
  # over the outer hyperparameter grid, x = (beta_psi, beta_p, beta_pos, field)
  # stacked. Carrying the field block (not just the betas) means downstream
  # derived quantities (delta_p, delta_cover) can marginalize the joint
  # betas+field posterior instead of a marginal-only diagonal. The arm-count-
  # agnostic computation lives in joint_postprocess_shared.R, shared with the
  # 2-arm `.occu_jc_postprocess()`; only the arm list and the field index differ.
  p_beta    <- p_psi + p_p + p_pos
  mcar      <- isTRUE(ctx$mcar)
  # One latent field block per coupled spatial field. The multi-block layout
  # reports every ICAR/structured field's 0-based offset in `field_starts`;
  # the single-block joint layout reports the one field's offset in `phi_start`.
  # The correlated MCAR field is a SINGLE block of p sub-fields laid out
  # contiguously (slot a*n_cells + cell for sub-field a), so the decode treats it
  # as p sub-fields over one contiguous run -- every per-field summary below
  # (demeaning, z-tables, block slicing) is shared with the independent path.
  field_starts0 <- layout$field_starts %||% layout$phi_start
  if (mcar) {
    n_fields  <- as.integer(ctx$n_fields_mcar)
    field_idx <- as.integer(field_starts0[[1L]] + seq_len(n_fields * n_cells))
  } else {
    n_fields  <- length(field_starts0)
    field_idx <- as.integer(unlist(lapply(field_starts0,
                                          function(s0) s0 + seq_len(n_cells))))
  }
  bfv <- .tobs_joint_beta_field_vcov(
    fit, modes, w, ok_cells,
    arms = list(list(idx = bpsi_idx, mean = beta_psi_m),
                list(idx = bp_idx,   mean = beta_p_m),
                list(idx = bpos_idx, mean = beta_pos_m)),
    field_idx = field_idx, n_cells = n_cells, n_fields = n_fields,
    n_threads = ctx$n_threads)
  sds_beta       <- bfv$sds_beta
  beta_block     <- bfv$beta_block
  field_demeaned <- bfv$field_demeaned
  field_sd       <- bfv$field_sd
  Vj             <- bfv$Vj

  # Per-group RE BLUPs (gcol33/tulpaObs#56, #102, #103). The RE blocks trail the
  # n_fields field blocks, so term i's blocks sit at layout positions
  # n_fields + block_start .. (+ n_blocks - 1); each block's latent is a
  # contiguous run in `modes`. The per-(group, coefficient) posterior-mean offset
  # is the grid-weighted mean of those columns (centred per coefficient), the SD
  # their grid-weighted posterior SD. A random intercept / uncorrelated slope has
  # one `iid` block per coefficient (`block c` -> column `c`); a correlated slope
  # is one `miid` block whose latent is coefficient-major (field a, group g at
  # (a-1)*n_groups + g), reshaped to [n_groups x n_coefs]. `latent_idx` is the
  # coefficient-major column run (the predict draws reshape it identically). The
  # per-coefficient variance / correlation is filled from the hyper axes below.
  re_descs <- ctx$re_descs %||% list()
  re_sig_names <- .occu_cover_re_sigma_names(re_descs)
  re_terms <- vector("list", length(re_descs))
  if (length(re_descs) > 0L) {
    bstart <- layout$block_start
    bsize  <- layout$block_size
    for (i in seq_along(re_descs)) {
      d   <- re_descs[[i]]
      nc  <- d$n_coefs; ng <- d$n_groups
      lay <- n_fields + d$block_start + seq_len(d$n_blocks) - 1L
      if (is.null(bstart) || length(bstart) < max(lay)) next
      B_mean <- matrix(0, ng, nc); B_sd <- matrix(0, ng, nc); lat <- integer(0)
      grid_moments <- function(cols) {
        u_mod <- modes[, cols, drop = FALSE]
        u_hat <- as.numeric(crossprod(w, u_mod))
        list(mean = u_hat, var = as.numeric(crossprod(w, u_mod^2)) - u_hat^2)
      }
      if (!isTRUE(d$correlated)) {
        for (cc in seq_len(nc)) {
          cols <- bstart[lay[cc]] + seq_len(bsize[lay[cc]])   # length n_groups
          mom  <- grid_moments(cols)
          B_mean[, cc] <- mom$mean - mean(mom$mean)
          B_sd[, cc]   <- sqrt(pmax(mom$var, 0))
          lat <- c(lat, cols)
        }
      } else {
        cols <- bstart[lay[1L]] + seq_len(bsize[lay[1L]])     # n_coefs*n_groups
        mom  <- grid_moments(cols)
        Hm <- matrix(mom$mean, ng, nc); Hs <- matrix(sqrt(pmax(mom$var, 0)), ng, nc)
        for (cc in seq_len(nc)) B_mean[, cc] <- Hm[, cc] - mean(Hm[, cc])
        B_sd[] <- Hs
        lat <- cols
      }
      # Back-transform a slope coefficient's BLUP from the standardized covariate
      # the fit ran on to its natural units (`b_raw = b_std / scale`); the
      # intercept's scale is 1. The per-coefficient SD is rescaled below the same
      # way; correlation is scale-free.
      sc <- d$coef_scales %||% rep(1, nc)
      if (any(sc != 1)) {
        B_mean <- sweep(B_mean, 2L, sc, "/")
        B_sd   <- sweep(B_sd,   2L, sc, "/")
      }
      re_terms[[i]] <- c(d, list(blup_mat = B_mean, blup_sd_mat = B_sd,
                                 prior_pos = lay, latent_idx = as.integer(lat)))
    }
    re_terms <- Filter(Negate(is.null), re_terms)
  }

  sd_psi <- sds_beta[seq_len(p_psi)]
  sd_p   <- sds_beta[p_psi + seq_len(p_p)]
  sd_pos <- sds_beta[p_psi + p_p + seq_len(p_pos)]

  means <- c(beta_psi_m, beta_p_m, beta_pos_m)
  sds   <- c(sd_psi,    sd_p,     sd_pos)

  par_names <- c(
    paste0("psi_", pi_list[[1L]]$coef_names),
    paste0("p_",   pi_list[[2L]]$coef_names),
    paste0("pos_", pi_list[[3L]]$coef_names)
  )

  # Hyperparams: surface sigma, alpha (and optionally sigma_pos when on the
  # outer grid) from the joint posterior moments. Recompute on the filtered
  # grid -- fit$theta_mean uses the engine's full weights, which are NaN
  # when any cell's log_marginal is non-finite.
  tg_full   <- fit$theta_grid
  tg_ok     <- tg_full[ok_cells, , drop = FALSE]
  tg_names  <- colnames(tg_full)
  hyper_means <- numeric(0)
  hyper_sds   <- numeric(0)
  hyper_vals  <- list()
  hyper_names <- character(0)
  pick <- function(name, public = name) {
    j <- match(name, tg_names)
    if (is.na(j)) return(invisible(NULL))
    vals <- as.numeric(tg_ok[, j])
    m <- sum(w * vals)
    v <- sum(w * vals^2) - m^2
    hyper_means[[public]] <<- m
    hyper_sds  [[public]] <<- sqrt(max(v, 0))
    hyper_vals [[public]] <<- vals
    hyper_names <<- c(hyper_names, public)
  }
  # On the latent path the pos arm's phi axis IS the cover-latent SD; surface it
  # as `sigma_u` rather than the engine's generic `phi_pos`.
  phi_pos_public <- if (is_latent) "sigma_u" else "phi_pos"
  # pick2 reads a multi-block axis column (`b<k>.<axis>`) under a public name.
  pick2 <- function(public, col) {
    j <- match(col, tg_names)
    if (is.na(j)) return(invisible(NULL))
    vals <- as.numeric(tg_ok[, j])
    m <- sum(w * vals)
    v <- sum(w * vals^2) - m^2
    hyper_means[[public]] <<- m
    hyper_sds  [[public]] <<- sqrt(max(v, 0))
    hyper_vals [[public]] <<- vals
    hyper_names <<- c(hyper_names, public)
  }
  # put_derived stores the grid-weighted moments of a DERIVED per-cell quantity
  # (e.g. a sigma / correlation reconstructed from log-Cholesky axes), so it is
  # marginalized over the joint posterior rather than plugged in at the mode.
  put_derived <- function(public, vals) {
    vals <- as.numeric(vals)
    mn <- sum(w * vals); vv <- sum(w * vals^2) - mn^2
    hyper_means[[public]] <<- mn
    hyper_sds  [[public]] <<- sqrt(max(vv, 0))
    hyper_vals [[public]] <<- vals
    hyper_names <<- c(hyper_names, public)
  }
  # Sorbye-Rue geo-mean marginal SD (gcol33/tulpaObs#221). `sigma` on this path
  # is the raw amplitude against the unscaled intrinsic precision Q = D - W --
  # the shared/trend occupancy field supports only icar here (bym2 is coerced
  # to icar with rho fixed to 1 upstream, `.occu_cover_spatial_fields()`), so
  # every field on this path shares the same Sorbye-Rue scale. That differs
  # from the #204 NUTS path's `field_sd`, which states the geo-mean marginal SD
  # directly, by sqrt(scale_q) -- a FIXED constant given the graph, not an
  # independent grid axis, so it is reported alongside `sigma` (`fit$field_sd`
  # / `field_sd_trend...`) rather than folded into `hyper_names`: doing that
  # would add a row/column to the joint parameter vcov perfectly collinear
  # with `sigma`'s (correlation exactly 1), singular by construction, which is
  # exactly what broke `chol(V)` in `test-occu-cover-joint.R` on the first cut
  # of this fix. Kept OUT of `means`/`sds`/`vcov`/`draws` for that reason.
  field_sd_summary <- list()
  put_field_sd <- function(sigma_name, out_name) {
    if (is.null(hyper_means[[sigma_name]])) return(invisible())
    scale_q <- sqrt(.occu_cover_icar_scale(adj))
    field_sd_summary[[out_name]] <<- list(
      mean = hyper_means[[sigma_name]] * scale_q,
      sd   = hyper_sds[[sigma_name]]   * scale_q)
  }
  # Clean a coefficient name for use in a hyperparameter name: `(Intercept)` ->
  # `intercept`, other punctuation collapsed to `_`.
  .re_coef_tag <- function(x) {
    x <- gsub("\\(Intercept\\)", "intercept", x)
    gsub("(^_|_$)", "", gsub("[^A-Za-z0-9]+", "_", x))
  }
  # A per-group RE block also forces the multi-block driver (its iid block is an
  # extra prior block), so the field axes carry the `b<k>.` names even without a
  # trend field. Each RE block's variance is its `b<n_fields+i>.sigma` axis.
  has_re     <- !is.null(ctx$re_spec)
  has_any_re <- length(re_descs) > 0L
  # Arm-specific cover field blocks (gcol33/tulpaObs#110) trail the occupancy
  # field blocks: fields 1..n_occ_fields are the shared occupancy intercept +
  # trends (copied to cover with alpha), the rest are the non-copied pos-arm
  # fields. Each reports its own sigma (b<k>.sigma) with no alpha copy axis.
  has_pos_armspec <- isTRUE(ctx$has_pos_armspec)
  n_occ_fields    <- ctx$n_occ_fields %||% n_fields
  pos_field_specs <- ctx$pos_field_specs %||% list()
  if (mcar) {
    # Free-Sigma MCAR hyperparameters (gcol33/tulpaObs#63). Reconstruct Sigma per
    # outer-grid cell from the log-Cholesky axes b1.L<i><j>, derive each field SD
    # (sigma_mcar<a>, a = 1 intercept .. p) and each cross-correlation
    # (rho_mcar_<a><b>), and carry the per-cell values so their grid-weighted
    # moments AND their between-cell covariance with the betas / other hypers are
    # marginalized over the joint posterior, not plugged in at the mode (the
    # marginalize-derived-quantities rule). The whole correlated field's copy
    # amplitude onto the cover arm is the alpha_mcar grid axis.
    p_f  <- n_fields
    m_lc <- p_f * (p_f + 1L) / 2L
    axis_nm <- character(m_lc); tt <- 1L
    for (j in seq_len(p_f)) for (i in j:p_f) {
      axis_nm[tt] <- sprintf("b1.L%d%d", i, j); tt <- tt + 1L
    }
    lc_cols <- match(axis_nm, tg_names)
    n_ok    <- length(ok_cells)
    sd_mat  <- matrix(NA_real_, n_ok, p_f)
    n_rho   <- p_f * (p_f - 1L) / 2L
    rho_mat <- matrix(NA_real_, n_ok, n_rho)
    for (k in seq_len(n_ok)) {
      L   <- .cover_mcar_logchol_to_L(as.numeric(tg_ok[k, lc_cols]), p_f)
      Sig <- L %*% t(L)
      sds_k <- sqrt(pmax(diag(Sig), 0))
      sd_mat[k, ] <- sds_k
      cc <- 1L
      for (a in seq_len(p_f - 1L)) for (b in (a + 1L):p_f) {
        rho_mat[k, cc] <- Sig[a, b] / max(sds_k[a] * sds_k[b], 1e-12); cc <- cc + 1L
      }
    }
    for (a in seq_len(p_f)) put_derived(sprintf("sigma_mcar%d", a), sd_mat[, a])
    cc <- 1L
    for (a in seq_len(p_f - 1L)) for (b in (a + 1L):p_f) {
      put_derived(sprintf("rho_mcar_%d%d", a, b), rho_mat[, cc]); cc <- cc + 1L
    }
    pick2("alpha_mcar", "b1.alpha")
    pick("phi_pos", phi_pos_public)
  } else if (has_trend || has_any_re || has_pos_armspec) {
    # Multi-block: block 1 is the intercept field, blocks 2.. the trend fields,
    # then the RE block(s). A single trend field keeps the bare
    # sigma_trend/alpha_trend names; several are indexed (sigma_trend1, ...).
    # Each RE block's SD is reported by `re_sig_names[i]` (sigma_re / sigma_re_p /
    # sigma_re_pos for a lone term on an arm; suffixed by grouping var for crossed
    # / nested terms sharing an arm).
    pick2("sigma", "b1.sigma")
    pick2("alpha", "b1.alpha")
    put_field_sd("sigma", "field_sd")
    for (j in seq_len(n_trend)) {
      suffix <- if (n_trend == 1L) "" else as.character(j)
      pick2(paste0("sigma_trend", suffix), sprintf("b%d.sigma", j + 1L))
      pick2(paste0("alpha_trend", suffix), sprintf("b%d.alpha", j + 1L))
      put_field_sd(paste0("sigma_trend", suffix), paste0("field_sd_trend", suffix))
    }
    # Arm-specific cover fields (gcol33/tulpaObs#110): blocks n_occ_fields+1 ..
    # n_fields, each a NON-copied ICAR with its own precision axis (b<k>.tau,
    # sigma = 1/sqrt(tau)) and NO alpha copy. Report the grid-weighted marginal SD
    # (marginalize-derived-quantities). A lone intercept field keeps the bare
    # `sigma_pos_field`; a covariate column is suffixed by its name.
    for (j in seq_along(pos_field_specs)) {
      spec  <- pos_field_specs[[j]]
      blk_k <- n_occ_fields + j
      # sigma_pos_field for the cover arm, sigma_p_field for the detection arm.
      base_nm <- paste0("sigma_", if (identical(spec$arm, "p")) "p" else "pos",
                        "_field")
      nm <- if (isTRUE(spec$is_intercept)) base_nm
            else paste0(base_nm, "_", .re_coef_tag(spec$column_name))
      tau_col <- sprintf("b%d.tau", blk_k)
      sig_col <- sprintf("b%d.sigma", blk_k)
      if (sig_col %in% tg_names) {
        pick2(nm, sig_col)
      } else if (tau_col %in% tg_names) {
        put_derived(nm, 1.0 / sqrt(as.numeric(tg_ok[, match(tau_col, tg_names)])))
      }
    }
    # Per-term RE variance components. An intercept / uncorrelated-slope term
    # has one `b<P>.sigma` axis per coefficient; a correlated-slope term has one
    # `miid` block whose log-Cholesky axes b<P>.L<ij> reconstruct a free Sigma,
    # marginalized to per-coefficient SDs + cross-correlations over the grid. The
    # per-coefficient SD (and correlation) are stored back on `re_terms` so the
    # fit summary keeps the structured covariance, not just the scalar names.
    # Back-transform a slope coefficient's reported SD from the standardized
    # covariate the fit ran on to its natural units (`sigma_raw = sigma_std /
    # scale`); the intercept scale is 1. Divides the stored hyper mean / SD / per-
    # cell values so the fit summary AND its vcov are on the natural scale.
    rescale_hyper <- function(nm, s) {
      if (s == 1 || is.null(hyper_means[[nm]])) return(invisible())
      hyper_means[[nm]] <<- hyper_means[[nm]] / s
      hyper_sds  [[nm]] <<- hyper_sds  [[nm]] / s
      hyper_vals [[nm]] <<- hyper_vals [[nm]] / s
    }
    for (i in seq_along(re_terms)) {
      trm  <- re_terms[[i]]
      base <- re_sig_names[i]
      nc   <- trm$n_coefs; cn <- trm$coef_names; PP <- trm$prior_pos
      sc   <- trm$coef_scales %||% rep(1, nc)
      tags <- vapply(cn, .re_coef_tag, character(1))
      sigma_vec <- stats::setNames(rep(NA_real_, nc), cn)
      cor_mat   <- NULL
      if (!isTRUE(trm$correlated)) {
        for (cc in seq_len(nc)) {
          nm <- if (nc == 1L) base else paste0(base, "_", tags[cc])
          pick2(nm, sprintf("b%d.sigma", PP[cc]))
          rescale_hyper(nm, sc[cc])
          sigma_vec[cc] <- hyper_means[[nm]] %||% NA_real_
        }
      } else {
        p_f <- nc; P1 <- PP[1L]
        axis_nm <- character(p_f * (p_f + 1L) / 2L); tt <- 1L
        for (jj in seq_len(p_f)) for (ii in jj:p_f) {
          axis_nm[tt] <- sprintf("b%d.L%d%d", P1, ii, jj); tt <- tt + 1L
        }
        lc_cols <- match(axis_nm, tg_names)
        n_ok2 <- length(ok_cells)
        sd_mat2 <- matrix(NA_real_, n_ok2, p_f)
        n_rho2  <- p_f * (p_f - 1L) / 2L
        rho_mat2 <- matrix(NA_real_, n_ok2, max(n_rho2, 1L))
        for (k in seq_len(n_ok2)) {
          L   <- .cover_mcar_logchol_to_L(as.numeric(tg_ok[k, lc_cols]), p_f)
          Sig <- L %*% t(L); sds_k <- sqrt(pmax(diag(Sig), 0))
          sd_mat2[k, ] <- sds_k; rr <- 1L
          for (a in seq_len(p_f - 1L)) for (b in (a + 1L):p_f) {
            rho_mat2[k, rr] <- Sig[a, b] / max(sds_k[a] * sds_k[b], 1e-12); rr <- rr + 1L
          }
        }
        cor_mat <- diag(nc)
        for (cc in seq_len(nc)) {
          nm <- paste0(base, "_", tags[cc])
          put_derived(nm, sd_mat2[, cc] / sc[cc])   # natural-scale slope SD
          sigma_vec[cc] <- hyper_means[[nm]]
        }
        cbase <- sub("^sigma", "cor", base); rr <- 1L
        for (a in seq_len(nc - 1L)) for (b in (a + 1L):nc) {
          nm <- paste0(cbase, "_", tags[a], "_", tags[b]); put_derived(nm, rho_mat2[, rr])
          cor_mat[a, b] <- cor_mat[b, a] <- hyper_means[[nm]]; rr <- rr + 1L
        }
      }
      re_terms[[i]]$sigma <- sigma_vec
      re_terms[[i]]$cor   <- cor_mat
    }
    pick("phi_pos", phi_pos_public)
  } else {
    pick("sigma"); pick("alpha"); pick("phi_pos", phi_pos_public)
    put_field_sd("sigma", "field_sd")
  }
  if (length(hyper_names) > 0L) {
    means <- c(means, unlist(hyper_means)[hyper_names])
    sds   <- c(sds,   unlist(hyper_sds)[hyper_names])
    par_names <- c(par_names, hyper_names)
  }

  names(means) <- par_names
  names(sds)   <- par_names

  # Parameter-surface vcov by the law of total covariance over the outer grid
  # (joint_postprocess_shared.R -- shared with .occu_jc_postprocess).
  V <- .tobs_joint_param_vcov(modes, w, bfv$beta_idx, beta_block, p_beta,
                              hyper_names, hyper_vals, hyper_means,
                              means, sds, par_names, Vj)

  n_draws <- 1000L
  draws <- .rmvn(n_draws, means, V)
  colnames(draws) <- par_names

  # Split the stacked per-field summaries into one block of n_cells per coupled
  # field. Occupancy fields are blocks 1..n_occ_fields (the intercept field --
  # the back-compat `spatial_field` -- then the coupled trends); the
  # arm-specific cover fields (gcol33/tulpaObs#110) are the trailing blocks.
  fs <- .tobs_joint_field_split(field_demeaned, field_sd, n_cells, n_fields,
                                n_occ_fields, coupled_trends)
  fblocks         <- fs$blocks
  field_z_table   <- fs$field_z_table
  field_intercept <- fs$intercept
  field_table     <- fs$field_table
  trend_means     <- fs$trend_means
  trend_tables    <- fs$trend_tables
  trend_labels    <- fs$trend_labels
  field_trend       <- fs$field_trend
  trend_field_table <- fs$trend_field_table

  # Arm-specific cover field posteriors (gcol33/tulpaObs#110): the per-cell z
  # tables for the independent cover-arm field(s), surfaced separately from the
  # occupancy fields so the user can inspect the cover trend map. The first is the
  # intercept field (bare `pos_field` / `pos_field_table`); a covariate column is
  # keyed by its name.
  pos_field_blocks <- if (has_pos_armspec && n_fields > n_occ_fields)
                        fblocks[(n_occ_fields + 1L):n_fields] else list()
  pos_field_means  <- lapply(pos_field_blocks, function(b) b$mean)
  pos_field_tables <- lapply(pos_field_blocks, field_z_table)
  if (length(pos_field_specs) == length(pos_field_means)) {
    pos_field_labels <- vapply(pos_field_specs, function(s)
      if (isTRUE(s$is_intercept)) "(Intercept)" else s$column_name, character(1))
    names(pos_field_means)  <- pos_field_labels
    names(pos_field_tables) <- pos_field_labels
  }
  pos_field       <- if (length(pos_field_means))  pos_field_means[[1L]]  else NULL
  pos_field_table <- if (length(pos_field_tables)) pos_field_tables[[1L]] else NULL

  # Joint betas+field posterior for downstream derived-quantity prediction
  # (delta_p / delta_cover marginalized over the full correlated posterior).
  # `joint_means` carries every field in the same demeaned convention as
  # `spatial_field`; `joint_vcov` is the law-of-total-covariance Vj (NULL on
  # the older-tulpa diagonal fallback). Fields are stacked in block order
  # (occupancy intercept, occupancy trends, then arm-specific cover fields).
  field_par_names <- .tobs_joint_field_par_names(
    n_cells, max(n_occ_fields - 1L, 0L), n_fields - n_occ_fields)
  joint_par_names <- c(
    paste0("psi_",   pi_list[[1L]]$coef_names),
    paste0("p_",     pi_list[[2L]]$coef_names),
    paste0("pos_",   pi_list[[3L]]$coef_names),
    field_par_names
  )
  joint_means <- c(beta_psi_m, beta_p_m, beta_pos_m, field_demeaned)
  names(joint_means) <- joint_par_names
  if (!is.null(Vj)) dimnames(Vj) <- list(joint_par_names, joint_par_names)

  log_lik_val <- sum(w * fit$log_marginal[ok_cells])

  # Record the fixed within-unit dispersion (sigma_eps / beta precision) so the
  # latent-path predict can reconstruct the marginal cover (it pairs with the
  # integrated sigma_u reported in `means`).
  if (is_latent) {
    model$cover_latent_disp2 <- disp2_fixed
    model$cover_latent_nquad <- n_quad_latent
  }
  # Record the pos-arm dispersion the spec held fixed (sigma_pos for non-latent;
  # the latent path integrates sigma_u on the grid instead). The pointwise
  # log-likelihood reads it to score the cover term at the fitted dispersion
  # rather than a bare unit default (gcol33/tulpaObs#34).
  if (!is_latent) model$cover_pos_disp <- sigma_pos_init

  # Spatial summary. The correlated MCAR field reports its per-field SDs
  # (sigma_mcar, intercept first) and cross-correlations (rho_mcar) alongside the
  # single copy amplitude (alpha_mcar); sigma_mean / alpha_mean stay populated
  # (intercept-field SD + copy) so the shared print / summary / predict layer
  # reads it without branching on the field structure.
  if (mcar) {
    rho_nms <- grep("^rho_mcar_", names(hyper_means), value = TRUE)
    spatial_summary <- list(
      type = "mcar", graph = adj,
      sigma_mean = unname(hyper_means["sigma_mcar1"]),
      alpha_mean = unname(hyper_means["alpha_mcar"]),
      sigma_trend_mean = if (n_fields >= 2L)
        unname(hyper_means["sigma_mcar2"]) else NULL,
      alpha_trend_mean = NULL,
      sigma_mcar = vapply(seq_len(n_fields), function(a)
        unname(hyper_means[[sprintf("sigma_mcar%d", a)]]), numeric(1)),
      rho_mcar   = if (length(rho_nms))
        unname(unlist(hyper_means[rho_nms])) else numeric(0),
      rho_mcar_names = rho_nms,
      alpha_mcar = unname(hyper_means["alpha_mcar"]))
  } else {
    fsd <- field_sd_summary[["field_sd"]]
    fsd_trend_nm <- if (n_trend == 1L) "field_sd_trend" else "field_sd_trend1"
    fsd_trend <- field_sd_summary[[fsd_trend_nm]]
    spatial_summary <- list(
      type = "icar", graph = adj,
      sigma_mean = unname(hyper_means["sigma"]),
      alpha_mean = unname(hyper_means["alpha"]),
      # Geo-mean marginal SD (Sorbye-Rue), the #204 NUTS `field_sd` convention
      # (gcol33/tulpaObs#221) -- `sigma_mean` above is the raw amplitude.
      field_sd_mean = fsd$mean, field_sd_sd = fsd$sd,
      sigma_trend_mean = if (has_trend)
        unname(hyper_means[if (n_trend == 1L) "sigma_trend"
                           else "sigma_trend1"]) else NULL,
      alpha_trend_mean = if (has_trend)
        unname(hyper_means[if (n_trend == 1L) "alpha_trend"
                           else "alpha_trend1"]) else NULL,
      field_sd_trend_mean = if (has_trend) fsd_trend$mean else NULL,
      field_sd_trend_sd   = if (has_trend) fsd_trend$sd   else NULL)
  }

  structure(c(list(
    draws        = draws,
    means        = means,
    sds          = sds,
    vcov         = V,
    n_samples    = n_draws,
    n_params     = length(means),
    log_prob     = rep(log_lik_val, n_draws),
    log_lik      = log_lik_val,
    N            = if (isTRUE(model$ragged)) model$n_visits_valid else sum(model$valid)),
    .tobs_na_nuts_diagnostics(n_draws),
    .tobs_promote_pareto_k(fit),
    .tobs_promote_outer_grid(fit),
    list(
    col_names    = par_names,
    param_names  = par_names,
    process_info = pi_list,
    model        = model,
    spatial      = spatial_summary,
    spatial_field = field_intercept,
    trend_field   = field_trend,
    trend_fields  = if (length(trend_means))  trend_means  else NULL,
    field_table  = field_table,
    trend_field_table  = trend_field_table,
    trend_field_tables = if (length(trend_tables)) trend_tables else NULL,
    trend_weight  = if (has_trend) trend_labels[[1L]]
                    else if (length(pos_field_specs)) {
                      cols <- Filter(Negate(is.null),
                                     lapply(pos_field_specs, `[[`, "column_name"))
                      nonint <- Filter(function(s) !isTRUE(s$is_intercept),
                                       pos_field_specs)
                      if (length(nonint)) nonint[[1L]]$column_name else NULL
                    } else NULL,
    trend_weights = if (has_trend) trend_labels        else NULL,
    # Arm-specific cover field (gcol33/tulpaObs#110): `field_specs` labels every
    # field block (shared occupancy vs pos arm) + its weight column so the draw
    # substrate maps each block to (occ, pos) amplitudes; the pos-field tables are
    # the independent cover field posterior for user inspection.
    field_specs      = ctx$field_specs,
    has_pos_armspec  = has_pos_armspec,
    pos_field        = pos_field,
    pos_field_table  = pos_field_table,
    pos_fields       = if (length(pos_field_means))  pos_field_means  else NULL,
    pos_field_tables = if (length(pos_field_tables)) pos_field_tables else NULL,
    joint_par_names = joint_par_names,
    joint_means     = joint_means,
    joint_vcov      = Vj,
    method       = "joint",
    positive     = model$positive,
    # Per-term RE summaries (gcol33/tulpaObs#56, #102, #103): a flat list, one
    # entry per RE block, each carrying its arm, grouping var + observed levels,
    # variance component, centred per-group BLUPs + SDs, group count, and latent
    # column indices (for the marginalized predict draws). A lone term on an arm
    # is keyed by the arm ("psi" / "p" / "pos"); crossed / nested terms sharing an
    # arm are keyed "<arm>:<var>". ranef() stacks them; predict() sums every term
    # on the predicted arm.
    re           = if (length(re_terms)) {
                     keys <- .occu_cover_re_keys(
                       vapply(re_terms, `[[`, character(1), "arm"),
                       vapply(re_terms, function(t) as.character(t$var),
                              character(1)))
                     stats::setNames(lapply(re_terms, function(t) {
                       # Intercept term -> per-group BLUP vector (back-compat);
                       # slope term -> [n_groups x n_coefs] matrix with coef-named
                       # columns, plus the slope covariate names predict() weights
                       # rows by. `sigma` is the per-coefficient SD; `cor` the free
                       # cross-coefficient correlation matrix (NULL when scalar).
                       nc <- t$n_coefs
                       blup <- if (nc == 1L) as.numeric(t$blup_mat[, 1L]) else {
                         m <- t$blup_mat; colnames(m) <- t$coef_names; m }
                       blup_sd <- if (nc == 1L) as.numeric(t$blup_sd_mat[, 1L]) else {
                         m <- t$blup_sd_mat; colnames(m) <- t$coef_names; m }
                       covnms <- if (isTRUE(t$has_intercept)) t$coef_names[-1L]
                                 else t$coef_names
                       list(arm        = t$arm,
                            var        = if (is.na(t$var)) NULL else t$var,
                            sigma      = t$sigma,
                            cor        = t$cor,
                            blup       = blup,
                            blup_sd    = blup_sd,
                            n_groups   = t$n_groups,
                            n_coefs    = nc,
                            coef_names = t$coef_names,
                            covariate_names = covnms,
                            coef_scales = t$coef_scales %||% rep(1, nc),
                            has_intercept   = isTRUE(t$has_intercept),
                            correlated = t$correlated,
                            levels     = t$levels,
                            latent_idx = t$latent_idx)
                     }), keys)
                   } else NULL,
    joint_fit    = fit,
    convergence  = list(converged = TRUE, n_iter = NA_integer_)
  )), class = c("tobs_fit", "tulpa_fit"))
}
