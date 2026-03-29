#' @export
print.tulpaOcc <- function(x, ...) {
  type_label <- switch(x$model_type,
    single = "Single-season occupancy model",
    dynamic = "Multi-season dynamic occupancy model",
    community = "Community occupancy model",
    integrated = sprintf("Integrated occupancy model (%d sources)", x$n_sources),
    jsdm = sprintf("Joint species distribution model (%d species)", x$n_species)
  )
  cat(type_label, "\n")

  if (x$model_type == "single") {
    cat(sprintf("  Sites: %d, Max visits: %d\n", x$n_sites, x$max_visits))
  } else if (x$model_type == "dynamic") {
    cat(sprintf("  Sites: %d, Seasons: %d, Max visits: %d\n",
                x$n_sites, x$n_seasons, x$max_visits))
  } else if (x$model_type == "community") {
    cat(sprintf("  Sites: %d, Species: %d, Max visits: %d\n",
                x$n_sites, x$n_species, x$max_visits))
  } else if (x$model_type == "integrated") {
    cat(sprintf("  Sites: %d, Sources: %d\n", x$n_sites, x$n_sources))
  } else if (x$model_type == "jsdm") {
    cat(sprintf("  Sites: %d, Species: %d\n", x$n_sites, x$n_species))
  }

  for (pi in x$process_info) {
    cat(sprintf("  %s covariates (%d): %s\n",
                pi$name, pi$p, paste(pi$coef_names, collapse = ", ")))
  }

  if (x$model_type == "single" && !is.null(x$naive_occ)) {
    cat(sprintf("  Naive occupancy: %.1f%%\n", 100 * x$naive_occ))
  }
  if (x$model_type == "community") {
    cat(sprintf("  Species RE: intercept on psi and p\n"))
  }

  invisible(x)
}

#' @export
print.tulpaOcc_fit <- function(x, ...) {
  model <- x$model
  type_label <- switch(model$model_type,
    single = "tulpaOcc fit (single-season occupancy, NUTS)",
    dynamic = "tulpaOcc fit (dynamic occupancy, NUTS)",
    community = "tulpaOcc fit (community occupancy, NUTS)"
  )
  cat(type_label, "\n")

  if (model$model_type == "single") {
    cat(sprintf("  Sites: %d, Max visits: %d\n", model$n_sites, model$max_visits))
  } else if (model$model_type == "dynamic") {
    cat(sprintf("  Sites: %d, Seasons: %d, Max visits: %d\n",
                model$n_sites, model$n_seasons, model$max_visits))
  } else if (model$model_type == "community") {
    cat(sprintf("  Sites: %d, Species: %d\n", model$n_sites, model$n_species))
  }

  cat(sprintf("  Samples: %d, Step size: %.4f\n", x$n_samples, x$epsilon))
  n_div <- sum(x$divergent)
  if (n_div > 0) cat(sprintf("  WARNING: %d divergent transitions\n", n_div))
  cat("\n")

  # Print intercept probabilities
  for (nm in names(x$intercepts)) {
    label <- switch(nm,
      psi  = "Mean occupancy (intercept)",
      psi1 = "Mean initial occupancy (intercept)",
      p    = "Mean detection (intercept)",
      gamma   = "Mean colonization (intercept)",
      epsilon = "Mean extinction (intercept)"
    )
    if (!is.null(label)) {
      cat(sprintf("%s: %.3f\n", label, x$intercepts[[nm]]))
    }
  }

  invisible(x)
}

# summary inherited from tulpa::summary.tulpa_fit
