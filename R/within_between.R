# =============================================================================
# within_between.R - Mundlak-style within/between decomposition (Phase 1d)
#
# Splits each numeric column into a per-group mean (between) and the deviation
# from that mean (within). Driven by the MOTIVATE resurvey workflow: a year
# variable can be decomposed into "mean year per plot" (between-plot
# heterogeneity in survey timing) and "year - mean(year per plot)" (within-
# plot temporal trend), so the within-slope identifies the actual temporal
# signal independent of when each plot was visited.
#
# Reference: Mundlak (1978); Bell & Jones (2015), Pol. Sci. Res. Methods.
# =============================================================================


#' Within / between (Mundlak) decomposition for longitudinal data
#'
#' For each variable in `vars`, adds two new columns to `data`:
#'
#' * `<var><suffix[1]>` - the per-group mean of `var` (the **between**-group
#'   component, constant within each group),
#' * `<var><suffix[2]>` - `var - <var><suffix[1]>` (the **within**-group
#'   deviation, summing to zero within each group).
#'
#' The decomposition is exact: `data[[var]] == data[[btw]] + data[[wtn]]`.
#' Including both columns in a regression (`~ x_btw + x_wtn`) separates the
#' cross-group association from the within-group association without imposing
#' that they be the same coefficient (the random-effects assumption).
#'
#' Typical use in the MOTIVATE cover-resurvey workflow:
#'
#' ```r
#' plots <- within_between(plots, group = "plot_id", vars = "year")
#' fit   <- tobs(~ year_btw + year_wtn + bym2(graph = adj), data = plots,
#'               family = cover("beta"), y = plots$cover,
#'               method = "nested_laplace")
#' ```
#'
#' `year_wtn` then carries the within-plot temporal trend; `year_btw` absorbs
#' across-plot timing heterogeneity (older versus newer baseline plots).
#'
#' @param data A data frame.
#' @param group Character name (or character vector of names) identifying the
#'   grouping column(s). With multiple names, grouping is by the interaction.
#' @param vars Character vector of numeric columns to decompose. Must all
#'   exist in `data`.
#' @param suffix Length-2 character vector giving the suffixes for the
#'   between- and within-group columns. Default `c("_btw", "_wtn")`.
#' @param na.rm Logical; passed to `mean()` when computing per-group means.
#'   Default `TRUE`.
#'
#' @return The input `data` with `2 * length(vars)` new columns appended.
#'
#' @references
#' Mundlak, Y. (1978). On the pooling of time series and cross section data.
#' *Econometrica*, 46(1), 69-85.
#'
#' Bell, A. & Jones, K. (2015). Explaining fixed effects: random effects
#' modeling of time-series cross-sectional and panel data.
#' *Political Science Research and Methods*, 3(1), 133-153.
#'
#' @export
#' @examples
#' set.seed(1)
#' plots <- data.frame(
#'   plot_id = rep(1:5, each = 4),
#'   year    = rep(2000:2003, times = 5) + rep(c(0, 5, 10, 15, 20), each = 4),
#'   cover   = stats::runif(20)
#' )
#' decomposed <- within_between(plots, group = "plot_id", vars = "year")
#' head(decomposed)
within_between <- function(data, group, vars,
                           suffix = c("_btw", "_wtn"),
                           na.rm  = TRUE) {
  if (!is.data.frame(data)) {
    stop("`data` must be a data frame.", call. = FALSE)
  }
  if (!is.character(group) || length(group) < 1L) {
    stop("`group` must be a character vector of column names.", call. = FALSE)
  }
  if (!is.character(vars) || length(vars) < 1L) {
    stop("`vars` must be a character vector of column names.", call. = FALSE)
  }
  if (!is.character(suffix) || length(suffix) != 2L ||
      identical(suffix[1], suffix[2])) {
    stop("`suffix` must be a length-2 character vector with distinct values.",
         call. = FALSE)
  }

  missing_grp <- setdiff(group, names(data))
  if (length(missing_grp)) {
    stop("Group column(s) not in `data`: ",
         paste(shQuote(missing_grp), collapse = ", "), ".", call. = FALSE)
  }
  missing_v <- setdiff(vars, names(data))
  if (length(missing_v)) {
    stop("Variable column(s) not in `data`: ",
         paste(shQuote(missing_v), collapse = ", "), ".", call. = FALSE)
  }
  not_num <- vars[!vapply(data[vars], is.numeric, logical(1))]
  if (length(not_num)) {
    stop("All `vars` must be numeric. Non-numeric: ",
         paste(shQuote(not_num), collapse = ", "), ".", call. = FALSE)
  }

  grp <- if (length(group) == 1L) {
    data[[group]]
  } else {
    interaction(data[group], drop = TRUE, sep = ".")
  }
  grp <- as.factor(grp)

  for (v in vars) {
    nm_btw <- paste0(v, suffix[1])
    nm_wtn <- paste0(v, suffix[2])
    if (nm_btw %in% names(data) || nm_wtn %in% names(data)) {
      stop("Decomposed column name(s) already exist in `data`: ",
           paste(shQuote(intersect(c(nm_btw, nm_wtn), names(data))),
                 collapse = ", "),
           ". Pass a different `suffix`.", call. = FALSE)
    }
    grp_mean <- stats::ave(data[[v]], grp,
                           FUN = function(z) mean(z, na.rm = na.rm))
    data[[nm_btw]] <- grp_mean
    data[[nm_wtn]] <- data[[v]] - grp_mean
  }

  data
}
