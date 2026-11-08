# =============================================================================
# One generator per family.
#
# A family's generative process is written ONCE, as
#
#   .tobs_gen_<key>(model, theta)   -> one replicate, in the layout `y` takes
#
# where `model` is the fitted model layout (the binder's slots the process
# reads: designs, site maps, the visit pattern of `y`) and `theta` is one
# parameter set: the coefficient vector in the fit's draw-column layout, plus,
# for a family carrying latent structure, its realised field / random-effect
# offsets. Three routes reach it and none carries a copy of the process:
#
#   simulate_<family>()   builds a layout + theta from its truth arguments
#                         (`.tobs_sim_layout_<key>()`) and calls the generator;
#   simulate(fit)         maps each selected posterior draw to theta and calls
#                         it (`.tobs_simulate_<key>()`);
#   sbc(fit)              the registry's `simulate` slot calls it at the
#                         candidate theta (`.tobs_sbc_sim_gen()`).
#
# So the data a user simulates, the posterior predictive replicates and the
# data SBC calibrates against come from one implementation, and a per-family
# test asserts the three agree draw for draw at a common theta and seed.
# =============================================================================


# A generator may return the latent state it realised (occupancy, abundance)
# on a "latent" attribute, which simulate_<family>() reports as truth. The
# response alone is what a replicate is.
.tobs_gen_strip <- function(rep) {
  attr(rep, "latent") <- NULL
  rep
}

# One replicate per selected posterior draw: `theta_of(i)` turns the i-th
# selected draw into the generator's theta.
.tobs_gen_replicates <- function(model, n, theta_of, gen) {
  lapply(seq_len(n), function(i) .tobs_gen_strip(gen(model, theta_of(i))))
}

# The generator of family `key`. Looked up by name at call time, so a table
# built when the namespace loads can name a generator defined in any file.
.tobs_generator <- function(key) {
  fn <- get0(paste0(".tobs_gen_", key), envir = asNamespace("tulpaObs"),
             mode = "function", inherits = FALSE)
  if (is.null(fn)) {
    stop("No generator `.tobs_gen_", key, "()` is defined.", call. = FALSE)
  }
  fn
}

# The SBC replicate at theta: seed the stream, build the model layout, map
# theta to the generator's form, draw, and wrap the response as the data set
# the shared pool and refit read. `key` names the family generator; `layout`
# maps spec -> the model layout (the fitted model by default); `theta` maps
# (spec, theta, layout) -> the generator's theta, and is evaluated before the
# generator draws, so a realisation it draws (a fresh field) comes first on
# the stream; `shape` maps (replicate, spec) -> the data set.
.tobs_sbc_sim_gen <- function(key, shape = .tobs_sbc_shape_y,
                              theta = function(spec, th, lt) th,
                              layout = function(spec) spec$model) {
  function(spec, theta_in, seed) {
    set.seed(seed)
    lt <- layout(spec)
    th <- theta(spec, theta_in, lt)
    shape(.tobs_gen_strip(.tobs_generator(key)(lt, th)), spec)
  }
}

# The data-set shape of a family whose replicate is the response alone, on the
# observed data set's cells and site labels.
.tobs_sbc_shape_y <- function(rep, spec) {
  obs <- spec$data_obs
  .tobs_sbc_dataset(obs$cells, rep, obs$site)
}

# A 1 x p draw matrix holding theta, the form a fit's `draws` slot takes.
.tobs_theta_row <- function(theta) {
  matrix(theta, nrow = 1L, dimnames = list(NULL, names(theta)))
}

# The realised offset of each of a model's first `n_arms` arms: the fitted
# field the arm's predictor carries, NULL where it carries none.
.tobs_arm_offsets <- function(model, n_arms) {
  lapply(seq_len(n_arms), function(k) .tobs_eta_offset(model, k))
}

# An arm offset as a C++ generator takes it: the offset, or an empty vector
# for an arm that carries none.
.tobs_off_arg <- function(off) if (is.null(off)) numeric(0) else as.numeric(off)

# The fit route of a family whose theta is one posterior draw row (selected
# per replicate, the draw sample.int(n, 1) picks), the fitted offsets of its
# first `n_arms` arms, and the fixed nuisance values in `extra`.
.tobs_simulate_draw_gen <- function(object, nsim, gen, n_arms, extra = list()) {
  model <- object$model
  draws <- object$draws
  off   <- .tobs_arm_offsets(model, n_arms)
  .tobs_gen_replicates(model, nsim, function(i)
    c(list(beta = draws[sample.int(nrow(draws), 1L), ], off = off), extra), gen)
}

# The coefficient count of each of a model's first `n_arms` arms.
.tobs_arm_p <- function(model, n_arms) {
  vapply(model$process_info[seq_len(n_arms)], function(pp) as.integer(pp$p),
         integer(1))
}
