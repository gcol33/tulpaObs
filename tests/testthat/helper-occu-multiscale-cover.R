# Bind a multiscale model from a simulate_occu_multiscale_cover() output,
# mirroring the dispatcher (.dispatch_occu_multiscale_cover) so a kernel check
# sees exactly the model the fitters do. `y_pos` overrides the simulated cover
# (used to inject missing-at-random cover).
.omc_bind_model <- function(sim, positive,
                            occ = ~ x_cell, theta = ~ x_plot,
                            det = ~ x_pdet, pos = ~ x_cov,
                            y_pos = sim$y_pos) {
  occ_f <- stats::as.formula(
    paste(deparse(occ), "+ icar(graph = sim$adj, group_var = \"cell\")"))
  si <- tulpaObs:::.occu_cover_spatial_fields(occ_f, sim$data)
  vd_det <- tulpaObs:::.normalize_visits(NULL, det, nrow(sim$y), ncol(sim$y))
  vd_pos <- tulpaObs:::.normalize_visits(NULL, pos, nrow(sim$y), ncol(sim$y))
  tulpaObs:::.tobs_build_occu_multiscale_cover(
    occ_formula = si$fe, theta_formula = theta,
    det_formula = vd_det$det_formula, pos_formula = vd_pos$det_formula,
    data = sim$data, y = sim$y, y_pos = y_pos,
    plot_cell = as.integer(sim$data[[si$group_var]]),
    n_cells = nrow(si$fields[[1L]]$graph), positive = positive,
    det_visit_formula = vd_det$det_visit_formula, det_visit_data = vd_det$visits,
    pos_visit_formula = vd_pos$det_visit_formula, pos_visit_data = vd_pos$visits)
}
