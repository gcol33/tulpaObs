// nmix_progress.h
// tulpaObs-side factory over tulpa's shared GridProgress reporter
// (<tulpa/nested_progress.h>, gcol33/tulpa#45). The N-mixture areal/SPDE
// spatial fitters run their own outer hyperparameter-grid loops (they do not
// route through tulpa's run_nested_laplace_grid), so they construct the same
// reporter directly. One factory keeps the construction logic single-source
// across the four nmix spatial translation units.

#ifndef TULPAOBS_NMIX_PROGRESS_H
#define TULPAOBS_NMIX_PROGRESS_H

#include <tulpa/nested_progress.h>
#include <memory>
#include <string>

namespace tulpaObs {

// Returns a live GridProgress when either channel is wanted, else nullptr
// (zero overhead). `on` gates the Rcout console line (the verbose/TTY channel);
// the heartbeat `file` is written whenever it is non-empty, independent of `on`,
// so a detached fit (on = false) with a file still gets its ETA -- the channel
// it exists for (gcol33/tulpaObs#43). The caller ticks it once per completed
// outer-grid cell and finishes after the loop; see the four cpp_nmix_*_spatial
// / cpp_nested_laplace_nmix_* entries.
inline std::unique_ptr<tulpa_progress::GridProgress> make_grid_progress(
        const char* label, int total, bool on,
        int every, double throttle, const std::string& file,
        const std::string& unit = "cells") {
    std::unique_ptr<tulpa_progress::GridProgress> gp;
    if (on || !file.empty()) {
        gp.reset(new tulpa_progress::GridProgress(label, total, every,
                                                  throttle, file,
                                                  /*emit_console=*/on, unit));
    }
    return gp;
}

// Same factory, but reads the four progress knobs from the scoped
// `tulpa.nl_progress` R option that tobs() sets (gcol33/tulpaObs#43). The
// in-tree EM / Newton fitters that run on the main R thread (every Laplace
// fitter does -- it is a .Call, no OpenMP at this level) call this directly
// instead of threading four arguments through their cpp signatures. MUST be
// called on the main thread (it reads an R option). `width` is the loop
// concurrency for the ETA (1 for a serial EM/Newton).
inline std::unique_ptr<tulpa_progress::GridProgress> make_grid_progress_from_option(
        const char* label, int total, int width = 1) {
    SEXP optS = Rf_GetOption1(Rf_install("tulpa.nl_progress"));
    if (optS == R_NilValue || TYPEOF(optS) != VECSXP) return nullptr;
    Rcpp::List opt(optS);
    bool on = opt.containsElementNamed("progress") &&
              Rcpp::as<bool>(opt["progress"]);
    std::string file = opt.containsElementNamed("progress_file")
                         ? Rcpp::as<std::string>(opt["progress_file"])
                         : std::string();
    int every = opt.containsElementNamed("progress_every")
                  ? Rcpp::as<int>(opt["progress_every"]) : 0;
    double throttle = opt.containsElementNamed("progress_throttle")
                        ? Rcpp::as<double>(opt["progress_throttle"]) : 2.0;
    // These in-tree fitters are EM / Newton iterations, not an outer grid.
    auto gp = make_grid_progress(label, total, on, every, throttle, file, "iter");
    if (gp) gp->set_width(width > 0 ? width : 1);
    return gp;
}

}  // namespace tulpaObs

#endif  // TULPAOBS_NMIX_PROGRESS_H
