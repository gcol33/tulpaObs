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

// Returns a live GridProgress when `on`, else nullptr (zero overhead). The
// caller ticks it once per completed outer-grid cell and finishes after the
// loop; see the four cpp_nmix_*_spatial / cpp_nested_laplace_nmix_* entries.
inline std::unique_ptr<tulpa_progress::GridProgress> make_grid_progress(
        const char* label, int total, bool on,
        int every, double throttle, const std::string& file) {
    std::unique_ptr<tulpa_progress::GridProgress> gp;
    if (on) {
        gp.reset(new tulpa_progress::GridProgress(label, total, every,
                                                  throttle, file));
    }
    return gp;
}

}  // namespace tulpaObs

#endif  // TULPAOBS_NMIX_PROGRESS_H
