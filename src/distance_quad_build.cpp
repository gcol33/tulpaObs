// distance_quad_build.cpp
// R-level factory for the per-fit detection quadrature (distance_quad.h): built
// ONCE per fit and threaded as an external pointer into every repeated .Call()
// the fit makes into the compiled distance kernel (cpp_distance_site_sweep,
// cpp_distance_total_log_lik, cpp_distance_laplace_fixed, cpp_distance_nuts,
// cpp_distance_grouped_oracle, cpp_distance_ploglik_batch), instead of paying
// the Gauss-Legendre Newton-Raphson root-find on every call (gcol33/tulpaObs#165).

#include "distance_quad.h"
#include <Rcpp.h>
#include <vector>

// [[Rcpp::export]]
SEXP cpp_distance_build_quad(Rcpp::NumericVector cutpoints, int transect,
                             int quad_order) {
    std::vector<double> cut(cutpoints.begin(), cutpoints.end());
    return tulpaObs::dist_quad_wrap(tulpaObs::dist_build_quad(cut, transect, quad_order));
}
