// re_chol.cpp
// The correlation Cholesky factor of a random-effect term from its raw
// parameters, through the engine's own map, for rebuilding group effects from
// sampler draws.

#include <Rcpp.h>
#include <vector>
#include <tulpa/lkj_chol.h>

// [[Rcpp::export]]
Rcpp::NumericMatrix cpp_re_chol_factor(Rcpp::NumericVector raw, int n) {
    if (raw.size() != n * (n - 1) / 2) {
        Rcpp::stop("cpp_re_chol_factor: %d raw values for a %d x %d factor.",
                   (int) raw.size(), n, n);
    }
    std::vector<double> r(raw.begin(), raw.end());
    std::vector<double> L_flat((std::size_t) n * n, 0.0);
    tulpa::build_L_from_raw(r.data(), n, L_flat.data());
    Rcpp::NumericMatrix L(n, n);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) L(i, j) = L_flat[(std::size_t) i * n + j];
    }
    return L;
}
