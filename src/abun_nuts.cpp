// abun_nuts.cpp
// NUTS target for the single-species N-mixture (abun()): the (single-line)
// instantiation of the shared count-marginal NUTS machinery
// (marginal_count_nuts.h) with the Royle (2004) per-site kernel. The R reference
// .tobs_abun_nuts_logpost (R/abun_nuts.R) is the oracle; cpp_abun_nuts_joint_logpost
// is its byte-for-byte cross-check before driving tulpa's NUTS engine.

#include <Rcpp.h>
#include "nmix_kernel.h"
#include "marginal_count_nuts.h"

using namespace Rcpp;

// Full-vector joint log-posterior + gradient, the Rcpp cross-check entry for the
// R oracle .tobs_abun_nuts_logpost.
// [[Rcpp::export]]
Rcpp::List cpp_abun_nuts_joint_logpost(Rcpp::List spec, Rcpp::NumericVector theta,
                                       double sigma_beta, double sigma_logr) {
    tulpaObs::CountNutsData d = tulpaObs::count_nuts_build_data(spec);
    if ((int) theta.size() != d.total)
        Rcpp::stop("theta length %d != expected %d", (int) theta.size(), d.total);
    Rcpp::NumericVector grad(d.total);
    const double lp = tulpaObs::count_nuts_eval(
        d, theta.begin(), sigma_beta, sigma_logr, grad.begin(),
        &tulpaObs::compute_nmix_site);
    return Rcpp::List::create(Rcpp::Named("lp") = lp, Rcpp::Named("grad") = grad);
}

// Run NUTS on the N-mixture target via tulpa's engine and the shared FullGradFn.
// [[Rcpp::export]]
Rcpp::List cpp_abun_nuts(Rcpp::List spec, Rcpp::NumericVector theta0,
                         double sigma_beta, double sigma_logr,
                         Rcpp::Nullable<Rcpp::NumericVector> inv_metric,
                         int n_iter, int n_warmup, int max_treedepth,
                         double adapt_delta, int seed, bool verbose) {
    return tulpaObs::count_nuts_run(
        spec, theta0, sigma_beta, sigma_logr, inv_metric, n_iter, n_warmup,
        max_treedepth, adapt_delta, seed, verbose, &tulpaObs::compute_nmix_site);
}
