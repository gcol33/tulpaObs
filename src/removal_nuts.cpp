// removal_nuts.cpp
// NUTS target for the removal-sampling family (removal()): the (single-line)
// instantiation of the shared count-marginal NUTS machinery
// (marginal_count_nuts.h) with the removal per-site kernel (removal_kernel.h).
// The R reference .tobs_removal_nuts_logpost (R/removal_nuts.R) is the oracle;
// cpp_removal_nuts_joint_logpost is its byte-for-byte cross-check.

#include <Rcpp.h>
#include "tobs_shape.h"
#include "removal_kernel.h"
#include "marginal_count_nuts.h"

using namespace Rcpp;

// Full-vector joint log-posterior + gradient, the cross-check for the R oracle.
// [[Rcpp::export]]
Rcpp::List cpp_removal_nuts_joint_logpost(Rcpp::List spec, Rcpp::NumericVector theta,
                                          double sigma_beta, double sigma_logr) {
    tulpaObs::CountNutsData d = tulpaObs::count_nuts_build_data(spec);
    if ((int) theta.size() != d.total)
        Rcpp::stop("theta length %d != expected %d", (int) theta.size(), d.total);
    Rcpp::NumericVector grad(d.total);
    const double lp = tulpaObs::count_nuts_eval(
        d, theta.begin(), sigma_beta, sigma_logr, grad.begin(),
        &tulpaObs::compute_removal_site);
    return Rcpp::List::create(Rcpp::Named("lp") = lp, Rcpp::Named("grad") = grad);
}

// Run NUTS on the removal target via tulpa's engine and the shared FullGradFn.
// [[Rcpp::export]]
Rcpp::List cpp_removal_nuts(Rcpp::List spec, Rcpp::NumericVector theta0,
                            double sigma_beta, double sigma_logr,
                            Rcpp::Nullable<Rcpp::NumericVector> inv_metric,
                            int n_iter, int n_warmup, int max_treedepth,
                            double adapt_delta, int seed, bool verbose) {
    return tulpaObs::count_nuts_run(
        spec, theta0, sigma_beta, sigma_logr,
        tulpaObs::shape::optional_numeric(inv_metric.get(), "inv_metric"),
        n_iter, n_warmup,
        max_treedepth, adapt_delta, seed, verbose, &tulpaObs::compute_removal_site);
}
