// newton_step.h
// The two blocks every Laplace / nested-Laplace fitter in this package closes
// its Newton iteration with:
//
//   * a Cholesky of the observed information that falls back to the
//     complete-data Fisher when the observed block is not positive definite;
//   * a backtracking line search that halves the step until the objective the
//     step is meant to increase comes back within slack of its current value.
//
// Both decide WHICH step is taken; neither touches what a step evaluates to, so
// they are family- and field-agnostic. Every route calls these, which is what
// keeps the convergence knobs -- halving depth, acceptance slack, fallback
// policy, and the message a fallback failure prints -- at one value each.

#ifndef TULPAOBS_NEWTON_STEP_H
#define TULPAOBS_NEWTON_STEP_H

#include <Rcpp.h>
#include <RcppEigen.h>
#include <Eigen/Cholesky>
#include <Eigen/Dense>
#include <cstdio>
#include <string>

namespace tulpaObs {

// Halvings a line search tries before it gives up on the iteration. Twelve
// halvings take the step to 2^-12 of the Newton direction, below which a step
// that still fails to raise the objective is at a mode rather than overshooting.
constexpr int kBacktrackMaxHalvings = 12;

// Halvings the negative-binomial dispersion profile tries. The log-r step is
// clamped to [-1.5, 1.5] before the search starts and the axis is only weakly
// identified from below, so a profile near a boundary needs a finer step than
// the coefficient blocks do before it moves at all.
constexpr int kDispersionMaxHalvings = 25;

// Halvings the shared community field solve tries. Each trial re-sweeps every
// species over the whole field, so the step is cheap to shrink relative to the
// sweep that scores it, and the field enters the objective through a sum over
// species that a single overshooting node can dominate.
constexpr int kFieldMaxHalvings = 20;

// Slack on the backtracking accept test. A Newton step at a converged mode moves
// the objective by less than the rounding error of the sweep that measures it, so
// a strict `obj_try >= obj_cur` would reject the final step and exhaust the
// halving loop; accepting a decrease this small lets the step through instead.
constexpr double kLineSearchSlack = 1e-10;

// Backtracking line search along a Newton direction. `trial(step)` evaluates the
// objective at the trial state and leaves that state in the caller's scratch;
// `commit(step)` moves the accepted trial into the fitter's state. The step
// starts at 1 and halves until the trial objective is finite and within `slack`
// of the current one. Returns whether a step was taken; false leaves the state
// untouched and ends the iteration.
template <class TrialFn, class CommitFn>
inline bool newton_backtrack(double obj_cur, TrialFn trial, CommitFn commit,
                             int max_halvings = kBacktrackMaxHalvings,
                             double slack = kLineSearchSlack) {
    double step = 1.0;
    for (int h = 0; h < max_halvings; ++h) {
        const double obj_try = trial(step);
        if (R_finite(obj_try) && obj_try >= obj_cur - slack) {
            commit(step);
            return true;
        }
        step *= 0.5;
    }
    return false;
}

// The grid coordinates a fit is conditioning on, formatted for the fallback
// warning at the four decimals the hyperparameter grids are specified to. An
// empty `n2` leaves the second coordinate out.
inline std::string newton_grid_label(const char* n1, double v1,
                                     const char* n2 = "", double v2 = 0.0) {
    char buf[96];
    if (n2 != nullptr && n2[0] != '\0')
        std::snprintf(buf, sizeof(buf), "%s %.4f, %s %.4f", n1, v1, n2, v2);
    else
        std::snprintf(buf, sizeof(buf), "%s %.4f", n1, v1);
    return std::string(buf);
}

// Solve H delta = grad by Cholesky, falling back to the complete-data Fisher
// when H is not positive definite. `rebuild_fisher()` returns that matrix
// carrying whatever prior block and ridge the route's H carries; it is built
// only on the fallback. Returns false when the fallback is not positive definite
// either, leaving `delta` untouched -- the caller ends its iteration at the
// current estimate. `where` names the grid point the failure happened at.
template <class RebuildFn>
inline bool solve_with_fisher_fallback(const Eigen::MatrixXd& H,
                                       const Eigen::VectorXd& grad,
                                       RebuildFn rebuild_fisher,
                                       const std::string& where, int iter,
                                       bool verbose, Eigen::VectorXd& delta) {
    Eigen::LLT<Eigen::MatrixXd> chol(H);
    if (chol.info() == Eigen::Success) {
        delta = chol.solve(grad);
        return true;
    }
    const Eigen::MatrixXd H_f = rebuild_fisher();
    Eigen::LLT<Eigen::MatrixXd> chol_f(H_f);
    if (chol_f.info() != Eigen::Success) {
        Rcpp::warning("Cholesky failure at " + where + ", iteration " +
                      std::to_string(iter) + ": the complete-data Fisher "
                      "fallback is not positive definite either. The Newton "
                      "iteration stops at the current estimate.");
        return false;
    }
    delta = chol_f.solve(grad);
    if (verbose) Rcpp::Rcout << "    (Fisher fallback)\n";
    return true;
}

}  // namespace tulpaObs

#endif  // TULPAOBS_NEWTON_STEP_H
