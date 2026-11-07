// graph_components.h
// Connected-component structure of an areal graph, as the in-tree areal drivers
// need it: an intrinsic field (icar, bym2's structured part) has one constant
// null direction per connected component, identified by one sum-to-zero each,
// and a BYM2 structured precision is scaled per component (diag(node_prec) Q,
// node_prec constant within a component; Freni-Sterrantino, Ventrucci & Rue
// 2018). A default-constructed value is the connected case: one component, no
// per-node multiplier.

#ifndef TULPAOBS_GRAPH_COMPONENTS_H
#define TULPAOBS_GRAPH_COMPONENTS_H

#include <Rcpp.h>
#include <vector>

namespace tulpaObs {

struct GraphComponents {
    std::vector<int> comp;          // 0-based component per node; empty = one component
    int n_comp = 1;
    std::vector<double> node_prec;  // per-node precision multiplier; empty = 1

    double prec(int s) const { return node_prec.empty() ? 1.0 : node_prec[s]; }
    int    of(int s)   const { return comp.empty() ? 0 : comp[s]; }
};

// From the R side: `node_comp` (0-based labels, length n or 0) and `node_prec`
// (length n or 0). Lengths are checked here so a kernel never indexes past them.
inline GraphComponents graph_components_from_R(const Rcpp::IntegerVector& node_comp,
                                               const Rcpp::NumericVector& node_prec,
                                               int n) {
    GraphComponents gc;
    if (node_comp.size() > 0) {
        if (node_comp.size() != n)
            Rcpp::stop("node_comp must have one label per spatial unit (%d).", n);
        gc.comp.assign(node_comp.begin(), node_comp.end());
        int mx = 0;
        for (int c : gc.comp) {
            if (c < 0) Rcpp::stop("node_comp labels must be 0-based and non-negative.");
            if (c > mx) mx = c;
        }
        gc.n_comp = mx + 1;
    }
    if (node_prec.size() > 0) {
        if (node_prec.size() != n)
            Rcpp::stop("node_prec must have one value per spatial unit (%d).", n);
        gc.node_prec.assign(node_prec.begin(), node_prec.end());
    }
    return gc;
}

}  // namespace tulpaObs

#endif  // TULPAOBS_GRAPH_COMPONENTS_H
