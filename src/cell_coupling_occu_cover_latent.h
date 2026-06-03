// cell_coupling_occu_cover_latent.h
// Stateful `CellCouplingSpec` for the latent-cover-per-unit occu_cover model
// (the principled counterpart of the mean / median cover aggregation): the
// cover arm carries a per-unit cover random effect u_i ~ N(0, sigma_u^2) shared
// across the unit's detected visits, integrated out per unit. See the math in
// occu_cover_latent.h.
//
// Arm layout matches the aggregated spec -- the pos arm carries ONE row per
// detected occupancy unit (cell) -- but the cover term is the marginal log M_i
// over u_i rather than a single log f_pos of the unit's mean cover:
//
//   det case:  log p_cell = occu_det_psi_p_block + sum_{pos rows j} log M_ij
//   nodet case: log p_cell = occu_nodet_block (pos arm contributes nothing)
//
// Because the cover predictor is unit-level (one eta per pos row), log M_ij is
// a scalar function of that single eta, so the kernel needs only a scalar score
// / observed information per pos row -- no within-arm off-diagonal Hessian, the
// same buffer shape the aggregated spec writes.
//
// Unlike the per-visit / aggregated specs this one is STATEFUL: the per-unit
// cover data (lognormal sufficient statistics or the raw beta cover values) and
// the fixed within-unit dispersion `disp2` (sigma_eps for lognormal, the beta
// precision phi for beta) are captured at construction, indexed by GLOBAL pos
// row. sigma_u is read per call from the pos arm's phi slot (the outer-grid
// axis). The probabilist Gauss-Hermite rule (beta only) is built once at
// construction.

#ifndef TULPAOBS_CELL_COUPLING_OCCU_COVER_LATENT_H
#define TULPAOBS_CELL_COUPLING_OCCU_COVER_LATENT_H

#include <tulpa/cell_coupling.h>
#include <tulpa/gauss_hermite.h>
#include "occu_coupling_shared.h"   // sigmoid_ / occu_det_psi_p_block / occu_nodet_block
#include "occu_cover_latent.h"      // LognormalLatent / BetaLatent / LatentMarginal
#include <string>
#include <vector>

namespace tulpaObs {

template <class PosLatent>
class OccuCoverLatentCoupling final : public tulpa::CellCouplingSpec {
public:
    // `site_data` is indexed by GLOBAL pos-arm row: each entry is the detected
    // cover values of that occupancy unit (already filtered to detected
    // visits). `disp2` is the fixed within-unit dispersion. `n_quad` is the GH
    // node count for the beta marginal (ignored by lognormal).
    OccuCoverLatentCoupling(const std::vector<std::vector<double>>& site_data,
                            double disp2, int n_quad)
        : disp2_(disp2),
          gh_(tulpa::gauss_hermite_prob(n_quad < 1 ? 1 : n_quad)) {
        site_.reserve(site_data.size());
        for (const auto& v : site_data) {
            site_.push_back(PosLatent::make_site(v.data(), (int) v.size()));
        }
    }

    std::vector<int> arm_ids() const override { return {0, 1, 2}; }

    double evaluate_cell(int                       /*cell_idx*/,
                         const tulpa::CellEtas&     etas,
                         const tulpa::CellResponse& y_cell,
                         tulpa::CellDerivs&         out) const override {
        const int    Jc      = etas.n_rows_in_arm(1);
        const double psi     = sigmoid_(etas.eta(0, 0));
        const double sigma_u = y_cell.phi(2);          // outer-grid cover-latent SD
        const bool   want_hess = !out.grad_only;

        bool any_det = false;
        for (int v = 0; v < Jc; v++) {
            if (y_cell.y(1, v) > 0.5) { any_det = true; break; }
        }

        if (any_det) {
            double cell_ll = occu_det_psi_p_block(psi, etas, y_cell, Jc,
                                                  want_hess, out);
            // One pos row per detected occupancy unit in this cell. Each row's
            // cover data is captured by GLOBAL row index; integrate its latent.
            const int n_pos = etas.n_rows_in_arm(2);
            for (int j = 0; j < n_pos; j++) {
                const int grow = y_cell.arm_rows[2][j];
                const double eta_pos = etas.eta(2, j);
                const LatentMarginal lm = PosLatent::marginal(
                    site_[grow], eta_pos, disp2_, sigma_u, gh_, out.curvature);
                cell_ll += lm.log_m;
                out.arm_grad[2][j] = lm.score;
                if (want_hess) out.arm_neg_hess_diag[2][j] = lm.neg_hess;
            }
            return cell_ll;
        }

        return occu_nodet_block(psi, etas, Jc, want_hess, out);
    }

    std::string name() const override {
        return std::string(PosLatent::spec_name());
    }

    bool thread_safe() const override { return true; }

private:
    std::vector<typename PosLatent::SiteData> site_;
    double                                    disp2_;
    tulpa::GaussHermite                       gh_;
};

typedef OccuCoverLatentCoupling<LognormalLatent> OccuCoverLognormalLatentCoupling;
typedef OccuCoverLatentCoupling<BetaLatent>      OccuCoverBetaLatentCoupling;

} // namespace tulpaObs

#endif // TULPAOBS_CELL_COUPLING_OCCU_COVER_LATENT_H
