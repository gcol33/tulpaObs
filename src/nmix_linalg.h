// nmix_linalg.h
// Small dense linear-algebra helpers shared by the N-mixture Laplace-EM drivers
// (nmix_community_em.cpp and nmix_community_spatial.cpp). Kept in one header so
// the symmetric-inverse and log-determinant primitives have a single source
// rather than a per-driver copy.

#ifndef TULPAOBS_NMIX_LINALG_H
#define TULPAOBS_NMIX_LINALG_H

#include <RcppEigen.h>
#include <limits>

namespace tulpaObs {

// Symmetric inverse of a small PD matrix via LLT, with a tiny ridge fallback
// when the Cholesky fails (an indefinite observed-information block away from
// the mode, or a barely-singular precision).
inline Eigen::MatrixXd nmix_safe_inverse(const Eigen::MatrixXd& M) {
    const int d = static_cast<int>(M.rows());
    Eigen::LLT<Eigen::MatrixXd> llt(M);
    if (llt.info() == Eigen::Success)
        return llt.solve(Eigen::MatrixXd::Identity(d, d));
    Eigen::MatrixXd Mr = M;
    double md = M.diagonal().cwiseAbs().mean();
    if (!(md > 0)) md = 1.0;
    Mr.diagonal().array() += 1e-8 * md;
    return Mr.ldlt().solve(Eigen::MatrixXd::Identity(d, d));
}

// log|M| for a symmetric PD matrix via its Cholesky factor; NaN if not PD.
inline double nmix_logdet_spd(const Eigen::MatrixXd& M) {
    Eigen::LLT<Eigen::MatrixXd> llt(M);
    if (llt.info() != Eigen::Success)
        return std::numeric_limits<double>::quiet_NaN();
    double ld = 0.0;
    const Eigen::MatrixXd& L = llt.matrixL();
    for (int i = 0; i < L.rows(); ++i) ld += std::log(L(i, i));
    return 2.0 * ld;
}

// Top-block (p_top x p_top) marginal covariance = the top-left block of M^{-1},
// solved as M X = [I; 0]. When `constrain`, an intrinsic (rank-deficient,
// sum-to-zero) field sub-block at [field_start, field_start + field_len) is
// pinned with a large quadratic penalty kappa (the penalty-method form of the
// sum-to-zero constraint) before the inverse, so the flat (intercept,
// field-mean) direction does not blow up the marginal covariance and the
// intercept variance is the constrained data-precision value. `M` is taken by
// value so the caller's matrix is untouched. NaN block if M is not PD.
inline Eigen::MatrixXd nmix_constrained_top_cov(Eigen::MatrixXd M, int n_x,
                                                int p_top, int field_start,
                                                int field_len, bool constrain) {
    if (constrain && field_len > 0) {
        double md = M.diagonal().head(p_top).cwiseAbs().mean();
        if (!(md > 0)) md = 1.0;
        const double kappa = 1e6 * md;
        for (int i = 0; i < field_len; ++i)
            for (int j = 0; j < field_len; ++j)
                M(field_start + i, field_start + j) += kappa;
    }
    Eigen::LLT<Eigen::MatrixXd> chol(M);
    if (chol.info() != Eigen::Success)
        return Eigen::MatrixXd::Constant(
            p_top, p_top, std::numeric_limits<double>::quiet_NaN());
    Eigen::MatrixXd E = Eigen::MatrixXd::Zero(n_x, p_top);
    E.topLeftCorner(p_top, p_top).setIdentity();
    Eigen::MatrixXd Hinv_cols = chol.solve(E);
    return Hinv_cols.topLeftCorner(p_top, p_top);
}

}  // namespace tulpaObs

#endif  // TULPAOBS_NMIX_LINALG_H
