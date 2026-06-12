# tulpaObs#94 - rank-1 all-undetected (p,p) cross-Hessian

## Problem

In `occu()` / `occu_cover()` / `occu_multiscale_cover()` joint_coupled, the
all-undetected occupancy-mixture branch (`occu_nodet_block` ->
`nodet_mixture_block`, `src/occu_coupling_shared.h`) wrote a **dense V x V
detection cross-Hessian** per site (V = visits at that site), where the block is
analytically the rank-1 `a * p p^T`. Occupancy data is sparse, so almost every
site hits the no-detection branch, and the per-iteration Gauss-Newton scatter was
**O(sum_s V_s^2)** -- the profiled bottleneck (~99.9% of inner-solve wall time at
EVA scale).

## Why it is rank-1, and why it collapses

`L = w P0 + (1 - w)`, `P0 = prod_v (1 - p_v)`: the density depends on the visit
etas only through the scalar `P0`, so every cross-visit second derivative factors
through it. The off-diagonal (p,p) block is exactly `a p p^T`,
`a = -w(1-w) P0 / L^2`.

Each detection row's eta is linear in the joint dofs via its chain
`c_v = chain(row_v)` (design row -> beta_p dofs). The dense scatter therefore
adds `sum_{v != w} a p_v p_w c_v c_w^T`, which collapses to a single rank-1 in
dof space:

    a u u^T,   u = sum_v p_v c_v.

For the targeted families the p arm chains to the (few) shared beta_p dofs, so
`a u u^T` lands inside the already-dense beta_p x beta_p block -- the same dof
pairs the dense path wrote, so no new Hessian nonzeros, and the per-site cost
drops from O(V^2 * p_beta^2) to O(V * p_beta) + O(p_beta^2).

Diagonal bookkeeping: the eta-space block is
`diag(nh_p) + a (p p^T - diag(p^2))`. To use the full rank-1, the block folds the
rank-1's own diagonal into the stored diagonal (`nh_p[v] -= a p_v^2`) and emits
the full `a p p^T`; the diagonal-row scatter + `a u u^T` then reproduce the exact
block. For V = 1 the fold and the rank-1 exactly cancel, recovering the dense
result.

## Implementation

- **tulpa (engine)**: `CellDerivs` gains an optional per-arm rank-1 self-cross
  descriptor (`arm_cross_rank1_coef`, `arm_cross_rank1_vec`, appended last;
  `inst/include/tulpa/cell_coupling.h`). `scatter_cell_coupling_branch_impl`
  (`src/nested_laplace_joint_multi.h`) provides the per-cell scratch and, for a
  self block with a non-zero coefficient, accumulates `u` over the arm's rows
  (`accumulate_self_rank1_u`, merged by dof) and scatters the symmetric
  `coef * u u^T` (`scatter_self_rank1_{dense,sparse}`) in place of the O(V^2)
  dense `arm_cross_hess[k][k]` loop. Cross-arm blocks and the dense path are
  unchanged; the descriptor is honoured only on the single-response path
  (`n_batch_ == 1`).
- **tulpaObs**: `nodet_mixture_block` emits `(a, p)` into the descriptor and
  folds `nh_p` when the engine supplies the buffers (else byte-identical dense
  `cross_p_p`); `occu_nodet_block` selects the rank-1 path on the single-response
  path. Shared by the `occu_only`, `occu_cover`, and `occu_cover_latent` specs.

## Scope

The single-mixture (p,p) block (occu_only / occu_cover / occu_cover_latent) is
done. `occu_multiscale_cover` keeps the dense path: its no-detection (p,p) block
is a nested mixture (two global rank-1 terms plus one per non-detected plot), not
a single rank-1, and it is not the profiled bottleneck. A multi-term descriptor
would be the clean follow-up if that family ever dominates.

## Validation

- Dense path unchanged: the FD coupling tests (`test-occu-only-coupling.R`,
  `test-occu-cover-coupling.R`, `test-occu-multiscale-cover-coupling.R`) all pass
  -- the rank-1 emission is opt-in (engine-supplied buffers), and the FD test
  harness constructs `CellDerivs` without them, so it stays on the dense path.
- Equivalence (the definitive check): same deterministic occu_cover joint_coupled
  fit under the OLD installed code (dense) vs the NEW code (rank-1) agrees to FP
  tolerance on coef means, SDs, the spatial field, and logLik (lognormal + beta
  arms). The Newton mode is gradient-driven and the gradients are untouched; only
  the Hessian summation order changes. Both engine branches are covered:
  **dense Newton** (`n_x < 200`, N=48) to <= 1.4e-14 (`_ab_94.R`), and **sparse
  Newton** (`n_x >= 200`, N=260) to <= 2.3e-13 (`_ab_sparse_94.R`).
- Whole-suite smoke (`TULPAOBS_FAST = 1`): 0 failures, 2396 assertions.
- occu() SVC joint_coupled recovery (`test-occu-svc-joint-recovery.R`, the
  issue's profiled family, via the `occu_only` spec): full recovery passes.

## Perf (`_perf_94.R`, occu_cover joint_coupled, N=250, max.iter=8)

tulpa phase profiler, `scatter` phase (ms), old (dense) vs new (rank-1):

    J     old scatter   new scatter   scatter speedup
    8        23,498        4,008          5.9x
    16       44,033        4,648          9.5x
    32      109,108        6,245         17.5x

Old scatter grows super-linearly in J (the per-site V^2 term); new scatter is
nearly flat (~4-6k ms, the residual non-(p,p) cost), confirming the
O(sum_s V_s^2) -> O(sum_s V_s) reduction. The speedup widens with visits per
site, so at EVA scale (J_cap up to ~40, where scatter was the profiled ~99.9%
bottleneck) the win is larger still. Wall-clock of the whole capped fit went
from 37.7s -> 11.3s (J=8) to 123.8s -> 13.5s (J=32).
