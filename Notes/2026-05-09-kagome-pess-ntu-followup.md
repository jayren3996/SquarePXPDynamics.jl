# Kagome PESS NTU Follow-Up Plan

Date: 2026-05-09. Author: agent-assisted planning session. Companion to `2026-05-09-kagome-pxp-pivot-plan.md`.

## Why this note exists

The kagome PESS PXP work landing in PR-1 (spec: `docs/superpowers/specs/2026-05-09-kagome-pess-pxp-dynamics-design.md`, plan: `docs/superpowers/plans/2026-05-09-kagome-pess-pxp-dynamics.md`) implements **Simple Update** as the truncation kernel. SU is the published canonical kagome ground-state algorithm but its track record for real-time dynamics is sparse. The honest read is that PESS+SU is a **prototype** for the eventual NTU-based dynamics implementation.

This note records the scope of the planned NTU follow-up PR so it isn't lost between PRs. Priority and timing depend on the **Layer-2 indicator outcome** in PR-1:

- Layer-2 PASS at `1e-3`: NTU PR queued at lower priority. Do scientific scarfinder runs with the SU prototype first; upgrade to NTU when ready for production-quality results.
- Layer-2 BORDERLINE: NTU PR is the next focused work after PR-1 lands. Disclaim accuracy in any pre-NTU scientific writeup.
- Layer-2 FAIL by orders of magnitude: NTU PR is dispatched immediately. Hold all scientific claims until NTU lands and re-runs the indicator.

## What stays unchanged

The infrastructure built in PR-1 is **algorithm-agnostic** with respect to SU vs NTU. NTU replaces only the truncation kernel.

Reusable directly:

- All of `KagomeGeometry.jl` — Bravais lattice, sublattice partition, 9-site UC, neighbor relations, color schedule, `kagome_star_sites`.
- All of `KagomePESS.jl` — state container, site/simplex tensor representation, λ bookkeeping, `product_pess`, `random_pess`, accessors.
- All of `KagomeModels.jl` — 32×32 PXP star Hamiltonian, blockade projector, cluster stabilizer Hamiltonian.
- All of `KagomeGates.jl` — dense and projected 5-site gates.
- All of `KagomeSchedules.jl` — 3-coloring, first/second-order schedule layers.
- All of `KagomeObservables.jl` — local expectations and blockade diagnostics.
- All of `KagomeEvolution.jl` — schedule driver, projected_pxp_step_kagome, run loop. The driver calls `apply_star_gate_*!`; we add a `:ntu` update mode alongside `:simple`.
- All of `KagomeScarFinder.jl` — search loop, seeds, ranking. The `update::Symbol` field in `KagomeScarFinderConfig` already accepts a symbol; extend it to accept `:ntu`.
- `test/util_kagome_finite_ed.jl` — finite-cluster vector helper, torus PXP Hamiltonian, torus Z observables. The Layer-2 integration test re-runs against NTU and is the primary success criterion for the NTU PR.
- `test/test_kagome_*.jl` — Layer-1 kernel tests (parametrize over `update` mode), Layer-3 stabilizer benchmark (re-run against NTU), most other tests are kernel-agnostic.

## What changes

New module `src/KagomeNTU.jl` implementing the NTU kernel. Public function:

```julia
apply_star_gate_ntu_pess!(state::KagomePESS, gate::AbstractMatrix, center::KagomeCoord;
                          cutoff::Real, maxdim::Integer,
                          neighborhood::Symbol = :star_shell,
                          als_iterations::Integer = 50,
                          als_tolerance::Real = 1e-9) -> NTUDiagnostics
```

`KagomeEvolution.jl` dispatches on the `update` symbol:

```julia
if update === :simple
    apply_star_gate_simple_update_pess!(...)
elseif update === :ntu
    apply_star_gate_ntu_pess!(...)
else
    throw(ArgumentError("update must be :simple or :ntu"))
end
```

`KagomeScarFinderConfig` already carries `update::Symbol`; just relax the validator to accept `:ntu`.

## NTU algorithm sketch (kagome-PESS specific)

NTU = Neighborhood Tensor Update (Dziarmaga PRB 104, 094411, 2021. arXiv: [2107.06635](https://arxiv.org/abs/2107.06635)). Core idea: instead of using the mean-field environment that SU relies on, use an **exactly contractible local neighborhood metric** to define the "best" replacement tensor for each site after the gate is applied.

For kagome PESS with the 5-site PXP star gate at center `c`:

1. **Choose the neighborhood `N`** around the star. The minimal "star shell" is the 5-site star plus the 4 outermost simplex tensors (the 4 simplices that the NN sites connect to outside the star). This neighborhood has 5 sites + 6 simplex tensors = 11 tensors. The neighborhood must be exactly contractible — no infinite-lattice approximation.

2. **Form the "exact target"**: contract the absorbed cluster + gate (same as SU step 3) but DON'T immediately decompose. The target is the post-gate cluster wavefunction restricted to the neighborhood.

3. **Define the local error metric**: for a candidate set of replacement site/simplex tensors `{T̃, S̃}`, the squared distance to the exact target, weighted by the neighborhood metric, is a quadratic-in-each-tensor function. This makes alternating least squares (ALS) the natural optimizer.

4. **ALS optimization**: starting from the SU result as the initial guess (warm start), iterate:
   - Fix all tensors except one site/simplex; solve the local linear system for that tensor's optimal entries given the others.
   - Cycle through all 7 tensors of the cluster (5 sites + 2 simplices).
   - Iterate until the squared-distance residual converges or `als_iterations` cap.

5. **Truncate to fixed D**: each updated tensor is truncated by SVD on the relevant bond, with the new λ extracted as in SU.

6. **Writeback**: same as SU writeback (Task 9.5 in the SU plan).

Cost: each ALS iteration is `O(D^k)` for the neighborhood contraction, where k depends on the neighborhood size. For the kagome star-shell neighborhood, k ≈ 5-6 (much better than the `D^18` we'd hit on triangular). NTU is more expensive per gate application than SU but produces qualitatively better dynamics.

### Subtleties to watch

- **Non-uniqueness of the ALS fixed point**: the local minimum of the ALS landscape isn't unique up to gauge; convergence checks should be on the squared-distance residual, not on tensor entries.
- **Initial guess matters**: warm-starting from the SU result (which we have!) makes ALS converge in 5-15 iterations typically. Cold-start from random can take 100+ or fail to converge.
- **PESS-specific: the simplex tensors have no physical leg**, so the ALS update for a simplex tensor is purely a geometric optimization over its 3-leg layout. This is structurally different from updating a site tensor (which has phys + 2 simplex legs).
- **Lambda absorption convention**: NTU typically works with absorbed-lambda tensors throughout the ALS loop, then extracts new λ at the truncation step. Same convention as SU.

## Validation strategy for the NTU PR

Same three-layer structure as SU, with the indicator promoted back to a gate:

1. **Layer-1 kernel ED**: re-run all SU kernel tests with `update=:ntu`. Expectation: identical or better than SU (NTU should be more accurate than SU when both can be computed).
2. **Layer-2 torus integration**: **promoted back to acceptance gate.** Pass within `1e-3` absolute on per-sublattice `<Z>` at D=4, dt=0.01, 3 steps. NTU's whole purpose is to make this pass.
3. **Layer-3 stabilizer benchmark**: re-run with `update=:ntu`. Should pass within the same tolerance as SU (both are exact for commuting Hamiltonians up to truncation, and the stabilizer evolution stays in a low-D manifold).

Additional NTU-specific tests:
- ALS convergence: residual decreases monotonically until convergence; if not, log a warning.
- Warm-start dependency: results don't depend on whether we warm-start from SU or from a perturbed SU result (gauge invariance).
- Comparison vs SU: at small D where both are tractable, NTU should produce smaller `discarded_weight` and better-converged local observables.

## Estimated implementation effort

The NTU PR is roughly the same size as the SU PR (PR-1) in lines of code, dominated by the ALS optimizer and neighborhood construction. Conservative estimate: 1-2 weeks of focused work after PR-1 lands, assuming the SU infrastructure is sound. Faster if the executing agent has prior NTU experience or can adapt published reference implementations.

## Files this PR will create / modify

Create:

- `src/KagomeNTU.jl` — the new kernel.
- `test/test_kagome_ntu.jl` — kernel-agnostic Layer-1/3 re-runs + NTU-specific tests.

Modify:

- `src/KagomeEvolution.jl` — add `:ntu` dispatch.
- `src/KagomeScarFinder.jl` — relax `update` validator to accept `:ntu`.
- `src/TriangularPEPSDynamics.jl` — re-export NTU symbols.
- `test/test_kagome_evolution.jl` — promote Layer-2 to gate-test, parametrize over `update` mode.
- `README.md` — update Kagome Status to reflect NTU as the production path.

No changes needed in: `KagomeGeometry`, `KagomePESS`, `KagomeModels`, `KagomeGates`, `KagomeSchedules`, `KagomeObservables`, `KagomeSimpleUpdate` (kept as the prototype kernel for cross-checks).

## Decision log

- **Why NTU instead of full update + CTMRG?** NTU is the lightest dynamics-grade upgrade. CTMRG adds environment-contraction infrastructure that we don't otherwise need; NTU's neighborhood metric is exactly contractible without CTMRG. Full update with CTMRG is the correct end-state for accuracy but not the immediate next step.
- **Why NTU on PESS instead of NTU on plain iPEPS?** PESS captures kagome's frustration cleanly via simplex tensors; switching to plain iPEPS just to use NTU would throw away that geometric advantage and require larger D for the same accuracy. NTU on PESS is less published than NTU on iPEPS but the algorithmic adaptation is straightforward.
- **Why not BP-SU as the production kernel?** BP-SU is validated only for 2-body bond gates in the published 2D dynamics literature ([Tindall et al. 2025](https://arxiv.org/abs/2503.05693)). Adapting to the 5-body kagome PXP gate is research-grade work; NTU has a more direct precedent.
