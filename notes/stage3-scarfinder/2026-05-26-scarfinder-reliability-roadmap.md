# ScarFinder reliability roadmap — status

Single point of entry for the 2026-05-26 reliability work. Five phases turn
ScarFinder from "infrastructure shipped" into "audited and publishable".

## Phase status

| # | Phase                                          | Status     | Artifacts                                                                              |
|---|------------------------------------------------|------------|----------------------------------------------------------------------------------------|
| 1 | CTM throughput timing matrix + recipe          | shipped    | `scripts/run_ctm_timing_matrix.jl`, [recipe](2026-05-26-ctm-throughput-recipe.md)      |
| 2 | CTM trust calibration framework + thresholds   | shipped    | `scripts/run_ctm_chi_sensitivity.jl`, [calibration](2026-05-26-ctm-trust-calibration.md) |
| 3 | Tensor snapshot persistence (JLD2)             | shipped    | `src/CandidateSnapshots.jl`, `JLD2CandidateStore`                                      |
| 4 | ScarFinder audit harness                       | shipped    | `src/ScarFinderAudit.jl`, `scripts/run_scarfinder_audit.jl`                            |
| 5 | First audit campaign + acceptance verdict      | **shipped** | [acceptance verdict](2026-05-26-scarfinder-acceptance.md) — all 6 criteria pass; v2 with two-D primitive shows compression cost ~3e-4 at D_target=2 |
| 6 | Two-bond-dimension primitive (`compress_to_target_maxdim!`, `target_maxdim` in ScarFinderParams, two-D audit) | **shipped** | `src/IPEPSCompression.jl`, `test/test_ipeps_compression.jl` (9 testsets) |

"Shipped" means: code + tests + smoke pass. "In flight" means: framework in
place, the actual production data run is launched or pending.

## What "reliable" requires

A ScarFinder candidate trajectory is **reliable** when:

1. CTM measurements at the chosen `chi` pass a trust policy whose numerical
   thresholds were calibrated on a representative state (Phase 2).
2. The candidate ranking is stable across `(dt, D, chi, cutoff)` variation
   per the Phase 4 stability summary thresholds in
   [scarfinder acceptance](2026-05-26-scarfinder-acceptance.md).
3. The candidate's tensor snapshot is saved (Phase 3) and the run is
   reproducible.
4. The compute throughput is sufficient to make 1-3 tractable (Phase 1).

## How a user reproduces a verdict

```bash
# 0. Threading is fixed for this server (see Phase 1 recipe).
export JULIA_NUM_THREADS=64
export SQUAREPXP_CTM_BLAS_THREADS=1
export SQUAREPXP_CTM_STRIDED_THREADS=64
export SQUAREPXP_CTM_STRIDED_THREADED_MUL=true
export SQUAREPXP_CTM_PEPSKIT_SCHEDULER=dynamic

# 1. Per-D sensitivity sweep at the audit's evolution time. Repeat for each D.
SQUAREPXP_SENS_D=2 SQUAREPXP_SENS_EVOLVE_TIME=0.5 \
SQUAREPXP_SENS_CHI_VALUES=8,16,32,64 \
SQUAREPXP_SENS_OUTPUT_JSON=artifacts/v1/ctm_sensitivity_d2_t05.json \
julia --project=. scripts/run_ctm_chi_sensitivity.jl

# 2. Audit campaign using the sensitivity output as the trust calibration.
SQUAREPXP_AUDIT_LABEL=v1-rev-density-d2 \
SQUAREPXP_AUDIT_INITIAL_STATE=down \
SQUAREPXP_AUDIT_PROJECTION_TIME=0.0625 SQUAREPXP_AUDIT_ITERATIONS=8 \
SQUAREPXP_AUDIT_DT=0.02,0.01 \
SQUAREPXP_AUDIT_TARGET_D=2 SQUAREPXP_AUDIT_EVOLVE_D=4 \
SQUAREPXP_AUDIT_CHI=16,32 \
SQUAREPXP_AUDIT_OBJECTIVE=revival_density \
SQUAREPXP_AUDIT_TRUST=calibrated \
SQUAREPXP_AUDIT_CALIBRATION_JSON=artifacts/v1/ctm_sensitivity_d2_t05.json \
julia --project=. scripts/run_scarfinder_audit.jl
```

3. Read `artifacts/scarfinder_audit_v1-rev-density-d2.json` and check the
   acceptance criteria in [scarfinder acceptance](2026-05-26-scarfinder-acceptance.md).
   If all criteria pass for at least one row, the row's best iteration is
   reliable. Reload the snapshot via
   `load_square_ipeps_snapshot("artifacts/scarfinder-candidates/candidate_NNNNNN.jld2")`.

## Things explicitly *not* in this roadmap

- **CTM-aware/full-update evolution.** Simple-update is the only evolver.
  This is the next-tier reliability concern (better truncation behavior at
  large `D`), but it doesn't block the v1 verdict.
- **`:z_up` and `:checkerboard_*` initial states.** Both crash the
  simple-update star projection at t=0 on this lattice. `:down` is the
  working starting point. Fix needs a redesign of the star-split spectrum
  guard, not in scope here.
- **Additional CTM observables.** `correlator_ctm` and
  `correlation_length_ctm` shipped earlier; return/fidelity proxy and
  energy-variance-quality remain future scope.
- **Multi-objective ranking.** Single-objective audit only.

## Open items after the v2 campaign

The v2 campaign (2026-05-27) introduced the two-bond-dimension primitive
(separate `evolve_maxdim` vs `target_maxdim`) that defines what a ScarFinder
"candidate" actually is — a state that compresses cheaply from the faithful
evolution manifold. Compression cost is now the scar-quality signal. At
`projection_time=0.08, iterations=6, D_evol=4`, target_D ∈ {2,3,4} all pass
all six acceptance criteria, with compression truncerr falling 33× per
unit of D_target. CTM ranking is independent of D_target to 1e-5.

**Highest-priority follow-ups:**
- Longer-time campaigns: hold `D_evol=4` and ramp `projection_time × iterations`
  from 0.48 to 1.0+ to find the time at which D_target=2 stops capturing the
  trajectory (compression truncerr will rise).
- Add `compression_cost_weight` to `CompositeObjective` so scar-quality enters
  the candidate ranking directly, not just as a diagnostic.
- D_evol scaling: at `D_evol=6` and `8`, does the D_target=2 candidate change?
  If yes, evolution is undersampled at D_evol=4. If no, D_evol=4 is sufficient.
- D=3 sensitivity calibration before claims about other initial states or
  longer evolution times.

**Lower-priority follow-ups:**
- Variational compression vs the current per-bond SVD (potentially tighter
  truncerr at the same D_target for highly entangled states).
- CTM-aware / full-update evolution (only if the D_target=2 candidate
  diverges from larger D_target at long times).
- Record the handoff in the relevant `notes/` file once the user explicitly requests it.
