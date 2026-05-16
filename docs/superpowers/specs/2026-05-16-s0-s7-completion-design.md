# S0-S7 Completion Design

Date: 2026-05-16

## Goal

Fulfill the original S0-S7 project trajectory without destabilizing the
verified custom-ITensors iPEPS stack. This means preserving the current
working architecture, closing implementation gaps that still matter for
ScarFinder and S7 full-update readiness, and explicitly superseding plan items
that no longer match the chosen architecture.

## Current Baseline

The current `main` branch has a clean S0-S6 prototype plus S7a CTM trust
readiness:

- dense square-star PXP source-of-truth gates and Hamiltonians,
- finite PEPS and periodic Gamma-lambda iPEPS containers,
- QR-reduced five-site star simple update,
- deterministic five-color and serial iPEPS evolution,
- simple/local PXP and TFIM diagnostics,
- experimental PEPSKit CTMRG measurement adapter,
- CTM finite-chi trust policy and read-only gauge diagnostics,
- S6-lite ScarFinder orchestration, ranking, and CSV/JSON logging.

The focused S0-S7 baseline command passed in an isolated worktree:

```bash
julia --project=. test/runtests.jl \
  test_square_unitcells.jl \
  test_square_ipeps.jl \
  test_square_ipeps_s2.jl \
  test_star_simple_update.jl \
  test_ipeps_evolution.jl \
  test_observables_evolved.jl \
  test_ctm_trust.jl \
  test_gauge_diagnostics.jl \
  test_scarfinder.jl \
  test_public_docs.jl
```

Result:

```text
SquarePXPDynamics | 63540  63540
```

## Completion Policy

The original multistage plan is a source of intent, not a requirement to add
obsolete abstractions. Remaining work is classified as:

- **Complete now:** directly improves reproducibility, diagnostics,
  ScarFinder control, or S7b readiness.
- **Defer:** required for full-update/gauge-fixed CTMRG, but depends on an
  intermediate S7b feasibility or norm-matrix layer.
- **Supersede:** contradicted by an active architectural decision and not
  needed for current ScarFinder or benchmark workflows.

All new production behavior must be test-first. Documentation-only
reconciliation can be committed after Markdown review and `git diff --check`.

## Milestone Decisions

### S0: Baseline Locked

Status: complete.

The baseline conventions are locked by source, tests, and current verification:

- `:up` / basis index `1` is excited/Rydberg.
- `:down` / basis index `2` is vacancy/unexcited.
- star order is `(center, right, up, left, down)`.
- direction order is `:right`, `:up`, `:left`, `:down`.
- finite `SquarePEPSState` remains separate from periodic iPEPS machinery.

Remaining action:

- Add a short milestone reconciliation note so future agents do not interpret
  stale S0 notes as unfinished work.

### S0.5: Backend Decision Spike

Status: complete by decision, not by the original facade design.

The project chose a custom ITensors iPEPS storage/update path plus PEPSKit for
experimental CTMRG measurement. The planned `AbstractProjectionBackend` facade
is superseded unless a second update backend is actually introduced.

Complete now:

- Document that `src/IPEPSBackends.jl`, `ProjectionBackend`, and
  `projection_backend` are superseded by the current concrete
  `SquareIPEPSState` architecture.

Do not implement now:

- A speculative backend abstraction with only one concrete backend.

### S1: Periodic Unit Cells And iPEPS State

Status: functionally complete with small public-helper gaps.

Complete now:

- Add documented convenience helpers that expose existing state facts without
  changing storage:
  - `unitcell_reps(psi::SquareIPEPSState)`,
  - `physical_dim(psi, c::SquareCoord)`,
  - `simple_weight_dim(psi, c::SquareCoord, dir::Symbol)`,
  - `copy_state(psi::SquareIPEPSState)`.

Rationale:

These helpers match the original plan, are cheap, and improve downstream code
clarity without reviving the abandoned backend facade.

Supersede:

- `projection_backend(psi) isa AbstractProjectionBackend`.
- Parametric `SquareIPEPSState{B,S,W}`.

### S2: Backend Gate Conversion And Link Weights

Status: mostly complete for the custom ITensors path.

Complete now:

- Add `normalize_link_weights!(psi)` as a public mutating invariant helper.
- Add tests proving it normalizes every stored lambda, rejects invalid
  all-zero spectra, increments the state version when it changes data, and
  leaves already normalized spectra stable.
- Add a PXP `:x_plus` five-site dense-reference observable regression.

Defer:

- Public PEPSKit/TensorKit `square_pxp_gate_tensormap` and
  `square_pxp_star_localoperator` helpers, unless a concrete CTM or external
  API need appears.

Supersede:

- Separate `src/LinkWeights.jl` and `src/SquarePXPGates.jl` files. The current
  local placement in `SquareIPEPS.jl` is acceptable until file size or reuse
  pressure justifies a split.

### S3: QR-Reduced Five-Site Star Simple Update

Status: complete with diagnostic tightening.

Complete now:

- Extend star diagnostics so lambda minima distinguish:
  - pre-update touched bond minima,
  - post-split new center-leaf minima.

Rationale:

The original plan asked for minimum lambda information for every touched bond.
The current `min_lambda` field records only new split spectra. Adding a new
field preserves compatibility while making diagnostics more faithful.

### S4: iPEPS Evolution Driver

Status: complete with reproducibility metadata gap.

Complete now:

- Add model/protocol metadata to `EvolutionLog` or a stable companion
  diagnostic so raw evolution logs can identify:
  - schedule,
  - model family,
  - PXP projected/unprojected choice when applicable.

Rationale:

The model-aware protocol layer intentionally removed `projected` from current
`TrotterParams`; logging should still make runs reproducible.

### S5: Observables And Diagnostics

Status: complete for explicit simple/CTM APIs; original unified facade remains
unimplemented.

Complete now:

- Add a narrow `measure(psi; backend = :simple | :ctm, params = ...)` wrapper
  only if it stays thin and explicit.
- Keep `measure_simple`, `measure_ctm`, and named local observables as the
  primary APIs.
- Add missing dense-reference tests for PXP `:x_plus` simple observables.

Supersede:

- `refresh_environment!`, `current_environment`, and broad backend-dispatched
  `expectation` APIs for S5. They blur CTM context freshness and are not needed
  for the current explicit CTM adapter.

### S6: ScarFinder Orchestration

Status: complete after Slice 3. Defaults preserve S6-lite behavior; guarded
simple-energy correction is opt-in and diagnostic-only.

Completed in Slice 3:

- Add opt-in energy-correction parameters and result fields without changing
  the default S6-lite behavior:
  - `target_energy::Union{Nothing,Float64}`,
  - `correction_time`,
  - `correction_attempts`,
  - `correction_accepted`,
  - `correction_energy_before`,
  - `correction_energy_after`.
- Initial correction uses existing simple/local PXP energy only and must be
  documented as diagnostic, not physics-quality.
- When correction is enabled, it must never silently worsen the selected
  diagnostic objective. A rejected correction keeps the pre-correction state or
  explicitly records rejection.

Defer:

- CTM-trusted energy-targeted ScarFinder ranking until S7b supplies trusted
  gauge/environment readiness.

### S7: CTMRG And Gauge-Fixed Full Update

Status: S7a complete. S7b Slice 4 now provides CTM local norm diagnostics,
readiness checks, and a transactional D=1 product/no-op gauge path. D>1
mutating gauge conditioning remains for Slice 5.

Complete as S7b:

1. Add CTM local bond norm-matrix diagnostics. (Slice 4 complete.)
2. Add a readiness predicate that combines fresh CTM context, CTM trust
   assessment, and local norm-matrix quality.
   (Slice 4 complete.)
3. Add a transactional `fix_bond_gauge!` no-op/product-state path first.
   (Slice 4 complete.)
4. Extend to D>1 gauge conditioning only after norm matrices are validated.
   (Slice 5 remaining.)

Required diagnostics:

- finite entries,
- Hermiticity residual,
- positive-semidefinite eigenvalue floor,
- condition number or reciprocal condition number,
- bond-direction coverage,
- stale-context rejection,
- no partial mutation on failure.

Tests must compare observables and CTM summaries, not raw tensor entries.

## Implementation Slices

### Slice 1: Reconciliation And Small Public Helpers

Deliver:

- milestone reconciliation doc,
- S1 helper APIs,
- `normalize_link_weights!`,
- PXP `:x_plus` observable regression,
- decision log updates.

This slice has low scientific risk and clears stale S0-S2 ambiguity.

### Slice 2: Diagnostics And Reproducibility

Deliver:

- expanded star update lambda diagnostics,
- evolution model/protocol metadata,
- tests and README notes.

This slice strengthens auditability without changing update math.

### Slice 3: S5/S6 User-Facing Workflow

Deliver:

- optional thin `measure` wrapper if tests show it improves ergonomics,
- guarded simple-energy correction in ScarFinder,
- logs for accepted/rejected correction decisions.

This slice completes the non-CTM ScarFinder plan while preserving S6-lite
defaults. Slice 3 was completed with guarded simple-energy correction fields
and CSV/JSON logging; the optional thin `measure` wrapper was not added because
the explicit `measure_simple`/`measure_ctm` boundary remains clearer.

### Slice 4: S7b Feasibility And Readiness

Deliver:

- CTM local bond norm-matrix design/probe,
- read-only norm diagnostics,
- `ctm_ready_for_gauge_updates`,
- product/no-op `fix_bond_gauge!` path.

This slice should not attempt full D>1 gauge mutation until the PEPSKit
environment path is proven. Slice 4 was completed with PEPSKit `bondenv_fu`
norm diagnostics, readiness results, and the D=1 product/no-op
`fix_bond_gauge!` path.

### Slice 5: S7b Gauge Conditioning

Deliver:

- D>1 gauge conditioning,
- transactional mutation,
- stale-context invalidation,
- gauge-invariant simple and CTM regression tests,
- extended CTM verification.

This is the first slice that should claim full original S7 completion.

## Verification Gates

Every slice must run:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Any S7b slice must also run:

```bash
SQUAREPXP_EXTENDED_TESTS=1 julia --project=. test/runtests.jl test_pepskit_measurements.jl
```

Before claiming full S0-S7 completion, run a completion audit that maps every
original S0-S7 item and every supersession/defer decision to source, tests, or
documentation evidence.
