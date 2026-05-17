# SquarePXPDynamics.jl

`SquarePXPDynamics` is a Julia package for PEPS-based dynamics on the 2D square-lattice PXP model.

## Status

The package now contains the S0-S7 prototype pipeline for square-lattice PXP
dynamics plus a v1 infinite TFIM benchmark runner: dense local model
definitions, finite and periodic PEPS/iPEPS state containers, QR-reduced
five-site star updates, deterministic Trotter evolution, simple/local
observables, reproducible TFIM/PXP validation records, and ScarFinder
orchestration with explicit objectives, optional trusted finite-`chi` CTM
measurement backends, candidate metadata persistence, CTM trust/readiness
diagnostics, and transactional CTM gauge conditioning.

Simple/local observables are useful diagnostics for development and regression
tests, but they are not final CTMRG-quality measurements. ScarFinder still
defaults to this fast simple/local path for development. Physics-facing
candidate ranking should instead use `TrustedCTMBackend`, an explicit
objective such as `RevivalObjective` or `CompositeObjective`, and
`require_trusted_ctm = true`. Do not make physics claims from simple
diagnostics alone.

This checkout also contains PEPSKit/TensorKit-facing measurement code in
`src/PEPSKitMeasurements.jl` and S7a CTM trust helpers in `src/CTMTrust.jl`.
The PEPSKit CTMRG measurement adapter is shipped as an experimental S5c-facing
API and can now be selected as ScarFinder's trusted measurement backend.
Within that adapter, density,
blockade diagnostics, and five-site square-star PXP energy density use PEPSKit
CTMRG. The dense square-star Hamiltonian remains the source of truth for the
PXP energy operator, with site order `(center, right, up, left, down)` and
basis order `1 = :up`, `2 = :down`. PEPSKit and TensorKit therefore remain main
dependencies while this exported measurement surface is present.

The original S0-S7 implementation plan has been reconciled against the current
architecture in `docs/superpowers/specs/2026-05-16-s0-s7-completion-design.md`.
The current S0.5/S1 backend-facade items are superseded by the concrete custom
ITensors iPEPS stack unless a second update backend is introduced. S7b now has
CTM local norm-matrix diagnostics, readiness checks, and transactional
`fix_bond_gauge!` paths for D=1 no-op/product bonds and D>1 PEPSKit
bond-environment gauge conditioning.

## Package Layout

- `Project.toml`: package metadata, dependencies, compatibility bounds, and the test workspace.
- `src/SquarePXPDynamics.jl`: package module entrypoint.
- `src/*.jl`: implementation modules included by the entrypoint.
- `test/runtests.jl`: package test runner.
- `test/Project.toml`: test-only environment for Julia's workspace-based test dependency workflow.

## Currently shipped

- Generic spin-1/2 operators (`src/SpinOps.jl`).
- Square-lattice geometry and 5-site star scheduling helpers (`src/SquareGeometry.jl`).
- Dense square-star PXP Hamiltonian, blockade projector, and projected real/imaginary gates (`src/SquarePXP.jl`).
- Finite ITensors-backed square PEPS product-state construction (`src/SquarePEPS.jl`).
- Periodic square iPEPS product and checkerboard states in Gamma-lambda simple-update form (`src/SquareIPEPS.jl`).
- Periodic iPEPS helper APIs, link-weight normalization, and bond-entropy diagnostics (`src/SquareIPEPS.jl`).
- ITensor wrappers for dense square-star PXP gates (`src/SquareIPEPS.jl`).
- QR-reduced five-site star update with pre-update touched-link minima diagnostics via `project_star!` (`src/StarSimpleUpdate.jl`).
- Deterministic five-color and serial Trotter evolution with model metadata and log-normalization ledger diagnostics via `evolve!` (`src/IPEPSEvolution.jl`).
- Simple/local density, blockade, energy-density, and entropy observables via `measure_simple` (`src/Observables.jl`).
- Simple/local TFIM observables and reproducible JSON/CSV benchmark records via `run_benchmark` (`src/Benchmarks.jl`).
- Experimental PEPSKit/TensorKit CTMRG density, blockade, and five-site PXP energy measurement adapter via `measure_ctm` (`src/PEPSKitMeasurements.jl`).
- CTM finite-`chi` trust assessment and audit CSV output via `assess_ctm_trust` and `write_ctm_trust_csv` (`src/CTMTrust.jl`).
- PXP validation and convergence reports that compare finite ED all-down trajectories against
  matched iPEPS trajectories and optionally attach trusted finite-`chi` CTM
  measurement sweeps via `validate_pxp_ed_ipeps` and
  `write_pxp_validation_json`, or sweep `dt`, `D`, cutoff, and CTM finite-`chi`
  settings via `validate_pxp_convergence` and
  `write_pxp_convergence_json` (`src/PXPValidation.jl`).
- Read-only local simple-gauge diagnostics via `gauge_diagnostic_simple` (`src/GaugeDiagnostics.jl`).
- CTM local bond norm diagnostics, `ctm_ready_for_gauge_updates`,
  `pepskit_private_full_update_available`, and transactional `fix_bond_gauge!`
  gauge conditioning (`src/CTMGaugeReadiness.jl`).
- `scarfinder!` orchestration, guarded simple-energy correction,
  objective-based candidate ranking, optional trusted CTM measurement backends,
  JSON candidate metadata persistence, and CSV/JSON diagnostic logging
  (`src/ScarFinder.jl`).

## Not Yet Shipped

- CTM-aware/full-update evolution.
- Full physics audit reports across `dt`, `D`, `chi`, cutoff, unit cell, and
  update scheme before publication-quality ScarFinder claims.

## Minimal Example

```julia
using SquarePXPDynamics

cell = PeriodicSquareUnitCell(10, 10)
psi = product_square_ipeps(cell; state = :down, maxdim = 1)
params = TrotterParams(0.01, 1, :real, true, 1, 1e-12)

evolve!(psi, 0.01; params = params)
summary = measure_simple(psi)
```

`summary` contains simple/local diagnostics only. These are useful for smoke
tests and regression checks, but they are not CTMRG-quality measurements.

### TFIM Benchmark Smoke Run

```julia
using SquarePXPDynamics

spec = BenchmarkSpec(
    "tfim-j0",
    StaticModel(TFIMStarModel(0.0, 1.0)),
    PeriodicSquareUnitCell(3, 3),
    :z_up,
    0.02,
    TrotterParams(0.01, 1, :real, 1, 1e-12; schedule = :serial),
    1,
)

result = run_benchmark(spec; run_label = "local-smoke")
write_benchmark_json(result, "tfim-j0.json")
write_benchmark_csv([result], "tfim-j0.csv")
```

The v1 TFIM benchmark uses simple-update diagnostics. Treat these outputs as
implementation regression records, not CTMRG-quality physics estimates.
For small cells such as `3 x 3`, `run_finite_tfim_reference` provides a dense
finite-Hilbert-space TFIM reference trajectory with matching conventions.
For larger finite references, `run_finite_mps_tfim_reference` runs an
open-boundary square-lattice TFIM benchmark with a snake-MPS mapping and
ITensorMPS TDVP; `scripts/finite_mps_tfim_6x6.jl` is the 6x6 smoke script.

### PXP ED Benchmark

The package also includes an EDKit-backed finite PBC PXP benchmark path for
short-time dynamics in the fully symmetric sector. The `7 x 7` runner is:

```bash
julia --project=. scripts/pxp_ed_7x7_benchmark.jl
```

The script writes JSON by default to `scripts/pxp-ed-7x7.json`. Runtime knobs
are environment variables, for example
`PXP_ED_TOTAL_TIME=0.05 PXP_ED_M_MAX=40 julia --project=. scripts/pxp_ed_7x7_benchmark.jl`.

An experimental PEPSKit CTMRG measurement adapter is present as `measure_ctm`,
with CTMRG density, blockade, sublattice imbalance, checkerboard structure
factor, and five-site PXP energy diagnostics. Check the raw CTMRG convergence
information and finite-chi sensitivity before treating these measurements as
physics-quality observables. `measure_ctm_trusted` and
`TrustedCTMBackend` package that finite-`chi` sweep, the final CTM summary,
and the `assess_ctm_trust` result for downstream validation and ScarFinder
ranking.

```julia
params_ctm = PEPSKitCTMRGParams(8, 1e-8, 100, 0)
ctx = pepskit_ctmrg_context(psi; params = params_ctm)
energy = pxp_energy_density_ctm(psi, ctx)
diagnostics = ctm_diagnostics(ctx)
```

Inspect `diagnostics` before using CTM values for ranking, and repeat CTM
measurements at multiple `chi` values before trusting energy comparisons. A
`PEPSKitMeasurementContext` belongs to the exact state used at creation; if
`psi` is mutated by `evolve!`, `project_star!`, or link-weight setters, the old
context is stale and measurement calls throw. `ScarFinderCandidateScore.score`
is the selected objective score. For simple/local default runs it is still only
a diagnostic sorting key. ScarFinder CSV/JSON logs include objective metadata,
CTM trust fields when available, `log_norm_before`, `log_norm_after`,
`log_norm_delta`, `correction_accepted`, `correction_energy_before`, and
`correction_energy_after` so long projection sweeps and guarded correction
attempts can be screened. The public mutators update the state version; direct
edits to `psi.tensors`, `psi.link_weights`, or `psi.link_indices` are internal
mutable implementation details and can bypass cache-staleness bookkeeping.

For finite-`chi` validation sweeps, use:

```julia
points = validate_ctm_sweep(
    psi;
    params = [
        PEPSKitCTMRGParams(4, 1e-6, 50, 0),
        PEPSKitCTMRGParams(8, 1e-8, 100, 0),
    ],
)
assessment = assess_ctm_trust(points)
write_ctm_validation_csv(points, "ctm-validation.csv")
write_ctm_trust_csv(points, "ctm-trust.csv")
```

Each `CTMValidationPoint` records the CTM summary, the simple/local reference,
observable deltas, and CTMRG diagnostics for one parameter setting.
`assess_ctm_trust` compares the final finite-`chi` CTM measurements against
each other; it does not use the simple/local reference deltas stored in
`CTMValidationPoint`. A trusted assessment is a measurement-validation signal,
not permission by itself to run gauge-changing updates. `ctm_ready_for_gauge_updates`
adds the separate S7b checks for fresh contexts and local CTM bond norm
diagnostics. `fix_bond_gauge!` is transactional: D=1 product bonds are a no-op,
and D>1 bonds are conditioned with PEPSKit bond-environment factorization before
the updated tensors are written back to the Gamma-lambda iPEPS state.

### PXP validation reports

The first production-facing validation path is
`validate_pxp_ed_ipeps(PXPValidationConfig(...))`. It runs a finite periodic
PXP ED trajectory, evolves a matched all-down iPEPS trajectory on the same
unit cell, and reports density differences at shared sample times. Without CTM,
the simple-density difference is a local-environment diagnostic; for D>1 loopy
periodic PEPS it is not an exact finite contraction. Passing
`ctm_params = (...)` attaches `measure_ctm_trusted` output at every sample:
the final CTM measurement, the finite-`chi` sweep points, and the
`assess_ctm_trust` result.

For a fast JSON artifact without CTMRG:

```julia
config = PXPValidationConfig(3; total_time = 0.02, dt = 0.01)
report = validate_pxp_ed_ipeps(config; ctm_params = nothing)
write_pxp_validation_json(report, "artifacts/pxp_validation_report.json")
```

For tiny periodic cells, validation can attach an exact finite contraction
density alongside simple/local and CTM fields:

```julia
config = PXPValidationConfig(
    3;
    total_time = 0.02,
    dt = 0.02,
    maxdim = 2,
    exact_finite_observables = true,
    exact_finite_max_sites = 12,
)
report = validate_pxp_ed_ipeps(config; ctm_params = nothing)
```

The exact finite path is intentionally size-limited and uses dense `2^N`
contractions of the supplied `SquareIPEPSState`. It is a debugging and
tiny-cell validation reference, not exact ED dynamics and not a replacement for
CTM-backed thermodynamic measurements.

or from the shell:

```bash
julia --project=. scripts/validate_pxp_ed_ipeps.jl
```

This is a validation harness, not a ScarFinder ranking change. ScarFinder still
uses its existing simple/local default for fast development runs. Use the
trusted backend path below for CTM-gated ranking.

### ScarFinder trusted ranking

ScarFinder now accepts explicit objectives and measurement backends. The simple
default remains useful for smoke tests:

```julia
result = scarfinder!(
    psi;
    iterations = 3,
    params = params,
    objective = CompositeObjective(; revival = RevivalObjective()),
)
```

For physics-facing candidate ranking, use trusted finite-`chi` CTM measurement
and a hard trust gate:

```julia
backend = TrustedCTMBackend(
    [
        PEPSKitCTMRGParams(8, 1e-7, 100, 0; seed = 11),
        PEPSKitCTMRGParams(12, 1e-8, 150, 0; seed = 11),
    ],
    CTMTrustPolicy(),
)

result = scarfinder!(
    psi;
    iterations = 5,
    params = params,
    measurement = backend,
    objective = CompositeObjective(; revival = RevivalObjective()),
    require_trusted_ctm = true,
    candidate_store = JSONCandidateStore("artifacts/scarfinder-candidates"),
)
```

`JSONCandidateStore` writes metadata and score records for auditability. It
does not yet persist full tensor snapshots.

For coarse error-budget artifacts, run:

```julia
config = PXPConvergenceConfig(
    PXPValidationConfig(3; total_time = 0.02, dt = 0.01);
    dt_values = [0.02, 0.01],
    D_values = [1, 2],
    cutoff_values = [1e-10],
    chi_values = [8, 12],
)
report = validate_pxp_convergence(config)
write_pxp_convergence_json(report, "artifacts/pxp_convergence_report.json")
```

For the M1 PXP audit campaign, use `run_pxp_audit_campaign`. It runs the
small all-down ED/iPEPS validation grid, adds a reversibility report per grid
point, and writes both nested JSON and a flat CSV summary:

```julia
config = PXPAuditConfig(;
    n_values = [3],
    total_time = 0.02,
    dt_values = [0.02, 0.01],
    D_values = [1, 2],
    cutoff_values = [1e-12],
    chi_values = Int[],
)
report = run_pxp_audit_campaign(config)
write_pxp_audit_json(report, "artifacts/pxp_audit_report.json")
write_pxp_audit_csv(report, "artifacts/pxp_audit_summary.csv")
```

or from the shell:

```bash
julia --project=. scripts/pxp_audit_campaign.jl
```

Set `SQUAREPXP_AUDIT_CHI=8,12` to attach trusted CTM finite-`chi` sweeps.
Other useful overrides are `SQUAREPXP_AUDIT_N=3,4`,
`SQUAREPXP_AUDIT_DT=0.02,0.01,0.005`, `SQUAREPXP_AUDIT_D=1,2`,
`SQUAREPXP_AUDIT_CUTOFF=1e-10,1e-12`, `SQUAREPXP_AUDIT_TOTAL_TIME=0.02`,
`SQUAREPXP_AUDIT_JSON=...`, and `SQUAREPXP_AUDIT_CSV=...`.

The CSV summary is for bottleneck triage. Large `max_abs_density_error_simple`
with small reversibility drift means the simple/local no-CTM diagnostic has
separated from the ED density; for D>1 this may be a local-environment effect,
not an exact finite PEPS update error. CTM density error or rejected
`ctm_trust_status` points at finite-`chi` drift; large `max_truncerr` points at
bond-dimension/truncation pressure; large `log_norm_delta_abs` or reversibility
drifts point at persistence or round-trip stability. These are audit signals
only, not physics-grade claims.

## Development

Instantiate the package environment:

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

Load the package from the repository root:

```bash
julia --project=. -e 'using SquarePXPDynamics'
```

Run the package tests:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Run the slower extended CTM tests locally:

```bash
SQUAREPXP_EXTENDED_TESTS=1 julia --project=. -e 'using Pkg; Pkg.test()'
```

The same extended CTM suite is available in GitHub Actions through the
`Extended CTM` manual workflow and its weekly scheduled run.

The test suite includes API docstring coverage for exported names and Aqua.jl
package-quality checks.
