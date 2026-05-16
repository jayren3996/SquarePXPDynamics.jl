# SquarePXPDynamics.jl

`SquarePXPDynamics` is a Julia package for PEPS-based dynamics on the 2D square-lattice PXP model.

## Status

The package now contains the S0-S6 prototype pipeline for square-lattice PXP
dynamics plus a v1 infinite TFIM benchmark runner: dense local model
definitions, finite and periodic PEPS/iPEPS state containers, QR-reduced
five-site star updates, deterministic Trotter evolution, simple/local
observables, reproducible TFIM benchmark records, and ScarFinder
orchestration with optional guarded simple-energy correction.

Simple/local observables are useful diagnostics for development and regression
tests, but they are not final CTMRG-quality measurements. ScarFinder-lite
uses these simple/local diagnostics by default, with optional scheduled CTM
diagnostics supplied by caller callbacks. Do not make physics claims from
simple diagnostics alone.

This checkout also contains PEPSKit/TensorKit-facing measurement code in
`src/PEPSKitMeasurements.jl` and S7a CTM trust helpers in `src/CTMTrust.jl`.
The PEPSKit CTMRG measurement adapter is shipped as an experimental S5c-facing
API, not production ScarFinder validation. Within that adapter, density,
blockade diagnostics, and five-site square-star PXP energy density use PEPSKit
CTMRG. The dense square-star Hamiltonian remains the source of truth for the
PXP energy operator, with site order `(center, right, up, left, down)` and
basis order `1 = :up`, `2 = :down`. PEPSKit and TensorKit therefore remain main
dependencies while this exported measurement surface is present.

The original S0-S7 implementation plan has been reconciled against the current
architecture in `docs/superpowers/specs/2026-05-16-s0-s7-completion-design.md`.
The current S0.5/S1 backend-facade items are superseded by the concrete custom
ITensors iPEPS stack unless a second update backend is introduced. S7b now has
CTM local norm-matrix diagnostics, readiness checks, and a transactional D=1
product/no-op `fix_bond_gauge!` path; D>1 mutating gauge conditioning remains
the outstanding S7 work.

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
- Deterministic five-color Trotter evolution with model metadata and log-normalization ledger diagnostics via `evolve!` (`src/IPEPSEvolution.jl`).
- Simple/local density, blockade, energy-density, and entropy observables via `measure_simple` (`src/Observables.jl`).
- Simple/local TFIM observables and reproducible JSON/CSV benchmark records via `run_benchmark` (`src/Benchmarks.jl`).
- Experimental PEPSKit/TensorKit CTMRG density, blockade, and five-site PXP energy measurement adapter via `measure_ctm` (`src/PEPSKitMeasurements.jl`).
- CTM finite-`chi` trust assessment and audit CSV output via `assess_ctm_trust` and `write_ctm_trust_csv` (`src/CTMTrust.jl`).
- Read-only local simple-gauge diagnostics via `gauge_diagnostic_simple` (`src/GaugeDiagnostics.jl`).
- CTM local bond norm diagnostics, `ctm_ready_for_gauge_updates`, and the D=1 product/no-op `fix_bond_gauge!` path (`src/CTMGaugeReadiness.jl`).
- `scarfinder!` orchestration, guarded simple-energy correction, candidate ranking, and CSV/JSON diagnostic logging using simple/local diagnostics by default (`src/ScarFinder.jl`).

## Not Yet Shipped

- D>1 mutating full-update gauge conditioning.
- Production ScarFinder validation.

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
with CTMRG density, blockade, and five-site PXP energy diagnostics. Check the
raw CTMRG convergence information and finite-chi sensitivity before treating
these measurements as physics-quality observables.
ScarFinder is currently an orchestration layer over `evolve!` and
`measure_simple`, with optional guarded simple-energy correction and
callback-supplied CTM diagnostics for selected iterations. Production
CTMRG-quality ScarFinder validation is not yet shipped.

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
is a diagnostic sorting key, not a physics-quality energy target. ScarFinder
CSV/JSON logs include `log_norm_before`, `log_norm_after`, `log_norm_delta`,
`correction_accepted`, `correction_energy_before`, and
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
diagnostics. The current `fix_bond_gauge!` path is transactional and no-op for
D=1 product bonds; D>1 mutating gauge conditioning is still deferred.

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
