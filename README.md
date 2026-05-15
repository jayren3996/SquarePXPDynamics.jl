# SquarePXPDynamics.jl

`SquarePXPDynamics` is a Julia package for PEPS-based dynamics on the 2D square-lattice PXP model.

## Status

The package now contains the S0-S6 prototype pipeline for square-lattice PXP
dynamics: dense local model definitions, finite and periodic PEPS/iPEPS state
containers, QR-reduced five-site star updates, deterministic Trotter evolution,
simple/local observables, and S6-lite ScarFinder orchestration.

Simple/local observables are useful diagnostics for development and regression
tests, but they are not final CTMRG-quality measurements. ScarFinder-lite
uses these simple/local diagnostics by default, with optional scheduled CTM
diagnostics supplied by caller callbacks. Do not make physics claims from
simple diagnostics alone.

This checkout also contains PEPSKit/TensorKit-facing measurement code in
`src/PEPSKitMeasurements.jl`. The PEPSKit CTMRG measurement adapter is shipped
as an experimental S5c-facing API, not production ScarFinder validation.
Within that adapter, density, blockade diagnostics, and five-site square-star
PXP energy density use PEPSKit CTMRG. The dense square-star Hamiltonian remains
the source of truth for the PXP energy operator, with site order `(center,
right, up, left, down)` and basis order `1 = :up`, `2 = :down`.
PEPSKit and TensorKit therefore remain main dependencies while this exported
measurement surface is present.

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
- Periodic link-weight helpers and bond-entropy diagnostics (`src/SquareIPEPS.jl`).
- ITensor wrappers for dense square-star PXP gates (`src/SquareIPEPS.jl`).
- QR-reduced five-site star update via `project_star!` (`src/StarSimpleUpdate.jl`).
- Deterministic five-color Trotter evolution via `evolve!` (`src/IPEPSEvolution.jl`).
- Simple/local density, blockade, energy-density, and entropy observables via `measure_simple` (`src/Observables.jl`).
- Experimental PEPSKit/TensorKit CTMRG density, blockade, and five-site PXP energy measurement adapter via `measure_ctm` (`src/PEPSKitMeasurements.jl`).
- S6-lite `scarfinder!` orchestration, candidate ranking, and CSV/JSON diagnostic logging using simple/local diagnostics by default (`src/ScarFinder.jl`).

## Not Yet Shipped

- Production CTMRG convergence policy and finite-chi validation suitable for physics claims.
- Full-update gauge fixing.
- Energy targeting or correction.
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
An experimental PEPSKit CTMRG measurement adapter is present as `measure_ctm`,
with CTMRG density, blockade, and five-site PXP energy diagnostics. Check the
raw CTMRG convergence information and finite-chi sensitivity before treating
these measurements as physics-quality observables.
ScarFinder-lite is currently a scaffold/orchestration layer over `evolve!` and
`measure_simple`, with optional callback-supplied CTM diagnostics for selected
iterations; production CTMRG-quality ScarFinder validation is not yet shipped.

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

The test suite includes API docstring coverage for exported names and Aqua.jl
package-quality checks.
