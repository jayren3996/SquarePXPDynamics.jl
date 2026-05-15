# SquarePXPDynamics.jl

`SquarePXPDynamics` is a Julia package for PEPS-based dynamics on the 2D square-lattice PXP model.

## Status

The package now contains the S0-S6 prototype pipeline for square-lattice PXP
dynamics: dense local model definitions, finite and periodic PEPS/iPEPS state
containers, QR-reduced five-site star updates, deterministic Trotter evolution,
simple/local observables, and S6-lite ScarFinder orchestration.

Simple/local observables are useful diagnostics for development and regression
tests, but they are not final CTMRG-quality measurements. ScarFinder-lite
currently uses these simple/local diagnostics. Do not make physics claims from
simple diagnostics alone.

This checkout also contains PEPSKit/TensorKit-facing measurement code in
`src/PEPSKitMeasurements.jl`. It is treated as experimental S5b-facing surface
in this cleanup audit, not production ScarFinder validation.

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
- S6-lite `scarfinder!` orchestration using simple/local diagnostics (`src/ScarFinder.jl`).

## Not Yet Shipped

- Production PEPSKit/TensorKit CTMRG measurement adapter.
- CTMRG-quality observables suitable for physics claims.
- Full-update gauge fixing.
- Energy targeting or correction.
- Candidate ranking.
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
ScarFinder-lite is currently a scaffold/orchestration layer over `evolve!` and
`measure_simple`; PEPSKit CTMRG-quality measurement integration remains planned
work before production ScarFinder validation.

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
