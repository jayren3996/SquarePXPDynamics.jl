# SquarePXPDynamics.jl

`SquarePXPDynamics` is a Julia package for PEPS-based dynamics on the 2D square-lattice PXP model.

## Status

The package is in an early square-lattice restart. The current code is a small, tested foundation: square geometry, dense 5-site projected PXP gates, and a minimal ITensors-backed square PEPS product-state container.

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
- Minimal square PEPS product-state construction with ITensors (`src/SquarePEPS.jl`).

## Planned Work

- Square PEPS update kernels for fixed-bond-dimension evolve-project loops.
- Blockade-violation diagnostics on square nearest-neighbor edges.
- ScarFinder orchestration and low-entanglement candidate ranking.
- ScarFinder driver.

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
