# Architecture

## Package Shape

- Confirmed: The package is named `SquarePXPDynamics` and targets Julia `1.12`.
- Confirmed: Main dependencies are `ITensors`, `LinearAlgebra`, `PEPSKit`, and
  `TensorKit`; tests additionally use `Aqua` and `Test`.
- Source: `Project.toml`
- Source: `test/Project.toml`

## Module Map On Current Main

- `src/SpinOps.jl`: dense spin-1/2 operators and Kronecker embedding helpers.
- `src/SquareGeometry.jl`: coordinates, directions, nearest neighbors,
  square-star sites, and five-color scheduling.
- `src/SquarePXP.jl`: dense PXP star Hamiltonian, blockade projector, and
  projected/unprojected real/imaginary gates.
- `src/SquarePEPS.jl`: finite ITensors-backed square PEPS product states.
- `src/SquareUnitCells.jl`: periodic rectangular unit cells, update centers,
  bond keys, and compatibility checks.
- `src/SquareIPEPS.jl`: periodic Gamma-lambda iPEPS states, link weights,
  mutation/version tracking, log normalization, and dense-gate ITensor wrappers.
- `src/StarSimpleUpdate.jl`: QR-reduced five-site star update through
  `project_star!`.
- `src/IPEPSEvolution.jl`: Trotter parameters, five-color evolution, and
  `EvolutionLog` diagnostics.
- `src/Observables.jl`: simple/local observables and `measure_simple`.
- `src/PEPSKitMeasurements.jl`: experimental PEPSKit/TensorKit CTMRG
  measurement adapter, diagnostics, validation sweeps, and CSV output.
- `src/ScarFinder.jl`: ScarFinder-lite orchestration, ranking, and CSV/JSON
  logs.
- Source: `src/SquarePXPDynamics.jl`
- Source: `src/*.jl`
- Source: `README.md`

## Data And Control Flow

- Confirmed: Dense local gates are constructed first, converted to ITensors,
  applied by `project_star!`, orchestrated by `evolve!`, then measured by
  simple/local or CTM-backed measurement functions.
- Confirmed: ScarFinder sits above this stack and should call evolution and
  measurement APIs rather than manipulating low-level tensor indices.
- Source: `src/SquarePXPDynamics.jl`
- Source: `src/StarSimpleUpdate.jl`
- Source: `src/IPEPSEvolution.jl`
- Source: `src/ScarFinder.jl`

## TFIM Feature Branch Architecture

- Confirmed: Branch `codex/infinite-tfim-benchmark` adds `src/StarModels.jl`
  and `src/Benchmarks.jl`, threads model protocols through star update and
  evolution, adds TFIM observables, and exports benchmark runner APIs.
- Confirmed: Current `main` does not yet include these branch files at the time
  this memory was created.
- Source: `docs/superpowers/plans/2026-05-15-infinite-tfim-benchmark.md`
- Source: `docs/superpowers/notes/2026-05-15-current-work-infinite-tfim-benchmark.md`
- Source: `git status`
