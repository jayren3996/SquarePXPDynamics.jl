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
  helper APIs, copy/normalization helpers, mutation/version tracking, log
  normalization, and dense-gate ITensor wrappers.
- `src/GaugeDiagnostics.jl`: read-only simple-gauge bond diagnostics.
- `src/StarModels.jl`: static PXP/TFIM star-model protocol, Hamiltonian, gate,
  and convention helpers.
- `src/StarSimpleUpdate.jl`: QR-reduced five-site star update through
  `project_star!`, including touched-link and split diagnostics.
- `src/IPEPSEvolution.jl`: Trotter parameters, five-color/serial evolution,
  model metadata, and `EvolutionLog` diagnostics.
- `src/Observables.jl`: simple/local PXP and TFIM observables.
- `src/PEPSKitMeasurements.jl`: experimental PEPSKit/TensorKit CTMRG
  measurement adapter, diagnostics, validation sweeps, and CSV output.
- `src/CTMTrust.jl`: finite-chi CTM trust policy and trust CSV output.
- `src/CTMGaugeReadiness.jl`: S7b CTM bond norm diagnostics, gauge-readiness
  checks, and transactional `fix_bond_gauge!`.
- `src/Benchmarks.jl`: simple-update benchmark runner and JSON/CSV writers.
- `src/FiniteTFIMReference.jl`: dense small-cell finite TFIM reference.
- `src/FiniteMPSTFIMReference.jl`: open-boundary finite MPS TFIM reference.
- `src/FinitePXPEEDBenchmark.jl`: EDKit-backed finite PBC PXP benchmark.
- `src/ScarFinder.jl`: ScarFinder orchestration, guarded simple-energy
  correction, ranking, and CSV/JSON logs.
- Source: `src/SquarePXPDynamics.jl`
- Source: `src/*.jl`
- Source: `README.md`

## Data And Control Flow

- Confirmed: Dense local gates are constructed first, converted to ITensors,
  applied by `project_star!`, orchestrated by `evolve!`, then measured by
  simple/local or CTM-backed measurement functions.
- Confirmed: ScarFinder sits above this stack and should call evolution and
  measurement APIs rather than manipulating low-level tensor indices.
- Confirmed: CTM-backed measurement contexts are version-guarded. S7a trust
  consumes CTM sweep records, while S7b gauge conditioning additionally
  requires fresh contexts, finite-chi trust, and local CTM bond norm diagnostics
  before mutating the Gamma-lambda state.
- Source: `src/SquarePXPDynamics.jl`
- Source: `src/StarSimpleUpdate.jl`
- Source: `src/IPEPSEvolution.jl`
- Source: `src/ScarFinder.jl`
- Source: `src/PEPSKitMeasurements.jl`
- Source: `src/CTMGaugeReadiness.jl`

## Benchmark And Reference Paths

- Confirmed: The v1 infinite TFIM benchmark uses the same square-star update
  machinery and records simple-update diagnostics as reproducible regression
  data, not CTMRG-quality physics estimates.
- Confirmed: Finite references now include dense periodic TFIM for small cells,
  open-boundary MPS TFIM for larger finite comparisons, and EDKit-backed PBC
  PXP dynamics in the fully symmetric sector.
- Source: `docs/superpowers/plans/2026-05-15-infinite-tfim-benchmark.md`
- Source: `docs/superpowers/notes/2026-05-15-current-work-infinite-tfim-benchmark.md`
- Source: `README.md`
- Source: `src/FiniteTFIMReference.jl`
- Source: `src/FiniteMPSTFIMReference.jl`
- Source: `src/FinitePXPEEDBenchmark.jl`
