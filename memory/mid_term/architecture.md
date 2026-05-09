# Architecture

- Confirmed: The project is a single root Julia package named `TriangularPEPSDynamics` using the root `Project.toml`. PEPS code should remain under `src/`, with tests under `test/`.
  - Source: `AGENTS.md`
  - Source: `Project.toml`
  - Source: `src/TriangularPEPSDynamics.jl`

- Confirmed: Current module layout:
  - `Geometry.jl`: axial coordinates, triangular directions, stars, and 7-coloring.
  - `SpinOps.jl`: dense spin-1/2 operators and embedding helpers.
  - `Models.jl`: PXP, blockade, cluster/stabilizer, diagonal star, and Ising helpers.
  - `Gates.jl`: dense real/imaginary gates and projected gates.
  - `Schedules.jl`: first- and second-order schedule layers.
  - `SolvableModels.jl`: exact benchmark helpers.
  - `States.jl`: iPEPS unit cells, tensors, bond indices/lambdas, seeds, and hard truncation.
  - `Observables.jl`: local expectations and blockade diagnostics.
  - `SimpleUpdate.jl`: identity/product-gate paths, `D=1` dense oracle, and current general star update path.
  - `Evolution.jl`: scheduled evolution APIs and projected PXP step diagnostics.
  - `ScarFinder.jl`: seed loop, evolve-project iteration, hard truncation, and candidate ranking.
  - Source: `src/TriangularPEPSDynamics.jl`

- Confirmed: `TriangularIPEPS` stores unit cell, physical indices, virtual bond indices, site tensors, and lambda spectra. Opposite logical bonds share lambda vectors, while one-site self-loop tensor indices remain distinct so a tensor never repeats an ITensors `Index`.
  - Source: `src/States.jl`
  - Source: `Notes/implemented_peps_algorithm_detail.md`

- Confirmed: Evolution supports both prebuilt-gate and Hamiltonian/gate-builder APIs. Hamiltonian-based second-order evolution builds half-step gates through schedule layer scales; prebuilt-gate second-order uses a matrix square root of the full-step gate.
  - Source: `src/Evolution.jl`
  - Source: `test/test_evolution.jl`

- Confirmed: `ScarFinderConfig` currently treats positional `maxdim` as `dynamics_maxdim`; `scar_maxdim` defaults to it and must not exceed it.
  - Source: `src/ScarFinder.jl`
  - Source: `test/test_scar_finder.jl`
