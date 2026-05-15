# Infinite TFIM Benchmark Framework Design

Date: 2026-05-15

## Goal

Build a benchmark framework that reuses the current square-star iPEPS evolution machinery to run infinite-system Transverse Field Ising Model (TFIM) benchmarks. The first implementation should keep the existing PXP behavior intact, add TFIM as a second star-gate model, and produce reproducible benchmark records from simple-update diagnostics. CTM-quality comparison is a later validation tier, not a first-version physics claim.

## Context

The repository currently targets square-lattice PXP dynamics with PEPS/iPEPS methods. The relevant existing pieces are:

- Dense five-site square-star PXP Hamiltonian and gates in `src/SquarePXP.jl`.
- Dense square-star gate to ITensor conversion in `src/SquareIPEPS.jl`.
- QR-reduced five-site simple update in `src/StarSimpleUpdate.jl`.
- Deterministic five-color Trotter evolution in `src/IPEPSEvolution.jl`.
- Simple/local observables and experimental PEPSKit CTM measurement hooks.

The current update path is structurally suitable for any dense five-site star gate, but `project_star!` still chooses PXP gates directly. The benchmark framework should separate "which local model supplies the gate" from "how a five-site star gate is applied to the iPEPS state."

## Literature Anchors

Primary benchmark model:

- 2D TFIM on the square lattice, because it has strong published infinite-iPEPS and dynamics references.

Useful reference points:

- Jordan, Orus, Vidal, Verstraete, and Cirac, "Classical Simulation of Infinite-Size Quantum Lattice Systems in Two Spatial Dimensions", Phys. Rev. Lett. 101, 250602 (2008), https://doi.org/10.1103/PhysRevLett.101.250602. This is the baseline infinite PEPS TFIM ground-state benchmark and critical-transition reference.
- Schmitt, Rams, Dziarmaga, Heyl, and Zurek, "Quantum phase transition dynamics in the two-dimensional transverse-field Ising model", Science Advances 8, eabl6850 (2022), https://arxiv.org/abs/2106.09046. This provides infinite/large-system TFIM critical dynamics context and iPEPS method expectations.
- Arias Espinoza and Corboz, "Spectral functions with infinite projected entangled-pair states", Phys. Rev. B 110, 094314 (2024), https://arxiv.org/abs/2405.10628. This shows TFIM as an iPEPS real-time local-quench benchmark and identifies unit-cell/time-window limitations.
- Vovrosh et al., "Simulating dynamics of the two-dimensional transverse-field Ising model: a comparative study of large-scale classical numerics", arXiv:2511.19340 (2025), https://arxiv.org/abs/2511.19340. Markus Schmitt's group page lists this work as accepted in Physical Review Research: https://www.computational-quantum.science/publications/33_classical_simulation_benchmark/. This is finite-system cross-method dynamics context and should be secondary for the infinite benchmark.

## Architecture

Add a narrow star-model abstraction. Julia containers that store models or
protocols must be parametric rather than storing `AbstractStarModel` directly.

```julia
abstract type AbstractStarModel end

struct PXPStarModel <: AbstractStarModel
    projected::Bool
end

struct TFIMStarModel{T<:Real} <: AbstractStarModel
    J::T
    h::T
end
```

The shared model interface should be:

```julia
star_hamiltonian(model)
star_gate(model, step; evolution = :real)
star_gate_itensor(model, site_indices, step; evolution = :real)
```

`PXPStarModel(true)` must reproduce the current projected PXP gate behavior. `PXPStarModel(false)` must reproduce the current unprojected PXP gate behavior.

For TFIM, use the five-site square-star convention `(center, right, up, left, down)`:

```text
h_c = -h X_c - (J/2) Z_c (Z_right + Z_up + Z_left + Z_down)
```

The factor `J/2` avoids double-counting nearest-neighbor bonds when the star term is summed over every lattice site. This convention is the source of truth for TFIM energy-density measurements in the first benchmark framework.

The implementation should expose convention helpers or constants that tests can
pin explicitly:

```julia
star_site_order() == (:center, :right, :up, :left, :down)
tfim_pauli_convention() == (:Z_up_is_plus_one, :X_field)
```

These names are less important than the requirement: basis ordering, Pauli
signs, star direction order, ITensor index order, priming convention, Trotter
sign, and field/bond normalization must be tested directly before any PEPS
evolution tests run.

## Evolution Changes

Extend the existing update path without duplicating it. The update kernel should
receive a concrete model:

```julia
project_star!(
    psi,
    center,
    step;
    model = PXPStarModel(true),
    evolution = :real,
    maxdim = psi.maxdim,
    cutoff = 1e-12,
    split_order = (:right, :up, :left, :down),
)
```

`project_star!` should obtain the gate through `star_gate_itensor(model, phys, step; evolution)` and then keep the current absorption, QR reduction, SVD splitting, reconstruction, and diagnostic behavior.

Keep Trotter algorithm parameters model-agnostic. `TrotterParams` should describe
the schedule and truncation choices only:

```julia
struct TrotterParams
    dt::Float64
    order::Int
    evolution::Symbol
    maxdim::Int
    cutoff::Float64
    split_order::NTuple{4,Symbol}
end
```

Model choice belongs in a protocol object:

```julia
abstract type AbstractModelProtocol end

struct StaticModel{M<:AbstractStarModel} <: AbstractModelProtocol
    model::M
end

model_at(protocol::StaticModel, time, step) = protocol.model
```

`evolve!` should ask `model_at(protocol, time, step)` for the model at each
full Trotter step and pass that concrete model into each `project_star!` call.
The first implementation only needs `StaticModel`, but this split leaves a
clean path for later ramps or annealing protocols without redesigning the
Trotter API.

The current constructor shape must remain supported for PXP compatibility. The
compatibility path should map the old `projected::Bool` argument and the current
`evolve!(...; projected = true)` keyword to `StaticModel(PXPStarModel(projected))`
internally. New benchmark code should use explicit protocols instead of storing
model choices inside `TrotterParams`. If exact source compatibility for
`TrotterParams(dt, order, evolution, projected, maxdim, cutoff)` conflicts with
strictly model-agnostic `TrotterParams`, preserve behavior with a small
compatibility wrapper accepted by `evolve!`; do not carry the legacy
`projected` field into the new benchmark API.

`evolve!` should continue to use the five-color square-star schedule. No arbitrary Hamiltonian framework, graph support, or new lattice support should be introduced.

## Observables

Add simple/local TFIM observables first:

```julia
local_x_simple(psi, c)::Float64
local_y_simple(psi, c)::Float64
local_z_simple(psi, c)::Float64
nearest_neighbor_zz_simple(psi, c, dir)::Float64
tfim_energy_density_star_simple(psi, model::TFIMStarModel)::Float64
tfim_energy_density_decomposed_simple(psi, model::TFIMStarModel)::Float64
measure_tfim_simple(psi, model::TFIMStarModel)::TFIMObservableSummary
```

`tfim_energy_density_star_simple` should average the same
`star_hamiltonian(model)` over all unit-cell representatives through the
existing star-patch expectation logic. `tfim_energy_density_decomposed_simple`
should compute:

```text
e = -h <X> - J (<ZZ>_right + <ZZ>_up)
```

The summary should include the discrepancy between the star-patch and
decomposed energy estimates. This is a direct diagnostic for direction mapping
and bond double-counting mistakes. `nearest_neighbor_zz_simple` should support
canonical directions `:right` and `:up`, with the same direction validation
conventions already used by PXP diagnostics.

The first `TFIMObservableSummary` should include:

- mean `<X>`
- mean `<Y>`
- mean `<Z>`
- optional even/odd `<Z>` for detecting symmetry breaking in larger unit cells
- mean nearest-neighbor `<ZZ>` over right/up bonds
- star-patch TFIM energy density
- decomposed TFIM energy density
- absolute star/decomposed energy discrepancy
- absolute discarded imaginary part for every Hermitian observable
- mean bond entropy
- max bond entropy

## Benchmark Cases

### Tier 0: Dense Star Gate Correctness

Validate the TFIM star model independently of PEPS:

- `star_hamiltonian(TFIMStarModel(J, h))` is Hermitian.
- `star_gate(model, dt; evolution = :real)` is unitary.
- `star_gate(model, dt; evolution = :imaginary)` is finite and non-unitary in the expected direction.
- Single-star real-time gates compose for one isolated star Hamiltonian.
- Product-basis diagonal energies match the `J/2` bond double-counting convention.
- Named convention tests pin site order, Pauli convention, basis-state energies,
  and real-time/imaginary-time sign conventions.
- `PXPStarModel` gates match the current `square_pxp_gate` and `projected_square_pxp_gate` exactly.

### Tier 1: Product-State Short-Time Dynamics

Use infinite product-state iPEPS with `D = 1` to catch sign, basis, and normalization bugs:

- `|up...up>` and `|down...down>` product states have expected initial `<Z>`, `<X>`, `<ZZ>`, and TFIM energy.
- Short real-time evolution produces finite deterministic diagnostics.
- For tiny `dt`, first nonzero short-time changes match dense local reference calculations within a conservative tolerance.
- Exact `J = 0` independent-spin dynamics from `|Z+>` matches `<Z(t)> = cos(2 h t)` at `D = 1`, up to the chosen Pauli/time-direction convention for `<Y(t)>`.
- Exact `h = 0` classical Ising dynamics keeps all-up/all-down product-state observables stationary up to phase.
- Short-time derivative checks include observables whose first nonzero change is second order in time, not only first-order changes.

### Tier 1.5: Finite Schedule Reference

Add a tiny finite periodic `L x L` dense or sparse reference that uses the same
star Hamiltonian convention, five-color schedule, and first-/second-order
Trotter coefficients as the iPEPS driver. This test is not a physics benchmark;
it exists to catch wrong color order, wrong half-steps, wrong periodic neighbor
mapping, wrong center selection, and incorrect `dt` scaling per sublayer before
the `10x10` infinite smoke runs are trusted.

### Tier 2: Infinite iPEPS Regression Runs

Use periodic cells compatible with the five-color schedule, starting with `10x10`.

Run fixed smoke configurations:

- exact single-spin case: `J = 0`, `h = 1.0`
- classical stationary case: `J = 1.0`, `h = 0`
- small field ratio, such as `h/J = 0.5`
- near critical field ratio, using the literature value around `h/J = 3.04`
- large field ratio, such as `h/J = 5.0` or approximately `2 h_c`

Use at least these initial states:

- `:z_up`
- `:z_down`
- `:x_plus`

Optional later initial states are `:checkerboard_z` and
`:random_product_seeded`.

For each case record:

- model parameters
- unit-cell dimensions
- initial state
- `D`, `dt`, Trotter order, total time, cutoff
- max truncation error
- mean/max bond entropy
- log-norm before/after/delta
- full time series of TFIM simple observable summaries and evolution diagnostics

Convergence sweeps should vary `D` and `dt` separately. The framework should report the data; it should not claim physics-quality agreement from simple-update diagnostics alone.

### Tier 3: CTM Validation Later

Once the PEPSKit adapter is stable enough, add CTM-backed TFIM observables and finite-`chi` sweeps:

- CTM energy density
- CTM `<X>` and `<Z>`
- CTM nearest-neighbor `<ZZ>`
- CTM convergence diagnostics and finite-`chi` sensitivity

Only this tier should be used for direct comparison to infinite-iPEPS or QMC literature values.

## Benchmark Run Records

Add a small benchmark runner API:

```julia
struct BenchmarkSpec{P<:AbstractModelProtocol,C}
    name::String
    protocol::P
    cell::C
    initial_state::Symbol
    total_time::Float64
    trotter::TrotterParams
    measure_every::Int
end

struct BenchmarkSample
    step::Int
    time::Float64
    observables::TFIMObservableSummary
    diagnostics::EvolutionDiagnostics
end

struct BenchmarkResult
    name::String
    run_label::String
    metadata::BenchmarkMetadata
    samples::Vector{BenchmarkSample}
    final_state_summary::BenchmarkSample
end

run_benchmark(spec)::BenchmarkResult
write_benchmark_json(result, path)
write_benchmark_csv(results, path)
```

`EvolutionDiagnostics` should be a lightweight per-sample summary derived from
`EvolutionLog` and `StarUpdateInfo`, not a second evolution implementation.

`BenchmarkResult` should be a reproducibility record and a trajectory, not just
a final table row. It should include:

- benchmark name
- run label supplied by caller
- model protocol type and parameters
- cell dimensions
- initial state
- evolution parameters
- sample vector with `step`, `time`, observables, and diagnostics
- final-state summary copied from the last sample for convenience
- package version if available

JSON is the canonical output for full metadata. CSV is a deterministic flattened
time-series table with one row per sample, including columns such as:

```text
name,run_label,step,time,J,h,D,dt,order,mean_x,mean_y,mean_z,
zz_right,zz_up,energy_star,energy_decomposed,energy_discrepancy,
max_truncerr,mean_bond_entropy,max_bond_entropy,lognorm_delta
```

Use `JSON3.jl` for JSON output rather than a hand-rolled serializer. Benchmark
serialization should reject nonfinite floats unless a field is explicitly
optional and represented as `nothing`.

## Scope Boundaries

Do not add:

- general Hamiltonian packaging
- arbitrary graph or lattice support
- spectral functions
- local-quench inserted-cell dynamics
- broad CTMRG infrastructure beyond the existing adapter
- GPU backends or symmetry machinery

This benchmark framework should answer one near-term question: whether the current PXP star-gate iPEPS machinery can reproduce controlled infinite-system dynamics on a model with stronger external reference data.

## Testing And Acceptance

The first implementation is acceptable when these gates pass:

### Gate A: Compatibility

- Existing PXP tests pass unchanged.
- `PXPStarModel(true)` and `PXPStarModel(false)` reproduce current projected and
  unprojected PXP gates exactly or within a stated floating tolerance.
- Current `project_star!` and `evolve!` PXP call shapes remain supported.

### Gate B: Dense TFIM Convention

- Hamiltonian is Hermitian.
- Real-time gate is unitary within tolerance.
- Imaginary-time gate is finite and non-unitary in the expected direction.
- Product-basis energies match the `J/2` star convention.
- Site order, Pauli convention, direction convention, and Trotter sign tests pass.

### Gate C: Exact Limits

- `J = 0` independent-spin dynamics matches analytic curves at `D = 1`.
- `h = 0` all-up/all-down states are stationary in observables.
- Tiny-`dt` short-time expansion matches dense local reference values.

### Gate D: Schedule

- Five-color schedule has no same-color star overlaps.
- First- and second-order Trotter coefficient sums are correct.
- Finite dense/sparse schedule reference passes, or is explicitly marked as an
  extended test if runtime is too high for default CI.

### Gate E: Infinite Regression

- Short `10x10` infinite TFIM runs complete deterministically.
- Result JSON contains full metadata and time-series samples.
- CSV flattening is deterministic.
- Documentation states that simple-update observables are regression diagnostics,
  not CTM-quality physics.

## Implementation Decisions

- Put `AbstractStarModel`, `PXPStarModel`, `TFIMStarModel`, and dense star-gate helpers in a new `src/StarModels.jl` module.
- Put `AbstractModelProtocol`, `StaticModel`, and `model_at` in the same module unless the implementation naturally needs a separate `src/Protocols.jl`.
- Export the model types and core model functions because callers need them to configure benchmark runs. Add docstrings for every exported symbol to satisfy the existing public-doc tests.
- Keep `TrotterParams` model-agnostic for new code. Preserve old PXP constructor and keyword call shapes through compatibility shims that build `StaticModel(PXPStarModel(projected))`.
- Put benchmark runner types and file writers in a new `src/Benchmarks.jl` module. Export `BenchmarkSpec`, `BenchmarkSample`, `BenchmarkResult`, `run_benchmark`, `write_benchmark_json`, and `write_benchmark_csv` only after their tests are in place.
- Add `JSON3.jl` for benchmark JSON serialization. Reject nonfinite floats during result construction or serialization so invalid benchmark records fail loudly.
