# Infinite TFIM Benchmark Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a v1 infinite TFIM benchmark framework that reuses the existing five-site square-star iPEPS update path while preserving PXP behavior.

**Architecture:** Introduce a narrow `StarModels` layer for dense five-site model gates and static model protocols. Thread concrete models through `project_star!` and `evolve!`, add TFIM simple observables, then build a small benchmark runner that records time-series samples to JSON and CSV.

**Tech Stack:** Julia 1.12, ITensors.jl, LinearAlgebra, JSON3.jl, existing `SquarePXPDynamics` modules, `Test`.

---

## File Structure

- Create `src/StarModels.jl`: `AbstractStarModel`, `PXPStarModel`, `TFIMStarModel`, convention helpers, dense Hamiltonians/gates, ITensor gate conversion, `StaticModel`, and `model_at`.
- Modify `src/SquarePXPDynamics.jl`: include/use/export `StarModels`, benchmark APIs, and TFIM observables.
- Modify `src/StarSimpleUpdate.jl`: accept `model = PXPStarModel(projected)` while keeping the legacy `projected` keyword.
- Modify `src/IPEPSEvolution.jl`: make new `TrotterParams` model-agnostic, add `LegacyPXPParams` compatibility, add `protocol` support, and pass models into `project_star!`.
- Modify `src/SquareIPEPS.jl`: support product states `:z_up`, `:z_down`, and `:x_plus` for benchmark initialization.
- Modify `src/SquarePEPS.jl`: keep finite PEPS product-state aliases aligned with iPEPS.
- Modify `src/Observables.jl`: add local Pauli expectations, TFIM summary, TFIM energy measurements, and imaginary-residual diagnostics.
- Create `src/Benchmarks.jl`: benchmark specs, samples, metadata, diagnostics summaries, runner, JSON writer, CSV writer.
- Modify `Project.toml` and `test/Project.toml`: add `JSON3`.
- Create `test/test_star_models.jl`: dense model, convention, PXP compatibility, and exact star-gate tests.
- Create `test/test_tfim_observables.jl`: local/product-state TFIM observables and exact-limit checks.
- Create `test/test_tfim_schedule_reference.jl`: finite periodic schedule reference for Trotter layer coefficients, full-sweep coverage, and center mapping.
- Create `test/test_benchmarks.jl`: benchmark result time-series, JSON/CSV serialization, deterministic short runs.
- Modify `test/test_ipeps_evolution.jl` and `test/test_star_simple_update.jl`: compatibility assertions for old and new model/protocol APIs.
- Modify `test/test_public_docs.jl`: public doc coverage for newly exported symbols.

## Task 1: Add JSON3 Dependency

**Files:**
- Modify: `Project.toml`
- Modify: `test/Project.toml`

- [ ] **Step 1: Add JSON3 with Julia's package manager**

Run:

```bash
julia --project=. -e 'using Pkg; Pkg.add("JSON3")'
julia --project=test -e 'using Pkg; Pkg.add("JSON3")'
```

Expected: both commands finish successfully and update the project manifests.

- [ ] **Step 2: Verify dependency loads**

Run:

```bash
julia --project=. -e 'using JSON3; println(JSON3.write(Dict("ok" => true)))'
julia --project=test -e 'using JSON3; println(JSON3.write(Dict("ok" => true)))'
```

Expected output contains:

```text
{"ok":true}
```

- [ ] **Step 3: Commit dependency changes**

Run:

```bash
git add Project.toml Manifest.toml test/Project.toml test/Manifest.toml
git commit -m "deps: add JSON3 for benchmark records"
```

Expected: commit succeeds and stages only dependency files.

## Task 2: Add Star Model Abstraction With Dense TFIM Gates

**Files:**
- Create: `src/StarModels.jl`
- Modify: `src/SquarePXPDynamics.jl`
- Test: `test/test_star_models.jl`
- Modify: `test/test_public_docs.jl`

- [ ] **Step 1: Write failing star-model tests**

Create `test/test_star_models.jl`:

```julia
using Test
using LinearAlgebra
using SquarePXPDynamics

@testset "star model conventions" begin
    @test star_site_order() == (:center, :right, :up, :left, :down)
    @test tfim_pauli_convention() == (:Z_up_is_plus_one, :X_field)
end

@testset "PXP star model reproduces existing gates" begin
    dt = 0.037
    @test star_hamiltonian(PXPStarModel(false)) ≈ square_pxp_star_hamiltonian()
    @test star_gate(PXPStarModel(false), dt; evolution = :real) ≈
          square_pxp_gate(dt; evolution = :real)
    @test star_gate(PXPStarModel(true), dt; evolution = :real) ≈
          projected_square_pxp_gate(dt; evolution = :real)
    @test star_gate(PXPStarModel(false), dt; evolution = :imaginary) ≈
          square_pxp_gate(dt; evolution = :imaginary)
    @test star_gate(PXPStarModel(true), dt; evolution = :imaginary) ≈
          projected_square_pxp_gate(dt; evolution = :imaginary)
end

@testset "TFIM dense star Hamiltonian convention" begin
    model = TFIMStarModel(2.0, 3.0)
    H = star_hamiltonian(model)
    @test size(H) == (32, 32)
    @test H ≈ H'
    @test tfim_product_basis_energy(model, (:up, :up, :up, :up, :up)) ≈ -4.0
    @test tfim_product_basis_energy(model, (:up, :down, :up, :down, :up)) ≈ 0.0
    @test tfim_product_basis_energy(model, (:down, :up, :up, :up, :up)) ≈ 4.0
end

@testset "TFIM dense gates" begin
    model = TFIMStarModel(1.0, 0.7)
    dt = 0.011
    U = star_gate(model, dt; evolution = :real)
    G = star_gate(model, dt; evolution = :imaginary)
    @test U' * U ≈ I atol = 1e-12 rtol = 1e-12
    @test all(isfinite, real.(G))
    @test norm(G' * G - I) > 1e-6
    @test star_gate(model, dt; evolution = :real) *
          star_gate(model, dt; evolution = :real) ≈
          star_gate(model, 2dt; evolution = :real) atol = 1e-12 rtol = 1e-12
    @test_throws ArgumentError star_gate(model, dt; evolution = :bad)
    @test_throws ArgumentError TFIMStarModel(Inf, 1.0)
    @test_throws ArgumentError TFIMStarModel(1.0, NaN)
end

@testset "static model protocol" begin
    model = TFIMStarModel(1.0, 2.0)
    protocol = StaticModel(model)
    @test model_at(protocol, 0.0, 1) === model
    @test model_at(protocol, 1.0, 17) === model
end
```

- [ ] **Step 2: Run tests to verify they fail**

Run:

```bash
julia --project=test -e 'using Pkg; Pkg.test("SquarePXPDynamics"; test_args=["test_star_models.jl"])'
```

Expected: FAIL with undefined names such as `PXPStarModel`, `TFIMStarModel`, or `star_site_order`.

- [ ] **Step 3: Implement `src/StarModels.jl`**

Create `src/StarModels.jl`:

```julia
module StarModels

using LinearAlgebra
using ITensors

using ..SpinOps: pauli_x, pauli_z, identity2, kron_all, embed_one_site
using ..SquarePXP:
    SQUARE_STAR_SITES,
    square_pxp_star_hamiltonian,
    square_pxp_gate,
    projected_square_pxp_gate

export AbstractStarModel,
    PXPStarModel,
    TFIMStarModel,
    AbstractModelProtocol,
    StaticModel,
    star_site_order,
    tfim_pauli_convention,
    tfim_product_basis_energy,
    star_hamiltonian,
    star_gate,
    star_gate_itensor,
    model_at

const STAR_SITE_ORDER = (:center, :right, :up, :left, :down)
const TFIM_PAULI_CONVENTION = (:Z_up_is_plus_one, :X_field)

abstract type AbstractStarModel end
abstract type AbstractModelProtocol end

"""
    PXPStarModel(projected)

Five-site square-star PXP model wrapper. `projected = true` selects the
blockade-projected dense gate and `projected = false` selects the unprojected
dense gate.
"""
struct PXPStarModel <: AbstractStarModel
    projected::Bool
end

"""
    TFIMStarModel(J, h)

Square-lattice transverse-field Ising model star term using
`-h X_center - (J/2) Z_center * sum(Z_neighbor)`.
"""
struct TFIMStarModel{T<:Real} <: AbstractStarModel
    J::T
    h::T
    function TFIMStarModel(J::T, h::T) where {T<:Real}
        isfinite(J) || throw(ArgumentError("J must be finite"))
        isfinite(h) || throw(ArgumentError("h must be finite"))
        return new{T}(J, h)
    end
end
TFIMStarModel(J::Real, h::Real) = TFIMStarModel(promote(J, h)...)

"""
    StaticModel(model)

Static model protocol returning the same star model for every time and step.
"""
struct StaticModel{M<:AbstractStarModel} <: AbstractModelProtocol
    model::M
end

"""Return the physical site order used by dense five-site star operators."""
star_site_order() = STAR_SITE_ORDER

"""Return the TFIM Pauli convention used by dense star operators."""
tfim_pauli_convention() = TFIM_PAULI_CONVENTION

"""Return the model active at `time` and integer Trotter `step`."""
model_at(protocol::StaticModel, time, step) = protocol.model

star_hamiltonian(model::PXPStarModel) = square_pxp_star_hamiltonian()

function _tfim_dense_type(model::TFIMStarModel)
    return promote_type(typeof(model.J), typeof(model.h), Float64)
end

function star_hamiltonian(model::TFIMStarModel)
    T = _tfim_dense_type(model)
    X = Matrix{T}(pauli_x())
    Z = Matrix{T}(pauli_z())
    I2 = Matrix{T}(identity2())
    H = zeros(T, 2^SQUARE_STAR_SITES, 2^SQUARE_STAR_SITES)
    H .-= model.h .* embed_one_site(X, 1, SQUARE_STAR_SITES)
    for site in 2:SQUARE_STAR_SITES
        ops = [I2 for _ in 1:SQUARE_STAR_SITES]
        ops[1] = Z
        ops[site] = Z
        H .-= (model.J / 2) .* kron_all(ops)
    end
    return ComplexF64.(H)
end

function star_gate(model::PXPStarModel, step::Real; evolution::Symbol = :real)
    return model.projected ?
           projected_square_pxp_gate(step; evolution) :
           square_pxp_gate(step; evolution)
end

function star_gate(model::TFIMStarModel, step::Real; evolution::Symbol = :real)
    finite_step = Float64(step)
    isfinite(finite_step) || throw(ArgumentError("step must be finite"))
    H = star_hamiltonian(model)
    if evolution === :real
        return exp((-im * finite_step) .* H)
    elseif evolution === :imaginary
        return exp((-finite_step) .* H)
    else
        throw(ArgumentError("evolution must be :real or :imaginary"))
    end
end

function _assign_gate_entries!(G::ITensor, dense, out_phys, in_phys)
    ranges = ntuple(_ -> 1:2, SQUARE_STAR_SITES)
    for out_values in Iterators.product(ranges...)
        row = LinearIndices(ntuple(_ -> 1:2, SQUARE_STAR_SITES))[out_values...]
        for in_values in Iterators.product(ranges...)
            col = LinearIndices(ntuple(_ -> 1:2, SQUARE_STAR_SITES))[in_values...]
            G[Tuple(vcat(Any[out_phys[i] => out_values[i] for i in 1:SQUARE_STAR_SITES],
                         Any[in_phys[i] => in_values[i] for i in 1:SQUARE_STAR_SITES]))...] =
                dense[row, col]
        end
    end
    return G
end

function star_gate_itensor(model::AbstractStarModel, site_indices, step::Real; evolution::Symbol = :real)
    length(site_indices) == SQUARE_STAR_SITES ||
        throw(ArgumentError("site_indices must contain five physical indices"))
    in_phys = collect(site_indices)
    out_phys = prime.(in_phys)
    G = ITensor(ComplexF64, out_phys..., in_phys...)
    return _assign_gate_entries!(G, star_gate(model, step; evolution), out_phys, in_phys)
end

function tfim_product_basis_energy(model::TFIMStarModel, states)
    length(states) == SQUARE_STAR_SITES ||
        throw(ArgumentError("states must contain five site labels"))
    zvalue(s) =
        s in (:up, :z_up, 1, +1) ? 1.0 :
        s in (:down, :z_down, 2, -1) ? -1.0 :
        throw(ArgumentError("basis state must be :up/:down, :z_up/:z_down, 1/2, or +/-1"))
    z = map(zvalue, states)
    return -(model.J / 2) * z[1] * sum(z[2:end])
end

end
```

- [ ] **Step 4: Include and export `StarModels` APIs**

Patch `src/SquarePXPDynamics.jl` so `StarModels` is included after `SquareIPEPS.jl` and before modules that will consume it:

```julia
include("StarModels.jl")
```

Add the imports:

```julia
using .StarModels:
    AbstractStarModel,
    PXPStarModel,
    TFIMStarModel,
    AbstractModelProtocol,
    StaticModel,
    star_site_order,
    tfim_pauli_convention,
    tfim_product_basis_energy,
    star_hamiltonian,
    star_gate,
    star_gate_itensor,
    model_at
```

Add the exports:

```julia
export AbstractStarModel, PXPStarModel, TFIMStarModel
export AbstractModelProtocol, StaticModel, model_at
export star_site_order, tfim_pauli_convention, tfim_product_basis_energy
export star_hamiltonian, star_gate, star_gate_itensor
```

- [ ] **Step 5: Run star-model tests**

Run:

```bash
julia --project=test -e 'using Pkg; Pkg.test("SquarePXPDynamics"; test_args=["test_star_models.jl"])'
```

Expected: PASS.

- [ ] **Step 6: Run public docs tests**

Run:

```bash
julia --project=test -e 'using Pkg; Pkg.test("SquarePXPDynamics"; test_args=["test_public_docs.jl"])'
```

Expected: PASS. If it fails with missing docs for newly exported symbols, add docstrings in `src/StarModels.jl` for the exact missing names and rerun.

- [ ] **Step 7: Commit star models**

Run:

```bash
git add src/StarModels.jl src/SquarePXPDynamics.jl test/test_star_models.jl test/test_public_docs.jl
git commit -m "feat: add square-star model abstraction"
```

Expected: commit succeeds.

## Task 3: Thread Star Models Through Star Updates And Evolution

**Files:**
- Modify: `src/StarSimpleUpdate.jl`
- Modify: `src/IPEPSEvolution.jl`
- Modify: `src/SquarePXPDynamics.jl`
- Test: `test/test_star_simple_update.jl`
- Test: `test/test_ipeps_evolution.jl`

- [ ] **Step 1: Add failing compatibility and new-model tests**

Append to `test/test_star_simple_update.jl`:

```julia
@testset "project_star accepts explicit star models" begin
    cell = PeriodicSquareUnitCell(10, 10)
    center = SquareCoord(5, 5)
    dt = 0.01

    legacy = product_square_ipeps(cell; state = :down, maxdim = 1)
    explicit = product_square_ipeps(cell; state = :down, maxdim = 1)
    project_star!(legacy, center, dt; projected = true, maxdim = 1)
    project_star!(explicit, center, dt; model = PXPStarModel(true), maxdim = 1)

    @test local_density_simple(explicit, center) ≈ local_density_simple(legacy, center) atol = 1e-12
    @test log_norm(explicit) ≈ log_norm(legacy) atol = 1e-12

    tfim = product_square_ipeps(cell; state = :up, maxdim = 1)
    info = project_star!(tfim, center, 0.0; model = TFIMStarModel(1.0, 0.0), maxdim = 1)
    @test info.max_truncerr ≥ 0
    @test isfinite(log_norm(tfim))
end
```

Append to `test/test_ipeps_evolution.jl`:

```julia
@testset "evolve accepts explicit static model protocol" begin
    cell = PeriodicSquareUnitCell(10, 10)
    params = TrotterParams(0.01, 1, :real, 1, 1e-12)

    legacy = product_square_ipeps(cell; state = :down, maxdim = 1)
    explicit = product_square_ipeps(cell; state = :down, maxdim = 1)

    legacy_log = evolve!(legacy, 0.01; dt = 0.01, order = 1, evolution = :real, projected = true)
    explicit_log = evolve!(
        explicit,
        0.01;
        params = params,
        protocol = StaticModel(PXPStarModel(true)),
    )

    @test explicit_log.nsteps == legacy_log.nsteps
    @test explicit_log.max_truncerr ≈ legacy_log.max_truncerr atol = 1e-12
    @test log_norm(explicit) ≈ log_norm(legacy) atol = 1e-12
end

@testset "legacy TrotterParams constructor remains accepted" begin
    old = TrotterParams(0.01, 1, :real, true, 1, 1e-12)
    cell = PeriodicSquareUnitCell(10, 10)
    psi = product_square_ipeps(cell; state = :down, maxdim = 1)
    log = evolve!(psi, 0.01; params = old)
    @test log.nsteps == 1
    @test isfinite(log.max_truncerr)
end
```

- [ ] **Step 2: Run tests to verify they fail**

Run:

```bash
julia --project=test -e 'using Pkg; Pkg.test("SquarePXPDynamics"; test_args=["test_star_simple_update.jl","test_ipeps_evolution.jl"])'
```

Expected: FAIL because `project_star!` does not accept `model`, `TrotterParams(0.01, 1, :real, 1, 1e-12)` does not exist, or `evolve!` does not accept `protocol`.

- [ ] **Step 3: Update `StarSimpleUpdate` imports and signature**

Patch `src/StarSimpleUpdate.jl` near imports:

```julia
using ..StarModels: AbstractStarModel, PXPStarModel, star_gate_itensor
```

Change `project_star!` signature to:

```julia
function project_star!(
    psi::SquareIPEPSState,
    center::SquareCoord,
    step::Real;
    model::Union{AbstractStarModel,Nothing} = nothing,
    evolution::Symbol = :real,
    projected::Bool = true,
    maxdim::Integer = psi.maxdim,
    cutoff::Real = 1e-12,
    split_order = _STAR_DIRECTIONS,
)::StarUpdateInfo
```

Replace the existing gate selection block with:

```julia
    actual_model = model === nothing ? PXPStarModel(projected) : model
    gate = star_gate_itensor(actual_model, phys, finite_step; evolution = evolution)
```

Leave all QR, SVD, reconstruction, and transaction code unchanged.

- [ ] **Step 4: Update `TrotterParams` and compatibility params**

Patch `src/IPEPSEvolution.jl` imports:

```julia
using ..StarModels: AbstractModelProtocol, StaticModel, PXPStarModel, model_at
```

Replace the current `TrotterParams` with:

```julia
struct TrotterParams
    dt::Float64
    order::Int
    evolution::Symbol
    maxdim::Int
    cutoff::Float64
    split_order::NTuple{4,Symbol}

    function TrotterParams(
        dt::Real,
        order::Integer,
        evolution::Symbol,
        maxdim::Integer,
        cutoff::Real,
        split_order = (:right, :up, :left, :down),
    )
        step = Float64(dt)
        isfinite(step) && step > 0 || throw(ArgumentError("dt must be finite and positive"))
        ord = Int(order)
        ord in (1, 2) || throw(ArgumentError("order must be 1 or 2"))
        evolution in (:real, :imaginary) ||
            throw(ArgumentError("evolution must be :real or :imaginary"))
        dim = Int(maxdim)
        dim >= 1 || throw(ArgumentError("maxdim must be at least 1"))
        trunc_cutoff = Float64(cutoff)
        isfinite(trunc_cutoff) && trunc_cutoff >= 0 ||
            throw(ArgumentError("cutoff must be finite and nonnegative"))
        dirs = Tuple(split_order)
        length(dirs) == 4 || throw(ArgumentError("split_order must contain four directions"))
        return new(step, ord, evolution, dim, trunc_cutoff, dirs)
    end
end

struct LegacyPXPParams
    trotter::TrotterParams
    protocol::StaticModel{PXPStarModel}
end

function TrotterParams(
    dt::Real,
    order::Integer,
    evolution::Symbol,
    projected::Bool,
    maxdim::Integer,
    cutoff::Real,
)
    return LegacyPXPParams(
        TrotterParams(dt, order, evolution, maxdim, cutoff),
        StaticModel(PXPStarModel(projected)),
    )
end
```

- [ ] **Step 5: Update evolution internals for protocols**

Change `EvolutionLog.params` to:

```julia
params::TrotterParams
```

Add helper:

```julia
_unwrap_params(params::TrotterParams, protocol) = params, protocol
_unwrap_params(params::LegacyPXPParams, protocol) =
    protocol === nothing ? (params.trotter, params.protocol) :
    throw(ArgumentError("do not pass protocol with legacy PXP TrotterParams"))
```

Update `_evolve_with_params!` signature:

```julia
function _evolve_with_params!(
    psi::SquareIPEPSState,
    total_time::Real,
    params::TrotterParams,
    protocol::AbstractModelProtocol,
)::EvolutionLog
```

Inside the step loop, use `step_index` and `time_before_step`:

```julia
    for step_index = 1:nsteps
        time_before_step = (step_index - 1) * params.dt
        step_model = model_at(protocol, time_before_step, step_index)
        for (color, layer_dt) in sequence
            centers = update_centers(psi.unitcell, color)
            stars_are_disjoint_mod_unitcell(psi.unitcell, centers) ||
                throw(ArgumentError("color layer $color has overlapping wrapped stars"))
            infos = StarUpdateInfo[]
            for center in centers
                push!(
                    infos,
                    project_star!(
                        psi,
                        center,
                        layer_dt;
                        model = step_model,
                        evolution = params.evolution,
                        maxdim = params.maxdim,
                        cutoff = params.cutoff,
                        split_order = params.split_order,
                    ),
                )
            end
            push!(layer_infos, infos)
        end
    end
```

Change `evolve!` keyword signature to include:

```julia
protocol::Union{AbstractModelProtocol,Nothing} = nothing,
```

Build actual params and protocol:

```julia
    raw_params = if params === nothing
        dt === nothing && throw(UndefKeywordError(:dt))
        TrotterParams(dt, order, evolution, maxdim, cutoff)
    else
        params
    end
    actual_params, params_protocol = _unwrap_params(raw_params, protocol)
    actual_protocol =
        params_protocol === nothing ? StaticModel(PXPStarModel(projected)) : params_protocol
    return _evolve_with_params!(psi, total_time, actual_params, actual_protocol)
```

- [ ] **Step 6: Export `LegacyPXPParams` only if public docs require it**

If `test/test_public_docs.jl` checks every exported type, do not export `LegacyPXPParams`. Keep it internal. Ensure `using .IPEPSEvolution:` in `src/SquarePXPDynamics.jl` still imports only `TrotterParams`, `EvolutionLog`, `trotter_sequence`, and `evolve!`.

- [ ] **Step 7: Run targeted evolution tests**

Run:

```bash
julia --project=test -e 'using Pkg; Pkg.test("SquarePXPDynamics"; test_args=["test_star_simple_update.jl","test_ipeps_evolution.jl"])'
```

Expected: PASS.

- [ ] **Step 8: Commit evolution integration**

Run:

```bash
git add src/StarSimpleUpdate.jl src/IPEPSEvolution.jl src/SquarePXPDynamics.jl test/test_star_simple_update.jl test/test_ipeps_evolution.jl
git commit -m "feat: thread star models through evolution"
```

Expected: commit succeeds.

## Task 4: Add Benchmark Product-State Aliases

**Files:**
- Modify: `src/SquareIPEPS.jl`
- Modify: `src/SquarePEPS.jl`
- Test: `test/test_square_ipeps.jl`
- Test: `test/test_square_peps.jl`

- [ ] **Step 1: Write failing product-state alias tests**

Append to `test/test_square_ipeps.jl`:

```julia
@testset "benchmark product-state aliases" begin
    cell = PeriodicSquareUnitCell(1, 1)
    z_up = product_square_ipeps(cell; state = :z_up, maxdim = 1)
    z_down = product_square_ipeps(cell; state = :z_down, maxdim = 1)
    x_plus = product_square_ipeps(cell; state = :x_plus, maxdim = 1)

    @test local_density_simple(z_up, SquareCoord(1, 1)) ≈ 1.0 atol = 1e-12
    @test local_density_simple(z_down, SquareCoord(1, 1)) ≈ 0.0 atol = 1e-12
    @test local_density_simple(x_plus, SquareCoord(1, 1)) ≈ 0.5 atol = 1e-12
end
```

Append to `test/test_square_peps.jl`:

```julia
@testset "finite PEPS benchmark product-state aliases" begin
    @test product_square_peps(1, 1; state = :z_up).tensors[SquareCoord(1, 1)][1] ≈ 1
    @test product_square_peps(1, 1; state = :z_down).tensors[SquareCoord(1, 1)][2] ≈ 1
    plus = product_square_peps(1, 1; state = :x_plus).tensors[SquareCoord(1, 1)]
    @test plus[1] ≈ inv(sqrt(2))
    @test plus[2] ≈ inv(sqrt(2))
end
```

- [ ] **Step 2: Run tests to verify they fail**

Run:

```bash
julia --project=test -e 'using Pkg; Pkg.test("SquarePXPDynamics"; test_args=["test_square_ipeps.jl","test_square_peps.jl"])'
```

Expected: FAIL with `ArgumentError` for unsupported product states.

- [ ] **Step 3: Add a shared state-vector helper in both constructors**

In `src/SquareIPEPS.jl`, replace the local basis selection in `product_square_ipeps` with:

```julia
function _product_state_vector(state::Symbol)
    if state in (:up, :z_up)
        return ComplexF64[1, 0]
    elseif state in (:down, :z_down)
        return ComplexF64[0, 1]
    elseif state === :x_plus
        return ComplexF64[inv(sqrt(2)), inv(sqrt(2))]
    else
        throw(ArgumentError("state must be :up, :down, :z_up, :z_down, or :x_plus"))
    end
end
```

Use `_product_state_vector(state)` wherever the constructor currently maps `:up`/`:down`.

In `src/SquarePEPS.jl`, add the same helper or import-free local equivalent:

```julia
function _product_state_vector(state::Symbol)
    if state in (:up, :z_up)
        return ComplexF64[1, 0]
    elseif state in (:down, :z_down)
        return ComplexF64[0, 1]
    elseif state === :x_plus
        return ComplexF64[inv(sqrt(2)), inv(sqrt(2))]
    else
        throw(ArgumentError("state must be :up, :down, :z_up, :z_down, or :x_plus"))
    end
end
```

- [ ] **Step 4: Run product-state tests**

Run:

```bash
julia --project=test -e 'using Pkg; Pkg.test("SquarePXPDynamics"; test_args=["test_square_ipeps.jl","test_square_peps.jl"])'
```

Expected: PASS.

- [ ] **Step 5: Commit product-state aliases**

Run:

```bash
git add src/SquareIPEPS.jl src/SquarePEPS.jl test/test_square_ipeps.jl test/test_square_peps.jl
git commit -m "feat: add benchmark product-state aliases"
```

Expected: commit succeeds.

## Task 5: Add TFIM Simple Observables And Exact-Limit Checks

**Files:**
- Modify: `src/Observables.jl`
- Modify: `src/SquarePXPDynamics.jl`
- Test: `test/test_tfim_observables.jl`

- [ ] **Step 1: Write failing TFIM observable tests**

Create `test/test_tfim_observables.jl`:

```julia
using Test
using SquarePXPDynamics

@testset "TFIM product-state observables" begin
    cell = PeriodicSquareUnitCell(10, 10)
    c = SquareCoord(5, 5)
    model = TFIMStarModel(1.0, 0.5)

    up = product_square_ipeps(cell; state = :z_up, maxdim = 1)
    down = product_square_ipeps(cell; state = :z_down, maxdim = 1)
    xplus = product_square_ipeps(cell; state = :x_plus, maxdim = 1)

    @test local_z_simple(up, c) ≈ 1.0 atol = 1e-12
    @test local_z_simple(down, c) ≈ -1.0 atol = 1e-12
    @test local_x_simple(xplus, c) ≈ 1.0 atol = 1e-12
    @test local_y_simple(xplus, c) ≈ 0.0 atol = 1e-12
    @test nearest_neighbor_zz_simple(up, c, :right) ≈ 1.0 atol = 1e-12
    @test nearest_neighbor_zz_simple(up, c, :up) ≈ 1.0 atol = 1e-12
    @test tfim_energy_density_decomposed_simple(up, model) ≈ -2.0 atol = 1e-12
    @test tfim_energy_density_star_simple(up, model) ≈ -2.0 atol = 1e-12

    summary = measure_tfim_simple(up, model)
    @test summary.mean_x ≈ 0.0 atol = 1e-12
    @test summary.mean_y ≈ 0.0 atol = 1e-12
    @test summary.mean_z ≈ 1.0 atol = 1e-12
    @test summary.zz_right ≈ 1.0 atol = 1e-12
    @test summary.zz_up ≈ 1.0 atol = 1e-12
    @test summary.energy_density_discrepancy ≈ 0.0 atol = 1e-12
    @test summary.max_imag_abs ≤ 1e-10
end

@testset "TFIM h=0 stationary product observables" begin
    cell = PeriodicSquareUnitCell(10, 10)
    model = TFIMStarModel(1.0, 0.0)
    psi = product_square_ipeps(cell; state = :z_up, maxdim = 1)
    before = measure_tfim_simple(psi, model)
    evolve!(psi, 0.01; params = TrotterParams(0.01, 1, :real, 1, 1e-12), protocol = StaticModel(model))
    after = measure_tfim_simple(psi, model)
    @test after.mean_z ≈ before.mean_z atol = 1e-10
    @test after.zz_right ≈ before.zz_right atol = 1e-10
    @test after.zz_up ≈ before.zz_up atol = 1e-10
    @test after.energy_density_star ≈ before.energy_density_star atol = 1e-10
end

@testset "TFIM J=0 independent spin dynamics" begin
    cell = PeriodicSquareUnitCell(10, 10)
    h = 1.0
    model = TFIMStarModel(0.0, h)
    psi = product_square_ipeps(cell; state = :z_up, maxdim = 1)
    total_time = 0.02
    evolve!(psi, total_time; params = TrotterParams(0.01, 1, :real, 1, 1e-12), protocol = StaticModel(model))
    summary = measure_tfim_simple(psi, model)
    @test summary.mean_z ≈ cos(2h * total_time) atol = 1e-6 rtol = 1e-6
    @test abs(abs(summary.mean_y) - abs(sin(2h * total_time))) ≤ 1e-6
end
```

- [ ] **Step 2: Run tests to verify they fail**

Run:

```bash
julia --project=test -e 'using Pkg; Pkg.test("SquarePXPDynamics"; test_args=["test_tfim_observables.jl"])'
```

Expected: FAIL with undefined names such as `local_x_simple`, `measure_tfim_simple`, or `TFIMObservableSummary`.

- [ ] **Step 3: Implement local Pauli observables**

In `src/Observables.jl`, import TFIM tools:

```julia
using ..SpinOps: pauli_x, pauli_y, pauli_z
using ..StarModels: TFIMStarModel, star_hamiltonian
```

Add exports:

```julia
export local_x_simple, local_y_simple, local_z_simple, nearest_neighbor_zz_simple
export tfim_energy_density_star_simple, tfim_energy_density_decomposed_simple
export TFIMObservableSummary, measure_tfim_simple
```

Add helpers:

```julia
function _local_operator_simple(psi::SquareIPEPSState, c::SquareCoord, O)::ComplexF64
    A = _site_patch_tensor(psi, c)
    p = only(siteinds(A))
    return _expectation_from_patch(A, O, (p,))
end

local_x_simple(psi::SquareIPEPSState, c::SquareCoord)::Float64 =
    _real_expectation(_local_operator_simple(psi, c, pauli_x()))

local_y_simple(psi::SquareIPEPSState, c::SquareCoord)::Float64 =
    _real_expectation(_local_operator_simple(psi, c, pauli_y()))

local_z_simple(psi::SquareIPEPSState, c::SquareCoord)::Float64 =
    _real_expectation(_local_operator_simple(psi, c, pauli_z()))
```

Use the existing two-site patch pattern from `nearest_neighbor_density_simple` to add:

```julia
function nearest_neighbor_zz_simple(
    psi::SquareIPEPSState,
    c::SquareCoord,
    dir::Symbol,
)::Float64
    dir in (:right, :up) ||
        throw(ArgumentError("dir must be :right or :up"))
    n = neighbor(psi.unitcell, c, dir)
    theta, phys = _two_site_patch_tensor(psi, c, n)
    value = _expectation_from_patch(theta, kron(pauli_z(), pauli_z()), phys)
    return _real_expectation(value)
end
```

- [ ] **Step 4: Implement TFIM summary and energy functions**

Add:

```julia
function _mean_over_reps(psi::SquareIPEPSState, f)
    reps = psi.unitcell.reps
    return sum(f(c) for c in reps) / length(reps)
end

function _mean_complex_over_reps(psi::SquareIPEPSState, f)
    reps = psi.unitcell.reps
    return sum(f(c) for c in reps) / length(reps)
end

function tfim_energy_density_star_simple(psi::SquareIPEPSState, model::TFIMStarModel)::Float64
    H = star_hamiltonian(model)
    value = _mean_complex_over_reps(psi, c -> star_expectation_simple(psi, c, H))
    return _real_expectation(value)
end

function tfim_energy_density_decomposed_simple(psi::SquareIPEPSState, model::TFIMStarModel)::Float64
    mean_x = _mean_over_reps(psi, c -> local_x_simple(psi, c))
    zz_right = _mean_over_reps(psi, c -> nearest_neighbor_zz_simple(psi, c, :right))
    zz_up = _mean_over_reps(psi, c -> nearest_neighbor_zz_simple(psi, c, :up))
    return -model.h * mean_x - model.J * (zz_right + zz_up)
end

struct TFIMObservableSummary
    mean_x::Float64
    mean_y::Float64
    mean_z::Float64
    z_even::Float64
    z_odd::Float64
    zz_right::Float64
    zz_up::Float64
    energy_density_star::Float64
    energy_density_decomposed::Float64
    energy_density_discrepancy::Float64
    x_imag_abs::Float64
    y_imag_abs::Float64
    z_imag_abs::Float64
    zz_imag_abs::Float64
    energy_imag_abs::Float64
    max_imag_abs::Float64
    mean_bond_entropy::Float64
    max_bond_entropy::Float64
end
```

Add `measure_tfim_simple`:

```julia
function measure_tfim_simple(psi::SquareIPEPSState, model::TFIMStarModel)::TFIMObservableSummary
    reps = psi.unitcell.reps
    raw_x = [ _local_operator_simple(psi, c, pauli_x()) for c in reps ]
    raw_y = [ _local_operator_simple(psi, c, pauli_y()) for c in reps ]
    raw_z = [ _local_operator_simple(psi, c, pauli_z()) for c in reps ]
    mean_x = _real_expectation(sum(raw_x) / length(raw_x))
    mean_y = _real_expectation(sum(raw_y) / length(raw_y))
    mean_z = _real_expectation(sum(raw_z) / length(raw_z))
    even = [raw_z[i] for (i, c) in enumerate(reps) if iseven(c.x + c.y)]
    odd = [raw_z[i] for (i, c) in enumerate(reps) if isodd(c.x + c.y)]
    z_even = isempty(even) ? mean_z : _real_expectation(sum(even) / length(even))
    z_odd = isempty(odd) ? mean_z : _real_expectation(sum(odd) / length(odd))
    zz_right = _mean_over_reps(psi, c -> nearest_neighbor_zz_simple(psi, c, :right))
    zz_up = _mean_over_reps(psi, c -> nearest_neighbor_zz_simple(psi, c, :up))
    H = star_hamiltonian(model)
    raw_energy = _mean_complex_over_reps(psi, c -> star_expectation_simple(psi, c, H))
    energy_star = _real_expectation(raw_energy)
    energy_decomposed = -model.h * mean_x - model.J * (zz_right + zz_up)
    x_imag = maximum(abs, imag.(raw_x))
    y_imag = maximum(abs, imag.(raw_y))
    z_imag = maximum(abs, imag.(raw_z))
    summary = TFIMObservableSummary(
        mean_x,
        mean_y,
        mean_z,
        z_even,
        z_odd,
        zz_right,
        zz_up,
        energy_star,
        energy_decomposed,
        abs(energy_star - energy_decomposed),
        x_imag,
        y_imag,
        z_imag,
        0.0,
        abs(imag(raw_energy)),
        maximum((x_imag, y_imag, z_imag, abs(imag(raw_energy)))),
        mean_bond_entropy(psi),
        max_bond_entropy(psi),
    )
    all(isfinite, Tuple(summary)) || throw(ArgumentError("TFIM observable summary must be finite"))
    return summary
end
```

If `Tuple(summary)` is not supported for the struct, replace that line with:

```julia
    all(isfinite, (getfield(summary, f) for f in fieldnames(TFIMObservableSummary))) ||
        throw(ArgumentError("TFIM observable summary must be finite"))
```

- [ ] **Step 5: Import and export TFIM observables at top level**

Patch `src/SquarePXPDynamics.jl`:

```julia
using .Observables:
    local_x_simple,
    local_y_simple,
    local_z_simple,
    nearest_neighbor_zz_simple,
    tfim_energy_density_star_simple,
    tfim_energy_density_decomposed_simple,
    TFIMObservableSummary,
    measure_tfim_simple

export local_x_simple, local_y_simple, local_z_simple, nearest_neighbor_zz_simple
export tfim_energy_density_star_simple, tfim_energy_density_decomposed_simple
export TFIMObservableSummary, measure_tfim_simple
```

- [ ] **Step 6: Run TFIM observable tests**

Run:

```bash
julia --project=test -e 'using Pkg; Pkg.test("SquarePXPDynamics"; test_args=["test_tfim_observables.jl"])'
```

Expected: PASS.

- [ ] **Step 7: Run public docs tests**

Run:

```bash
julia --project=test -e 'using Pkg; Pkg.test("SquarePXPDynamics"; test_args=["test_public_docs.jl"])'
```

Expected: PASS after docstrings are added for exported observable names.

- [ ] **Step 8: Commit TFIM observables**

Run:

```bash
git add src/Observables.jl src/SquarePXPDynamics.jl test/test_tfim_observables.jl test/test_public_docs.jl
git commit -m "feat: add TFIM simple observables"
```

Expected: commit succeeds.

## Task 6: Add Finite Schedule Reference Tests

**Files:**
- Create: `test/test_tfim_schedule_reference.jl`
- Modify: `src/IPEPSEvolution.jl` only if tests expose a schedule defect

- [ ] **Step 1: Write finite schedule reference tests**

Create `test/test_tfim_schedule_reference.jl`:

```julia
using Test
using LinearAlgebra
using SquarePXPDynamics

function _finite_site_index(x, y, L)
    return mod1(x, L) + (mod1(y, L) - 1) * L
end

function _finite_star_sites(center, L)
    x, y = center
    return (
        _finite_site_index(x, y, L),
        _finite_site_index(x + 1, y, L),
        _finite_site_index(x, y + 1, L),
        _finite_site_index(x - 1, y, L),
        _finite_site_index(x, y - 1, L),
    )
end

@testset "trotter coefficient sums" begin
    p1 = TrotterParams(0.1, 1, :real, 2, 1e-12)
    p2 = TrotterParams(0.1, 2, :real, 2, 1e-12)
    @test sum(step for (_, step) in trotter_sequence(p1)) ≈ 5 * p1.dt
    @test [color for (color, _) in trotter_sequence(p1)] == collect(1:5)
    @test [color for (color, _) in trotter_sequence(p2)] == [1, 2, 3, 4, 5, 4, 3, 2, 1]
    @test sum(step for (color, step) in trotter_sequence(p2) if color == 1) ≈ p2.dt
    @test sum(step for (color, step) in trotter_sequence(p2) if color == 5) ≈ p2.dt
end

@testset "five-color schedule has disjoint finite periodic stars" begin
    cell = PeriodicSquareUnitCell(10, 10)
    covered_centers = Set{Tuple{Int,Int}}()
    for color in 1:5
        centers = update_centers(cell, color)
        @test stars_are_disjoint_mod_unitcell(cell, centers)
        occupied = Set{Tuple{Int,Int}}()
        for c in centers
            push!(covered_centers, (wrap(cell, c).x, wrap(cell, c).y))
            for s in square_star_sites(c)
                wrapped = wrap(cell, s)
                key = (wrapped.x, wrapped.y)
                @test !(key in occupied)
                push!(occupied, key)
            end
        end
    end
    @test length(covered_centers) == length(cell.reps)
    @test covered_centers == Set((c.x, c.y) for c in cell.reps)
end

@testset "finite reference center mapping" begin
    L = 5
    center = (3, 3)
    @test _finite_star_sites(center, L) == (13, 14, 18, 12, 8)
    @test _finite_star_sites((1, 1), L) == (1, 2, 6, 5, 21)
end
```

- [ ] **Step 2: Run finite schedule tests**

Run:

```bash
julia --project=test -e 'using Pkg; Pkg.test("SquarePXPDynamics"; test_args=["test_tfim_schedule_reference.jl"])'
```

Expected: PASS. If a test fails, inspect whether the test assumption or implementation differs from the existing five-color convention before editing production code.

- [ ] **Step 3: Commit schedule reference**

Run:

```bash
git add test/test_tfim_schedule_reference.jl src/IPEPSEvolution.jl
git commit -m "test: add finite TFIM schedule reference"
```

Expected: commit succeeds. If `src/IPEPSEvolution.jl` was not modified, omit it from `git add`.

## Task 7: Add Benchmark Runner And Time-Series Serialization

**Files:**
- Create: `src/Benchmarks.jl`
- Modify: `src/SquarePXPDynamics.jl`
- Test: `test/test_benchmarks.jl`
- Modify: `test/test_public_docs.jl`

- [ ] **Step 1: Write failing benchmark tests**

Create `test/test_benchmarks.jl`:

```julia
using Test
using JSON3
using SquarePXPDynamics

@testset "benchmark result records time series" begin
    spec = BenchmarkSpec(
        "tfim-j0",
        StaticModel(TFIMStarModel(0.0, 1.0)),
        PeriodicSquareUnitCell(10, 10),
        :z_up,
        0.02,
        TrotterParams(0.01, 1, :real, 1, 1e-12),
        1,
    )
    result = run_benchmark(spec; run_label = "unit-test")
    @test result.name == "tfim-j0"
    @test result.run_label == "unit-test"
    @test length(result.samples) == 3
    @test [s.step for s in result.samples] == [0, 1, 2]
    @test result.samples[end].time ≈ 0.02 atol = 1e-12
    @test result.final_state_summary === result.samples[end]
    @test result.samples[end].observables.mean_z ≈ cos(0.04) atol = 1e-6 rtol = 1e-6
end

@testset "benchmark JSON and CSV writers are deterministic" begin
    spec = BenchmarkSpec(
        "tfim-static",
        StaticModel(TFIMStarModel(1.0, 0.0)),
        PeriodicSquareUnitCell(10, 10),
        :z_up,
        0.01,
        TrotterParams(0.01, 1, :real, 1, 1e-12),
        1,
    )
    result = run_benchmark(spec; run_label = "serialize-test")
    dir = mktempdir()
    json_path = joinpath(dir, "result.json")
    csv_path = joinpath(dir, "result.csv")
    write_benchmark_json(result, json_path)
    write_benchmark_csv([result], csv_path)

    parsed = JSON3.read(read(json_path, String))
    @test parsed[:name] == "tfim-static"
    @test parsed[:run_label] == "serialize-test"
    @test parsed[:metadata][:observable_source] == "simple"
    @test length(parsed[:samples]) == 2

    csv = read(csv_path, String)
    @test startswith(csv, "name,run_label,observable_source,step,time,J,h,D,dt,order")
    @test occursin("tfim-static,serialize-test,simple,0,0.0,1.0,0.0", csv)
    @test occursin("tfim-static,serialize-test,simple,1,0.01,1.0,0.0", csv)
end

@testset "benchmark validation rejects unsupported inputs" begin
    @test_throws ArgumentError BenchmarkSpec(
        "bad-measure",
        StaticModel(TFIMStarModel(1.0, 0.0)),
        PeriodicSquareUnitCell(10, 10),
        :z_up,
        0.01,
        TrotterParams(0.01, 1, :real, 1, 1e-12),
        0,
    )
    @test_throws ArgumentError BenchmarkSpec(
        "bad-state",
        StaticModel(TFIMStarModel(1.0, 0.0)),
        PeriodicSquareUnitCell(10, 10),
        :bad,
        0.01,
        TrotterParams(0.01, 1, :real, 1, 1e-12),
        1,
    )
end
```

- [ ] **Step 2: Run tests to verify they fail**

Run:

```bash
julia --project=test -e 'using Pkg; Pkg.test("SquarePXPDynamics"; test_args=["test_benchmarks.jl"])'
```

Expected: FAIL with undefined names such as `BenchmarkSpec` or `run_benchmark`.

- [ ] **Step 3: Implement `src/Benchmarks.jl`**

Create `src/Benchmarks.jl`:

```julia
module Benchmarks

using JSON3

using ..SquareUnitCells: PeriodicSquareUnitCell
using ..SquareIPEPS: product_square_ipeps
using ..StarModels: AbstractModelProtocol, StaticModel, TFIMStarModel, model_at
using ..IPEPSEvolution: TrotterParams, EvolutionLog, evolve!
using ..Observables: TFIMObservableSummary, measure_tfim_simple

export BenchmarkSpec,
    BenchmarkMetadata,
    EvolutionDiagnostics,
    BenchmarkSample,
    BenchmarkResult,
    run_benchmark,
    write_benchmark_json,
    write_benchmark_csv

const BENCHMARK_INITIAL_STATES = (:z_up, :z_down, :x_plus)

struct BenchmarkSpec{P<:AbstractModelProtocol,C}
    name::String
    protocol::P
    cell::C
    initial_state::Symbol
    total_time::Float64
    trotter::TrotterParams
    measure_every::Int
    function BenchmarkSpec(
        name::AbstractString,
        protocol::P,
        cell::C,
        initial_state::Symbol,
        total_time::Real,
        trotter::TrotterParams,
        measure_every::Integer,
    ) where {P<:AbstractModelProtocol,C}
        isempty(name) && throw(ArgumentError("benchmark name must be nonempty"))
        initial_state in BENCHMARK_INITIAL_STATES ||
            throw(ArgumentError("initial_state must be :z_up, :z_down, or :x_plus"))
        total = Float64(total_time)
        isfinite(total) && total >= 0 ||
            throw(ArgumentError("total_time must be finite and nonnegative"))
        cadence = Int(measure_every)
        cadence >= 1 || throw(ArgumentError("measure_every must be at least 1"))
        return new{P,C}(String(name), protocol, cell, initial_state, total, trotter, cadence)
    end
end

struct BenchmarkMetadata
    observable_source::Symbol
    protocol_type::String
    model_type::String
    J::Union{Nothing,Float64}
    h::Union{Nothing,Float64}
    cell_width::Int
    cell_height::Int
    initial_state::Symbol
    total_time::Float64
    dt::Float64
    order::Int
    evolution::Symbol
    maxdim::Int
    cutoff::Float64
    measure_every::Int
end

struct EvolutionDiagnostics
    max_truncerr::Float64
    mean_bond_entropy::Float64
    max_bond_entropy::Float64
    log_norm_before::Float64
    log_norm_after::Float64
    log_norm_delta::Float64
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

function _finite_float(x, name)
    value = Float64(x)
    isfinite(value) || throw(ArgumentError("$name must be finite"))
    return value
end

function _diagnostics_from_log(log::EvolutionLog)
    return EvolutionDiagnostics(
        _finite_float(log.max_truncerr, "max_truncerr"),
        _finite_float(log.mean_bond_entropy, "mean_bond_entropy"),
        _finite_float(log.max_bond_entropy, "max_bond_entropy"),
        _finite_float(log.log_norm_before, "log_norm_before"),
        _finite_float(log.log_norm_after, "log_norm_after"),
        _finite_float(log.log_norm_delta, "log_norm_delta"),
    )
end

function _zero_diagnostics()
    return EvolutionDiagnostics(0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
end

function _metadata(spec::BenchmarkSpec)
    model = model_at(spec.protocol, 0.0, 0)
    J = model isa TFIMStarModel ? Float64(model.J) : nothing
    h = model isa TFIMStarModel ? Float64(model.h) : nothing
    return BenchmarkMetadata(
        :simple,
        string(typeof(spec.protocol)),
        string(typeof(model)),
        J,
        h,
        spec.cell.width,
        spec.cell.height,
        spec.initial_state,
        spec.total_time,
        spec.trotter.dt,
        spec.trotter.order,
        spec.trotter.evolution,
        spec.trotter.maxdim,
        spec.trotter.cutoff,
        spec.measure_every,
    )
end

function _tfim_model(spec::BenchmarkSpec, time, step)
    model = model_at(spec.protocol, time, step)
    model isa TFIMStarModel ||
        throw(ArgumentError("v1 benchmark runner supports TFIMStarModel protocols only"))
    return model
end

function run_benchmark(spec::BenchmarkSpec; run_label::AbstractString = "manual")
    psi = product_square_ipeps(spec.cell; state = spec.initial_state, maxdim = spec.trotter.maxdim)
    samples = BenchmarkSample[]
    model0 = _tfim_model(spec, 0.0, 0)
    push!(samples, BenchmarkSample(0, 0.0, measure_tfim_simple(psi, model0), _zero_diagnostics()))
    nsteps = round(Int, spec.total_time / spec.trotter.dt)
    isapprox(nsteps * spec.trotter.dt, spec.total_time; atol = 1e-12, rtol = 1e-10) ||
        throw(ArgumentError("total_time must be an integer multiple of trotter.dt"))
    for step in 1:nsteps
        log = evolve!(psi, spec.trotter.dt; params = spec.trotter, protocol = spec.protocol)
        if step % spec.measure_every == 0 || step == nsteps
            time = step * spec.trotter.dt
            model = _tfim_model(spec, time, step)
            push!(samples, BenchmarkSample(step, time, measure_tfim_simple(psi, model), _diagnostics_from_log(log)))
        end
    end
    isempty(samples) && throw(ArgumentError("benchmark produced no samples"))
    return BenchmarkResult(spec.name, String(run_label), _metadata(spec), samples, samples[end])
end

_json_scalar(x::Symbol) = String(x)
_json_scalar(x) = x

function _summary_dict(obs::TFIMObservableSummary)
    return Dict(String(f) => _json_scalar(getfield(obs, f)) for f in fieldnames(TFIMObservableSummary))
end

function _diag_dict(diag::EvolutionDiagnostics)
    return Dict(String(f) => getfield(diag, f) for f in fieldnames(EvolutionDiagnostics))
end

function _metadata_dict(meta::BenchmarkMetadata)
    return Dict(String(f) => _json_scalar(getfield(meta, f)) for f in fieldnames(BenchmarkMetadata))
end

function _sample_dict(sample::BenchmarkSample)
    return Dict(
        "step" => sample.step,
        "time" => sample.time,
        "observables" => _summary_dict(sample.observables),
        "diagnostics" => _diag_dict(sample.diagnostics),
    )
end

function _result_dict(result::BenchmarkResult)
    return Dict(
        "name" => result.name,
        "run_label" => result.run_label,
        "metadata" => _metadata_dict(result.metadata),
        "samples" => [_sample_dict(s) for s in result.samples],
        "final_state_summary" => _sample_dict(result.final_state_summary),
    )
end

function write_benchmark_json(result::BenchmarkResult, path::AbstractString)
    open(path, "w") do io
        JSON3.write(io, _result_dict(result))
        println(io)
    end
    return path
end

const CSV_HEADER = [
    "name", "run_label", "observable_source", "step", "time", "J", "h", "D", "dt", "order",
    "mean_x", "mean_y", "mean_z", "zz_right", "zz_up",
    "energy_star", "energy_decomposed", "energy_discrepancy",
    "max_truncerr", "mean_bond_entropy", "max_bond_entropy", "lognorm_delta",
]

function _csv_row(result::BenchmarkResult, sample::BenchmarkSample)
    meta = result.metadata
    obs = sample.observables
    diag = sample.diagnostics
    return join((
        result.name,
        result.run_label,
        meta.observable_source,
        sample.step,
        sample.time,
        meta.J,
        meta.h,
        meta.maxdim,
        meta.dt,
        meta.order,
        obs.mean_x,
        obs.mean_y,
        obs.mean_z,
        obs.zz_right,
        obs.zz_up,
        obs.energy_density_star,
        obs.energy_density_decomposed,
        obs.energy_density_discrepancy,
        diag.max_truncerr,
        diag.mean_bond_entropy,
        diag.max_bond_entropy,
        diag.log_norm_delta,
    ), ",")
end

function write_benchmark_csv(results, path::AbstractString)
    open(path, "w") do io
        println(io, join(CSV_HEADER, ","))
        for result in results
            for sample in result.samples
                println(io, _csv_row(result, sample))
            end
        end
    end
    return path
end

end
```

- [ ] **Step 4: Include and export benchmark APIs**

Patch `src/SquarePXPDynamics.jl` after `include("IPEPSEvolution.jl")`:

```julia
include("Benchmarks.jl")
```

Add imports:

```julia
using .Benchmarks:
    BenchmarkSpec,
    BenchmarkMetadata,
    EvolutionDiagnostics,
    BenchmarkSample,
    BenchmarkResult,
    run_benchmark,
    write_benchmark_json,
    write_benchmark_csv
```

Add exports:

```julia
export BenchmarkSpec, BenchmarkMetadata, EvolutionDiagnostics, BenchmarkSample, BenchmarkResult
export run_benchmark, write_benchmark_json, write_benchmark_csv
```

- [ ] **Step 5: Run benchmark tests**

Run:

```bash
julia --project=test -e 'using Pkg; Pkg.test("SquarePXPDynamics"; test_args=["test_benchmarks.jl"])'
```

Expected: PASS.

- [ ] **Step 6: Run public docs tests**

Run:

```bash
julia --project=test -e 'using Pkg; Pkg.test("SquarePXPDynamics"; test_args=["test_public_docs.jl"])'
```

Expected: PASS after docstrings are added for exported benchmark names.

- [ ] **Step 7: Commit benchmark runner**

Run:

```bash
git add src/Benchmarks.jl src/SquarePXPDynamics.jl test/test_benchmarks.jl test/test_public_docs.jl
git commit -m "feat: add TFIM benchmark runner"
```

Expected: commit succeeds.

## Task 8: Full Regression And Documentation Update

**Files:**
- Modify: `docs/superpowers/specs/2026-05-15-infinite-tfim-benchmark-design.md` only if implementation decisions diverged
- Modify: `README.md` if it already contains package usage examples

- [ ] **Step 1: Run focused new test files**

Run:

```bash
julia --project=test -e 'using Pkg; Pkg.test("SquarePXPDynamics"; test_args=["test_star_models.jl","test_tfim_observables.jl","test_tfim_schedule_reference.jl","test_benchmarks.jl"])'
```

Expected: PASS.

- [ ] **Step 2: Run full test suite**

Run:

```bash
julia --project=test -e 'using Pkg; Pkg.test("SquarePXPDynamics")'
```

Expected: PASS.

- [ ] **Step 3: Run a manual smoke benchmark**

Run:

```bash
julia --project=. -e 'using SquarePXPDynamics; spec = BenchmarkSpec("manual-tfim-j0", StaticModel(TFIMStarModel(0.0, 1.0)), PeriodicSquareUnitCell(10, 10), :z_up, 0.02, TrotterParams(0.01, 1, :real, 1, 1e-12), 1); result = run_benchmark(spec; run_label="manual-smoke"); write_benchmark_json(result, "tfim-smoke.json"); write_benchmark_csv([result], "tfim-smoke.csv"); println(length(result.samples)); println(result.samples[end].observables.mean_z)'
```

Expected output includes:

```text
3
```

and the second printed line is close to `cos(0.04) = 0.9992001066609779`.

- [ ] **Step 4: Remove manual smoke artifacts**

Run:

```bash
rm -f tfim-smoke.json tfim-smoke.csv
git status --short
```

Expected: no `tfim-smoke` files appear in git status.

- [ ] **Step 5: Update user-facing docs if README has usage examples**

If `README.md` already contains package usage examples, add this compact TFIM benchmark example:

````markdown
### TFIM Benchmark Smoke Run

```julia
using SquarePXPDynamics

spec = BenchmarkSpec(
    "tfim-j0",
    StaticModel(TFIMStarModel(0.0, 1.0)),
    PeriodicSquareUnitCell(10, 10),
    :z_up,
    0.02,
    TrotterParams(0.01, 1, :real, 1, 1e-12),
    1,
)

result = run_benchmark(spec; run_label = "local-smoke")
write_benchmark_json(result, "tfim-j0.json")
write_benchmark_csv([result], "tfim-j0.csv")
```

The v1 TFIM benchmark uses simple-update diagnostics. Treat these outputs as
implementation regression records, not CTM-quality physics estimates.
````

If `README.md` has no usage section, skip this edit and leave documentation in the design spec plus public docstrings.

- [ ] **Step 6: Run docs/public API tests after docs edits**

Run:

```bash
julia --project=test -e 'using Pkg; Pkg.test("SquarePXPDynamics"; test_args=["test_public_docs.jl"])'
```

Expected: PASS.

- [ ] **Step 7: Run git verification**

Run:

```bash
git diff --check
git status --short
```

Expected: `git diff --check` exits 0. `git status --short` lists only intentional implementation and documentation files.

- [ ] **Step 8: Commit final docs or verification-only changes**

If Step 5 changed docs, run:

```bash
git add README.md docs/superpowers/specs/2026-05-15-infinite-tfim-benchmark-design.md
git commit -m "docs: document TFIM benchmark smoke run"
```

If Step 5 did not change docs, do not create an empty commit.

## Self-Review Checklist

- Spec coverage: star-model abstraction is covered by Task 2; protocol/Trotter separation by Task 3; product-state benchmark starts by Task 4; TFIM observables, exact limits, imaginary residuals, and decomposed energy by Task 5; finite schedule reference with non-overlap and full-sweep coverage by Task 6; time-series benchmark records, `observable_source = :simple` metadata, and JSON3 output by Task 7; full regression and docs by Task 8.
- Red-flag scan: no unresolved implementation markers remain, and every code-changing step names exact files, commands, and expected outcomes.
- Type consistency: `TFIMStarModel`, `PXPStarModel`, `StaticModel`, `TrotterParams`, `BenchmarkSpec`, `BenchmarkSample`, `BenchmarkResult`, `TFIMObservableSummary`, and `EvolutionDiagnostics` are introduced before later tasks use them.
