# General Star Simple Update Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make projected 7-site PXP gates usable in the PEPS evolution path, with tests that pin down scheduling, unit-cell wrapping, truncation diagnostics, and the first fixed-bond-dimension projection contract.

**Architecture:** Start by hardening current tests around ordering and schedule semantics, then add a clear gate-builder evolution API so first- and second-order Trotter steps construct correct full/half-step gates. Implement the first non-product star-update path in small layers: exact `D=1` dense-star projection for product iPEPS, then the general Simple Update contraction/SVD path for `D>1`. Keep all work inside the root Julia project under `src/` and `test/`.

**Tech Stack:** Julia 1.12, ITensors.jl, LinearAlgebra, Test.

---

## Review Summary

- Dense model conventions are sound: `|up> = |0>`, `|down> = |1>`, center-first 7-site ordering, triangular neighbor order, `projected_gate = P_blockade * U`.
- The main blocker is `apply_star_gate_simple_update!`: it rejects non-product 7-site gates and does not use `maxdim`, `cutoff`, SVD truncation, or lambda updates.
- The second-order evolution API is ambiguous because it applies the same prebuilt gate on every color layer instead of constructing half-step gates.
- Current blockade diagnostics are dense-vector diagnostics only, not PEPS-level nearest-neighbor diagnostics.

## File Structure

- Modify: `src/SolvableModels.jl` to remove the misleading `:plus` alias.
- Modify: `src/Schedules.jl` to expose schedule layers with step multipliers.
- Modify: `src/Evolution.jl` to add Hamiltonian/gate-builder based evolution APIs.
- Modify: `src/SimpleUpdate.jl` to add non-product dense-star update support incrementally.
- Modify: `src/Observables.jl` to add PEPS-local blockade diagnostics.
- Modify: `src/TriangularPEPSDynamics.jl` to export new public helpers.
- Modify: `test/test_solvable_models.jl` for alias cleanup.
- Modify: `test/test_evolution.jl` for schedule application-count tests.
- Modify: `test/test_simple_update.jl` for non-product PXP and truncation tests.
- Modify: `test/test_observables.jl` for PEPS-level blockade diagnostics.
- Modify: `README.md` to document what is ScarFinder-facing and what remains approximate.

## Task 1: Fix Small Benchmark Ambiguity And Schedule Tests

**Files:**
- Modify: `src/SolvableModels.jl`
- Modify: `test/test_solvable_models.jl`
- Modify: `test/test_evolution.jl`

- [ ] **Step 1: Write failing tests for the benchmark alias and canonical colors**

Add to `test/test_solvable_models.jl`:

```julia
@test_throws ArgumentError cluster_center_z_expectation_exact(0.1; initial = :plus)
```

Add to `test/test_evolution.jl`:

```julia
@testset "color canonical centers match requested color" begin
    for color in 1:7
        @test star_color(color_canonical_center(color)) == color
    end
end
```

- [ ] **Step 2: Run targeted tests and verify failure**

Run:

```bash
julia --project=. -e 'using Pkg; Pkg.test(test_args=["solvable models", "evolution"])'
```

Expected: FAIL because `:plus` is still accepted.

- [ ] **Step 3: Remove the ambiguous alias**

Change `src/SolvableModels.jl`:

```julia
function cluster_center_z_expectation_exact(t::Real; initial::Symbol = :z_plus)
    if initial == :z_plus
        return cos(2t)
    else
        throw(ArgumentError("supported initial states: :z_plus"))
    end
end
```

- [ ] **Step 4: Run targeted tests and verify pass**

Run:

```bash
julia --project=. -e 'using Pkg; Pkg.test(test_args=["solvable models", "evolution"])'
```

Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add src/SolvableModels.jl test/test_solvable_models.jl test/test_evolution.jl
git commit -m "test: pin benchmark aliases and color centers"
```

## Task 2: Make Trotter Step Weights Explicit

**Files:**
- Modify: `src/Schedules.jl`
- Modify: `src/Evolution.jl`
- Modify: `src/TriangularPEPSDynamics.jl`
- Modify: `test/test_schedules.jl`
- Modify: `test/test_evolution.jl`

- [ ] **Step 1: Write failing schedule-layer tests**

Add to `test/test_schedules.jl`:

```julia
@test schedule_layers(:first) == [(color = c, scale = 1.0) for c in 1:7]
@test schedule_layers(:second) == vcat(
    [(color = c, scale = 0.5) for c in 1:7],
    [(color = c, scale = 0.5) for c in 7:-1:1],
)
@test_throws ArgumentError schedule_layers(:bad)
```

- [ ] **Step 2: Implement `schedule_layers`**

Add to `src/Schedules.jl`:

```julia
export first_order_colors, second_order_colors, schedule_layers

function schedule_layers(order::Symbol)
    if order === :first
        return [(color = c, scale = 1.0) for c in first_order_colors()]
    elseif order === :second
        return vcat(
            [(color = c, scale = 0.5) for c in first_order_colors()],
            [(color = c, scale = 0.5) for c in reverse(first_order_colors())],
        )
    else
        throw(ArgumentError("order must be :first or :second"))
    end
end
```

Update `src/TriangularPEPSDynamics.jl` to import and export `schedule_layers`.

- [ ] **Step 3: Add a gate-builder evolution API**

Add this method to `src/Evolution.jl` while keeping the existing `evolve_step!(state, gate; ...)` method for compatibility:

```julia
using ..Gates: dense_gate, projected_gate
using ..Schedules: schedule_layers

function evolve_step!(state::TriangularIPEPS,
                      H::AbstractMatrix,
                      dt::Real;
                      order::Symbol = :second,
                      update::Symbol = :simple,
                      evolution::Symbol = :real,
                      projected::Bool = false,
                      projector::Union{Nothing,AbstractMatrix} = nothing)
    update === :simple || throw(ArgumentError("update must be :simple"))
    for layer in schedule_layers(order)
        step = dt * layer.scale
        gate = if projected
            projector === nothing ?
                projected_gate(H, step; evolution) :
                projected_gate(H, step; evolution, projector)
        else
            dense_gate(H, step; evolution)
        end
        apply_star_gate_simple_update!(state, gate, color_canonical_center(layer.color))
    end
    return state
end
```

- [ ] **Step 4: Write evolution count tests using X rotations**

Add to `test/test_evolution.jl`:

```julia
@testset "Hamiltonian evolution uses first/second order step weights" begin
    α = 0.03
    Xsum = sum(embed_one_site(pauli_x(), site, 7) for site in 1:7)

    first = product_ipeps(OneSiteUnitCell(), :down; D = 1)
    evolve_step!(first, Xsum, α; order = :first, update = :simple)
    @test real(local_expectation(first, Coord(0, 0), pauli_z())) ≈ -cos(14 * α) atol = 1e-10

    second = product_ipeps(OneSiteUnitCell(), :down; D = 1)
    evolve_step!(second, Xsum, α; order = :second, update = :simple)
    @test real(local_expectation(second, Coord(0, 0), pauli_z())) ≈ -cos(14 * α) atol = 1e-10
end
```

These expected values reflect the current translational one-site schedule: each color applies one representative update after all seven star factors wrap to the same tensor. Second order has twice as many layers, but each layer uses a half step, so the total single-site rotation matches first order for this commuting test Hamiltonian.

- [ ] **Step 5: Run targeted tests**

Run:

```bash
julia --project=. -e 'using Pkg; Pkg.test(test_args=["schedules", "evolution"])'
```

Expected: PASS.

- [ ] **Step 6: Commit**

```bash
git add src/Schedules.jl src/Evolution.jl src/TriangularPEPSDynamics.jl test/test_schedules.jl test/test_evolution.jl
git commit -m "feat: add weighted trotter schedule layers"
```

## Task 3: Add Exact `D=1` Non-Product Star Update Oracle

**Files:**
- Modify: `src/SimpleUpdate.jl`
- Modify: `test/test_simple_update.jl`

- [ ] **Step 1: Write failing projected PXP `D=1` test**

Replace the current non-product throw assertion in `test/test_simple_update.jl` with:

```julia
@testset "projected PXP gate updates D=1 three-site product exactly" begin
    t = 0.17
    H = pxp_star_hamiltonian(projector_down(), pauli_x())
    Uproj = projected_gate(H, t; evolution = :real)
    state = product_ipeps(ThreeSiteUnitCell(), :down; D = 1)

    diag = apply_star_gate_simple_update!(state, Uproj, Coord(0, 0); maxdim = 1)

    @test diag isa SimpleUpdateDiagnostics
    @test diag.discarded_weight ≈ 0 atol = 1e-12
    @test real(local_expectation(state, Coord(0, 0), pauli_z())) ≈ -cos(2t) atol = 1e-10
    @test real(local_expectation(state, Coord(1, 0), pauli_z())) ≈ -1 atol = 1e-10
    @test real(local_expectation(state, Coord(2, 0), pauli_z())) ≈ -1 atol = 1e-10
end
```

- [ ] **Step 2: Run targeted test and verify failure**

Run:

```bash
julia --project=. -e 'using Pkg; Pkg.test(test_args=["simple update"])'
```

Expected: FAIL with the current non-product gate error.

- [ ] **Step 3: Implement exact `D=1` dense-star update path**

In `src/SimpleUpdate.jl`, before the non-product error, add a branch for all affected representative tensors having only dimension-1 bond legs:

```julia
if all(dim(idx) == 1 for rep in keys(state.tensors) for idx in inds(state.tensors[rep]) if idx != state.phys_inds[rep])
    local_vectors = Vector{Vector{ComplexF64}}()
    for sc in star
        rep = wrap_coord(state.unitcell, sc)
        T = state.tensors[rep]
        ph = state.phys_inds[rep]
        push!(local_vectors, ComplexF64[T[ph => 1, (idx => 1 for idx in inds(T) if idx != ph)...],
                                        T[ph => 2, (idx => 1 for idx in inds(T) if idx != ph)...]])
    end
    psi = local_vectors[1]
    for v in local_vectors[2:end]
        psi = kron(psi, v)
    end
    phi = G * psi
    tensor_phi = reshape(phi, ntuple(_ -> 2, _STAR_NSITES)...)
    for (rep, positions) in _star_positions_by_rep(state, star)
        rho = _one_site_density_from_star_tensor(tensor_phi, positions[1])
        vals, vecs = eigen(Hermitian(rho))
        v = vecs[:, argmax(vals)]
        _set_d1_site_vector!(state, rep, ComplexF64.(v))
    end
    dims = [dim(state.bond_inds[b]) for b in affected]
    return SimpleUpdateDiagnostics(0.0, affected, dims)
end
```

Add private helpers in the same file:

```julia
function _star_positions_by_rep(state::TriangularIPEPS, star)
    grouped = Dict{Coord,Vector{Int}}()
    for (i, sc) in enumerate(star)
        rep = wrap_coord(state.unitcell, sc)
        push!(get!(grouped, rep, Int[]), i)
    end
    return grouped
end

function _one_site_density_from_star_tensor(tensor_phi, position::Int)
    perm = (position, (i for i in 1:_STAR_NSITES if i != position)...)
    psi = reshape(permutedims(tensor_phi, perm), 2, :)
    rho = psi * psi'
    tr = real(sum(diag(rho)))
    tr == 0 && return Matrix{ComplexF64}(I, 2, 2) / 2
    return rho / tr
end

function _set_d1_site_vector!(state::TriangularIPEPS, rep::Coord, v::Vector{ComplexF64})
    ph = state.phys_inds[rep]
    binds = Tuple(idx for idx in inds(state.tensors[rep]) if idx != ph)
    T = ITensor(ComplexF64, ph, binds...)
    nrm = norm(v)
    nrm == 0 && throw(ArgumentError("cannot set zero local vector"))
    v = v / nrm
    for k in 1:2
        T[ph => k, (binds[d] => 1 for d in eachindex(binds))...] = v[k]
    end
    state.tensors[rep] = T
    return nothing
end
```

- [ ] **Step 4: Run targeted tests**

Run:

```bash
julia --project=. -e 'using Pkg; Pkg.test(test_args=["simple update"])'
```

Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add src/SimpleUpdate.jl test/test_simple_update.jl
git commit -m "feat: support exact d1 nonproduct star updates"
```

## Task 4: Define And Test Simple Update Diagnostics Contract

**Files:**
- Modify: `src/SimpleUpdate.jl`
- Modify: `test/test_simple_update.jl`

- [ ] **Step 1: Write diagnostics tests**

Add to `test/test_simple_update.jl`:

```julia
@testset "simple update diagnostics report affected bonds and dimensions" begin
    state = random_ipeps(OneSiteUnitCell(), 2; seed = 7)
    I128 = Matrix{ComplexF64}(I, 128, 128)
    diag = apply_star_gate_simple_update!(state, I128, Coord(0, 0); maxdim = 2)
    @test Set(diag.affected_bonds) == Set((Coord(0, 0), d) for d in 1:6)
    @test diag.output_bond_dims == fill(2, 6)
    for d in 1:6
        λ = bond_lambda(state, Coord(0, 0), d)
        @test all(λ .>= 0)
        @test norm(λ) ≈ sqrt(length(λ))
    end
end
```

- [ ] **Step 2: Run diagnostics test**

Run:

```bash
julia --project=. -e 'using Pkg; Pkg.test(test_args=["simple update"])'
```

Expected: PASS before general SVD implementation; this pins current identity behavior.

- [ ] **Step 3: Commit**

```bash
git add test/test_simple_update.jl
git commit -m "test: pin simple update diagnostics contract"
```

## Task 5: Implement General Dense-Star Simple Update Skeleton For `D>1`

**Files:**
- Modify: `src/SimpleUpdate.jl`
- Modify: `test/test_simple_update.jl`

- [ ] **Step 1: Add a failing fixed-`D` non-product smoke test**

Add to `test/test_simple_update.jl`:

```julia
@testset "projected PXP simple update preserves fixed D on random D=2 state" begin
    state = random_ipeps(OneSiteUnitCell(), 2; seed = 11)
    H = pxp_star_hamiltonian(projector_down(), pauli_x())
    Uproj = projected_gate(H, 0.02; evolution = :real)

    diag = apply_star_gate_simple_update!(state, Uproj, Coord(0, 0); maxdim = 2, cutoff = 1e-12)

    @test diag isa SimpleUpdateDiagnostics
    @test diag.discarded_weight >= 0
    @test all(dim(bond_index(state, Coord(0, 0), d)) <= 2 for d in 1:6)
    @test all(length(bond_lambda(state, Coord(0, 0), d)) <= 2 for d in 1:6)
    @test all(all(bond_lambda(state, Coord(0, 0), d) .>= 0) for d in 1:6)
end
```

- [ ] **Step 2: Run targeted test and verify failure**

Run:

```bash
julia --project=. -e 'using Pkg; Pkg.test(test_args=["simple update"])'
```

Expected: FAIL because the `D>1` non-product path is not implemented.

- [ ] **Step 3: Implement the conservative `D>1` skeleton**

Implement a private `_apply_general_star_gate_simple_update!` in `src/SimpleUpdate.jl` with this initial contract:

```julia
function _apply_general_star_gate_simple_update!(state::TriangularIPEPS,
                                                G::Matrix{ComplexF64},
                                                center::Coord;
                                                cutoff::Real,
                                                maxdim::Union{Nothing,Integer})
    maxdim === nothing && throw(ArgumentError("maxdim is required for general star updates"))
    maxdim >= 1 || throw(ArgumentError("maxdim must be >= 1"))

    # Implementation checkpoint A: reject one-site wrapped stars until repeated
    # representative writeback is designed explicitly.
    if state.unitcell isa OneSiteUnitCell
        throw(ArgumentError("general star updates for OneSiteUnitCell are not yet supported"))
    end

    star = star_sites(center)
    occurrence_tensors = [_copy_occurrence_tensor_with_lambdas(state, sc) for sc in star]
    updated_cluster = _apply_dense_gate_to_occurrence_cluster(occurrence_tensors, G)
    new_tensors, new_lambdas, discarded = _split_star_cluster(
        updated_cluster, state, center; cutoff = cutoff, maxdim = Int(maxdim)
    )

    for (rep, T) in new_tensors
        state.tensors[rep] = T
    end
    for (bond, λ) in new_lambdas
        state.lambdas[bond] = λ
    end

    affected = _affected_star_bonds(state, center)
    dims = [dim(state.bond_inds[b]) for b in affected]
    return SimpleUpdateDiagnostics(discarded, affected, dims)
end
```

Add these private helper signatures in the same task, with tests covering each one through the public `apply_star_gate_simple_update!` path:

```julia
_affected_star_bonds(state::TriangularIPEPS, center::Coord) -> Vector{Tuple{Coord,Int}}
_copy_occurrence_tensor_with_lambdas(state::TriangularIPEPS, c::Coord) -> ITensor
_apply_dense_gate_to_occurrence_cluster(occurrences::Vector{ITensor}, G::Matrix{ComplexF64}) -> ITensor
_split_star_cluster(cluster::ITensor, state::TriangularIPEPS, center::Coord; cutoff::Real, maxdim::Int) -> Tuple{Dict{Coord,ITensor},Dict{Tuple{Coord,Int},Vector{Float64}},Float64}
```

The first green implementation may support only `ThreeSiteUnitCell`; add an explicit test that `OneSiteUnitCell` raises the `ArgumentError` above for non-product `D>1` gates. Do not silently fall back to identity behavior.

- [ ] **Step 4: Run targeted tests**

Run:

```bash
julia --project=. -e 'using Pkg; Pkg.test(test_args=["simple update"])'
```

Expected: PASS for supported unit cells and clear `ArgumentError` tests for unsupported wrapped cases.

- [ ] **Step 5: Commit**

```bash
git add src/SimpleUpdate.jl test/test_simple_update.jl
git commit -m "feat: add general star simple update skeleton"
```

## Task 6: Add PEPS-Local Blockade Diagnostics

**Files:**
- Modify: `src/Observables.jl`
- Modify: `src/TriangularPEPSDynamics.jl`
- Modify: `test/test_observables.jl`

- [ ] **Step 1: Write diagnostics tests**

Add to `test/test_observables.jl`:

```julia
@testset "nearest-neighbor blockade diagnostics on product states" begin
    down = product_ipeps(OneSiteUnitCell(), :down; D = 1)
    @test local_blockade_violation(down, Coord(0, 0), 1) ≈ 0
    @test mean_blockade_violation(down, [Coord(0, 0)]) ≈ 0

    up = product_ipeps(OneSiteUnitCell(), :up; D = 1)
    @test local_blockade_violation(up, Coord(0, 0), 1) ≈ 1
    @test mean_blockade_violation(up, [Coord(0, 0)]) ≈ 1
end
```

- [ ] **Step 2: Implement product/local-environment diagnostics**

Add to `src/Observables.jl`:

```julia
export local_blockade_violation, mean_blockade_violation

function local_blockade_violation(state::TriangularIPEPS, c::Coord, d::Integer)
    pup_c = real(local_expectation(state, c, projector_up()))
    pup_n = real(local_expectation(state, neighbor(c, d), projector_up()))
    return pup_c * pup_n
end

function mean_blockade_violation(state::TriangularIPEPS, centers)
    vals = Float64[]
    for c in centers, d in 1:6
        push!(vals, local_blockade_violation(state, c, d))
    end
    return isempty(vals) ? 0.0 : sum(vals) / length(vals)
end
```

Import `neighbor` and `projector_up` in `src/Observables.jl`; export the new functions in `src/TriangularPEPSDynamics.jl`.

- [ ] **Step 3: Run targeted tests**

Run:

```bash
julia --project=. -e 'using Pkg; Pkg.test(test_args=["observables"])'
```

Expected: PASS.

- [ ] **Step 4: Commit**

```bash
git add src/Observables.jl src/TriangularPEPSDynamics.jl test/test_observables.jl
git commit -m "feat: add local blockade diagnostics"
```

## Task 7: End-To-End Projected PXP Smoke Test And Docs

**Files:**
- Modify: `test/test_evolution.jl`
- Modify: `README.md`

- [ ] **Step 1: Add projected PXP evolution smoke test**

Add to `test/test_evolution.jl`:

```julia
@testset "projected PXP evolution smoke test" begin
    state = product_ipeps(ThreeSiteUnitCell(), :down; D = 1)
    H = pxp_star_hamiltonian(projector_down(), pauli_x())
    evolve_step!(state, H, 0.01; order = :first, update = :simple, evolution = :real, projected = true)
    @test all(isfinite(real(local_expectation(state, c, pauli_z()))) for c in unit_cell_representatives(ThreeSiteUnitCell()))
    @test mean_blockade_violation(state, collect(unit_cell_representatives(ThreeSiteUnitCell()))) < 1e-10
end
```

- [ ] **Step 2: Update README boundary text**

Add a short section to `README.md`:

```markdown
## Current ScarFinder-Facing Boundary

The code supports dense 7-site projected PXP gates and an initial PEPS evolution path. The Simple Update backend is being built in stages: exact `D=1` non-product star updates first, then fixed-`D` SVD truncation with lambda updates. PEPS blockade diagnostics are currently local-environment diagnostics and should be treated as screening metrics until a stronger environment contraction is added.
```

- [ ] **Step 3: Run full suite**

Run:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: PASS.

- [ ] **Step 4: Commit**

```bash
git add test/test_evolution.jl README.md
git commit -m "test: add projected pxp evolution smoke test"
```

## Self-Review Checklist

- The plan keeps all Julia code under the root `src/` tree and tests under the root `test/` tree.
- The plan does not introduce a nested package or a nested Julia environment.
- The first implementation target is projected PXP evolution, not broad PEPS package generalization.
- Tests are added before behavior changes for each task.
- The highest-risk implementation, general `D>1` Simple Update, is isolated after schedule/API and exact `D=1` oracle tests.
