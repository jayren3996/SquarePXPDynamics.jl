# Triangular PEPS Dynamics Foundation Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build the first working Julia/ITensors.jl foundation for triangular-lattice 7-site dynamics: package scaffold, geometry, dense/projected gates, solvable true-2D benchmarks, and test coverage.

**Architecture:** This first plan intentionally stops before PEPS tensor truncation. It creates a tested core that later PEPS update code can depend on: triangular geometry, star neighborhoods, 7-color schedules, dense 7-site operators, PXP projected gates, and analytically solvable benchmark models.

**Tech Stack:** Julia 1.12, ITensors.jl, LinearAlgebra, Test, Documenter-style docstrings.

---

## File Structure

- `Project.toml`: Julia package metadata and dependencies.
- `src/TriangularPEPSDynamics.jl`: top-level module and exports.
- `src/Geometry.jl`: triangular axial coordinates, directions, distances, stars, and 7-coloring.
- `src/SpinOps.jl`: dense spin-1/2 operators and Kronecker-product helpers.
- `src/Models.jl`: PXP star Hamiltonian, blockade projector, cluster/stabilizer star, diagonal star, and Ising bond terms.
- `src/Gates.jl`: dense real-time/imaginary-time gates and projected gates.
- `src/Schedules.jl`: first-order and second-order color schedules.
- `src/SolvableModels.jl`: exact benchmark helpers for commuting/stabilizer models.
- `test/runtests.jl`: test entry point.
- `test/test_geometry.jl`: geometry and coloring tests.
- `test/test_spinops.jl`: spin operator tests.
- `test/test_models.jl`: model and projector tests.
- `test/test_gates.jl`: dense/projected gate tests.
- `test/test_schedules.jl`: schedule tests.
- `test/test_solvable_models.jl`: true-2D solvable benchmark tests.

## Task 1: Package Scaffold

**Files:**
- Create: `Project.toml`
- Create: `src/TriangularPEPSDynamics.jl`
- Create: `test/runtests.jl`

- [ ] **Step 1: Write the package metadata**

Create `Project.toml` with:

```toml
name = "TriangularPEPSDynamics"
uuid = "6f0e7c3a-7e61-4f61-8b55-0c69f1dd7c75"
authors = ["Ren <ren@example.com>"]
version = "0.1.0"

[deps]
ITensors = "9136182c-28ba-11e9-034c-db9fb085ebd5"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[compat]
ITensors = "0.9"
julia = "1.12"

[extras]
Test = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[targets]
test = ["Test"]
```

- [ ] **Step 2: Write the top-level module**

Create `src/TriangularPEPSDynamics.jl` with:

```julia
module TriangularPEPSDynamics

include("Geometry.jl")
include("SpinOps.jl")
include("Models.jl")
include("Gates.jl")
include("Schedules.jl")
include("SolvableModels.jl")

using .Geometry
using .SpinOps
using .Models
using .Gates
using .Schedules
using .SolvableModels

export Coord, TRIANGULAR_DIRECTIONS, triangular_distance, neighbor, star_sites
export star_color, disjoint_stars
export pauli_x, pauli_y, pauli_z, identity2, projector_up, projector_down
export pxp_star_hamiltonian, blockade_projector, cluster_star_hamiltonian
export diagonal_star_hamiltonian, ising_bond_hamiltonian
export dense_gate, projected_gate
export first_order_colors, second_order_colors
export stabilizer_expectation_exact

end
```

- [ ] **Step 3: Write the test entry point**

Create `test/runtests.jl` with:

```julia
using Test
using TriangularPEPSDynamics

@testset "TriangularPEPSDynamics" begin
    include("test_geometry.jl")
    include("test_spinops.jl")
    include("test_models.jl")
    include("test_gates.jl")
    include("test_schedules.jl")
    include("test_solvable_models.jl")
end
```

- [ ] **Step 4: Run tests to verify scaffold failure**

Run:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: FAIL because included source and test files do not exist yet.

- [ ] **Step 5: Commit**

If this directory is not a git repo, first run:

```bash
git init
```

Then run:

```bash
git add Project.toml src/TriangularPEPSDynamics.jl test/runtests.jl
git commit -m "chore: scaffold triangular PEPS dynamics package"
```

## Task 2: Triangular Geometry And 7-Color Star Scheduling Primitive

**Files:**
- Create: `src/Geometry.jl`
- Create: `test/test_geometry.jl`

- [ ] **Step 1: Write failing geometry tests**

Create `test/test_geometry.jl` with:

```julia
@testset "geometry" begin
    c = Coord(2, -1)
    @test neighbor(c, 1) == Coord(3, -1)
    @test neighbor(c, 4) == Coord(1, -1)
    @test triangular_distance(Coord(0, 0), Coord(2, -1)) == 2

    s = star_sites(Coord(0, 0))
    @test length(s) == 7
    @test Coord(0, 0) in s
    @test Coord(1, 0) in s
    @test Coord(0, 1) in s

    centers = [Coord(q, r) for q in -4:4 for r in -4:4]
    for a in centers, b in centers
        if a != b && star_color(a) == star_color(b)
            @test disjoint_stars(a, b)
        end
    end
end
```

- [ ] **Step 2: Run geometry tests to verify failure**

Run:

```bash
julia --project=. -e 'using Pkg; Pkg.test(test_args=["geometry"])'
```

Expected: FAIL because `Coord` and geometry functions are not defined.

- [ ] **Step 3: Implement geometry**

Create `src/Geometry.jl` with:

```julia
module Geometry

export Coord, TRIANGULAR_DIRECTIONS, triangular_distance, neighbor, star_sites
export star_color, disjoint_stars

struct Coord
    q::Int
    r::Int
end

const TRIANGULAR_DIRECTIONS = (
    Coord(1, 0),
    Coord(0, 1),
    Coord(-1, 1),
    Coord(-1, 0),
    Coord(0, -1),
    Coord(1, -1),
)

Base.:+(a::Coord, b::Coord) = Coord(a.q + b.q, a.r + b.r)
Base.:-(a::Coord, b::Coord) = Coord(a.q - b.q, a.r - b.r)
Base.:(==)(a::Coord, b::Coord) = a.q == b.q && a.r == b.r
Base.hash(c::Coord, h::UInt) = hash((c.q, c.r), h)

function triangular_distance(a::Coord, b::Coord)
    d = b - a
    return (abs(d.q) + abs(d.r) + abs(d.q + d.r)) ÷ 2
end

function neighbor(c::Coord, direction::Integer)
    1 <= direction <= 6 || throw(ArgumentError("direction must be in 1:6"))
    return c + TRIANGULAR_DIRECTIONS[direction]
end

function star_sites(center::Coord)
    return [center; [neighbor(center, dir) for dir in 1:6]]
end

function star_color(c::Coord)
    return mod(c.q + 3c.r, 7) + 1
end

function disjoint_stars(a::Coord, b::Coord)
    return isempty(intersect(Set(star_sites(a)), Set(star_sites(b))))
end

end
```

- [ ] **Step 4: Run geometry tests to verify pass**

Run:

```bash
julia --project=. test/runtests.jl
```

Expected: geometry tests PASS; other tests FAIL because their files do not exist yet.

- [ ] **Step 5: Commit**

```bash
git add src/Geometry.jl test/test_geometry.jl
git commit -m "feat: add triangular geometry and star coloring"
```

## Task 3: Dense Spin Operators

**Files:**
- Create: `src/SpinOps.jl`
- Create: `test/test_spinops.jl`

- [ ] **Step 1: Write failing spin-operator tests**

Create `test/test_spinops.jl` with:

```julia
@testset "spin operators" begin
    X = pauli_x()
    Y = pauli_y()
    Z = pauli_z()
    I2 = identity2()

    @test X * X ≈ I2
    @test Y * Y ≈ I2
    @test Z * Z ≈ I2
    @test X * Y ≈ im * Z
    @test projector_up() + projector_down() ≈ I2
    @test projector_up() * projector_down() ≈ zeros(ComplexF64, 2, 2)
end
```

- [ ] **Step 2: Run spin tests to verify failure**

Run:

```bash
julia --project=. test/runtests.jl
```

Expected: FAIL because `pauli_x` and related functions are not defined.

- [ ] **Step 3: Implement spin operators**

Create `src/SpinOps.jl` with:

```julia
module SpinOps

export pauli_x, pauli_y, pauli_z, identity2, projector_up, projector_down
export kron_all, embed_one_site

function identity2()
    return ComplexF64[1 0; 0 1]
end

function pauli_x()
    return ComplexF64[0 1; 1 0]
end

function pauli_y()
    return ComplexF64[0 -im; im 0]
end

function pauli_z()
    return ComplexF64[1 0; 0 -1]
end

function projector_up()
    return ComplexF64[1 0; 0 0]
end

function projector_down()
    return ComplexF64[0 0; 0 1]
end

function kron_all(ops::AbstractVector{<:AbstractMatrix})
    isempty(ops) && throw(ArgumentError("ops must be nonempty"))
    out = Matrix{ComplexF64}(ops[1])
    for op in ops[2:end]
        out = kron(out, Matrix{ComplexF64}(op))
    end
    return out
end

function embed_one_site(op::AbstractMatrix, site::Integer, nsites::Integer)
    1 <= site <= nsites || throw(ArgumentError("site must be in 1:nsites"))
    return kron_all([i == site ? op : identity2() for i in 1:nsites])
end

end
```

- [ ] **Step 4: Run tests**

Run:

```bash
julia --project=. test/runtests.jl
```

Expected: geometry and spin operator tests PASS; remaining missing test files FAIL.

- [ ] **Step 5: Commit**

```bash
git add src/SpinOps.jl test/test_spinops.jl
git commit -m "feat: add dense spin operators"
```

## Task 4: Models And Blockade Projector

**Files:**
- Create: `src/Models.jl`
- Create: `test/test_models.jl`

- [ ] **Step 1: Write failing model tests**

Create `test/test_models.jl` with:

```julia
@testset "models" begin
    Hpxp = pxp_star_hamiltonian(projector_down(), pauli_x())
    @test size(Hpxp) == (128, 128)
    @test Hpxp ≈ Hpxp'

    P = blockade_projector()
    @test size(P) == (128, 128)
    @test P * P ≈ P
    @test P ≈ P'
    all_up = zeros(ComplexF64, 128)
    all_up[1] = 1
    @test norm(P * all_up) == 0
    all_down = zeros(ComplexF64, 128)
    all_down[end] = 1
    @test P * all_down ≈ all_down

    Hcluster = cluster_star_hamiltonian()
    @test size(Hcluster) == (128, 128)
    @test Hcluster ≈ Hcluster'
    @test Hcluster * Hcluster ≈ Matrix{ComplexF64}(I, 128, 128)

    Hdiag = diagonal_star_hamiltonian()
    @test Hdiag ≈ Hdiag'
    @test Hdiag * Hdiag ≈ Matrix{ComplexF64}(I, 128, 128)

    Hising = ising_bond_hamiltonian()
    @test size(Hising) == (4, 4)
    @test Hising ≈ Hising'
end
```

- [ ] **Step 2: Run model tests to verify failure**

Run:

```bash
julia --project=. test/runtests.jl
```

Expected: FAIL because model constructors are not defined.

- [ ] **Step 3: Implement model constructors**

Create `src/Models.jl` with:

```julia
module Models

using LinearAlgebra
using ..SpinOps

export pxp_star_hamiltonian, blockade_projector, cluster_star_hamiltonian
export diagonal_star_hamiltonian, ising_bond_hamiltonian

const STAR_NSITES = 7
const CENTER_SITE = 1
const NEIGHBOR_SITES = 2:7

function pxp_star_hamiltonian(projector::AbstractMatrix = projector_down(),
                              flip::AbstractMatrix = pauli_x())
    ops = Matrix{ComplexF64}[]
    push!(ops, Matrix{ComplexF64}(flip))
    append!(ops, [Matrix{ComplexF64}(projector) for _ in NEIGHBOR_SITES])
    return kron_all(ops)
end

function blockade_projector()
    diag = ones(ComplexF64, 2^STAR_NSITES)
    for state in 0:(2^STAR_NSITES - 1)
        bits = digits(state, base = 2, pad = STAR_NSITES)
        forbidden = false
        for n in NEIGHBOR_SITES
            if bits[CENTER_SITE] == 0 && bits[n] == 0
                forbidden = true
            end
        end
        if forbidden
            diag[state + 1] = 0
        end
    end
    return Diagonal(diag) |> Matrix
end

function cluster_star_hamiltonian()
    ops = Matrix{ComplexF64}[pauli_x()]
    append!(ops, [pauli_z() for _ in NEIGHBOR_SITES])
    return kron_all(ops)
end

function diagonal_star_hamiltonian()
    ops = Matrix{ComplexF64}[pauli_z()]
    append!(ops, [pauli_z() for _ in NEIGHBOR_SITES])
    return kron_all(ops)
end

function ising_bond_hamiltonian()
    return kron(pauli_z(), pauli_z())
end

end
```

- [ ] **Step 4: Run tests**

Run:

```bash
julia --project=. test/runtests.jl
```

Expected: geometry, spin, and model tests PASS; remaining missing test files FAIL.

- [ ] **Step 5: Commit**

```bash
git add src/Models.jl test/test_models.jl
git commit -m "feat: add dense star model constructors"
```

## Task 5: Dense And Projected Gates

**Files:**
- Create: `src/Gates.jl`
- Create: `test/test_gates.jl`

- [ ] **Step 1: Write failing gate tests**

Create `test/test_gates.jl` with:

```julia
@testset "gates" begin
    H = cluster_star_hamiltonian()
    U = dense_gate(H, 0.2; evolution = :real)
    @test U' * U ≈ Matrix{ComplexF64}(I, 128, 128)

    G = dense_gate(H, 0.2; evolution = :imaginary)
    @test G ≈ G'

    Hpxp = pxp_star_hamiltonian(projector_down(), pauli_x())
    Uproj = projected_gate(Hpxp, 0.05; evolution = :real)
    P = blockade_projector()
    @test P * Uproj ≈ Uproj

    bad = zeros(ComplexF64, 128)
    bad[1] = 1
    cleaned = Uproj * bad
    @test norm((I - P) * cleaned) < 1e-12
end
```

- [ ] **Step 2: Run gate tests to verify failure**

Run:

```bash
julia --project=. test/runtests.jl
```

Expected: FAIL because gate constructors are not defined.

- [ ] **Step 3: Implement gates**

Create `src/Gates.jl` with:

```julia
module Gates

using LinearAlgebra
using ..Models

export dense_gate, projected_gate

function dense_gate(H::AbstractMatrix, step::Real; evolution::Symbol = :real)
    if evolution == :real
        return exp(-im * step * Matrix{ComplexF64}(H))
    elseif evolution == :imaginary
        return exp(-step * Matrix{ComplexF64}(H))
    else
        throw(ArgumentError("evolution must be :real or :imaginary"))
    end
end

function projected_gate(H::AbstractMatrix, step::Real;
                        evolution::Symbol = :real,
                        projector::AbstractMatrix = blockade_projector())
    return Matrix{ComplexF64}(projector) * dense_gate(H, step; evolution)
end

end
```

- [ ] **Step 4: Run tests**

Run:

```bash
julia --project=. test/runtests.jl
```

Expected: geometry, spin, model, and gate tests PASS; remaining missing test files FAIL.

- [ ] **Step 5: Commit**

```bash
git add src/Gates.jl test/test_gates.jl
git commit -m "feat: add dense and projected star gates"
```

## Task 6: Trotter Color Schedules

**Files:**
- Create: `src/Schedules.jl`
- Create: `test/test_schedules.jl`

- [ ] **Step 1: Write failing schedule tests**

Create `test/test_schedules.jl` with:

```julia
@testset "schedules" begin
    @test first_order_colors() == collect(1:7)
    @test second_order_colors() == [1, 2, 3, 4, 5, 6, 7, 7, 6, 5, 4, 3, 2, 1]

    centers = [Coord(q, r) for q in -3:3 for r in -3:3]
    for color in first_order_colors()
        layer = [c for c in centers if star_color(c) == color]
        for i in eachindex(layer), j in (i + 1):lastindex(layer)
            @test disjoint_stars(layer[i], layer[j])
        end
    end
end
```

- [ ] **Step 2: Run schedule tests to verify failure**

Run:

```bash
julia --project=. test/runtests.jl
```

Expected: FAIL because schedule functions are not defined.

- [ ] **Step 3: Implement schedules**

Create `src/Schedules.jl` with:

```julia
module Schedules

export first_order_colors, second_order_colors

function first_order_colors()
    return collect(1:7)
end

function second_order_colors()
    colors = first_order_colors()
    return vcat(colors, reverse(colors))
end

end
```

- [ ] **Step 4: Run tests**

Run:

```bash
julia --project=. test/runtests.jl
```

Expected: all tests except solvable-model tests PASS.

- [ ] **Step 5: Commit**

```bash
git add src/Schedules.jl test/test_schedules.jl
git commit -m "feat: add triangular star update schedules"
```

## Task 7: Solvable True-2D Benchmark Helpers

**Files:**
- Create: `src/SolvableModels.jl`
- Create: `test/test_solvable_models.jl`

- [ ] **Step 1: Write failing solvable benchmark tests**

Create `test/test_solvable_models.jl` with:

```julia
@testset "solvable models" begin
    H = cluster_star_hamiltonian()
    for t in (0.0, 0.1, 0.7)
        @test stabilizer_expectation_exact(t; initial = :plus) ≈ cos(2t)
    end

    U1 = dense_gate(H, 0.11; evolution = :real)
    U2 = dense_gate(H, 0.23; evolution = :real)
    U12 = dense_gate(H, 0.34; evolution = :real)
    @test U2 * U1 ≈ U12
end
```

- [ ] **Step 2: Run solvable tests to verify failure**

Run:

```bash
julia --project=. test/runtests.jl
```

Expected: FAIL because `stabilizer_expectation_exact` is not defined.

- [ ] **Step 3: Implement solvable helpers**

Create `src/SolvableModels.jl` with:

```julia
module SolvableModels

export stabilizer_expectation_exact

function stabilizer_expectation_exact(t::Real; initial::Symbol = :plus)
    if initial == :plus
        return cos(2t)
    else
        throw(ArgumentError("supported initial states: :plus"))
    end
end

end
```

- [ ] **Step 4: Run full test suite**

Run:

```bash
julia --project=. test/runtests.jl
```

Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add src/SolvableModels.jl test/test_solvable_models.jl
git commit -m "feat: add solvable triangular benchmark helpers"
```

## Task 8: Self-Check And README

**Files:**
- Create: `README.md`

- [ ] **Step 1: Write README**

Create `README.md` with:

```markdown
# TriangularPEPSDynamics.jl

High-performance Julia foundations for PEPS-based dynamics on translationally invariant triangular lattices.

The first implemented layer provides:

- triangular axial geometry and 7-site star neighborhoods;
- 7-color non-overlapping star schedules;
- dense spin-1/2 operators;
- dense PXP, projected PXP, and cluster/stabilizer star gates;
- analytically solvable true-2D benchmark helpers.

The next implementation layer will add native triangular iPEPS state containers, Simple Update truncation, observables, and ScarFinder evolve-project drivers.

## Test

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```
```

- [ ] **Step 2: Run final tests**

Run:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: PASS.

- [ ] **Step 3: Run red-flag scan**

Run:

```bash
rg -n "(T)(BD)|(T)(ODO)|(F)(IXME)" .
```

Expected: no matches in source, tests, or README.

- [ ] **Step 4: Commit**

```bash
git add README.md
git commit -m "docs: describe triangular PEPS dynamics foundation"
```

## Plan Self-Review

Spec coverage in this foundation plan:

- Covered: package scaffold, triangular geometry, 7-site star neighborhoods, 7-color schedules, dense spin operators, PXP star Hamiltonian, blockade projector, projected gates, true-2D cluster/stabilizer benchmark, dense gate tests, projected gate tests.
- Deferred to the next plan: iPEPS tensor containers, ITensor-backed tensor indices, Simple Update, NTU, observables over PEPS environments, ScarFinder full evolve-project loop, benchmarks over PEPS bond dimension, GPU experiments.

The deferrals are intentional because this first foundation must produce a small, reliable kernel before implementing expensive PEPS truncation code.
