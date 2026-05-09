# Kagome PESS PXP Dynamics Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build the kagome PESS PXP dynamics pipeline from scratch (geometry → state container → models → gates → schedules → simple update → observables → evolution → scarfinder), with three-layer validation, in the existing `KagomePXPDynamics` Julia project, after removing the triangular predecessor baseline.

**Architecture:** New `Kagome*.jl` modules under `src/`, mirroring the triangular module layout but specialized for the PESS ansatz on kagome with a 9-site enlarged unit cell. The Simple Update kernel applies the 5-site projected PXP gate via cluster contraction through the PESS network and HOSVD-based decomposition back to PESS form. No silent fallbacks — paths not yet supported throw clear `ArgumentError`. Tests follow strict TDD: failing test, run-to-confirm-failure, minimal implementation, run-to-confirm-pass, commit.

**Tech Stack:** Julia 1.12, ITensors.jl, LinearAlgebra, Test. No new dependencies. Root `Project.toml` is unchanged.

**Reference:** [docs/superpowers/specs/2026-05-09-kagome-pess-pxp-dynamics-design.md](../specs/2026-05-09-kagome-pess-pxp-dynamics-design.md) and [Notes/2026-05-09-kagome-pxp-pivot-plan.md](../../../Notes/2026-05-09-kagome-pxp-pivot-plan.md).

---

## Task 1: KagomeGeometry — Coordinates, Sublattices, Neighbors, 9-Site UC, Coloring

**Files:**
- Create: `src/KagomeGeometry.jl`
- Create: `test/test_kagome_geometry.jl`
- Modify: `src/KagomePXPDynamics.jl`
- Modify: `test/runtests.jl`

- [ ] **Step 1: Write failing geometry tests**

Create `test/test_kagome_geometry.jl`:

```julia
@testset "kagome geometry" begin
    @testset "coord and sublattice" begin
        c = KagomeCoord(0, 0, :A)
        @test c.n1 == 0 && c.n2 == 0 && c.sublat === :A
        @test_throws ArgumentError KagomeCoord(0, 0, :X)  # invalid sublattice
    end

    @testset "coordination 4: every site has exactly 4 neighbors" begin
        for sublat in (:A, :B, :C)
            c = KagomeCoord(0, 0, sublat)
            nbrs = [kagome_neighbor(c, d) for d in 1:4]
            @test length(unique(nbrs)) == 4
        end
    end

    @testset "two-triangles-per-site invariant" begin
        # Each kagome site belongs to one up-triangle and one down-triangle.
        for sublat in (:A, :B, :C)
            c = KagomeCoord(1, 2, sublat)
            up = up_triangle_of(c)
            down = down_triangle_of(c)
            @test up !== down
            @test c in triangle_sites(up)
            @test c in triangle_sites(down)
        end
    end

    @testset "5-site star occupies 5 distinct positions" begin
        c = KagomeCoord(0, 0, :A)
        star = kagome_star_sites(c)
        @test length(star) == 5
        @test length(unique(star)) == 5
        @test c == star[1]   # center first
    end

    @testset "9-site UC: 5-star occupies 5 distinct reps" begin
        uc = NineSiteKagomeUC()
        reps = unit_cell_representatives(uc)
        @test length(reps) == 9
        for c in reps
            star = kagome_star_sites(c)
            wrapped = [wrap_kagome_coord(uc, sc) for sc in star]
            @test length(unique(wrapped)) == 5
        end
    end

    @testset "3-coloring: same-color stars are vertex-disjoint" begin
        # Stars centered on the same sublattice (A, B, or C) form one color.
        # Disjointness is checked over a small window of centers.
        centers = [KagomeCoord(n1, n2, s) for n1 in -2:2 for n2 in -2:2 for s in (:A, :B, :C)]
        for color in 1:3
            same = filter(c -> kagome_star_color(c) == color, centers)
            for i in eachindex(same), j in (i+1):lastindex(same)
                @test disjoint_kagome_stars(same[i], same[j])
            end
        end
    end
end
```

- [ ] **Step 2: Run failing tests**

```bash
julia --project=/Users/ren/Codex/PEPs -e 'using Pkg; Pkg.test(test_args=["kagome geometry"])' 2>&1 | tail -20
```

Expected: FAIL on undefined `KagomeCoord`, etc.

- [ ] **Step 3: Implement KagomeGeometry.jl**

Create `src/KagomeGeometry.jl`:

```julia
module KagomeGeometry

export KagomeCoord, KagomeTriangleCoord
export NineSiteKagomeUC, KagomeUnitCell
export kagome_neighbor, kagome_star_sites, kagome_star_color, disjoint_kagome_stars
export up_triangle_of, down_triangle_of, triangle_sites
export unit_cell_representatives, wrap_kagome_coord

const _SUBLATTICES = (:A, :B, :C)

struct KagomeCoord
    n1::Int
    n2::Int
    sublat::Symbol
    function KagomeCoord(n1::Integer, n2::Integer, sublat::Symbol)
        sublat in _SUBLATTICES || throw(ArgumentError("sublat must be :A, :B, or :C"))
        return new(Int(n1), Int(n2), sublat)
    end
end

Base.:(==)(a::KagomeCoord, b::KagomeCoord) = a.n1 == b.n1 && a.n2 == b.n2 && a.sublat === b.sublat
Base.hash(c::KagomeCoord, h::UInt) = hash((c.n1, c.n2, c.sublat), h)

struct KagomeTriangleCoord
    n1::Int
    n2::Int
    orientation::Symbol  # :up or :down
    function KagomeTriangleCoord(n1::Integer, n2::Integer, orientation::Symbol)
        orientation in (:up, :down) || throw(ArgumentError("orientation must be :up or :down"))
        return new(Int(n1), Int(n2), orientation)
    end
end

Base.:(==)(a::KagomeTriangleCoord, b::KagomeTriangleCoord) =
    a.n1 == b.n1 && a.n2 == b.n2 && a.orientation === b.orientation
Base.hash(t::KagomeTriangleCoord, h::UInt) = hash((t.n1, t.n2, t.orientation), h)

# Within an up-triangle anchored at Bravais (n1, n2), the three sites are:
#   A at (n1, n2, :A), B at (n1, n2, :B), C at (n1, n2, :C).
# Within a down-triangle anchored at (n1, n2):
#   A at (n1+1, n2, :A), B at (n1, n2+1, :B), C at (n1, n2, :C)  -- one possible convention.
# The exact convention is chosen here; tests pin the choice.

up_triangle_of(c::KagomeCoord) = KagomeTriangleCoord(c.n1, c.n2, :up)

function down_triangle_of(c::KagomeCoord)
    if c.sublat === :A
        return KagomeTriangleCoord(c.n1 - 1, c.n2, :down)
    elseif c.sublat === :B
        return KagomeTriangleCoord(c.n1, c.n2 - 1, :down)
    else  # :C
        return KagomeTriangleCoord(c.n1, c.n2, :down)
    end
end

function triangle_sites(t::KagomeTriangleCoord)
    if t.orientation === :up
        return (KagomeCoord(t.n1, t.n2, :A),
                KagomeCoord(t.n1, t.n2, :B),
                KagomeCoord(t.n1, t.n2, :C))
    else
        return (KagomeCoord(t.n1 + 1, t.n2, :A),
                KagomeCoord(t.n1, t.n2 + 1, :B),
                KagomeCoord(t.n1, t.n2, :C))
    end
end

# 4 NNs of c: 2 in up-triangle, 2 in down-triangle (excluding c itself).
function kagome_neighbor(c::KagomeCoord, d::Integer)
    1 <= d <= 4 || throw(ArgumentError("direction must be in 1:4"))
    up = collect(triangle_sites(up_triangle_of(c)))
    down = collect(triangle_sites(down_triangle_of(c)))
    up_others = filter(s -> s != c, up)
    down_others = filter(s -> s != c, down)
    nbrs = vcat(up_others, down_others)
    length(nbrs) == 4 || error("kagome_neighbor: expected 4 NNs, got $(length(nbrs))")
    return nbrs[d]
end

# Star: center + 4 NNs.
kagome_star_sites(c::KagomeCoord) = vcat([c], [kagome_neighbor(c, d) for d in 1:4])

# 3-coloring: by sublattice. Stars centered on different sublattices may overlap;
# stars centered on the SAME sublattice within the natural lattice are vertex-disjoint
# only if their Bravais separation is large enough. Our color scheme: color = sublattice index.
kagome_star_color(c::KagomeCoord) = findfirst(==(c.sublat), _SUBLATTICES)

function disjoint_kagome_stars(a::KagomeCoord, b::KagomeCoord)
    return isempty(intersect(Set(kagome_star_sites(a)), Set(kagome_star_sites(b))))
end

abstract type KagomeUnitCell end

struct NineSiteKagomeUC <: KagomeUnitCell end

function unit_cell_representatives(::NineSiteKagomeUC)
    return [KagomeCoord(p, q, s) for p in 0:2 for q in 0:2 for s in _SUBLATTICES]
end

function wrap_kagome_coord(::NineSiteKagomeUC, c::KagomeCoord)
    return KagomeCoord(mod(c.n1, 3), mod(c.n2, 3), c.sublat)
end

end  # module
```

- [ ] **Step 4: Wire into the package**

Modify `src/KagomePXPDynamics.jl` to add:

```julia
include("KagomeGeometry.jl")
using .KagomeGeometry: KagomeCoord, KagomeTriangleCoord, NineSiteKagomeUC,
                       kagome_neighbor, kagome_star_sites, kagome_star_color,
                       disjoint_kagome_stars, up_triangle_of, down_triangle_of,
                       triangle_sites, unit_cell_representatives, wrap_kagome_coord

# Note: unit_cell_representatives is now exported by both States (triangular)
# and KagomeGeometry. Resolve via explicit method dispatch — no ambiguity since
# they take different unit-cell types.
export KagomeCoord, KagomeTriangleCoord, NineSiteKagomeUC
export kagome_neighbor, kagome_star_sites, kagome_star_color, disjoint_kagome_stars
export up_triangle_of, down_triangle_of, triangle_sites, wrap_kagome_coord
```

Modify `test/runtests.jl` to add `include("test_kagome_geometry.jl")` to the testset.

- [ ] **Step 5: Run targeted tests, debug to green**

```bash
julia --project=/Users/ren/Codex/PEPs -e 'using Pkg; Pkg.test(test_args=["kagome geometry"])' 2>&1 | tail -30
```

Expected: PASS. If the "9-site UC: 5-star occupies 5 distinct reps" test fails, the down_triangle_of convention needs a different offset — adjust until the enumeration test passes.

- [ ] **Step 6: Run full suite to confirm no regression**

```bash
julia --project=/Users/ren/Codex/PEPs -e 'using Pkg; Pkg.test()' 2>&1 | tail -10
```

Expected: PASS.

- [ ] **Step 7: Commit**

```bash
git add src/KagomeGeometry.jl src/KagomePXPDynamics.jl test/test_kagome_geometry.jl test/runtests.jl
git commit -m "feat: add kagome geometry with 9-site unit cell and 3-coloring"
```

---

## Task 2: KagomeModels — 32×32 PXP Hamiltonian, Blockade Projector, Stabilizer

**Files:**
- Create: `src/KagomeModels.jl`
- Create: `test/test_kagome_models.jl`
- Modify: `src/KagomePXPDynamics.jl`
- Modify: `test/runtests.jl`

- [ ] **Step 1: Write failing model tests**

Create `test/test_kagome_models.jl`:

```julia
using LinearAlgebra

basis_index_5(bits) = 1 + sum(bits[i] << (5 - i) for i in 1:5)

@testset "kagome models" begin
    @testset "PXP star Hamiltonian shape and hermiticity" begin
        H = pxp_kagome_star_hamiltonian()
        @test size(H) == (32, 32)
        @test H ≈ H'
    end

    @testset "blockade projector idempotence" begin
        P = kagome_blockade_projector()
        @test size(P) == (32, 32)
        @test P * P ≈ P
        @test P ≈ P'
    end

    @testset "blockade projector kills all-up; preserves all-down" begin
        all_up = zeros(ComplexF64, 32); all_up[1] = 1
        all_down = zeros(ComplexF64, 32); all_down[end] = 1
        P = kagome_blockade_projector()
        @test norm(P * all_up) == 0
        @test P * all_down ≈ all_down
    end

    @testset "blockade projector edge enumeration: 6 internal edges" begin
        # Center-NN edges (4) + within-up-triangle NN-NN (1) + within-down-triangle NN-NN (1) = 6.
        # Test specific 5-site bit configurations against the projector.
        # Convention: site 1 = center, sites 2,3 = up-triangle NNs, sites 4,5 = down-triangle NNs.
        P = kagome_blockade_projector()

        # Center up + any single NN up: forbidden (center-NN edge).
        for nn in 2:5
            v = zeros(ComplexF64, 32)
            bits = [0, 1, 1, 1, 1]; bits[nn] = 0  # center=up, this NN=up, others=down
            v[basis_index_5(bits)] = 1
            @test norm(P * v) == 0
        end

        # NN_2 up + NN_3 up (within up-triangle): forbidden.
        v = zeros(ComplexF64, 32)
        v[basis_index_5([1, 0, 0, 1, 1])] = 1  # center=down, NN2=up, NN3=up, NN4=down, NN5=down
        @test norm(P * v) == 0

        # NN_4 up + NN_5 up (within down-triangle): forbidden.
        v = zeros(ComplexF64, 32)
        v[basis_index_5([1, 1, 1, 0, 0])] = 1  # center=down, NN2,3=down, NN4=up, NN5=up
        @test norm(P * v) == 0

        # NN_2 up + NN_4 up (across triangles, NOT adjacent on kagome): allowed.
        v = zeros(ComplexF64, 32)
        v[basis_index_5([1, 0, 1, 0, 1])] = 1
        @test P * v ≈ v
    end

    @testset "kagome cluster star Hamiltonian: hermitian and squares to identity" begin
        K = kagome_cluster_star_hamiltonian()
        @test size(K) == (32, 32)
        @test K ≈ K'
        @test K * K ≈ Matrix{ComplexF64}(I, 32, 32)
    end
end
```

- [ ] **Step 2: Run failing tests**

```bash
julia --project=/Users/ren/Codex/PEPs -e 'using Pkg; Pkg.test(test_args=["kagome models"])' 2>&1 | tail -20
```

Expected: FAIL on undefined functions.

- [ ] **Step 3: Implement KagomeModels.jl**

Create `src/KagomeModels.jl`:

```julia
module KagomeModels

using LinearAlgebra
using ..SpinOps

export pxp_kagome_star_hamiltonian, kagome_blockade_projector, kagome_cluster_star_hamiltonian

const STAR_NSITES_KAGOME = 5
const CENTER_SITE = 1
const NEIGHBOR_SITES = 2:5

# Edge list for the kagome 5-site star blockade:
# (1,2), (1,3), (1,4), (1,5)  -- center-NN edges
# (2,3)                        -- within-up-triangle NN-NN edge
# (4,5)                        -- within-down-triangle NN-NN edge
const KAGOME_STAR_EDGES = (
    (1, 2), (1, 3), (1, 4), (1, 5),
    (2, 3),
    (4, 5),
)

site_bit(state::Integer, site::Integer) = (state >> (STAR_NSITES_KAGOME - site)) & 1

function pxp_kagome_star_hamiltonian(projector::AbstractMatrix = projector_down(),
                                     flip::AbstractMatrix = pauli_x())
    size(projector) == (2, 2) || throw(ArgumentError("projector must be 2x2"))
    size(flip) == (2, 2) || throw(ArgumentError("flip must be 2x2"))

    ops = Matrix{ComplexF64}[]
    push!(ops, Matrix{ComplexF64}(flip))
    append!(ops, [Matrix{ComplexF64}(projector) for _ in NEIGHBOR_SITES])
    return kron_all(ops)
end

function kagome_blockade_projector()
    diag = ones(ComplexF64, 2^STAR_NSITES_KAGOME)
    for state in 0:(2^STAR_NSITES_KAGOME - 1)
        forbidden = false
        for (i, j) in KAGOME_STAR_EDGES
            # excited == 0 (up); forbidden if both endpoints excited.
            if site_bit(state, i) == 0 && site_bit(state, j) == 0
                forbidden = true
                break
            end
        end
        if forbidden
            diag[state + 1] = 0
        end
    end
    return Matrix(Diagonal(diag))
end

function kagome_cluster_star_hamiltonian()
    ops = Matrix{ComplexF64}[pauli_x()]
    append!(ops, [pauli_z() for _ in NEIGHBOR_SITES])
    return kron_all(ops)
end

end  # module
```

- [ ] **Step 4: Wire into package**

Modify `src/KagomePXPDynamics.jl`:

```julia
include("KagomeModels.jl")
using .KagomeModels: pxp_kagome_star_hamiltonian, kagome_blockade_projector,
                     kagome_cluster_star_hamiltonian
export pxp_kagome_star_hamiltonian, kagome_blockade_projector, kagome_cluster_star_hamiltonian
```

Add `include("test_kagome_models.jl")` to `test/runtests.jl`.

- [ ] **Step 5: Run targeted then full suite**

```bash
julia --project=/Users/ren/Codex/PEPs -e 'using Pkg; Pkg.test(test_args=["kagome models"])' 2>&1 | tail -20
julia --project=/Users/ren/Codex/PEPs -e 'using Pkg; Pkg.test()' 2>&1 | tail -10
```

Expected: PASS.

- [ ] **Step 6: Commit**

```bash
git add src/KagomeModels.jl src/KagomePXPDynamics.jl test/test_kagome_models.jl test/runtests.jl
git commit -m "feat: add 5-site kagome PXP Hamiltonian and blockade projector"
```

---

## Task 3: KagomeGates — Dense and Projected 5-Site Gates

**Files:**
- Create: `src/KagomeGates.jl`
- Create: `test/test_kagome_gates.jl`
- Modify: `src/KagomePXPDynamics.jl`
- Modify: `test/runtests.jl`

- [ ] **Step 1: Write failing gate tests**

Create `test/test_kagome_gates.jl`:

```julia
using LinearAlgebra

@testset "kagome gates" begin
    H = kagome_cluster_star_hamiltonian()

    @testset "dense gate: real-time is unitary, imaginary-time is hermitian" begin
        U = dense_kagome_gate(H, 0.2; evolution = :real)
        @test U' * U ≈ Matrix{ComplexF64}(I, 32, 32)
        G = dense_kagome_gate(H, 0.2; evolution = :imaginary)
        @test G ≈ G'
    end

    @testset "projected gate matches P * U" begin
        Hpxp = pxp_kagome_star_hamiltonian()
        U = dense_kagome_gate(Hpxp, 0.05; evolution = :real)
        Uproj = projected_kagome_gate(Hpxp, 0.05; evolution = :real)
        P = kagome_blockade_projector()
        @test Uproj ≈ P * U
        @test P * Uproj ≈ Uproj
    end

    @testset "projected gate kills forbidden output" begin
        Hpxp = pxp_kagome_star_hamiltonian()
        Uproj = projected_kagome_gate(Hpxp, 0.05; evolution = :real)
        bad = zeros(ComplexF64, 32); bad[1] = 1   # all-up
        cleaned = Uproj * bad
        P = kagome_blockade_projector()
        @test norm((I - P) * cleaned) < 1e-12
    end

    @testset "size validation" begin
        @test_throws ArgumentError dense_kagome_gate(zeros(ComplexF64, 4, 4), 0.1)
        @test_throws ArgumentError dense_kagome_gate(H, 0.1; evolution = :bad)
        @test_throws ArgumentError projected_kagome_gate(zeros(ComplexF64, 4, 4), 0.1)
    end
end
```

- [ ] **Step 2: Run failing tests**

```bash
julia --project=/Users/ren/Codex/PEPs -e 'using Pkg; Pkg.test(test_args=["kagome gates"])' 2>&1 | tail -15
```

Expected: FAIL on undefined functions.

- [ ] **Step 3: Implement KagomeGates.jl**

Create `src/KagomeGates.jl`:

```julia
module KagomeGates

using LinearAlgebra
using ..KagomeModels: kagome_blockade_projector

export dense_kagome_gate, projected_kagome_gate

function _validate_kagome_star_hamiltonian(H::AbstractMatrix)
    size(H) == (32, 32) || throw(ArgumentError("H must be 32x32 (5-site kagome star)"))
    return nothing
end

function dense_kagome_gate(H::AbstractMatrix, step::Real; evolution::Symbol = :real)
    _validate_kagome_star_hamiltonian(H)
    if evolution === :real
        return exp(-im * step * Matrix{ComplexF64}(H))
    elseif evolution === :imaginary
        return exp(-step * Matrix{ComplexF64}(H))
    else
        throw(ArgumentError("evolution must be :real or :imaginary"))
    end
end

function projected_kagome_gate(H::AbstractMatrix, step::Real;
                               evolution::Symbol = :real,
                               projector::AbstractMatrix = kagome_blockade_projector())
    _validate_kagome_star_hamiltonian(H)
    size(projector) == size(H) || throw(ArgumentError("projector must match H dimensions"))
    return Matrix{ComplexF64}(projector) * dense_kagome_gate(H, step; evolution)
end

end  # module
```

- [ ] **Step 4: Wire and test**

Modify `src/KagomePXPDynamics.jl`:

```julia
include("KagomeGates.jl")
using .KagomeGates: dense_kagome_gate, projected_kagome_gate
export dense_kagome_gate, projected_kagome_gate
```

Add `include("test_kagome_gates.jl")` to `test/runtests.jl`.

```bash
julia --project=/Users/ren/Codex/PEPs -e 'using Pkg; Pkg.test()' 2>&1 | tail -10
```

Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add src/KagomeGates.jl src/KagomePXPDynamics.jl test/test_kagome_gates.jl test/runtests.jl
git commit -m "feat: add dense and projected 5-site kagome gates"
```

---

## Task 4: KagomeSchedules — 3-Color Trotter Schedule

**Files:**
- Create: `src/KagomeSchedules.jl`
- Create: `test/test_kagome_schedules.jl`
- Modify: `src/KagomePXPDynamics.jl`
- Modify: `test/runtests.jl`

- [ ] **Step 1: Failing tests**

Create `test/test_kagome_schedules.jl`:

```julia
@testset "kagome schedules" begin
    @test first_order_kagome_colors() == [1, 2, 3]
    @test second_order_kagome_colors() == [1, 2, 3, 3, 2, 1]
    @test kagome_schedule_layers(:first) == [(color = c, scale = 1.0) for c in 1:3]
    @test kagome_schedule_layers(:second) == vcat(
        [(color = c, scale = 0.5) for c in 1:3],
        [(color = c, scale = 0.5) for c in 3:-1:1],
    )
    @test_throws ArgumentError kagome_schedule_layers(:bad)

    # canonical center per color matches the sublattice
    @test kagome_color_canonical_center(1).sublat === :A
    @test kagome_color_canonical_center(2).sublat === :B
    @test kagome_color_canonical_center(3).sublat === :C
end
```

- [ ] **Step 2: Implement**

Create `src/KagomeSchedules.jl`:

```julia
module KagomeSchedules

using ..KagomeGeometry: KagomeCoord

export first_order_kagome_colors, second_order_kagome_colors
export kagome_schedule_layers, kagome_color_canonical_center

first_order_kagome_colors() = [1, 2, 3]
second_order_kagome_colors() = [1, 2, 3, 3, 2, 1]

function kagome_schedule_layers(order::Symbol)
    if order === :first
        return [(color = c, scale = 1.0) for c in first_order_kagome_colors()]
    elseif order === :second
        return vcat(
            [(color = c, scale = 0.5) for c in 1:3],
            [(color = c, scale = 0.5) for c in 3:-1:1],
        )
    else
        throw(ArgumentError("order must be :first or :second"))
    end
end

function kagome_color_canonical_center(color::Integer)
    1 <= color <= 3 || throw(ArgumentError("color must be in 1:3"))
    sublat = (:A, :B, :C)[color]
    return KagomeCoord(0, 0, sublat)
end

end  # module
```

- [ ] **Step 3: Wire, test, commit**

Modify `KagomePXPDynamics.jl` to include and re-export. Add to `runtests.jl`. Run full suite.

```bash
git add src/KagomeSchedules.jl src/KagomePXPDynamics.jl test/test_kagome_schedules.jl test/runtests.jl
git commit -m "feat: add kagome 3-color trotter schedule"
```

---

## Task 5: KagomePESS — State Container, Initializers

**Files:**
- Create: `src/KagomePESS.jl`
- Create: `test/test_kagome_pess.jl`
- Modify: `src/KagomePXPDynamics.jl`
- Modify: `test/runtests.jl`

- [ ] **Step 1: Failing tests**

Create `test/test_kagome_pess.jl`:

```julia
using ITensors
using Random

@testset "kagome PESS" begin
    uc = NineSiteKagomeUC()

    @testset "product :down state has all sites in |down>" begin
        state = product_pess(uc, :down; D = 1)
        @test state isa KagomePESS
        for c in unit_cell_representatives(uc)
            ph = phys_index_pess(state, c)
            v_down = ITensor(ComplexF64[0, 1], ph)
            T = site_tensor_pess(state, c)
            # Project T onto v_down; remaining is simplex-bond legs only (all dim 1).
            proj = T * dag(v_down)
            @test abs(scalar(proj) - 1) < 1e-12
        end
    end

    @testset "simplex tensors exist for all 18 triangles in 9-site UC" begin
        state = product_pess(uc, :down; D = 1)
        triangles = collect(keys(state.simplex_tensors))
        @test length(triangles) == 18  # 9 up + 9 down
        for t in triangles
            S = state.simplex_tensors[t]
            @test ndims(S) == 3  # 3 simplex legs, no phys
        end
    end

    @testset "lambdas: each shared bond shared by reference" begin
        state = product_pess(uc, :down; D = 1)
        for c in unit_cell_representatives(uc)
            for which in (:up, :down)
                key_here = (c, which)
                λ_here = simplex_bond_lambda(state, c, which)
                @test all(λ_here .>= 0)
                @test norm(λ_here) ≈ sqrt(length(λ_here))
            end
        end
    end

    @testset "random_pess: shape and reproducibility" begin
        s1 = random_pess(uc, 2; seed = 42)
        s2 = random_pess(uc, 2; seed = 42)
        c0 = KagomeCoord(0, 0, :A)
        T1 = site_tensor_pess(s1, c0)
        T2 = site_tensor_pess(s2, c0)
        @test array(T1) ≈ array(T2)

        s3 = random_pess(uc, 2; seed = 43)
        T3 = site_tensor_pess(s3, c0)
        @test !(array(T1) ≈ array(T3))
    end
end
```

- [ ] **Step 2: Implement**

Create `src/KagomePESS.jl`:

```julia
module KagomePESS

using ITensors
using Random
using ..KagomeGeometry: KagomeCoord, KagomeTriangleCoord, KagomeUnitCell, NineSiteKagomeUC,
                       up_triangle_of, down_triangle_of, triangle_sites,
                       unit_cell_representatives, wrap_kagome_coord

export KagomePESS, product_pess, random_pess
export site_tensor_pess, simplex_tensor_pess, phys_index_pess
export simplex_bond_index, simplex_bond_lambda

# A site-simplex bond is identified by (site_coord, :up | :down).
const KagomeBondKey = Tuple{KagomeCoord, Symbol}

struct KagomePESS{UC<:KagomeUnitCell}
    unitcell::UC
    site_phys_inds::Dict{KagomeCoord, Index}
    site_simplex_inds::Dict{KagomeBondKey, Index}
    simplex_tensors::Dict{KagomeTriangleCoord, ITensor}
    site_tensors::Dict{KagomeCoord, ITensor}
    lambdas::Dict{KagomeBondKey, Vector{Float64}}
end

site_tensor_pess(state::KagomePESS, c::KagomeCoord) =
    state.site_tensors[wrap_kagome_coord(state.unitcell, c)]
simplex_tensor_pess(state::KagomePESS, t::KagomeTriangleCoord) =
    state.simplex_tensors[_wrap_triangle(state.unitcell, t)]
phys_index_pess(state::KagomePESS, c::KagomeCoord) =
    state.site_phys_inds[wrap_kagome_coord(state.unitcell, c)]
simplex_bond_index(state::KagomePESS, c::KagomeCoord, which::Symbol) =
    state.site_simplex_inds[(wrap_kagome_coord(state.unitcell, c), which)]
simplex_bond_lambda(state::KagomePESS, c::KagomeCoord, which::Symbol) =
    state.lambdas[(wrap_kagome_coord(state.unitcell, c), which)]

# Wrap a triangle coord into the 9-site UC. Both n1, n2 mod 3.
function _wrap_triangle(::NineSiteKagomeUC, t::KagomeTriangleCoord)
    return KagomeTriangleCoord(mod(t.n1, 3), mod(t.n2, 3), t.orientation)
end

function _build_pess_indices(uc::NineSiteKagomeUC, D::Integer)
    site_phys_inds = Dict{KagomeCoord, Index}()
    site_simplex_inds = Dict{KagomeBondKey, Index}()
    lambdas = Dict{KagomeBondKey, Vector{Float64}}()

    for c in unit_cell_representatives(uc)
        site_phys_inds[c] = Index(2, "phys,$(c.n1),$(c.n2),$(c.sublat)")
        for which in (:up, :down)
            key = (c, which)
            site_simplex_inds[key] = Index(D, "ssbond,$(c.n1),$(c.n2),$(c.sublat),$(which)")
            lambdas[key] = ones(Float64, D)
        end
    end
    return site_phys_inds, site_simplex_inds, lambdas
end

# Build all simplex tensors for a 9-site UC. There are 9 up- and 9 down-triangles.
function _build_simplex_tensors(uc::NineSiteKagomeUC,
                                site_simplex_inds::Dict{KagomeBondKey, Index};
                                rng::Union{Nothing,AbstractRNG} = nothing,
                                identity_init::Bool = false)
    simplex_tensors = Dict{KagomeTriangleCoord, ITensor}()
    for n1 in 0:2, n2 in 0:2, orientation in (:up, :down)
        t = KagomeTriangleCoord(n1, n2, orientation)
        sites = triangle_sites(t)
        # Each site's "which" leg points to this triangle's orientation.
        legs = Tuple(site_simplex_inds[(wrap_kagome_coord(uc, sites[i]), orientation)] for i in 1:3)
        D = dim(legs[1])
        if identity_init
            data = zeros(ComplexF64, D, D, D)
            for k in 1:D
                data[k, k, k] = 1.0
            end
            S = ITensor(data, legs...)
        else
            S = ITensor(randn(rng, ComplexF64, D, D, D), legs...)
        end
        simplex_tensors[t] = S
    end
    return simplex_tensors
end

# Note: in PESS the same site-simplex Index is shared between the site tensor
# and the simplex tensor. The site tensor for site c has 3 legs: phys + ssbond_up + ssbond_down.
# The simplex tensor for triangle t (e.g., :up) has 3 legs: ssbond_up of each of its 3 sites.
# Since the indices are shared by construction, contraction works automatically.

function product_pess(uc::NineSiteKagomeUC, state_symbol::Symbol; D::Integer = 1)
    D >= 1 || throw(ArgumentError("D must be >= 1"))
    site_phys_inds, site_simplex_inds, lambdas = _build_pess_indices(uc, D)

    local_vec = if state_symbol === :up
        ComplexF64[1, 0]
    elseif state_symbol === :down
        ComplexF64[0, 1]
    else
        throw(ArgumentError("state_symbol must be :up or :down"))
    end

    site_tensors = Dict{KagomeCoord, ITensor}()
    for c in unit_cell_representatives(uc)
        ph = site_phys_inds[c]
        b_up = site_simplex_inds[(c, :up)]
        b_down = site_simplex_inds[(c, :down)]
        T = ITensor(ComplexF64, ph, b_up, b_down)
        for k in 1:2
            T[ph => k, b_up => 1, b_down => 1] = local_vec[k]
        end
        site_tensors[c] = T
    end

    simplex_tensors = _build_simplex_tensors(uc, site_simplex_inds; identity_init = true)
    return KagomePESS(uc, site_phys_inds, site_simplex_inds, simplex_tensors, site_tensors, lambdas)
end

function random_pess(uc::NineSiteKagomeUC, D::Integer; seed::Union{Nothing,Integer} = nothing)
    D >= 1 || throw(ArgumentError("D must be >= 1"))
    site_phys_inds, site_simplex_inds, lambdas = _build_pess_indices(uc, D)
    rng = seed === nothing ? Random.default_rng() : Random.MersenneTwister(seed)
    site_tensors = Dict{KagomeCoord, ITensor}()
    for c in unit_cell_representatives(uc)
        ph = site_phys_inds[c]
        b_up = site_simplex_inds[(c, :up)]
        b_down = site_simplex_inds[(c, :down)]
        data = randn(rng, ComplexF64, 2, D, D)
        site_tensors[c] = ITensor(data, ph, b_up, b_down)
    end
    simplex_tensors = _build_simplex_tensors(uc, site_simplex_inds; rng = rng)
    return KagomePESS(uc, site_phys_inds, site_simplex_inds, simplex_tensors, site_tensors, lambdas)
end

end  # module
```

- [ ] **Step 3: Wire, test, commit**

Wire into `KagomePXPDynamics.jl`. Add `include("test_kagome_pess.jl")` to `runtests.jl`. Run.

```bash
git add src/KagomePESS.jl src/KagomePXPDynamics.jl test/test_kagome_pess.jl test/runtests.jl
git commit -m "feat: add kagome PESS state container with site and simplex tensors"
```

---

## Task 6: KagomeObservables — Local Expectations, Blockade Diagnostics

**Files:**
- Create: `src/KagomeObservables.jl`
- Create: `test/test_kagome_observables.jl`
- Modify: `src/KagomePXPDynamics.jl`
- Modify: `test/runtests.jl`

- [ ] **Step 1: Failing tests**

Create `test/test_kagome_observables.jl`:

```julia
@testset "kagome observables" begin
    uc = NineSiteKagomeUC()

    @testset "local Z on product :down" begin
        state = product_pess(uc, :down; D = 1)
        for c in unit_cell_representatives(uc)
            @test local_expectation_kagome(state, c, pauli_z()) ≈ -1
            @test local_expectation_kagome(state, c, projector_up()) ≈ 0
        end
    end

    @testset "local Z on product :up" begin
        state = product_pess(uc, :up; D = 1)
        for c in unit_cell_representatives(uc)
            @test local_expectation_kagome(state, c, pauli_z()) ≈ 1
        end
    end

    @testset "blockade diagnostics on product states" begin
        down = product_pess(uc, :down; D = 1)
        @test mean_blockade_violation_kagome(down, unit_cell_representatives(uc)) ≈ 0

        up = product_pess(uc, :up; D = 1)
        @test mean_blockade_violation_kagome(up, unit_cell_representatives(uc)) ≈ 1
    end
end
```

- [ ] **Step 2: Implement**

Create `src/KagomeObservables.jl`:

```julia
module KagomeObservables

using LinearAlgebra
using ITensors
using ..KagomeGeometry: KagomeCoord, kagome_neighbor, wrap_kagome_coord
using ..KagomePESS: KagomePESS, site_tensor_pess, phys_index_pess
using ..SpinOps: projector_up

export local_expectation_kagome, mean_blockade_violation_kagome

function local_expectation_kagome(state::KagomePESS, c::KagomeCoord, op::AbstractMatrix)
    size(op) == (2, 2) || throw(ArgumentError("op must be 2x2"))
    rep = wrap_kagome_coord(state.unitcell, c)
    T = state.site_tensors[rep]
    ph = state.site_phys_inds[rep]
    ph_prime = prime(ph)
    op_T = ITensor(Matrix{ComplexF64}(op), ph_prime, ph)
    Tdag = prime(dag(T), ph)
    num = scalar(Tdag * op_T * T)
    denom = scalar(dag(T) * T)
    return num / denom
end

function mean_blockade_violation_kagome(state::KagomePESS, centers)
    vals = Float64[]
    for c in centers, d in 1:4
        pup_c = real(local_expectation_kagome(state, c, projector_up()))
        pup_n = real(local_expectation_kagome(state, kagome_neighbor(c, d), projector_up()))
        push!(vals, clamp(pup_c * pup_n, 0.0, 1.0))
    end
    return isempty(vals) ? 0.0 : sum(vals) / length(vals)
end

end  # module
```

- [ ] **Step 3: Wire, test, commit**

```bash
git add src/KagomeObservables.jl src/KagomePXPDynamics.jl test/test_kagome_observables.jl test/runtests.jl
git commit -m "feat: add kagome PESS observables and blockade diagnostics"
```

---

## Task 7: KagomeSimpleUpdate — Identity And Site-Product Paths

**Goal:** Get the easy paths working first. Identity gate is a no-op; site-product gates apply per-rep single-site operators. Sets up the dispatch skeleton inside `apply_star_gate_simple_update_pess!`.

**Files:**
- Create: `src/KagomeSimpleUpdate.jl`
- Create: `test/test_kagome_simple_update.jl`
- Modify: `src/KagomePXPDynamics.jl`
- Modify: `test/runtests.jl`

- [ ] **Step 1: Failing tests for identity and site-product**

Create `test/test_kagome_simple_update.jl`:

```julia
using LinearAlgebra
using ITensors

@testset "kagome simple update: easy paths" begin
    uc = NineSiteKagomeUC()
    c0 = KagomeCoord(0, 0, :A)

    @testset "identity gate is a no-op" begin
        state = product_pess(uc, :down; D = 1)
        T_before = copy(site_tensor_pess(state, c0))
        I32 = Matrix{ComplexF64}(I, 32, 32)
        diag = apply_star_gate_simple_update_pess!(state, I32, c0)
        T_after = site_tensor_pess(state, c0)
        @test array(T_before) ≈ array(T_after)
        @test diag.discarded_weight == 0
        @test local_expectation_kagome(state, c0, pauli_z()) ≈ -1
    end

    @testset "site-symmetric u^⊗5 on D=1 product :down" begin
        theta = 0.37
        u = cos(theta) * Matrix{ComplexF64}(I, 2, 2) - im * sin(theta) * pauli_x()
        G = u
        for _ in 2:5
            G = kron(G, u)
        end
        state = product_pess(uc, :down; D = 1)
        diag = apply_star_gate_simple_update_pess!(state, G, c0)
        @test diag isa SimpleUpdateDiagnostics
        @test local_expectation_kagome(state, c0, pauli_z()) ≈ -cos(2 * theta) atol = 1e-10
    end

    @testset "general non-product gate at D>1 throws (placeholder until Task 9)" begin
        state = random_pess(uc, 2; seed = 11)
        H = pxp_kagome_star_hamiltonian()
        Uproj = projected_kagome_gate(H, 0.02; evolution = :real)
        @test_throws ArgumentError apply_star_gate_simple_update_pess!(
            state, Uproj, c0; maxdim = 2)
    end
end
```

- [ ] **Step 2: Implement skeleton**

Create `src/KagomeSimpleUpdate.jl`:

```julia
module KagomeSimpleUpdate

using LinearAlgebra
using ITensors
using ..KagomeGeometry: KagomeCoord, kagome_star_sites, wrap_kagome_coord
using ..KagomePESS: KagomePESS, site_tensor_pess, phys_index_pess

# Reuse the diagnostics struct from triangular SimpleUpdate.
using ..SimpleUpdate: SimpleUpdateDiagnostics

export apply_star_gate_simple_update_pess!

const _STAR_NSITES_KAGOME = 5

function apply_star_gate_simple_update_pess!(state::KagomePESS,
                                             gate::AbstractMatrix,
                                             center::KagomeCoord;
                                             cutoff::Real = 1e-12,
                                             maxdim::Union{Nothing,Integer} = nothing)
    size(gate) == (32, 32) || throw(ArgumentError("gate must be 32x32 (5-site kagome star)"))

    G = Matrix{ComplexF64}(gate)
    Iref = Matrix{ComplexF64}(I, 32, 32)
    star = kagome_star_sites(center)

    # Path 1: identity gate.
    if norm(G - Iref) <= 1e-12 * max(norm(Iref), 1.0)
        affected = _affected_kagome_bonds(state, center)
        dims = [dim(state.site_simplex_inds[(c, w)]) for (c, w) in affected]
        return SimpleUpdateDiagnostics(0.0, _affected_to_triangular_keys(affected), dims)
    end

    # Path 2: site-product factorization.
    factors = _try_factorize_product_gate_n(G, _STAR_NSITES_KAGOME)
    if factors !== nothing
        # Group by rep, require translational consistency, average, apply per rep.
        rep_factors = Dict{KagomeCoord, Vector{Matrix{ComplexF64}}}()
        for (i, sc) in enumerate(star)
            rep = wrap_kagome_coord(state.unitcell, sc)
            push!(get!(rep_factors, rep, Matrix{ComplexF64}[]), factors[i])
        end
        for (rep, fs) in rep_factors
            ref = fs[1]
            ref_norm2 = sum(abs2, ref)
            ref_norm2 < 1e-28 && error("vanishing factor at rep $rep")
            for f in fs
                c_scalar = sum(conj.(ref) .* f) / ref_norm2
                if norm(f - c_scalar * ref) > 1e-8 * max(norm(ref), 1.0)
                    error("kagome simple update: gate factors at star positions sharing rep $rep " *
                          "are not scalar multiples of a common single-site operator")
                end
            end
            u = ref / sqrt(ref_norm2 / 2)
            T = state.site_tensors[rep]
            ph = state.site_phys_inds[rep]
            u_T = ITensor(u, prime(ph), ph)
            state.site_tensors[rep] = noprime(u_T * T)
        end
        affected = _affected_kagome_bonds(state, center)
        dims = [dim(state.site_simplex_inds[(c, w)]) for (c, w) in affected]
        return SimpleUpdateDiagnostics(0.0, _affected_to_triangular_keys(affected), dims)
    end

    # Path 3: general dense 5-site gate. Implemented in Tasks 8-9.
    throw(ArgumentError(
        "general dense 5-site Simple Update for kagome PESS is not implemented yet; " *
        "only identity gates and site-product gates are currently supported"))
end

function _affected_kagome_bonds(state::KagomePESS, center::KagomeCoord)
    star = kagome_star_sites(center)
    bonds = Tuple{KagomeCoord, Symbol}[]
    for (i, sc) in enumerate(star)
        rep = wrap_kagome_coord(state.unitcell, sc)
        # Center: both up and down simplex bonds are within the cluster.
        # Each NN: one of (up, down) is within the cluster (the shared one with center).
        if i == 1
            push!(bonds, (rep, :up))
            push!(bonds, (rep, :down))
        else
            # Determine whether this NN shares the up-triangle or the down-triangle with the center.
            # Sites in the same up-triangle as `center`: have up_triangle_of(nn) == up_triangle_of(center).
            # Implementation will pick the right one.
            nn_up = KagomePXPDynamics.KagomeGeometry.up_triangle_of(sc)
            center_up = KagomePXPDynamics.KagomeGeometry.up_triangle_of(center)
            push!(bonds, (rep, nn_up == center_up ? :up : :down))
        end
    end
    return bonds
end

# The triangular SimpleUpdateDiagnostics expects affected_bonds as Vector{Tuple{Coord,Int}}.
# For kagome we use a different bond key shape; convert to a placeholder int representation
# for now, or extend SimpleUpdateDiagnostics. Easiest: extend the struct OR define a kagome
# variant. Decision: keep the existing diagnostics struct, encode kagome bonds via a stable
# string-int hash. Cleaner long-term: introduce a generic Diagnostics struct.
# For first cut, define a kagome-local Diagnostics:
function _affected_to_triangular_keys(bonds)
    # Placeholder: shoehorn (KagomeCoord, Symbol) into (Coord, Int) by indexing.
    # The actual SimpleUpdateDiagnostics fields stay typed against triangular Coord;
    # for kagome we use a separate KagomeSimpleUpdateDiagnostics defined below.
    return [(KagomePXPDynamics.Geometry.Coord(b[1].n1, b[1].n2), b[2] === :up ? 1 : 2)
            for b in bonds]
end

# Generalized site-product factorization for n sites (currently n=5 for kagome, n=7 for triangular).
function _try_factorize_product_gate_n(G::AbstractMatrix, n::Int; tol::Real = 1e-10)
    @assert size(G, 1) == size(G, 2) == 2^n
    peeled = Matrix{ComplexF64}[]
    rest = Matrix{ComplexF64}(G)
    for _ in 1:(n - 1)
        m = size(rest, 1) ÷ 2
        rest_4d = reshape(rest, 2, m, 2, m)
        rest_perm = permutedims(rest_4d, (1, 3, 2, 4))
        rest_mat = reshape(rest_perm, 4, m * m)
        F = svd(rest_mat)
        if length(F.S) > 1 && F.S[2] / max(F.S[1], eps()) > tol
            return nothing
        end
        s = F.S[1]
        u_right = reshape(F.U[:, 1] * sqrt(s), 2, 2)
        rest = reshape(F.Vt[1, :] * sqrt(s), m, m)
        push!(peeled, u_right)
    end
    push!(peeled, rest)
    factors = reverse(peeled)
    reconstructed = factors[1]
    for factor in factors[2:end]
        reconstructed = kron(reconstructed, factor)
    end
    if norm(reconstructed - G) > 1e-8 * max(norm(G), 1.0)
        return nothing
    end
    return factors
end

end  # module
```

**Important note for the executing agent:** The diagnostics-key mismatch above is a placeholder. The clean fix is to introduce a `KagomeSimpleUpdateDiagnostics` struct in `KagomeSimpleUpdate.jl` mirroring `SimpleUpdateDiagnostics` but with `Vector{KagomeBondKey}` for `affected_bonds`. Replace `SimpleUpdateDiagnostics` returns above with the kagome variant, and update tests accordingly. This is a 5-line change but critical to keeping the bond bookkeeping correct.

- [ ] **Step 3: Run targeted tests**

```bash
julia --project=/Users/ren/Codex/PEPs -e 'using Pkg; Pkg.test(test_args=["kagome simple update: easy paths"])' 2>&1 | tail -25
```

Expected: PASS for identity and site-product, FAIL (correctly throwing) for general gate.

- [ ] **Step 4: Wire, run full suite, commit**

```bash
git add src/KagomeSimpleUpdate.jl src/KagomePXPDynamics.jl test/test_kagome_simple_update.jl test/runtests.jl
git commit -m "feat: kagome PESS simple update for identity and site-product gates"
```

---

## Task 8: KagomeSimpleUpdate — D=1 Dense Product Oracle

**Goal:** Same correctness pattern as the triangular D=1 oracle, adapted for kagome's 5-site cluster: build the 32-dim cluster vector from per-site vectors, apply G, decompose result back to per-site vectors via dominant-eigenvector projection.

**Files:**
- Modify: `src/KagomeSimpleUpdate.jl`
- Modify: `test/test_kagome_simple_update.jl`

- [ ] **Step 1: Failing test — projected kagome PXP at D=1 from all-down**

Append to `test/test_kagome_simple_update.jl`:

```julia
@testset "projected kagome PXP gate updates D=1 product exactly" begin
    t = 0.17
    H = pxp_kagome_star_hamiltonian()
    Uproj = projected_kagome_gate(H, t; evolution = :real)
    state = product_pess(NineSiteKagomeUC(), :down; D = 1)
    c0 = KagomeCoord(0, 0, :A)

    diag = apply_star_gate_simple_update_pess!(state, Uproj, c0; maxdim = 1)

    @test diag isa SimpleUpdateDiagnostics  # or KagomeSimpleUpdateDiagnostics if introduced
    @test diag.discarded_weight ≈ 0 atol = 1e-12
    @test real(local_expectation_kagome(state, c0, pauli_z())) ≈ -cos(2t) atol = 1e-10
    # NN sites still in :down (PXP from all-down only flips center).
    for d in 1:4
        @test real(local_expectation_kagome(state, kagome_neighbor(c0, d), pauli_z())) ≈ -1 atol = 1e-10
    end
end
```

- [ ] **Step 2: Implement the D=1 oracle path in KagomeSimpleUpdate.jl**

Add to `src/KagomeSimpleUpdate.jl`, before the throwing branch in `apply_star_gate_simple_update_pess!`:

```julia
    # Path 3a: D=1 dense product oracle (when all simplex bonds have dim 1).
    if _is_d1_pess(state)
        local_vectors = [_d1_pess_site_vector(state, wrap_kagome_coord(state.unitcell, sc)) for sc in star]
        psi = local_vectors[1]
        for v in local_vectors[2:end]
            psi = kron(psi, v)
        end
        phi = G * psi

        # Try product-state factorization first; otherwise fall back to dominant eigenvectors.
        product_factors = _try_factorize_product_state_n(phi, _STAR_NSITES_KAGOME)
        if product_factors !== nothing
            for (rep, positions) in _kagome_star_positions_by_rep(state, star)
                v = _common_factor_for_positions_n(product_factors, positions)
                _set_d1_pess_site_vector!(state, rep, v)
            end
        else
            tensor_phi = reshape(phi, ntuple(_ -> 2, _STAR_NSITES_KAGOME)...)
            for (rep, positions) in _kagome_star_positions_by_rep(state, star)
                rho = _one_site_density_from_kagome_star_tensor(tensor_phi, positions[1])
                vals, vecs = eigen(Hermitian(rho))
                v = ComplexF64.(vecs[:, argmax(vals)])
                _set_d1_pess_site_vector!(state, rep, v)
            end
        end

        affected = _affected_kagome_bonds(state, center)
        dims = [dim(state.site_simplex_inds[(c, w)]) for (c, w) in affected]
        return SimpleUpdateDiagnostics(0.0, _affected_to_triangular_keys(affected), dims)
    end
```

Add the helpers:

```julia
function _is_d1_pess(state::KagomePESS)
    for key in keys(state.site_simplex_inds)
        dim(state.site_simplex_inds[key]) == 1 || return false
    end
    return true
end

function _d1_pess_site_vector(state::KagomePESS, rep::KagomeCoord)
    T = state.site_tensors[rep]
    ph = state.site_phys_inds[rep]
    bup = state.site_simplex_inds[(rep, :up)]
    bdn = state.site_simplex_inds[(rep, :down)]
    return ComplexF64[T[ph => 1, bup => 1, bdn => 1], T[ph => 2, bup => 1, bdn => 1]]
end

function _set_d1_pess_site_vector!(state::KagomePESS, rep::KagomeCoord, v::Vector{ComplexF64})
    nrm = norm(v)
    nrm == 0 && throw(ArgumentError("cannot set zero local vector"))
    v = v / nrm
    ph = state.site_phys_inds[rep]
    bup = state.site_simplex_inds[(rep, :up)]
    bdn = state.site_simplex_inds[(rep, :down)]
    T = ITensor(ComplexF64, ph, bup, bdn)
    for k in 1:2
        T[ph => k, bup => 1, bdn => 1] = v[k]
    end
    state.site_tensors[rep] = T
end

function _kagome_star_positions_by_rep(state::KagomePESS, star)
    grouped = Dict{KagomeCoord, Vector{Int}}()
    for (i, sc) in enumerate(star)
        rep = wrap_kagome_coord(state.unitcell, sc)
        push!(get!(grouped, rep, Int[]), i)
    end
    return grouped
end

# These two are isomorphic to the existing triangular helpers; copy and adapt.
function _try_factorize_product_state_n(vec::AbstractVector, n::Int; tol::Real = 1e-10)
    length(vec) == 2^n || throw(ArgumentError("vec length must be 2^n"))
    factors = Vector{ComplexF64}[]
    rest = ComplexF64.(vec)
    for _ in 1:(n - 1)
        block = length(rest) ÷ 2
        rest_mat = Matrix{ComplexF64}(undef, 2, block)
        rest_mat[1, :] .= view(rest, 1:block)
        rest_mat[2, :] .= view(rest, (block + 1):(2 * block))
        F = svd(rest_mat)
        if length(F.S) > 1 && F.S[2] / max(F.S[1], eps()) > tol
            return nothing
        end
        s = F.S[1]
        push!(factors, ComplexF64.(F.U[:, 1] * sqrt(s)))
        rest = ComplexF64.(F.Vt[1, :] * sqrt(s))
    end
    push!(factors, rest)
    reconstructed = factors[1]
    for factor in factors[2:end]
        reconstructed = kron(reconstructed, factor)
    end
    if norm(reconstructed - vec) > 1e-8 * max(norm(vec), 1.0)
        return nothing
    end
    return factors
end

function _common_factor_for_positions_n(factors, positions)
    ref = factors[positions[1]]
    ref_norm2 = sum(abs2, ref)
    ref_norm2 < 1e-28 && throw(ArgumentError("cannot use a zero product factor"))
    aligned = zeros(ComplexF64, 2)
    for pos in positions
        f = factors[pos]
        c = sum(conj.(ref) .* f) / ref_norm2
        if norm(f - c * ref) > 1e-8 * max(norm(ref), 1.0)
            return ref
        end
        aligned .+= f / c
    end
    return aligned / length(positions)
end

function _one_site_density_from_kagome_star_tensor(tensor_phi, position::Int)
    others = [i for i in 1:_STAR_NSITES_KAGOME if i != position]
    perm = (position, others...)
    psi = reshape(permutedims(tensor_phi, perm), 2, :)
    rho = psi * psi'
    tr = real(sum(diag(rho)))
    tr == 0 && return Matrix{ComplexF64}(I, 2, 2) / 2
    return rho / tr
end
```

- [ ] **Step 3: Run targeted, then full suite, commit**

```bash
julia --project=/Users/ren/Codex/PEPs -e 'using Pkg; Pkg.test(test_args=["projected kagome PXP gate updates D=1"])'
julia --project=/Users/ren/Codex/PEPs -e 'using Pkg; Pkg.test()' 2>&1 | tail -10
```

```bash
git add src/KagomeSimpleUpdate.jl test/test_kagome_simple_update.jl
git commit -m "feat: kagome simple update D=1 dense product oracle"
```

---

## Task 9: KagomeSimpleUpdate — General 5-Site Cluster + HOSVD Decomposition

**Goal:** The hard task. Replace the explicit-throw branch with the cluster-and-HOSVD-split path. Layered as five sub-steps with their own tests, mirroring the triangular plan's structure but adapted for PESS.

This task is intentionally designed to land in **multiple commits**, one per sub-step.

### Task 9.1: Failing kernel test (random unitary at D=2, no truncation)

**Files:**
- Modify: `test/test_kagome_simple_update.jl`
- Create: `test/util_kagome_finite_ed.jl` (partial — just `cluster_vector_from_pess`)
- Modify: `test/runtests.jl`

- [ ] **Step 1: Add the cluster-vector helper**

Create `test/util_kagome_finite_ed.jl`:

```julia
using ITensors
using LinearAlgebra
using KagomePXPDynamics

"""
    cluster_vector_from_pess(state, center) -> Vector{ComplexF64}

Contract the 5 site tensors and 2 simplex tensors of the star at `center`
into a single 32-dim vector. External simplex bonds (those of NNs that
exit the cluster) are summed against all-ones environments so the result
is a consistent recipe for before/after gate comparisons.
"""
function cluster_vector_from_pess(state::KagomePESS, center::KagomeCoord)
    star = kagome_star_sites(center)
    reps = [wrap_kagome_coord(state.unitcell, sc) for sc in star]

    # Get the two simplex tensors of triangles containing the center.
    up_t = up_triangle_of(center)
    down_t = down_triangle_of(center)
    # Wrap to the 9-site UC.
    Sup = state.simplex_tensors[KagomeTriangleCoord(mod(up_t.n1, 3), mod(up_t.n2, 3), :up)]
    Sdn = state.simplex_tensors[KagomeTriangleCoord(mod(down_t.n1, 3), mod(down_t.n2, 3), :down)]

    # Rename phys indices per position so they're distinct.
    fresh_phys = [Index(2, "phys_pos_$(i)") for i in 1:5]
    site_tensors = [state.site_tensors[reps[i]] for i in 1:5]
    renamed = [replaceind(site_tensors[i], state.site_phys_inds[reps[i]], fresh_phys[i]) for i in 1:5]

    # Identify external simplex bonds: NNs (positions 2..5) each have one external bond
    # (to a triangle outside the star). Center has zero external bonds.
    cluster = Sup * Sdn
    for i in 1:5
        cluster = cluster * renamed[i]
    end

    # Trace remaining external bonds (those that didn't get contracted in the cluster) with all-ones.
    for idx in inds(cluster)
        if !(idx in fresh_phys)
            cluster = cluster * ITensor(ones(ComplexF64, dim(idx)), idx)
        end
    end

    return ComplexF64.(reshape(array(cluster, fresh_phys...), 32))
end
```

Add `include("util_kagome_finite_ed.jl")` to `test/runtests.jl` (before the testset block).

- [ ] **Step 2: Add the failing kernel test**

Append to `test/test_kagome_simple_update.jl`:

```julia
@testset "kagome general star kernel: random unitary with maxdim large enough is exact" begin
    using Random
    rng = MersenneTwister(2026_05_09_2)
    state = random_pess(NineSiteKagomeUC(), 2; seed = 11)
    psi_before = cluster_vector_from_pess(state, KagomeCoord(0, 0, :A))

    A = randn(rng, ComplexF64, 32, 32)
    Q, _ = qr(A); G = Matrix{ComplexF64}(Q)

    diag = apply_star_gate_simple_update_pess!(state, G, KagomeCoord(0, 0, :A);
                                               maxdim = 64, cutoff = 0.0)

    psi_after = cluster_vector_from_pess(state, KagomeCoord(0, 0, :A))
    expected = G * psi_before
    overlap = abs(dot(psi_after, expected))
    @test overlap ≈ norm(psi_after) * norm(expected) atol = 1e-8
    @test diag.discarded_weight ≈ 0 atol = 1e-10
end
```

- [ ] **Step 3: Run, confirm failure, commit**

```bash
julia --project=/Users/ren/Codex/PEPs -e 'using Pkg; Pkg.test(test_args=["kagome general star kernel"])' 2>&1 | tail -25
```

Expected: FAIL (current code throws on general path).

```bash
git add test/util_kagome_finite_ed.jl test/test_kagome_simple_update.jl test/runtests.jl
git commit -m "test: failing kernel test for general kagome 5-site simple update"
```

### Task 9.2: λ Absorption Helper

- [ ] **Step 1: Add helper + test**

In `src/KagomeSimpleUpdate.jl`:

```julia
"""
    _absorb_lambda_into_kagome_star(state, center) -> NamedTuple

Returns absorbed copies of the 5 star site tensors and the 2 simplex tensors
spanning the star, with sqrt(lambda) absorbed on every simplex bond leg of
the site tensors. Pure: does not mutate state.
"""
function _absorb_lambda_into_kagome_star(state::KagomePESS, center::KagomeCoord)
    star = kagome_star_sites(center)
    reps = [wrap_kagome_coord(state.unitcell, sc) for sc in star]

    absorbed_sites = Vector{ITensor}(undef, 5)
    for i in 1:5
        T = state.site_tensors[reps[i]]
        for which in (:up, :down)
            bind = state.site_simplex_inds[(reps[i], which)]
            λ = state.lambdas[(reps[i], which)]
            T = T * ITensor(sqrt.(λ), bind)
        end
        absorbed_sites[i] = T
    end

    up_t = up_triangle_of(center)
    down_t = down_triangle_of(center)
    Sup = state.simplex_tensors[KagomeTriangleCoord(mod(up_t.n1, 3), mod(up_t.n2, 3), :up)]
    Sdn = state.simplex_tensors[KagomeTriangleCoord(mod(down_t.n1, 3), mod(down_t.n2, 3), :down)]

    return (sites = absorbed_sites, Sup = Sup, Sdn = Sdn, reps = reps)
end
```

Add a kernel test for round-trip with λ = 1 (a no-op) in `test/test_kagome_simple_update.jl`. Run, commit.

```bash
git add src/KagomeSimpleUpdate.jl test/test_kagome_simple_update.jl
git commit -m "feat: kagome lambda absorption helper for star simple update"
```

### Task 9.3: Cluster Build With Gate

- [ ] **Step 1: Add helper + test**

```julia
function _build_kagome_cluster_with_gate(absorbed::NamedTuple, G::Matrix{ComplexF64})
    size(G) == (32, 32) || throw(ArgumentError("gate must be 32x32"))

    # Build per-position out-physical indices; remap absorbed sites' phys indices
    # to fresh "in" indices that the gate will consume.
    in_phys = Vector{Index}(undef, 5)
    out_phys = Vector{Index}(undef, 5)
    for i in 1:5
        ph_old = nothing
        for idx in inds(absorbed.sites[i])
            if hastags(idx, "phys")
                ph_old = idx; break
            end
        end
        in_phys[i] = Index(2, "in_phys_pos_$(i)")
        out_phys[i] = Index(2, "out_phys_pos_$(i)")
        absorbed.sites[i] = replaceind(absorbed.sites[i], ph_old, in_phys[i])
    end

    G_data = reshape(G, ntuple(_ -> 2, 10)...)
    G_tensor = ITensor(G_data, out_phys..., in_phys...)

    cluster = absorbed.Sup * absorbed.Sdn
    for i in 1:5
        cluster = cluster * absorbed.sites[i]
    end
    cluster = cluster * G_tensor

    return cluster, out_phys
end
```

Add a test that identity gate keeps the cluster unchanged (up to phys index relabeling). Commit.

```bash
git add src/KagomeSimpleUpdate.jl test/test_kagome_simple_update.jl
git commit -m "feat: kagome cluster contraction with gate helper"
```

### Task 9.4: HOSVD Decomposition Back To PESS

- [ ] **Step 1: Add helper + test**

The decomposition is the most algorithmically dense step. Approach:

1. **Refit each site tensor**: for each of the 5 site positions, the new site tensor is determined by SVD-cutting the cluster across `(out_phys[i], external_bonds_of_position_i) | rest`. The left side becomes the new site tensor; the right side accumulates into the rest.
2. **Refit simplex tensors via HOSVD**: after extracting the 5 site tensors, the residual core has 6 internal "simplex-leg" indices (3 per simplex) plus 2 simplex tensors' worth of structure. Apply HOSVD to each simplex's 3-leg core to obtain the new simplex tensors.

This is non-trivial to write correctly on first attempt. **The executing agent should expect this sub-step to require 2-3 iterations**, with intermediate failing tests guiding the fix. Rough skeleton:

```julia
function _hosvd_decompose_kagome_cluster(cluster::ITensor,
                                          out_phys::Vector{Index},
                                          absorbed::NamedTuple;
                                          cutoff::Real, maxdim::Int)
    # Step A: peel each site position via SVD across (out_phys[i], external bonds of pos i) | rest.
    # Step B: residual cluster has the 6 internal simplex-leg indices.
    # Step C: HOSVD each simplex tensor: SVD-truncate each of its 3 legs.
    # Step D: return (new_site_tensors, new_simplex_tensors, new_lambdas, total_discarded).

    # Implementation involves careful index tracking. Recommend writing on
    # paper first with explicit index labels for a small example (D=2).
    error("Task 9.4 implementation pending")
end
```

Add a test that with `maxdim = 64, cutoff = 0` and identity gate input, the decomposed PESS reconstructs the original (within tolerance). Commit when passing.

```bash
git add src/KagomeSimpleUpdate.jl test/test_kagome_simple_update.jl
git commit -m "feat: kagome HOSVD decomposition back to PESS form"
```

### Task 9.5: Writeback And Wire Together

- [ ] **Step 1: Implement writeback**

```julia
function _extract_and_writeback_kagome!(state::KagomePESS,
                                         center::KagomeCoord,
                                         new_sites::Vector{ITensor},
                                         new_simplices::Dict{Symbol, ITensor},
                                         new_lambdas::Dict{KagomeBondKey, Vector{Float64}},
                                         absorbed::NamedTuple)
    # Step A: divide out sqrt(lambda) on external bonds of each NN site.
    # Step B: relabel internal bond indices to fresh canonical Index objects.
    # Step C: write back state.site_tensors, state.simplex_tensors, state.bond_inds, state.lambdas.
    # Step D: enforce lambda normalization: nonneg descending, norm(lambda) == sqrt(length(lambda)).
    # Step E: shared-bond invariant: opposite ends of each updated bond share lambda by reference.
    error("Task 9.5 implementation pending")
end
```

- [ ] **Step 2: Wire the four helpers together in `apply_star_gate_simple_update_pess!`**

Replace the `throw(ArgumentError(...))` line in Path 3 with calls to the four helpers in sequence. Update the diagnostics return.

- [ ] **Step 3: Run Task 9.1's failing kernel test → should now PASS**

```bash
julia --project=/Users/ren/Codex/PEPs -e 'using Pkg; Pkg.test(test_args=["kagome general star kernel"])'
```

- [ ] **Step 4: Full suite**

```bash
julia --project=/Users/ren/Codex/PEPs -e 'using Pkg; Pkg.test()' 2>&1 | tail -10
```

Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add src/KagomeSimpleUpdate.jl test/test_kagome_simple_update.jl
git commit -m "feat: kagome general 5-site simple update with HOSVD decomposition"
```

---

## Task 10: Layer-1 Coverage — Truncation, Hermiticity, Site-Product Equivalence

**Files:**
- Modify: `test/test_kagome_simple_update.jl`

- [ ] **Step 1: Add coverage tests** (mirroring the triangular Task 8 tests, adapted to 32-dim cluster):

```julia
@testset "kagome general star kernel: truncation activates" begin
    using Random
    rng = MersenneTwister(2026_05_09_3)
    state = random_pess(NineSiteKagomeUC(), 2; seed = 31)
    A = randn(rng, ComplexF64, 32, 32); Q, _ = qr(A); G = Matrix{ComplexF64}(Q)

    diag = apply_star_gate_simple_update_pess!(state, G, KagomeCoord(0, 0, :A);
                                               maxdim = 2, cutoff = 0.0)
    @test diag.discarded_weight > 0
    @test all(d -> d <= 2, diag.output_bond_dims)
end

@testset "kagome general star kernel: imaginary-time gate keeps tensors finite" begin
    state = random_pess(NineSiteKagomeUC(), 2; seed = 37)
    H = pxp_kagome_star_hamiltonian()
    G = exp(-0.01 * Matrix{ComplexF64}(H))
    @test G ≈ G' atol = 1e-12
    apply_star_gate_simple_update_pess!(state, G, KagomeCoord(0, 0, :A);
                                        maxdim = 4, cutoff = 1e-12)
    for c in unit_cell_representatives(NineSiteKagomeUC())
        @test all(isfinite, array(site_tensor_pess(state, c)))
    end
end
```

- [ ] **Step 2: Run, fix, commit**

```bash
git add test/test_kagome_simple_update.jl
git commit -m "test: layer-1 coverage for kagome truncation and hermiticity"
```

---

## Task 11: KagomeEvolution — Step Driver, Run Loop

**Files:**
- Create: `src/KagomeEvolution.jl`
- Create: `test/test_kagome_evolution.jl`
- Modify: `src/KagomePXPDynamics.jl`
- Modify: `test/runtests.jl`

- [ ] **Step 1: Implement following the triangular Evolution.jl pattern**

Adapt the triangular `Evolution.jl` to call `apply_star_gate_simple_update_pess!` over the kagome 3-color schedule. Add `projected_pxp_step_kagome!`, `imaginary_projected_pxp_step_kagome!`, `run_projected_pxp_kagome!`. Implement `KagomeProjectedPXPStepDiagnostics` mirroring the triangular variant.

- [ ] **Step 2: Add tests**: identity-evolution preserves product state; projected PXP smoke test on D=1 product :down 3-site UC; step helpers report correct schedule lengths (3 for first-order, 6 for second-order).

- [ ] **Step 3: Wire, run full suite, commit**

```bash
git add src/KagomeEvolution.jl src/KagomePXPDynamics.jl test/test_kagome_evolution.jl test/runtests.jl
git commit -m "feat: kagome projected PXP evolution driver"
```

---

## Task 12: Layer-2 Torus Integration Test (Indicator, Not Gate)

**Status note for the executing agent:** This task is **demoted to dynamics-fidelity indicator** per the spec's revised acceptance criteria. The test is added and run, the result is recorded, but **a failing Layer-2 test does NOT block this PR**. The result interprets whether the SU prototype is acceptable for first scientific scarfinder use or whether NTU should be prioritized. See spec section "Acceptance Criteria → Indicator" for the outcome interpretation matrix.

The test code is identical to what a gate-test would be; only the orchestrator's reaction to a failure differs. Execute Steps 1-3 normally; in Step 4 record the outcome rather than treat it as a hard failure.

**Files:**
- Modify: `test/util_kagome_finite_ed.jl`
- Modify: `test/test_kagome_evolution.jl`

- [ ] **Step 1: Add torus helpers**

Append to `test/util_kagome_finite_ed.jl`:

```julia
"""
    build_kagome_torus_pxp_hamiltonian(L1, L2) -> (H::Matrix{ComplexF64}, sites::Vector{KagomeCoord})

Build the L1 x L2 supercell kagome torus PXP Hamiltonian (3 * L1 * L2 sites).
"""
function build_kagome_torus_pxp_hamiltonian(L1::Integer, L2::Integer)
    # ... iterate over Bravais sites, construct PXP terms with periodic wrap.
    # 12-site (L1=L2=2) is the recommended starting size; Hilbert space 4096.
    error("torus PXP builder pending")
end

function kagome_torus_local_z_per_sublattice(vec::Vector{ComplexF64}, L1, L2)
    # Compute mean <Z_i> for each sublattice A, B, C.
    error("kagome torus Z helper pending")
end

function kagome_torus_initial_all_down(L1, L2)
    n = 3 * L1 * L2
    v = zeros(ComplexF64, 2^n); v[end] = 1; return v
end
```

- [ ] **Step 2: Add the integration test**

```julia
@testset "kagome projected PXP iPEPS matches small torus ED at short times" begin
    L1, L2 = 2, 2  # 12-site torus
    H_torus, _ = build_kagome_torus_pxp_hamiltonian(L1, L2)
    psi_torus = kagome_torus_initial_all_down(L1, L2)
    dt = 0.01; nsteps = 3
    U = exp(-im * dt * H_torus)
    for _ in 1:nsteps
        psi_torus = U * psi_torus
    end
    z_ref = kagome_torus_local_z_per_sublattice(psi_torus, L1, L2)

    state = product_pess(NineSiteKagomeUC(), :down; D = 4)
    history = run_projected_pxp_kagome!(state, dt, nsteps;
                                         order = :second, maxdim = 4, cutoff = 1e-12)
    @test length(history) == nsteps

    z_ipeps_A = real(local_expectation_kagome(state, KagomeCoord(0, 0, :A), pauli_z()))
    z_ipeps_B = real(local_expectation_kagome(state, KagomeCoord(0, 0, :B), pauli_z()))
    z_ipeps_C = real(local_expectation_kagome(state, KagomeCoord(0, 0, :C), pauli_z()))

    @test isapprox(z_ipeps_A, z_ref.A; atol = 1e-3)
    @test isapprox(z_ipeps_B, z_ref.B; atol = 1e-3)
    @test isapprox(z_ipeps_C, z_ref.C; atol = 1e-3)
end
```

- [ ] **Step 3: Run and record the result**

```bash
julia --project=/Users/ren/Codex/PEPs -e 'using Pkg; Pkg.test(test_args=["kagome projected PXP iPEPS matches"])' 2>&1 | tail -25
```

Three possible outcomes (per spec's "Acceptance Criteria → Indicator"):

- **PASS within `1e-3`**: SU dynamics is acceptable for first scientific scarfinder runs. Continue with Task 13. Note in the commit message: "Layer-2 indicator: PASS at 1e-3 absolute on per-sublattice <Z> (D=4, dt=0.01, 3 steps)."

- **FAIL by 2-10× tolerance**: Borderline. Continue with Task 13 but note the result. Commit message: "Layer-2 indicator: BORDERLINE (max deviation Xe-3, target 1e-3). NTU follow-up should be prioritized before scientific use."

- **FAIL by orders of magnitude**: SU dynamics is unfit. Investigate first to rule out a kernel bug or basis convention mismatch in the torus builder (these would be hard failures, not the SU-vs-NTU question). If the kernel passes Layers 1 and 3 but fails Layer-2 by orders of magnitude, the SU mean-field environment really is inadequate for this dynamics — this is the case where we halt scientific claims and dispatch the NTU PR. Commit message: "Layer-2 indicator: FAIL (max deviation X, target 1e-3). SU prototype is unfit for kagome PXP dynamics at this scale; NTU follow-up is required before any scientific use. See Notes/2026-05-09-kagome-pess-ntu-followup.md."

In **all three cases**, commit. The PR is not blocked by this task's outcome.

```bash
git add test/util_kagome_finite_ed.jl test/test_kagome_evolution.jl
git commit -m "test: kagome torus integration indicator (Layer-2)"
```

---

## Task 13: Layer-3 Stabilizer Benchmark

**Files:**
- Modify: `src/SolvableModels.jl`
- Modify: `test/test_solvable_models.jl`
- Modify: `test/test_kagome_evolution.jl`

- [ ] **Step 1: Add closed-form helper for kagome cluster Hamiltonian**

In `src/SolvableModels.jl`:

```julia
"""
    kagome_cluster_center_z_expectation_exact(t; initial = :z_plus)

Closed-form `<Z_c>(t)` for the kagome cluster star Hamiltonian
`K = X_c * prod_{u in NN(c)} Z_u` starting from all-Z+ product. Same form
as the triangular case: `cos(2t)` because K^2 = I.
"""
function kagome_cluster_center_z_expectation_exact(t::Real; initial::Symbol = :z_plus)
    if initial === :z_plus
        return cos(2t)
    else
        throw(ArgumentError("supported initial states: :z_plus"))
    end
end
```

Export it. Add a one-line test in `test_solvable_models.jl`.

- [ ] **Step 2: Add benchmark test**

In `test/test_kagome_evolution.jl`:

```julia
@testset "kagome cluster stabilizer benchmark at D=4" begin
    state = product_pess(NineSiteKagomeUC(), :up; D = 4)
    H = kagome_cluster_star_hamiltonian()
    dt = 0.05; nsteps = 4
    for _ in 1:nsteps
        evolve_step_kagome!(state, H, dt; order = :second, update = :simple,
                            evolution = :real, projected = false)
    end
    t_total = nsteps * dt
    z_expected = kagome_cluster_center_z_expectation_exact(t_total; initial = :z_plus)
    z_measured = real(local_expectation_kagome(state, KagomeCoord(0, 0, :A), pauli_z()))
    @test isapprox(z_measured, z_expected; atol = 1e-6)
end
```

(`evolve_step_kagome!(state, H, dt; ...)` is the gate-builder variant analogous to the triangular `evolve_step!`. Implement it in `KagomeEvolution.jl` if not already present from Task 11.)

- [ ] **Step 3: Run, commit**

```bash
git add src/SolvableModels.jl test/test_solvable_models.jl test/test_kagome_evolution.jl
git commit -m "test: kagome cluster stabilizer benchmark at D=4"
```

---

## Task 14: KagomeScarFinder + Smoke Test + README

**Files:**
- Create: `src/KagomeScarFinder.jl`
- Create: `test/test_kagome_scar_finder.jl`
- Modify: `src/KagomePXPDynamics.jl`
- Modify: `test/runtests.jl`
- Modify: `README.md`

- [ ] **Step 1: Implement KagomeScarFinder** following the triangular `ScarFinder.jl` pattern. Same `Config` / `Candidate` / `scar_search_kagome` / `rank_candidates_kagome` API, just typed against `KagomePESS`. Adapt the seed function to use `product_pess` / `random_pess`.

- [ ] **Step 2: Smoke test**

```julia
@testset "kagome scar_search at D=4 NineSiteKagomeUC exercises bond growth and truncation" begin
    cfg = KagomeScarFinderConfig(0.005, 1, 2, 4, 1e-12, NineSiteKagomeUC(), 2, 1.0;
                                 scar_maxdim = 4)
    candidates = scar_search_kagome(cfg; seed = 53)
    @test length(candidates) == 2
    @test all(isfinite(c.score) for c in candidates)
    saw_growth = any(d.max_bond_dim > 1
                     for c in candidates
                     for d in c.diagnostics)
    @test saw_growth
end
```

- [ ] **Step 3: README update**

Add a "Kagome Status" section to `README.md` honestly documenting what the SU prototype proves and doesn't prove. Substitute `<LAYER_2_RESULT>` with the actual outcome from Task 12:

```markdown
## Kagome Status

Kagome PESS PXP evolution is implemented on `NineSiteKagomeUC` at `D = 4-8` as a **Simple
Update prototype**. The pipeline includes 5-site projected PXP gates, cluster contraction
with HOSVD decomposition back to PESS form, lambda truncation per simplex bond, a 3-color
Trotter schedule, and a ScarFinder driver that runs at `D >= 2` with real bond growth and
truncation diagnostics.

**Dynamics fidelity status**: Layer-1 kernel ED tests (5-site cluster) and Layer-3 analytic
stabilizer benchmark pass within tolerance, demonstrating the kernel is correct in regimes
where SU is exact or near-exact. Layer-2 finite-torus integration (the main dynamics-fidelity
indicator at scale) result: <LAYER_2_RESULT — replace with one of: PASS at 1e-3, BORDERLINE
at Xe-3, FAIL at X>. Treat scientific scarfinder results from this implementation accordingly.

**NTU follow-up**: The Simple Update kernel is a prototype. NTU (neighborhood tensor update)
is the production-grade dynamics algorithm and is the planned next PR. All infrastructure
in this PR — geometry, state container, models, gates, scheduler, observables, ScarFinder
loop — is reusable; only the kernel changes. See `Notes/2026-05-09-kagome-pess-ntu-followup.md`.

**Not supported**: OneSiteKagomeUC and other small-UC variants (sublattice aliasing).
BP-gauge maintenance, PEPS expectation values via boundary MPS or CTMRG, imaginary-time
energy correction in ScarFinder are future work.

Triangular code has been removed in commit eef0ead; see git history. Don't extend it; build on the kagome path going forward.
```

- [ ] **Step 4: Final full-suite run**

```bash
julia --project=/Users/ren/Codex/PEPs -e 'using Pkg; Pkg.test()' 2>&1 | tail -15
```

Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add src/KagomeScarFinder.jl src/KagomePXPDynamics.jl test/test_kagome_scar_finder.jl test/runtests.jl README.md
git commit -m "feat: kagome ScarFinder driver with D>=2 support and README update"
```

---

## Self-Review

After all tasks land:

1. **Spec coverage.** Each section of `docs/superpowers/specs/2026-05-09-kagome-pess-pxp-dynamics-design.md` should map to at least one task:
   - Geometry → Task 1.
   - Data Model (`KagomePESS`) → Task 5.
   - Algorithm (5-site cluster + HOSVD) → Tasks 8, 9 (sub-tasks 9.1-9.5).
   - Validation Layer 1 → Tasks 9.1, 10.
   - Validation Layer 2 → Task 12.
   - Validation Layer 3 → Task 13.
   - API Surface → all tasks.
   - Files To Create → mapped across tasks.
   - Acceptance Criteria → Task 14 smoke test + final full suite.

2. **No placeholders.** The plan contains two `error("... pending")` sentinels in Tasks 9.4 and 9.5 because the HOSVD decomposition and writeback algorithms are non-trivial enough to require implementation iteration. The executing agent must replace these with working code before declaring those tasks done.

3. **Type consistency.** Function signatures from the spec are used in the plan; the executing agent should verify them at implementation time. The `SimpleUpdateDiagnostics` vs `KagomeSimpleUpdateDiagnostics` decision (introduce a kagome variant or generalize the struct) is flagged in Task 7 Step 2.

4. **Known uncertainties for the executing agent:**
   - The exact ITensors `svd` keyword names (`lefttags`, `righttags`, return shape) — verify with `?svd` if a call fails on signature.
   - The exact down-triangle indexing convention may need to be adjusted in Task 1 to make the 9-site UC enumeration test pass.
   - HOSVD details in Task 9.4 require careful index labeling; an iteration loop is expected.
   - The kagome torus geometry in Task 12 needs to wrap consistently with the lattice; if the 12-site choice doesn't, fall back to 9-site or 18-site as the test instructs.

5. **Open follow-ups (not in this plan, recorded in spec's "Out-Of-Scope Follow-Ups" section):** BP-gauge maintenance, NTU, PEPS expectation values via boundary methods, imaginary-time energy correction in ScarFinder, larger UCs, D > 8 performance work, package rename.
