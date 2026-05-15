# Square iPEPS Evolution Engine Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use `superpowers:subagent-driven-development` (recommended) or `superpowers:executing-plans` to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build a reliable fixed-bond-dimension square-lattice iPEPS evolve-project engine with an explicit backend decision layer, diagnostics, and tests before adding ScarFinder orchestration.

**Architecture:** Preserve the existing finite `SquarePEPSState` as a small-system test oracle. First test whether PEPSKit.jl can own the infinite PEPS container, simple-update weights/interface, CTMRG environments, local measurements, and generic unit cells. If PEPSKit supports the required five-site star projection cleanly, this repo stays focused on square PXP conventions, dense five-site gates, star scheduling, ScarFinder orchestration, and finite exact oracle tests; otherwise implement only a custom simple-update projection backend with an abstract interface that can later swap in PEPSKit/CTMRG.

**Tech Stack:** Julia 1.12, LinearAlgebra, Test, Aqua, StableRNGs for randomized tests once random iPEPS constructors are introduced. Backend spike candidates: PEPSKit.jl 0.7, TensorKit.jl 0.15, and ITensors.jl 0.9 only for the custom fallback path.

---

## File Map

- Modify: `src/SquarePXPDynamics.jl`
  Root include/import/export entrypoint. Every exported symbol needs a docstring because `test/test_public_docs.jl` checks `Docs.hasdoc`.
- Keep unchanged except for possible shared helper use: `src/SquarePEPS.jl`
  Existing finite PEPS container and product constructor. Do not turn this into iPEPS.
- Modify: `src/SquareGeometry.jl`
  Keep `SquareCoord`, `square_neighbor`, `square_star_sites`, and `square_star_color`. Add direction constants only if they are public and documented.
- Create: `src/IPEPSBackends.jl`
  Backend-neutral interfaces for state construction, projection, evolution, measurements, and environment refresh. This file owns abstract types and small dispatch wrappers, not tensor algorithms.
- Create if PEPSKit route passes: `src/PEPSKitBackend.jl`
  Adapter from square PXP conventions into PEPSKit `InfinitePEPS`, `SUWeight`, `LocalOperator`, `CTMRGEnv`, `leading_boundary`, `reduced_densitymatrix`, and `expectation_value`.
- Create only if PEPSKit star projection fails cleanly: `src/CustomSimpleUpdateBackend.jl`
  ITensors-based five-site simple-update projection backend. Keep this limited to projection and cheap diagnostics; do not add CTMRG here.
- Create: `src/SquareUnitCells.jl`
  Periodic rectangular unit cells, wrapping, periodic neighbors, color centers, and unit-cell disjointness checks.
- Create: `src/SquareIPEPS.jl`
  Public repo-facing iPEPS facade over the selected backend. With PEPSKit, this should wrap PEPSKit state plus metadata; with the fallback backend, this owns custom ITensors storage.
- Create if custom fallback is selected: `src/LinkWeights.jl`
  Link-weight normalization, safe absorb/deabsorb, entropy, basic gauge diagnostics for the ITensors backend. With PEPSKit, prefer `SUWeight`.
- Create: `src/SquarePXPGates.jl`
  Backend-neutral dense five-site PXP gate/operator construction. Add TensorKit/PEPSKit conversion if PEPSKit is selected; add ITensor conversion only for the custom fallback.
- Create: `src/Observables.jl`
  Measurement facade. Prefer PEPSKit CTMRG-backed `expectation_value`/`reduced_densitymatrix`; keep cheap product-state checks for tests.
- Create later: `src/StarSimpleUpdate.jl`
  Backend-neutral star projection dispatch plus PEPSKit or custom implementation selected by the decision spike.
- Create later: `src/IPEPSEvolution.jl`
  Trotter scheduling, color sweeps, evolution logs.
- Create later: `src/ScarFinder.jl`
  Thin orchestration layer over evolution, energy correction, observables, and ranking.
- Do not create until PEPSKit CTMRG has been tested and rejected: `src/CTMRG.jl`
  Custom environment contraction and full-update/gauge-improvement infrastructure. This is not part of the first implementation route.

## State S0: Baseline Locked

**Purpose:** Confirm the current repo conventions before adding new iPEPS code.

**Files:**
- Read: `README.md`
- Read: `notes/README.md`
- Read: `src/SpinOps.jl`
- Read: `src/SquareGeometry.jl`
- Read: `src/SquarePXP.jl`
- Read: `src/SquarePEPS.jl`
- Read: `test/runtests.jl`

- [ ] **Step 1: Run existing tests**

Run:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: existing tests pass before feature work starts.

- [ ] **Step 2: Record invariants in the PR description**

Include these invariants:

```text
- `:up` / basis index 1 is the Rydberg or excited state.
- `:down` / basis index 2 is the unexcited/vacancy state.
- Dense star order is `(center, right, up, left, down)`.
- Direction order is `:right`, `:up`, `:left`, `:down`.
- Existing finite `SquarePEPSState` remains unchanged except for non-breaking helper reuse.
```

- [ ] **Step 3: Commit baseline note only if needed**

If the branch needs a documentation-only checkpoint:

```bash
git add notes/2026-05-15-chatgpt-pro-ipeps-review-plan.md notes/2026-05-15-ipeps-superpowers-multistage-plan.md
git commit -m "docs: add square ipeps implementation plan"
```

## State S0.5: Backend Decision Spike

**Purpose:** Decide the iPEPS backend before implementing any repo-native iPEPS tensor storage or CTMRG code.

**Current PEPSKit investigation result:** PEPSKit.jl 0.7 has the core infrastructure this project wants to avoid reimplementing. Its public docs and source provide `InfinitePEPS` storage with periodic unit cells, `SUWeight` bond weights, `CTMRGEnv` plus `leading_boundary`, `LocalOperator` with arbitrary-length local terms, `reduced_densitymatrix`, `expectation_value`, and generic unit-cell support. Its public simple-update evolution supports two-site nearest-neighbor updates and a three-site cluster path; a custom five-site square-star update is the unresolved integration point.

**PEPSKit references checked on 2026-05-15:**
- PEPSKit home/docs: https://quantumkithub.github.io/PEPSKit.jl/stable/
- PEPSKit library API: https://quantumkithub.github.io/PEPSKit.jl/dev/lib/lib/
- Heisenberg simple-update example: https://quantumkithub.github.io/PEPSKit.jl/stable/examples/heisenberg_su/
- Three-site simple-update example: https://quantumkithub.github.io/PEPSKit.jl/stable/examples/j1j2_su/
- PEPSKit v0.7.0 release notes: https://github.com/QuantumKitHub/PEPSKit.jl/releases/tag/v0.7.0

**Files:**
- Create: `src/IPEPSBackends.jl`
- Create if the PEPSKit route passes: `src/PEPSKitBackend.jl`
- Create only if the PEPSKit star-projection route fails: `src/CustomSimpleUpdateBackend.jl`
- Modify after decision: `Project.toml`
- Create: `test/test_backend_decision.jl`
- Modify: `test/runtests.jl`

- [ ] **Step 1: Record PEPSKit capability checklist**

Record this checklist in the PR description and keep it updated with exact PEPSKit version and Julia version:

```text
PEPSKit capability check:
- InfinitePEPS state storage: yes, via `InfinitePEPS(...; unitcell=(Nr, Nc))`.
- Simple-update weights: yes, via `SUWeight(peps)`.
- Built-in simple update: yes for two-site nearest-neighbor gates and PEPSKit's three-site cluster path.
- Five-site square-star update: must be tested; do not assume PEPSKit's public `time_evolve` handles it.
- CTMRG environments: yes, via `CTMRGEnv` and `leading_boundary`.
- Local observables: yes, via `LocalOperator`, `expectation_value`, and `reduced_densitymatrix`.
- Generic unit cells: yes, via matrix/unitcell constructors and periodic indexing.
```

- [ ] **Step 2: Add PEPSKit compatibility spike environment**

Before adding PEPSKit as a main dependency, test it in a temporary Julia environment:

```bash
julia --project=/tmp/ipeps-pepskit-spike -e 'using Pkg; Pkg.activate("/tmp/ipeps-pepskit-spike"); Pkg.add(["PEPSKit", "TensorKit"]); using PEPSKit, TensorKit; @show pkgversion(PEPSKit)'
```

Expected: PEPSKit loads on the current Julia version. If it does not, document the resolver error and choose the custom simple-update backend for the first implementation while preserving the abstract backend facade.

- [ ] **Step 3: Test PEPSKit state and weight construction**

In the spike, construct a `10 x 10` trivial-symmetry spin-1/2 iPEPS and its simple-update weights:

```julia
using PEPSKit, TensorKit
Nr, Nc = 10, 10
Dbond = 1
peps = InfinitePEPS(rand, Float64, ℂ^2, ℂ^Dbond; unitcell = (Nr, Nc))
wts = SUWeight(peps)
@assert size(peps) == (Nr, Nc)
@assert size(wts) == (2, Nr, Nc)
```

Expected: construction succeeds and dimensions match the square unit cell and horizontal/vertical weight layout.

- [ ] **Step 4: Test PEPSKit local observable path**

Construct at least one one-site density operator and one five-site square-star `LocalOperator` using repo star order `(center, right, up, left, down)` mapped to PEPSKit `CartesianIndex` sites:

```julia
lattice = fill(ℂ^2, 10, 10)
center = CartesianIndex(5, 5)
sites = (
    center,
    CartesianIndex(5, 6),
    CartesianIndex(4, 5),
    CartesianIndex(5, 4),
    CartesianIndex(6, 5),
)
Ostar = LocalOperator(lattice, sites => pxp_star_tensor_map)
```

Expected: PEPSKit accepts the five-site term. If it rejects the term shape or site layout, keep PEPSKit for container/CTMRG only if a measurement adapter is still practical; otherwise choose the custom fallback and defer CTMRG.

- [ ] **Step 5: Test PEPSKit CTMRG measurement route**

Use a small environment dimension on a product or random state:

```julia
χ = 2
env0 = CTMRGEnv(peps, ℂ^χ)
env, info = leading_boundary(env0, peps; tol = 1e-8, maxiter = 20, verbosity = 1)
value = expectation_value(peps, Ostar, env)
@assert isfinite(real(value))
```

Expected: CTMRG converges or returns finite diagnostics for the small smoke case. Do not implement repo-native CTMRG if this path works.

- [ ] **Step 6: Test five-site projection feasibility**

Try the PEPSKit-first projection route in this order:

```text
1. Prefer a public PEPSKit API that can evolve or project an arbitrary five-site `LocalOperator`.
2. If no public API exists, test whether a small adapter can reuse PEPSKit `InfinitePEPS`, `SUWeight`, `absorb_weight`, TensorKit SVD/truncation, and only stable exported APIs.
3. Do not build on private PEPSKit internals unless the adapter is isolated and version-pinned behind `PEPSKitBackend`.
```

Decision rule:

```text
Choose `PEPSKitBackend` if the five-site star projection can be implemented with PEPSKit state/weights and a small, tested adapter.
Choose `CustomSimpleUpdateBackend` if the five-site star projection requires broad private-API dependency, fights PEPSKit tensor layout, or cannot produce clear truncation diagnostics.
In both cases keep `ProjectionBackend` abstract so ScarFinder only calls `project_star!`, `evolve!`, `measure`, and `refresh_environment!`.
```

- [ ] **Step 7: Commit backend decision**

Commit the decision note and any spike tests that remain in the repo:

```bash
git add notes/2026-05-15-ipeps-superpowers-multistage-plan.md test/test_backend_decision.jl src/IPEPSBackends.jl Project.toml
git commit -m "docs: add ipeps backend decision layer"
```

## State S1: Periodic Unit Cells And iPEPS State

**Purpose:** Add the periodic rectangular iPEPS state model without time evolution, using the backend selected in S0.5.

**Files:**
- Create: `src/SquareUnitCells.jl`
- Create: `src/SquareIPEPS.jl`
- Create or modify: `src/IPEPSBackends.jl`
- Create if selected: `src/PEPSKitBackend.jl`
- Create only if fallback selected: `src/CustomSimpleUpdateBackend.jl`
- Modify: `src/SquarePXPDynamics.jl`
- Create: `test/test_square_unitcells.jl`
- Create: `test/test_square_ipeps.jl`
- Modify: `test/runtests.jl`

- [ ] **Step 1: Write periodic unit-cell tests**

Add tests covering:

```julia
cell = PeriodicSquareUnitCell(10, 10)
@test wrap(cell, SquareCoord(11, 1)) == SquareCoord(1, 1)
@test wrap(cell, SquareCoord(0, 1)) == SquareCoord(10, 1)
@test neighbor(cell, SquareCoord(10, 1), :right) == SquareCoord(1, 1)
@test neighbor(cell, SquareCoord(1, 1), :left) == SquareCoord(10, 1)
@test length(cell.reps) == 100
@test_throws ArgumentError PeriodicSquareUnitCell(0, 10)
@test_throws ArgumentError assert_five_color_compatible(PeriodicSquareUnitCell(4, 4))
```

Also test `update_centers(cell, color)` and `stars_are_disjoint_mod_unitcell(cell, centers)` for all five colors on a `10 x 10` cell.

- [ ] **Step 2: Implement `PeriodicSquareUnitCell`**

Define:

```julia
struct PeriodicSquareUnitCell <: SquareUnitCell
    Lx::Int
    Ly::Int
    reps::Vector{SquareCoord}
end
```

Constructor rules:

```text
Lx >= 1
Ly >= 1
reps = [SquareCoord(x, y) for y in 1:Ly for x in 1:Lx]
```

Use one-based coordinates to match the existing finite PEPS constructor.

- [ ] **Step 3: Implement periodic helpers**

Public documented helpers:

```julia
wrap(cell::PeriodicSquareUnitCell, c::SquareCoord)::SquareCoord
neighbor(cell::PeriodicSquareUnitCell, c::SquareCoord, dir::Symbol)::SquareCoord
update_centers(cell::PeriodicSquareUnitCell, color::Integer)::Vector{SquareCoord}
assert_five_color_compatible(cell::PeriodicSquareUnitCell)
stars_are_disjoint_mod_unitcell(cell::PeriodicSquareUnitCell, centers)
```

Keep `square_star_color` in `SquareGeometry.jl`; treat it as a scheduler, not as the physical unit cell.

- [ ] **Step 4: Write iPEPS constructor tests**

Cover:

```julia
cell = PeriodicSquareUnitCell(10, 10)
psi = product_square_ipeps(cell; state = :down, maxdim = 1)

@test psi.maxdim == 1
@test projection_backend(psi) isa AbstractProjectionBackend
@test length(unitcell_reps(psi)) == 100
@test all(physical_dim(psi, c) == 2 for c in cell.reps)
@test all(simple_weight_dim(psi, c, dir) == 1 for c in cell.reps for dir in (:right, :up))
@test density_simple(psi) ≈ 0 atol=1e-14
@test blockade_violation_simple(psi) ≈ 0 atol=1e-14
```

Checkerboard tests:

```julia
psi = checkerboard_square_ipeps(cell; excited_on = :even, maxdim = 1)
dens = sublattice_densities(psi)
@test dens.even ≈ 1 atol=1e-14
@test dens.odd ≈ 0 atol=1e-14
@test blockade_violation_simple(psi) ≈ 0 atol=1e-14
```

- [ ] **Step 5: Implement backend-neutral `SquareIPEPSState`**

Define the public facade so downstream code does not care whether PEPSKit or the custom fallback owns tensor storage:

```julia
abstract type AbstractProjectionBackend end

struct SquareIPEPSState{B<:AbstractProjectionBackend,S,W}
    unitcell::PeriodicSquareUnitCell
    backend::B
    state::S
    weights::W
    maxdim::Int
    metadata::Dict{Symbol,Any}
end
```

Required public helpers:

```julia
projection_backend(psi::SquareIPEPSState)
unitcell_reps(psi::SquareIPEPSState)
physical_dim(psi::SquareIPEPSState, c::SquareCoord)::Int
simple_weight_dim(psi::SquareIPEPSState, c::SquareCoord, dir::Symbol)::Int
copy_state(psi::SquareIPEPSState)
```

With `PEPSKitBackend`, `state` should be a PEPSKit `InfinitePEPS` and `weights` should be a PEPSKit `SUWeight`. With `CustomSimpleUpdateBackend`, `state` may be a custom ITensors storage object and `weights` may be canonical bond-weight dictionaries.

- [ ] **Step 6: Implement custom fallback bond keys only if needed**

Use a canonical periodic bond key so every undirected bond is stored once:

```julia
struct BondKey
    site::SquareCoord
    dir::Symbol
end
```

Only canonical directions should be stored in `BondKey`: `:right` and `:up`. `bondkey(cell, c, :left)` maps to the left neighbor's `:right` bond; `bondkey(cell, c, :down)` maps to the down neighbor's `:up` bond.

Define:

```julia
struct CustomIPEPSStorage
    tensors::Dict{SquareCoord,ITensor}
    physical_indices::Dict{SquareCoord,Index}
    link_indices::Dict{Tuple{SquareCoord,Symbol},Index}
    link_weights::Dict{BondKey,Vector{Float64}}
end
```

- [ ] **Step 7: Implement product and checkerboard constructors**

Rules:

```text
- `maxdim >= 1`.
- `state` accepts only `:up` or `:down`.
- `excited_on` accepts only `:even` or `:odd`.
- PEPSKit backend: construct TensorKit product tensors inside `InfinitePEPS` and initialize `SUWeight(peps)`.
- Custom fallback: one physical index per representative, four virtual indices per representative direction, shared neighbor endpoint indices, and one normalized `lambda` vector per canonical bond.
```

For `maxdim = 1`, product tensors should have only one nonzero physical amplitude and all link weights equal `[1.0]`.

- [ ] **Step 8: Add minimal simple observables needed by constructor tests**

In `src/Observables.jl`, implement product-state-compatible helpers that work before CTMRG is introduced:

```julia
density_simple(psi::SquareIPEPSState; sublattice = nothing)
sublattice_densities(psi::SquareIPEPSState)
blockade_violation_simple(psi::SquareIPEPSState)
```

At this stage these may be exact for `maxdim = 1` product/checkerboard states and may throw `ArgumentError` for unsupported higher-bond states if the full simple environment is not implemented yet. Do not return misleading values.

- [ ] **Step 9: Run focused and full tests**

Run:

```bash
julia --project=. -e 'using Pkg; Pkg.test(test_args=["square_unitcells"])'
julia --project=. -e 'using Pkg; Pkg.test()'
```

If `test_args` is not wired in the test runner, run full `Pkg.test()` and note that focused test routing is unavailable.

## State S2: Backend Gate Conversion And Link Weights

**Purpose:** Keep dense PXP gates as the source of truth, then convert them into the selected backend representation and make simple-update weight algebra safe.

**Files:**
- Create if custom fallback is selected: `src/LinkWeights.jl`
- Create: `src/SquarePXPGates.jl`
- Modify if PEPSKit route is selected: `src/PEPSKitBackend.jl`
- Modify if custom fallback is selected: `src/CustomSimpleUpdateBackend.jl`
- Modify: `src/SquarePXPDynamics.jl`
- Create if custom fallback is selected: `test/test_link_weights.jl`
- Create: `test/test_square_pxp_gates_backend.jl`
- Modify: `test/runtests.jl`

- [ ] **Step 1: Test backend gate conversion against dense source of truth**

For computational basis inputs in star order `(center, right, up, left, down)`, compare backend gate/operator output coefficients against `projected_square_pxp_gate(dt; evolution = :real)`.

PEPSKit route:

```julia
gate = square_pxp_gate_tensormap(dt; evolution = :real, projected = true)
op = square_pxp_star_localoperator(lattice, center, gate)
@test op isa PEPSKit.LocalOperator
```

Custom fallback route:

```julia
gate = square_pxp_gate_itensor(dt, phys; evolution = :real, projected = true)
@test hasinds(gate, prime.(phys)..., phys...)
```

- [ ] **Step 2: Implement `square_pxp_gate_tensormap` for PEPSKit if selected**

Convert the existing dense `32 x 32` matrix into a TensorKit `TensorMap` with five output physical spaces and five input physical spaces. Preserve repo star order exactly:

```julia
square_pxp_gate_tensormap(
    step::Real;
    evolution::Symbol = :real,
    projected::Bool = true,
)
```

- [ ] **Step 3: Implement `square_pxp_star_localoperator` for PEPSKit if selected**

Map `SquareCoord` center/right/up/left/down to PEPSKit `CartesianIndex` sites and return:

```julia
LocalOperator(lattice, sites => gate)
```

The adapter must document the coordinate convention and must use periodic wrapping from `PeriodicSquareUnitCell`.

- [ ] **Step 4: Test lambda absorption round trip for custom fallback only**

Use a random small ITensor with a known index `i`, absorb `[0.8, 0.6]`, deabsorb it, and check the original tensor is recovered within tolerance when all lambda values are above `atol`.

- [ ] **Step 5: Test safe deabsorb on tiny lambda for custom fallback only**

Use `lambda = [1.0, 0.0]` and assert all output tensor elements are finite. The inverse for zero or tiny entries must be `0.0`, not `Inf`.

- [ ] **Step 6: Implement `absorb_lambda` and `deabsorb_lambda` for custom fallback only**

Use a diagonal ITensor with a temporary `sim(i)` index, contract, and replace the temporary index back to `i`.

- [ ] **Step 7: Implement weight entropy helpers**

Public documented helpers:

```julia
weight_entropy(lambda::AbstractVector{<:Real})::Float64
bond_entropy(psi::SquareIPEPSState, c::SquareCoord, dir::Symbol)::Float64
all_bond_entropies(psi::SquareIPEPSState)::Dict{Tuple{SquareCoord,Symbol},Float64}
normalize_link_weights!(psi::SquareIPEPSState)
```

Test:

```julia
@test weight_entropy([1.0]) ≈ 0 atol=1e-14
@test weight_entropy([1 / sqrt(2), 1 / sqrt(2)]) ≈ log(2) atol=1e-14
```

- [ ] **Step 8: Implement `square_pxp_gate_itensor` for custom fallback only**

Signature:

```julia
square_pxp_gate_itensor(
    step::Real,
    phys::NTuple{5,Index};
    evolution::Symbol = :real,
    projected::Bool = true,
)
```

The ITensor has primed output physical indices followed by unprimed input physical indices.

- [ ] **Step 9: Run tests**

Run:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

## State S3: QR-Reduced Five-Site Star Simple Update

**Purpose:** Implement the core backend-selected projection operation that makes ScarFinder possible.

**Files:**
- Create: `src/StarSimpleUpdate.jl`
- Modify: `src/SquareIPEPS.jl`
- Modify if PEPSKit route is selected: `src/PEPSKitBackend.jl`
- Modify if custom fallback is selected: `src/CustomSimpleUpdateBackend.jl`
- Modify if custom fallback is selected: `src/LinkWeights.jl`
- Modify: `src/Observables.jl`
- Modify: `src/SquarePXPDynamics.jl`
- Create: `test/test_star_simple_update.jl`
- Modify: `test/runtests.jl`

- [ ] **Step 1: Add identity-update test**

Construct a `10 x 10` checkerboard iPEPS with `maxdim = 1`, apply `project_star!(psi, center, 0.0; evolution = :real, projected = true)`, and compare product observables plus simple-update weights before and after.

- [ ] **Step 2: Add flippable-center product test**

From all-down product state, apply one isolated star update at small real `dt`. Check center density against `sin(dt)^2` using the local star observable.

- [ ] **Step 3: Define diagnostics**

Create:

```julia
struct StarUpdateInfo
    center::SquareCoord
    max_truncerr::Float64
    truncerrs::Dict{Symbol,Float64}
    keptdims::Dict{Symbol,Int}
    min_lambda::Dict{Symbol,Float64}
    norm_factors::Dict{Symbol,Float64}
end
```

- [ ] **Step 4: Implement star validation**

For a wrapped center, gather `center/right/up/left/down` representatives and throw `ArgumentError` if the five representatives are not distinct. Repeated representatives require index-copy handling and must stay out of this first implementation.

- [ ] **Step 5: Implement PEPSKit star projection first if selected**

The PEPSKit implementation should keep PEPSKit objects as the source of state and weight truth:

```text
- Read/write PEPSKit `InfinitePEPS` tensors through the facade.
- Read/write PEPSKit `SUWeight` entries for the four center-leaf bonds.
- Prefer exported PEPSKit/TensorKit operations for weight absorption, QR/SVD, and truncation.
- Return `StarUpdateInfo` with the same fields as the custom backend.
- Keep any use of PEPSKit internals isolated in `PEPSKitBackend.jl` with a version guard.
```

If this step requires broad private PEPSKit internals or cannot preserve star-order diagnostics, stop and switch S3 implementation to `CustomSimpleUpdateBackend` without changing the public `project_star!` interface.

- [ ] **Step 6: Implement custom fallback lambda absorption plan only if needed**

For the five-site patch:

```text
- Absorb each external lambda into the tensor carrying that external leg.
- Absorb each center-leaf internal lambda exactly once, preferably into the leaf tensor.
- Record the minimum lambda for every touched bond.
```

- [ ] **Step 7: Implement custom fallback QR leaf reduction only if needed**

For each leaf, factor external virtual legs away from the active physical-plus-center-bond legs using ITensors `factorize(...; ortho="left", which_decomp="qr")`.

- [ ] **Step 8: Apply the five-site gate to the reduced core only if custom fallback is used**

Contract:

```julia
theta = gate * Acenter * Rright * Rup * Rleft * Rdown
theta = noprime(theta)
```

The resulting tensor should contain the five output physical indices and four QR reduced indices, not the twelve external virtual legs.

- [ ] **Step 9: Sequentially split leaves by SVD only if custom fallback is used**

For each direction in `(:right, :up, :left, :down)`, SVD across `(leaf_physical_index, leaf_qr_reduced_index)`, truncate by `maxdim` and `cutoff`, normalize singular values into the corresponding `BondKey`, and pass the removed norm into the remaining core.

- [ ] **Step 10: Reattach QR factors and deabsorb external weights only if custom fallback is used**

Reconstruct each leaf tensor as `Qleaf * leaf_active`, deabsorb external lambdas safely, normalize tensors, and store updated Gamma tensors.

- [ ] **Step 11: Run tests and inspect diagnostics**

Run:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected:

```text
- Zero-step update has near-zero truncation error.
- Flippable-center density matches direct two-level formula for D=1.
- Truncation diagnostics are finite for forced `maxdim = 1`.
```

## State S4: iPEPS Evolution Driver

**Purpose:** Wrap star updates in deterministic real/imaginary Trotter sweeps.

**Files:**
- Create: `src/IPEPSEvolution.jl`
- Modify: `src/SquarePXPGates.jl`
- Modify: `src/SquarePXPDynamics.jl`
- Create: `test/test_ipeps_evolution.jl`
- Modify: `test/runtests.jl`

- [ ] **Step 1: Test Trotter schedules**

Expected first order:

```julia
[(1, dt), (2, dt), (3, dt), (4, dt), (5, dt)]
```

Expected second order:

```julia
[(1, dt/2), (2, dt/2), (3, dt/2), (4, dt/2), (5, dt),
 (4, dt/2), (3, dt/2), (2, dt/2), (1, dt/2)]
```

- [ ] **Step 2: Implement `TrotterParams` and `EvolutionLog`**

Include `dt`, `order`, `evolution`, `projected`, `maxdim`, `cutoff`, total applied time, and collected `StarUpdateInfo` entries.

- [ ] **Step 3: Implement `evolve!`**

For each Trotter layer, collect centers by color, assert stars are disjoint modulo unit cell, build the selected backend gate/operator for each star, call `project_star!`, and collect diagnostics. `evolve!` must not branch on PEPSKit vs custom tensor details except through backend dispatch.

- [ ] **Step 4: Test zero-time and short-time behavior**

`total_time = 0` should leave simple observables unchanged. A short real-time sweep from a blockaded product state should keep blockade violation finite and small enough to diagnose regressions.

- [ ] **Step 5: Run tests**

Run:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

## State S5: Observables And Diagnostics

**Purpose:** Provide measurements for development and early ScarFinder ranking, preferring PEPSKit CTMRG-backed observables whenever the PEPSKit route is active.

**Files:**
- Modify: `src/Observables.jl`
- Modify if PEPSKit route is selected: `src/PEPSKitBackend.jl`
- Modify if custom fallback is selected: `src/LinkWeights.jl`
- Modify: `src/SquarePXPDynamics.jl`
- Create: `test/test_observables.jl`
- Modify: `test/runtests.jl`

- [ ] **Step 1: Implement PEPSKit CTMRG-backed expectations if selected**

Use PEPSKit `LocalOperator`, `CTMRGEnv`, `leading_boundary`, `expectation_value`, and `reduced_densitymatrix` through backend dispatch:

```julia
refresh_environment!(psi; chi, tol, maxiter)
expectation(psi, observable; environment = current_environment(psi))
```

This route should support one-site density, nearest-neighbor blockade checks, and five-site PXP energy terms. Do not add repo-native CTMRG here.

- [ ] **Step 2: Implement one-site and two-site simple expectations only for custom fallback**

Use Gamma-lambda absorbed local patches and exact bra-ket contraction of the local simple-update environment.

- [ ] **Step 3: Implement star expectation**

PEPSKit route: use the five-site `LocalOperator` and `expectation_value`. Custom fallback route: use the same star patch construction as `project_star!`, replacing the evolution gate with the local operator under test.

- [ ] **Step 4: Implement public observables**

Public documented functions:

```julia
rydberg_density(psi)
sublattice_densities(psi)
blockade_violation_simple(psi)
pxp_energy_density(psi)
mean_bond_entropy(psi)
max_bond_entropy(psi)
```

- [ ] **Step 5: Add product and dense-local tests**

Test all-up, all-down, checkerboard, and a `|+>`-like one-site product state against dense five-site calculations. For PEPSKit, include at least one CTMRG smoke test at small `chi` and compare product-state results against exact dense expectations.

- [ ] **Step 6: Run tests**

Run:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

## State S6: ScarFinder Orchestration

**Purpose:** Add ScarFinder as a thin driver over the engine, not as another tensor algorithm implementation.

**Files:**
- Create: `src/ScarFinder.jl`
- Modify: `src/IPEPSEvolution.jl`
- Modify: `src/Observables.jl`
- Modify: `src/SquarePXPDynamics.jl`
- Create: `test/test_scarfinder.jl`
- Modify: `test/runtests.jl`

- [ ] **Step 1: Define parameter and result structs**

Add `ScarFinderParams`, `ScarFinderIteration`, and `ScarFinderResult` with fields for energy, truncation error, bond entropy, blockade violation, densities, and acceptance.

- [ ] **Step 2: Implement `scarfinder!`**

Loop:

```text
evolve! -> normalize_backend_state! -> optional energy_correct! -> measure -> append log
```

Do not put low-level tensor index logic in `ScarFinder.jl`; it should only call backend-neutral projection, normalization, environment, and measurement functions.

- [ ] **Step 3: Implement guarded energy correction**

Use short imaginary-time attempts to reduce `abs(E - Etarget)`, keep the best state, and return an explicit accepted/rejected status.

- [ ] **Step 4: Test algorithmic invariants**

Cover:

```text
- `iterations = 0` returns the input state.
- Fixed RNG seed gives identical logs.
- Every log row has finite diagnostics.
- `target_energy = nothing` skips correction.
- Correction enabled never silently worsens the result without marking rejection.
```

- [ ] **Step 5: Run tests**

Run:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

## State S7: CTMRG And Gauge-Fixed Full Update

**Purpose:** Add environment-quality measurements and later full-update/gauge-fixing infrastructure after the simple-update engine is stable, using PEPSKit CTMRG before considering any repo-native CTMRG.

**Files:**
- Modify if PEPSKit route is selected: `src/PEPSKitBackend.jl`
- Create only after a documented PEPSKit CTMRG blocker: `src/CTMRG.jl`
- Modify: `src/Observables.jl`
- Modify if custom fallback is selected: `src/LinkWeights.jl`
- Modify: `src/SquarePXPDynamics.jl`
- Create: `test/test_ctmrg.jl`
- Modify: `test/runtests.jl`

- [ ] **Step 1: Test PEPSKit CTMRG route on repo PXP states**

For PEPSKit-backed states, initialize `CTMRGEnv`, run `leading_boundary`, and measure one-site density, nearest-neighbor blockade, and five-site PXP energy density through `expectation_value`.

- [ ] **Step 2: Wrap PEPSKit CTMRG diagnostics**

Expose backend-neutral diagnostics from PEPSKit's returned `info`, including `chi`, iteration status when available, convergence error, truncation error, and condition number. Do not hide non-convergence.

- [ ] **Step 3: Implement environment-backed public expectations**

Add backend-neutral `expectation`, `correlation`, and `log_fidelity_density` only after the PEPSKit CTMRG smoke tests pass. If the custom simple-update backend is active, this step should either convert to a PEPSKit-compatible state or remain unsupported with a clear `ArgumentError`; do not write CTMRG from scratch.

- [ ] **Step 4: Document custom CTMRG blocker before creating `src/CTMRG.jl`**

Only create repo-native CTMRG if all of the following are true:

```text
- PEPSKit CTMRG cannot contract the selected state representation.
- A PEPSKit conversion adapter is impractical or fails tested PXP observables.
- The blocker is recorded with a minimal reproducer and exact PEPSKit version.
- A separate implementation plan is written for custom CTMRG.
```

- [ ] **Step 5: Implement gauge fixing only after environment backend is stable**

Add:

```julia
fix_bond_gauge!(psi, env, c::SquareCoord, dir::Symbol; rcond = 1e-12)
```

This should condition local norm matrices before ALS/full-update truncation, not replace simple-update link weights.

- [ ] **Step 6: Benchmark against known models**

Use PEPSKit.jl Heisenberg simple-update/CTMRG examples as external sanity checks for implementation shape and expected diagnostic scales before applying CTMRG to PXP ScarFinder runs.

## Non-Negotiable Rules For Implementers

- [ ] Every exported symbol has a docstring.
- [ ] Every mutating function ends in `!`.
- [ ] Never silently truncate without a diagnostics struct.
- [ ] Never compare raw tensors after gauge-changing operations.
- [ ] Never invert tiny lambda values without regularization.
- [ ] Keep `:up` as Rydberg/excited and `:down` as vacancy/unexcited.
- [ ] Keep dense 32x32 PXP gates as local source of truth.
- [ ] Prefer PEPSKit for iPEPS storage, `SUWeight`, CTMRG, measurements, and generic unit cells when the S0.5 spike passes.
- [ ] Keep the five-site projection interface backend-neutral; ScarFinder must not depend on PEPSKit or ITensors internals.
- [ ] Do not implement custom CTMRG until the PEPSKit route has been tested and rejected with a reproducer.
- [ ] Keep CTMRG and ScarFinder in separate PRs from the star-update milestone.
- [ ] Run `julia --project=. -e 'using Pkg; Pkg.test()'` before claiming a state is complete.
