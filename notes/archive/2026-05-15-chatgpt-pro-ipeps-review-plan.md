# ChatGPT Pro iPEPS Review And Implementation Plan

Date saved: 2026-05-15

This note preserves the review comment and implementation plan supplied for the square-lattice PXP/iPEPS restart.

## 1. Current repo state

At commit `2adff05`, the repo is an early square-lattice restart. The README says the current foundation is square geometry, dense five-site projected PXP gates, and a minimal ITensors-backed finite square PEPS product-state container. It also lists PEPS update kernels, blockade diagnostics, and ScarFinder as planned work.

The current package has only `ITensors` and `LinearAlgebra` as main dependencies, with Julia `1.12` compatibility. The current module entrypoint exports `SpinOps`, `SquareGeometry`, `SquarePXP`, and `SquarePEPS` functionality.

The current PXP convention is important: `projector_up()` is `|0><0|`, `projector_down()` is `|1><1|`, and the existing square-star PXP term is `X_center * P_down_right * P_down_up * P_down_left * P_down_down`. In other words, the repo currently treats `:up`/basis index `1` as the Rydberg/excited state and `:down`/basis index `2` as the unexcited state allowed around a flippable center.

The current `SquarePEPSState` is finite, not iPEPS: it stores one ITensor per finite site, physical indices, nearest-neighbor link indices, and open-boundary links of dimension 1. Keep this finite PEPS code for small exact tests, but do not try to stretch it directly into the infinite algorithm.

ScarFinder itself is an evolve-project iteration: initialize a low-entanglement variational state, evolve under the Hamiltonian for a projection time, project/truncate back onto the variational manifold, optionally enforce energy/symmetry constraints, and repeat. The paper emphasizes that the projection/truncation step is subtle because naive truncation can drift in energy or collapse to trivial low-entanglement states; they use energy correction by imaginary-time evolution when needed.

## 2. Target architecture

Add the following modules, in this order:

```text
src/
  SquareUnitCells.jl      # true periodic/iPEPS unit cells, separate from scheduler colors
  SquareIPEPS.jl          # infinite PEPS state, tensors, link weights, constructors
  LinkWeights.jl          # λ-vectors, absorb/deabsorb helpers, gauge diagnostics
  SquarePXPGates.jl       # ITensor versions of dense PXP gates
  StarSimpleUpdate.jl     # five-site star update with QR reduction + SVD truncation
  IPEPSEvolution.jl       # real/imaginary-time sweeps and Trotter schedules
  Observables.jl          # density, blockade, energy, local ops, bond entropies
  CTMRG.jl                # later: accurate environments
  ScarFinder.jl           # orchestration, logs, candidate ranking
```

Keep the existing finite `SquarePEPSState` and tests. Add a new type instead of changing the existing one:

```julia
struct SquareIPEPSState
    unitcell::PeriodicSquareUnitCell
    tensors::Dict{SquareCoord,ITensor}          # Γ tensors
    physical_indices::Dict{SquareCoord,Index}
    link_indices::Dict{Tuple{SquareCoord,Symbol},Index}
    link_weights::Dict{BondKey,Vector{Float64}} # λ on each virtual bond
    maxdim::Int
    gauge::Symbol                               # initially :simple
end
```

Use Gamma-lambda simple-update gauge as the first serious implementation target. Store PEPS tensors as bare Gamma tensors and store one positive vector `lambda_b` per virtual bond. The represented network is Gamma tensors connected by diagonal lambda matrices. This is the minimum viable proper gauging needed for meaningful truncation.

A true global canonical form does not exist for generic PEPS in the same way it does for MPS, so the first implementation should be honest: call it `:simple` or `:simple_update` gauge, not `:canonical`. Later add local CTMRG/full-update gauge fixing. Gauge fixing is known to improve iPEPS stability and convergence in full-update-style algorithms.

## 3. Unit-cell design: fix this before coding dynamics

Do not use the current `FiveSiteSquareUC` as the general iPEPS unit cell. It is useful as a five-color scheduler, but it is not enough as a physical iPEPS unit cell for general 2D scar states.

Add:

```julia
struct PeriodicSquareUnitCell
    Lx::Int
    Ly::Int
    reps::Vector{SquareCoord}
end
```

with:

```julia
wrap(cell, c)::SquareCoord
neighbor(cell, c, dir)::SquareCoord
bondkey(cell, c, dir)::BondKey
```

Also keep:

```julia
square_star_color(c) = mod(c.x + 2c.y, 5) + 1
```

but treat it as a Trotter schedule color, not as the physical unit cell.

For exact five-color star sweeps on a rectangular periodic unit cell, require `Lx % 5 == 0` and `Ly % 5 == 0`. If you also want a checkerboard/CDW ansatz in the same rectangular cell, require even `Lx` and `Ly`; the simple default is therefore `10 x 10`. That is large but unambiguous. Once the update engine works, add oblique Bravais unit cells to reduce this to a 10-site color-plus-parity cell.

Recommended initial constructors:

```julia
product_square_ipeps(cell; state = :down, maxdim = 1)
checkerboard_square_ipeps(cell; excited_on = :even, maxdim = 1)
random_square_ipeps(cell; maxdim, rng, scale = 1e-2, around = nothing)
```

Use the repo's current convention consistently:

```julia
rydberg_projector() = projector_up()
vacancy_projector() = projector_down()
rydberg_density_op() = projector_up()
```

This avoids future confusion where one agent assumes `|1>` is the Rydberg state while the current code uses `|0>`.

## 4. Gate representation

Keep the dense 32x32 PXP gate as the source of truth, but add an ITensor wrapper.

```julia
square_pxp_gate_itensor(
    step::Real,
    phys::NTuple{5,Index};
    evolution::Symbol = :real,
    projected::Bool = true,
)
```

Use the star order already in the repo:

```julia
(center, right, up, left, down)
```

The ITensor gate should have output physical indices primed and input physical indices unprimed:

```julia
G[p_center', p_right', p_up', p_left', p_down',
  p_center, p_right, p_up, p_left, p_down]
```

Then the local application pattern is:

```julia
theta = G * Acenter * Aright * Aup * Aleft * Adown
theta = noprime(theta)
```

ITensors supports tensor construction, contraction by `*`, and index-based element setting; its SVD treats chosen ITensor indices collectively as matrix row indices and supports `maxdim`/`cutoff` truncation.

## 5. The most important algorithm: five-site star simple update

This is the core of the project. Implement this before ScarFinder.

### 5.1 Why not contract the full five-site patch directly?

A direct five-site star patch has five physical legs and twelve external virtual legs. Its raw dimension scales like:

```text
2^5 D^12
```

which is already too large at modest `D`. The correct implementation should use QR/factorization to strip external legs before applying the gate.

ITensors has `factorize(A, Linds...; ortho="left", which_decomp="qr")`, which can be used as a QR-like factorization where the left factor is orthogonal and the right factor carries the active indices.

### 5.2 Star update outline

Implement:

```julia
apply_star_gate!(
    psi::SquareIPEPSState,
    center::SquareCoord,
    gate::ITensor;
    maxdim::Int = psi.maxdim,
    cutoff::Float64 = 1e-12,
    split_order = (:right, :up, :left, :down),
    svd_alg::String = "divide_and_conquer",
)::StarUpdateInfo
```

`StarUpdateInfo` should include:

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

### 5.3 Detailed steps

For a star centered at `c`, define:

```julia
coords = (
    center = c,
    right  = square_neighbor(c, :right),
    up     = square_neighbor(c, :up),
    left   = square_neighbor(c, :left),
    down   = square_neighbor(c, :down),
)
```

After wrapping through the iPEPS unit cell, assert that these five representatives are distinct. If not, throw an error for now. Supporting repeated representatives requires tensor-copy index renaming and should be a later feature.

#### Step A: copy local tensors

```julia
A0 = copy(psi.tensors[center])
Ar = copy(psi.tensors[right])
Au = copy(psi.tensors[up])
Al = copy(psi.tensors[left])
Ad = copy(psi.tensors[down])
```

#### Step B: absorb link weights

For every incident bond of every star tensor:

1. Absorb external lambda into the corresponding tensor.
2. Absorb each internal center-leaf lambda exactly once, preferably on the leaf side.
3. Never absorb the same lambda twice.

Add helpers:

```julia
absorb_lambda(T::ITensor, lambda::Vector, i::Index)::ITensor
deabsorb_lambda(T::ITensor, lambda::Vector, i::Index; atol = 1e-14)::ITensor
```

Implementation detail: create a diagonal ITensor with old index `i` and temporary index `j`, contract, then replace `j => i`.

```julia
function absorb_lambda(T, lambda, i)
    j = sim(i)
    Lambda = ITensor(eltype(T), i, j)
    for n in 1:dim(i)
        Lambda[i => n, j => n] = lambda[n]
    end
    return replaceind(T * Lambda, j, i)
end
```

For inverse lambda, use a safe inverse:

```julia
invlambda[n] = lambda[n] > atol ? inv(lambda[n]) : 0.0
```

#### Step C: QR-reduce the four leaves

For each leaf, separate external legs from the active physical-plus-center-bond legs.

Example for the `:right` leaf:

```julia
external = (right_leg, up_leg, down_leg) # not the left leg to center
Qright, Rright = factorize(Ar, external...; ortho = "left", which_decomp = "qr")
```

`Qright` carries the external legs plus a new reduced index. `Rright` carries that reduced index, the leaf physical index, and the internal bond to the center.

Do this for all four leaves. The center tensor has no external legs in a five-site star, so leave it unreduced.

#### Step D: apply the five-site gate only to the small core

Contract:

```julia
theta =
    gate *
    Acenter *
    Rright *
    Rup *
    Rleft *
    Rdown

theta = noprime(theta)
```

Now `theta` should contain:

```text
five output physical indices
four QR reduced indices
```

and no twelve external virtual legs. This is the key scaling improvement.

#### Step E: sequentially split leaves from the center core

For each leaf in `split_order`, SVD the current `theta/rest` across:

```julia
left_indices = (leaf_physical_index, leaf_qr_reduced_index)
```

Use ITensors SVD:

```julia
U, S, V = svd(rest, left_indices...; maxdim, cutoff, alg = svd_alg)
```

Extract singular values from `S`, normalize them safely, and store them as the new link weight for the center-leaf bond. For exact reconstruction while keeping normalized lambda, multiply the removed norm into `V` before continuing:

```julia
sraw = diag_values(S)
scale = norm(sraw)
lambdanew = scale > 0 ? sraw / scale : fill(1 / sqrt(length(sraw)), length(sraw))
V *= scale
```

Then:

```julia
leaf_active = U
rest = V
psi.link_weights[bondkey(center, leafdir)] = lambdanew
```

The new common link index must replace the SVD-generated indices in both `leaf_active` and `rest`, so all tensors use the canonical link index stored in `psi.link_indices`.

#### Step F: reattach QR factors and deabsorb external lambda

For each leaf:

```julia
leaf_full = Qleaf * leaf_active
leaf_full = deabsorb_external_lambdas(leaf_full, psi, leaf_coord, excluding_internal_bond)
psi.tensors[leaf_coord] = normalize_tensor(leaf_full)
```

For the center:

```julia
psi.tensors[center] = normalize_tensor(rest)
```

For real-time evolution, normalization is mostly gauge management. For imaginary-time correction, normalization is essential.

#### Step G: return diagnostics

Return `StarUpdateInfo`. Do not silently drop truncation errors. ScarFinder will need them.

## 6. Evolution layer

Add:

```julia
struct TrotterParams
    dt::Float64
    order::Int          # 1 or 2
    evolution::Symbol   # :real or :imaginary
    projected::Bool
    maxdim::Int
    cutoff::Float64
end
```

Implement:

```julia
trotter_sequence(params)::Vector{Tuple{Int,Float64}}
```

First order:

```julia
[(1, dt), (2, dt), (3, dt), (4, dt), (5, dt)]
```

Second order:

```julia
[(1, dt/2), (2, dt/2), (3, dt/2), (4, dt/2),
 (5, dt),
 (4, dt/2), (3, dt/2), (2, dt/2), (1, dt/2)]
```

Then:

```julia
evolve!(
    psi::SquareIPEPSState,
    model::SquarePXPModel,
    total_time::Real;
    trotter::TrotterParams,
)::EvolutionLog
```

Each color update should:

1. collect all unit-cell centers with that color,
2. assert their stars are disjoint modulo the unit cell,
3. apply `apply_star_gate!` to each center,
4. collect truncation diagnostics.

Use projected gates by default. The existing ScarFinder paper emphasizes that the method relies on real-time evolution followed by projection back to a low-entanglement manifold; in this implementation the projection is exactly the SVD truncation back to fixed iPEPS bond dimension.

## 7. Gauge-fixing roadmap

### Phase 1 gauge: simple-update Gamma-lambda gauge

This is mandatory.

Implement:

```julia
bond_entropy(psi, bond)::Float64
all_bond_entropies(psi)::Dict{BondKey,Float64}
gauge_deviation(psi, bond)::Float64
normalize_link_weights!(psi)
normalize_tensors!(psi)
```

`bond_entropy` should use:

```text
p_i = lambda_i^2 / sum_j lambda_j^2
S = -sum_i p_i log(p_i)
```

`gauge_deviation` should initially be a local diagnostic: absorb surrounding lambda, contract a one-bond reduced norm matrix, and check how diagonal it is in the lambda basis. It does not have to be perfect, but it should catch severe gauge corruption.

### Phase 2 gauge: local QR regauging

Before every truncating SVD, the QR-reduced patch already provides a local gauge. Add optional:

```julia
local_regauge_star!(psi, center)
```

This should only change gauges, not observables. Test this by finite contraction on small PEPS.

### Phase 3 gauge: CTMRG/full-update gauge

Only after simple update works, add environment-based local gauge fixing:

```julia
fix_bond_gauge!(
    psi,
    env::CTMEnvironment,
    bond::BondKey;
    rcond = 1e-12,
)
```

The goal is to make the local norm matrix near identity before an ALS/full-update truncation. This is the more proper PEPS gauge, but it is not the first milestone. iPEPS gauge fixing is a known stabilizer for full-update algorithms, but it requires an environment backend.

## 8. Observable implementation

Implement observables in two layers.

### 8.1 Simple-update observables

These are cheap and good enough for early ScarFinder iteration.

Add:

```julia
expect_onesite_simple(psi, c, op)::ComplexF64
expect_twosite_simple(psi, c, dir, op1, op2)::ComplexF64
expect_star_simple(psi, c, opstar)::ComplexF64
density_simple(psi; sublattice = nothing)
blockade_violation_simple(psi)
energy_density_simple(psi, model::SquarePXPModel)
```

For `energy_density_simple`, compute:

```text
<X_c prod_{j in nn(c)} P_vac_j>
```

using the same five-site lambda-absorbed star patch as the update, but with bra-ket contraction instead of gate application.

Minimum observables for ScarFinder:

```julia
rydberg_density(psi)
sublattice_densities(psi)
blockade_violation_density(psi)
pxp_energy_density(psi)
max_bond_entropy(psi)
mean_bond_entropy(psi)
```

### 8.2 CTMRG observables

Add later:

```julia
struct CTMEnvironment
    chi::Int
    corners::Dict
    edges::Dict
end

ctmrg(psi; chi, maxiter = 10_000, tol = 1e-10)
expectation(psi, env, obs)
correlation(psi, env, op1, op2, r)
log_fidelity_density(psi, phi, env)
```

For serious 2D iPEPS measurements, CTMRG is eventually necessary. The ScarFinder paper notes that extending the method to isotropic 2D PEPS is natural but computationally expensive, so expect this environment layer to be a major later milestone rather than a quick helper.

## 9. ScarFinder orchestration

Once evolution and observables exist, ScarFinder is straightforward.

```julia
struct ScarFinderParams
    projection_time::Float64
    trotter_dt::Float64
    iterations::Int
    maxdim::Int
    cutoff::Float64
    target_energy::Union{Nothing,Float64}
    energy_tol::Float64
    imaginary_dt::Float64
    max_energy_correction_steps::Int
    rng_seed::Int
end

struct ScarFinderIteration
    iter::Int
    energy::Float64
    energy_error::Union{Nothing,Float64}
    max_truncerr::Float64
    mean_bond_entropy::Float64
    max_bond_entropy::Float64
    blockade_violation::Float64
    density_A::Float64
    density_B::Float64
    accepted::Bool
end

struct ScarFinderResult
    state::SquareIPEPSState
    log::Vector{ScarFinderIteration}
end
```

Main function:

```julia
function scarfinder!(
    psi::SquareIPEPSState,
    model::SquarePXPModel,
    params::ScarFinderParams;
    observables = default_scarfinder_observables(),
)
    for n in 1:params.iterations
        evolve!(psi, model, params.projection_time;
                trotter = TrotterParams(
                    params.trotter_dt,
                    2,
                    :real,
                    true,
                    params.maxdim,
                    params.cutoff,
                ))

        normalize_link_weights!(psi)
        normalize_tensors!(psi)

        if params.target_energy !== nothing
            energy_correct!(
                psi, model, params.target_energy;
                dt = params.imaginary_dt,
                maxsteps = params.max_energy_correction_steps,
                tol = params.energy_tol,
            )
        end

        push!(log, measure_scarfinder_iteration(...))
    end
    return ScarFinderResult(psi, log)
end
```

Energy correction:

```julia
energy_correct!(psi, model, Etarget; dt, maxsteps, tol)
```

Algorithm:

1. Measure current energy `E0`.
2. If `abs(E0 - Etarget) < tol`, return success.
3. Try short imaginary-time steps with sign chosen to move energy toward target.
4. Keep the state with smallest `abs(E - Etarget)`.
5. Abort if the energy error increases for several consecutive attempts.

This matches the ScarFinder warning that truncation can break energy conservation and may require imaginary-time correction after projection.

Candidate ranking:

```julia
score =
    w_entropy * mean_bond_entropy +
    w_blockade * blockade_violation +
    w_energy * abs(energy - target_energy) -
    w_revival * revival_quality
```

For early development, do not require revival quality inside the ScarFinder loop. First rank by energy stability, blockade violation, and entanglement. Then run a separate real-time validation trajectory from the final state.

## 10. Tests: concrete design

Keep the existing tests. Add new tests in this order.

### 10.1 Unit-cell tests

File:

```text
test/test_square_unitcells.jl
```

Tests:

```julia
@testset "periodic square unit cell" begin
    cell = PeriodicSquareUnitCell(10, 10)
    @test wrap(cell, SquareCoord(11, 1)) == SquareCoord(1, 1)
    @test neighbor(cell, SquareCoord(10, 1), :right) == SquareCoord(1, 1)
    @test neighbor(cell, SquareCoord(1, 1), :left) == SquareCoord(10, 1)
end
```

Color compatibility:

```julia
@testset "five-color schedule compatibility" begin
    cell = PeriodicSquareUnitCell(10, 10)
    for color in 1:5
        centers = update_centers(cell, color)
        @test stars_are_disjoint_mod_unitcell(cell, centers)
    end
end
```

Negative test:

```julia
@test_throws ArgumentError assert_five_color_compatible(PeriodicSquareUnitCell(4, 4))
```

### 10.2 iPEPS construction tests

File:

```text
test/test_square_ipeps.jl
```

Tests:

```julia
@testset "product square iPEPS" begin
    cell = PeriodicSquareUnitCell(10, 10)
    psi = product_square_ipeps(cell; state = :down, maxdim = 1)

    @test length(psi.tensors) == 100
    @test all(dim(physical_index(psi, c)) == 2 for c in cell.reps)
    @test all(length(lambda) == 1 for lambda in values(psi.link_weights))
    @test density_simple(psi) approx 0 atol=1e-14
    @test blockade_violation_simple(psi) approx 0 atol=1e-14
end
```

Checkerboard:

```julia
psi = checkerboard_square_ipeps(cell; excited_on = :even, maxdim = 1)
@test blockade_violation_simple(psi) approx 0 atol=1e-14
@test sublattice_densities(psi).even approx 1 atol=1e-14
@test sublattice_densities(psi).odd  approx 0 atol=1e-14
```

### 10.3 Gate tensor tests

File:

```text
test/test_square_pxp_gates_itensor.jl
```

Tests:

1. ITensor gate equals dense gate on every computational basis vector.
2. Real-time gate is unitary.
3. Projected gate annihilates locally forbidden star basis states.
4. `step = 0` gives identity on allowed states.

Dense comparison:

```julia
Gdense = projected_square_pxp_gate(dt; evolution = :real)
Git = square_pxp_gate_itensor(dt, phys; evolution = :real, projected = true)

for basis_in in Iterators.product(ntuple(_ -> 1:2, 5)...)
    # apply ITensor gate to basis ket
    # compare coefficients to Gdense[:, col]
end
```

### 10.4 Link-weight and gauge helper tests

File:

```text
test/test_link_weights.jl
```

Tests:

1. `absorb_lambda` followed by `deabsorb_lambda` returns the original tensor.
2. Near-zero lambda does not produce `NaN` or `Inf`.
3. `bond_entropy([1.0]) == 0`.
4. `bond_entropy([1/sqrt(2), 1/sqrt(2)]) approx log(2)`.

### 10.5 Star-update tests

File:

```text
test/test_star_simple_update.jl
```

Core tests:

```julia
@testset "zero-step star update is identity" begin
    cell = PeriodicSquareUnitCell(10, 10)
    psi = checkerboard_square_ipeps(cell; maxdim = 1)
    psi0 = deepcopy(psi)

    G = square_pxp_gate_itensor(0.0, star_physinds(psi, SquareCoord(1,1)))
    info = apply_star_gate!(psi, SquareCoord(1,1), G; maxdim = 1)

    @test info.max_truncerr approx 0 atol=1e-14
    @test gauge_invariant_distance_simple(psi, psi0) approx 0 atol=1e-12
end
```

Known simple product test:

From all-down product state, a flippable center under a small real-time PXP gate should rotate locally. Check the center density after one isolated star update against the exact two-level formula:

```text
|down> -> cos(dt)|down> - i sin(dt)|up>
<n> = sin^2(dt)
```

Use a finite PEPS patch or local star observable to avoid relying on CTMRG.

Truncation test:

```julia
psi = random_square_ipeps(cell; maxdim = 2, rng = StableRNG(1))
info = apply_star_gate!(psi, c, G; maxdim = 1, cutoff = 0.0)
@test maximum(values(info.keptdims)) <= 1
@test info.max_truncerr >= 0
@test all(isfinite, values(info.truncerrs))
```

### 10.6 Evolution tests

File:

```text
test/test_ipeps_evolution.jl
```

Tests:

1. `total_time = 0` leaves observables unchanged.
2. First- and second-order Trotter schedules have the expected color/step pattern.
3. One full real-time sweep from a blockaded product state keeps blockade violation small.
4. Imaginary-time step changes norm but final normalized observables remain finite.

### 10.7 Observable tests

File:

```text
test/test_observables.jl
```

Tests:

```julia
@test density_simple(product_down) approx 0
@test density_simple(product_up) approx 1
@test blockade_violation_simple(product_down) approx 0
@test blockade_violation_simple(checkerboard) approx 0
@test blockade_violation_simple(product_up) > 0
```

Energy tests:

1. Any computational-basis product state has zero PXP energy because `X` has zero diagonal expectation.
2. A simple product `|+>`-like state has energy matching a direct five-site dense calculation in the simple-update environment.

### 10.8 ScarFinder orchestration tests

File:

```text
test/test_scarfinder.jl
```

Do not test discovering a new 2D scar in CI. Test algorithmic invariants.

Tests:

1. `iterations = 0` returns the input state.
2. Fixed RNG seed gives identical logs.
3. Every log row has finite energy, entropy, truncation error, and blockade violation.
4. If `target_energy = nothing`, no energy correction runs.
5. If energy correction is enabled on a small mock system, final `abs(E - Etarget)` is no worse than before correction, or the result is explicitly marked `accepted = false`.

### 10.9 Regression tests against finite exact contraction

Use the existing finite `SquarePEPSState` for exact small tests.

Add test utilities, not exported package functions:

```julia
dense_state_from_finite_peps(psi::SquarePEPSState)
dense_square_pxp_hamiltonian(width, height)
```

For `2 x 2` or `3 x 3` open-boundary systems:

1. Build a finite product PEPS.
2. Apply a single dense star gate exactly.
3. Apply the PEPS star update with sufficiently large `maxdim`.
4. Contract to dense state and compare up to gauge/norm.

This is the strongest protection against broken index ordering.

## 11. CI and coding rules for future agents

1. Every exported function must have a docstring, because the repo already has a public-docstring test.
2. Every mutating function must end in `!`.
3. Do not compare raw tensors after gauge-changing operations. Compare observables, dense contractions, bond spectra, or gauge-invariant distances.
4. Never truncate a bond without returning truncation diagnostics.
5. Never silently invert a tiny lambda. Use a regularized inverse and record the minimum lambda.
6. Keep physical convention explicit: `:up` is Rydberg/excited in the current code.
7. Keep dense 32x32 gates as truth for local tests.
8. Add `StableRNGs` to `test/Project.toml` once random iPEPS tests are introduced.
9. Avoid adding CTMRG and ScarFinder in the same PR. CTMRG is large enough to deserve its own milestone.
10. ScarFinder should call `evolve!`, `energy_correct!`, and `measure`; it should not contain low-level tensor index logic.

## 12. Recommended milestone order

### Milestone 1: unit cells and iPEPS state

Deliver:

```julia
PeriodicSquareUnitCell
SquareIPEPSState
product_square_ipeps
checkerboard_square_ipeps
link/bond indexing
```

Pass:

```text
test_square_unitcells.jl
test_square_ipeps.jl
```

### Milestone 2: ITensor gates and link-weight helpers

Deliver:

```julia
square_pxp_gate_itensor
absorb_lambda
deabsorb_lambda
bond_entropy
```

Pass:

```text
test_square_pxp_gates_itensor.jl
test_link_weights.jl
```

### Milestone 3: QR-reduced five-site star update

Deliver:

```julia
apply_star_gate!
StarUpdateInfo
```

Pass:

```text
test_star_simple_update.jl
```

This is the hardest milestone.

### Milestone 4: iPEPS evolution

Deliver:

```julia
TrotterParams
trotter_sequence
evolve!
EvolutionLog
```

Pass:

```text
test_ipeps_evolution.jl
```

### Milestone 5: simple observables

Deliver:

```julia
density_simple
sublattice_densities
blockade_violation_simple
energy_density_simple
all_bond_entropies
```

Pass:

```text
test_observables.jl
```

### Milestone 6: ScarFinder driver

Deliver:

```julia
ScarFinderParams
ScarFinderResult
scarfinder!
energy_correct!
candidate ranking utilities
```

Pass:

```text
test_scarfinder.jl
```

### Milestone 7: CTMRG and environment-gauge/full-update improvements

Deliver:

```julia
CTMEnvironment
ctmrg
expectation with CTMRG
log_fidelity_density
fix_bond_gauge!
```

This is where the code becomes a serious 2D iPEPS research tool rather than a simple-update prototype.

## 13. Minimal first PR instruction to give an agent

Start with only this:

```text
Implement `PeriodicSquareUnitCell`, `SquareIPEPSState`,
`product_square_ipeps`, `checkerboard_square_ipeps`,
and link-weight storage. Do not implement time evolution yet.

Requirements:
- Preserve the current spin convention: `:up` is the Rydberg/excited state.
- Keep existing finite `SquarePEPSState` unchanged.
- Add tests for periodic wrapping, neighbor maps, bond keys, product iPEPS,
  checkerboard iPEPS, density, and blockade violation.
- Every exported symbol must have a docstring.
- `Pkg.test()` must pass.
```

The second PR should implement ITensor gate conversion and link-weight absorb/deabsorb. The third PR should implement the QR-reduced star update. That separation is important: if the star update fails, you want to know whether the bug is in the gate ordering, lambda absorption, QR split, SVD reconstruction, or observable measurement.

## Source Links From Original Review

- GitHub tree at commit `2adff05717469fc81dc48cbc0382a90ad08a38bc`: https://github.com/jayren3996/SquarePXPDynamics.jl/tree/2adff05717469fc81dc48cbc0382a90ad08a38bc
- Project.toml at that commit: https://raw.githubusercontent.com/jayren3996/SquarePXPDynamics.jl/2adff05717469fc81dc48cbc0382a90ad08a38bc/Project.toml
- Module entrypoint at that commit: https://raw.githubusercontent.com/jayren3996/SquarePXPDynamics.jl/2adff05717469fc81dc48cbc0382a90ad08a38bc/src/SquarePXPDynamics.jl
- SpinOps at that commit: https://raw.githubusercontent.com/jayren3996/SquarePXPDynamics.jl/2adff05717469fc81dc48cbc0382a90ad08a38bc/src/SpinOps.jl
- SquarePEPS at that commit: https://raw.githubusercontent.com/jayren3996/SquarePXPDynamics.jl/2adff05717469fc81dc48cbc0382a90ad08a38bc/src/SquarePEPS.jl
- ScarFinder paper: https://arxiv.org/html/2504.12383v2
- iPEPS gauge-fixing paper: https://arxiv.org/abs/1503.05345
- ITensors docs: https://itensor.github.io/ITensors.jl/stable/
- ITensors ITensor type docs: https://itensor.github.io/ITensors.jl/stable/ITensorType.html
- Public docs test at that commit: https://raw.githubusercontent.com/jayren3996/SquarePXPDynamics.jl/2adff05717469fc81dc48cbc0382a90ad08a38bc/test/test_public_docs.jl
