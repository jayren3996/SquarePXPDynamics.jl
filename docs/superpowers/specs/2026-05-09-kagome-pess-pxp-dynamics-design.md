# Kagome PESS PXP Dynamics Design

## Purpose

Build a Julia tensor-network tooling layer for fixed-bond-dimension projected PXP evolution on the **2D Kagome lattice**, using the **PESS (Projected Entangled Simplex State)** ansatz with a **Simple Update** kernel, **as a prototype for the eventual NTU-based dynamics implementation**. This is the kagome counterpart of the now-parked triangular work; it pivots away from triangular because the 7-site star + standard iPEPS combination has a structural sublattice-aliasing problem in cluster-and-split Simple Update and infeasible cluster scaling above D ≈ 4.

**Honest framing of SU vs NTU.** The PESS *ansatz* is a representation choice and is dynamics-agnostic. The PESS *Simple Update kernel* (apply gate → HOSVD-truncate to fixed D) is the published canonical kagome ground-state algorithm; its track record for real-time dynamics is sparse. NTU (Dziarmaga 2021) is the production-grade dynamics algorithm but substantially more complex to implement. This spec implements PESS+SU as a prototype because:

1. ScarFinder's evolve-project loop projects to a low-D manifold every step; this is the regime where SU truncation is most defensible (the projection is what makes SU acceptable for dynamics here, not its general accuracy on out-of-equilibrium states).
2. The geometry, state container, models, gates, scheduler, observables, and ScarFinder-loop infrastructure are **algorithm-agnostic** and will be reused unchanged when NTU lands.
3. Validation Layer-1 (kernel ED) and Layer-3 (analytic stabilizer benchmark) test the SU kernel in regimes where it should be exact or near-exact; passing these is meaningful even if dynamics fidelity at scale is not yet established.

The bar for "done" on the SU prototype is **passing Layers 1 and 3, plus running the ScarFinder loop end-to-end at D ≥ 2 with real bond growth and truncation observed**. Layer-2 (finite-torus integration) is **demoted from acceptance gate to dynamics-fidelity indicator**: we run it and record the result; passing within tolerance means SU dynamics is acceptable for first scientific work; failure by orders of magnitude means we halt scientific claims and queue the NTU PR. There is no published kagome PXP dynamics work to cross-check against, so finite-cluster ED is the only oracle.

The design context, literature scan, and decision rationale live in `Notes/2026-05-09-kagome-pxp-pivot-plan.md`. The NTU follow-up PR scope is sketched in `Notes/2026-05-09-kagome-pess-ntu-followup.md`. This spec assumes the four locked decisions (PESS ansatz, 9-site UC, vanilla λ environment, D=4/8/12 targets) hold; the SU-vs-NTU framing here is a clarification of what the SU prototype proves and doesn't prove.

## Locked Design Decisions

1. **Ansatz: PESS** with one rank-3 simplex tensor `S` per kagome triangle (no physical leg) and one rank-3 site tensor `T` per kagome site (1 physical + 2 simplex bonds). Plain iPEPS on kagome is rejected because (a) it wastes representation capacity on triangle correlations and (b) PESS is the published canonical kagome ansatz with state-of-the-art ground-state benchmarks (Xie et al. 2014).

2. **Unit cell: 9-site enlarged.** The natural 3-site UC has the same sibling-aliasing problem we hit on triangular (multiple star positions wrap to the same rep). A 3×3 enlargement of the natural UC gives 9 distinct reps, and the 5-site PXP star then occupies 5 distinct reps with no aliasing. Cost: 3× more reps to track than the natural UC; per-step work scales accordingly. Sublattice-aliased small UCs stay out of scope for the first PR.

3. **Gauge: vanilla λ-environment SU first.** Each simplex bond carries a non-negative descending λ vector with `‖λ‖ = √dim(λ)`, the same convention already shipped in `States.jl`. BP-gauge maintenance (Tindall-Fishman 2023) is a follow-on PR layered on top once the SU kernel and validation tests are passing. This avoids debugging two layers at once.

4. **Bond dimension scope.** D = 4 in CI tests. D = 8 in benchmarks. D ≥ 12 documented as out-of-scope (server-class memory; future work). The PESS cluster scaling is `2^5 × D^4` for the post-gate cluster — far better than the `2^7 × D^18` triangular wall — but per-iteration HOSVD costs grow as `D^k` and become uncomfortable past D ≈ 8 on a workstation.

## Scope

In scope:

- Kagome lattice geometry: Bravais coordinates, sublattice partition (A, B, C), 9-site enlarged unit cell (3×3 of the natural up-triangle).
- PESS state container with site tensors (1 phys + 2 simplex bonds) and simplex tensors (3 simplex bonds, no phys), with vanilla Vidal-form λ on every simplex bond.
- Dense 32-dim spin-1/2 5-site PXP star Hamiltonian builder, blockade projector for the 5-site star (5 internal star edges to check), real-time and imaginary-time projected gates.
- 5-site projected gate application via cluster contraction through PESS, HOSVD-based decomposition back to PESS form, λ truncation per simplex bond.
- Schedule for parallel star-gate application on kagome (3-coloring suffices because kagome stars centered on the same sublattice are vertex-disjoint).
- Local observables and blockade-violation diagnostics.
- ScarFinder driver adapted from triangular: deterministic seeds, evolve-project iterations, candidate ranking.
- Validation: kernel ED tests on the 5-site cluster, finite-torus integration test against full ED, kagome stabilizer benchmark with closed-form `<Z_c>(t)`.

Out of scope (explicit):

- Plain iPEPS on kagome.
- 3-site (or any aliased) UC for non-product gates.
- BP-gauge maintenance, loop-corrected BP, cluster-BP corrections.
- NTU, Full Update, CTMRG.
- Real PEPS expectation values via boundary MPS or environment contraction.
- Imaginary-time energy correction inside ScarFinder (depends on the previous item).
- Fermionic, symmetric, or GPU tensor backends.
- Larger unit cells beyond the 9-site enlargement.
- Reviving the triangular code path. Triangular spec, plan, and partial implementation stay in the repo as historical artifacts.

## Architecture

The kagome work lives alongside (not replacing) the triangular code in the same root Julia project. New modules under `src/`:

- `KagomeGeometry.jl` — Bravais lattice, sublattice partition, 9-site UC, neighbor relations, color schedule.
- `KagomePESS.jl` — PESS state container, site/simplex tensor builders, λ bookkeeping, product/random initializers.
- `KagomeModels.jl` — 5-site PXP star Hamiltonian (32×32), 5-site blockade projector, kagome cluster stabilizer Hamiltonian for the analytic benchmark.
- `KagomeGates.jl` — dense and projected 5-site gates over the 32-dim cluster Hilbert space.
- `KagomeSchedules.jl` — color partition for parallel star centers; first-order and second-order Trotter schedules.
- `KagomeSimpleUpdate.jl` — gate application via cluster contraction + HOSVD decomposition back to PESS, λ truncation, diagnostics.
- `KagomeObservables.jl` — local single-site expectations, blockade violation diagnostics adapted to the PESS structure.
- `KagomeEvolution.jl` — schedule driver, projected PXP step helpers, run loop.
- `KagomeScarFinder.jl` — search loop, seeds, ranking. Mirrors the triangular `ScarFinder.jl` API but typed against `KagomePESS`.

Existing code that stays untouched:

- `SpinOps.jl`, `Schedules.jl` (the abstract first/second-order weights are reused, but kagome has its own concrete color list), `SolvableModels.jl` (extended with kagome variants).

The shared `TriangularPEPSDynamics.jl` module re-exports kagome symbols alongside triangular ones. Or (alternative) a new top-level module `KagomePXPDynamics.jl` exists as a sibling. **Lean: same module, namespaced public symbols.** The package name is currently `TriangularPEPSDynamics`; that name is acknowledged as inaccurate post-pivot but renaming is out of scope for this PR (would touch every test file and the Project.toml).

## Geometry

### Kagome Bravais lattice and sublattice partition

Bravais vectors (in lattice units): `a1 = (2, 0)`, `a2 = (1, √3)`. The natural 3-site unit cell contains:
- Sublattice A at offset `(0, 0)`.
- Sublattice B at offset `(1, 0)`.
- Sublattice C at offset `(1/2, √3/2)`.

Site coordinates are stored as `(n1, n2, sublat)` where `n1, n2 ∈ ℤ` are Bravais indices and `sublat ∈ {A, B, C}`. Equivalently a `KagomeCoord` struct with three integer fields. The sublattice index makes wrap_coord cheap — no quotient computation per site.

### Coordination and neighbor relations

Each kagome site has **coordination 4**. The 4 NNs of a site at `(n1, n2, A)` are:
- 2 NNs in the up-triangle containing it: `(n1, n2, B)` and `(n1, n2, C)`.
- 2 NNs in the down-triangle containing it: `(n1-1, n2, B)` and `(n1, n2-1, C)`.

The exact NN offsets for B and C centers (and the indexing convention for down-triangles) are pinned down in the implementation plan via explicit enumeration tests rather than restated here. A site belongs to exactly **2 triangles** (one up, one down). Each triangle contains 3 sites of distinct sublattices.

### 9-site enlarged unit cell

Reps are `(p, q, sublat)` with `p, q ∈ {0, 1, 2}` and `sublat ∈ {A, B, C}`. Wrapping: `wrap_coord(c) = (mod(c.n1, 3), mod(c.n2, 3), c.sublat)`.

In this UC, the 5-site star centered at any rep occupies 5 distinct reps (verified explicitly by enumeration in tests). This is the minimal UC where no aliasing occurs for a 5-site PXP star.

### Color schedule

Kagome admits a clean **3-coloring** of stars centered on each sublattice: stars centered on A-sites form one color, B-sites another, C-sites another. Stars of the same color are vertex-disjoint (their 5-site footprints don't overlap) within the 9-site UC. So one Trotter sweep is 3 layers.

First-order schedule: `[A, B, C]`. Second-order schedule: `[A, B, C, C, B, A]` with each layer at half the time step.

## Data Model

### KagomePESS struct

```julia
struct KagomePESS{UC<:KagomeUnitCell}
    unitcell::UC
    site_phys_inds::Dict{KagomeCoord, Index}              # phys index per rep
    site_simplex_inds::Dict{Tuple{KagomeCoord, Symbol}, Index}  # one per (rep, :up | :down)
    simplex_inds::Dict{KagomeTriangleCoord, NTuple{3, Index}}   # one Index per leg of each simplex
    site_tensors::Dict{KagomeCoord, ITensor}              # 3-leg: phys + 2 simplex
    simplex_tensors::Dict{KagomeTriangleCoord, ITensor}   # 3-leg: 3 simplex, no phys
    lambdas::Dict{KagomeBondKey, Vector{Float64}}         # one per (site, :up | :down) bond
end
```

`KagomeBondKey` identifies a site-simplex bond. Opposite ends share the same `Index` and the same λ vector by reference, mirroring the existing triangular `States.jl` invariant.

`KagomeTriangleCoord` identifies a triangle by its sublattice-A corner Bravais position plus an `:up` or `:down` orientation. A 9-site UC has 9 up-triangles and 9 down-triangles = 18 triangle reps total; each has its own simplex tensor.

### Product-state and random initializers

```julia
product_pess(uc, state_symbol; D = 1) -> KagomePESS
random_pess(uc, D; seed = nothing) -> KagomePESS
```

`product_pess(uc, :down; D=1)` builds a state where every site is `|down⟩` and every simplex tensor is the identity-like rank-3 tensor (entries 1 on the all-ones diagonal of bond indices, 0 elsewhere). λ vectors initialized to all-ones with the right normalization.

## Algorithm

### 5-site PXP star gate

```text
H_star = X_c * P_↓1 * P_↓2 * P_↓3 * P_↓4
```

Acting on the 5-site Hilbert space (32-dim) where site 1 is the center and sites 2-5 are the 4 NNs in some fixed direction order. The projected real-time gate is:

```text
G_eff = P_blockade^(5) · exp(-i dt H_star)
```

where `P_blockade^(5)` is the 32×32 diagonal projector onto blockade-allowed configurations on the 5-site star. Blockade-forbidden configurations: any of the 4 center-NN bonds (4 edges) plus the 2 within-triangle NN-NN bonds (one per triangle the star spans) where both endpoints are `|up⟩`. Total: 6 edges to check per 5-site configuration.

### Star-gate application via PESS cluster contraction + HOSVD decomposition

Given a kagome PESS state and a star center `c`:

**Step 1: Identify the cluster.** The 5-site star at `c` plus the 2 simplex tensors of the triangles containing `c`. With the 9-site UC, all 5 sites and 2 simplices map to distinct reps.

**Step 2: Lambda absorption.** For each of the 5 site tensors, multiply on each of its 2 simplex bond legs by `√λ`. The 4 NNs each have one bond inside the cluster and one outside; both get `√λ` absorbed (the outside one gets divided out at writeback).

**Step 3: Cluster contraction with the gate.**
- Build the 7-tensor network: 5 absorbed site tensors + 2 simplex tensors, connected by 6 internal site-simplex bonds.
- Apply the dense 5-site gate `G_eff` (a 32×32 ITensor with 5 in-phys + 5 out-phys indices) to the 5 phys legs.
- Result: a single ITensor with 5 fresh out-phys legs (dim 2) + 4 external simplex bond legs (dim D each, one per NN's outgoing bond).
- Worst-case cluster size: `2^5 × D^4 = 32 × D^4` complex floats. D=4 → 8K entries. D=8 → 130K entries. D=16 → 2M entries. Trivially in-memory through D = 16 and beyond.

**Step 4: HOSVD-based decomposition back to PESS form.**

The decomposition target is: 5 new site tensors (each with phys + 2 simplex bonds, with the cluster-internal simplex bonds carrying truncated dimensions) + 2 new simplex tensors (each with 3 simplex bonds).

Algorithm:

1. **HOSVD pass over the 6 internal site-simplex bonds.** For each internal bond, compute the mode unfolding of the cluster across that bond, SVD-truncate to `maxdim` with `cutoff`, take the singular values as the new λ for that bond. The left/right singular vectors define the new "absorbed" site and simplex tensors on either end of that bond.
2. **Group truncated factors per simplex tensor.** Each simplex `S` has 3 internal bonds; its 3 truncated mode-unfoldings together define the new `S` via HOSVD recombination (Tucker-style).
3. **Group truncated factors per site tensor.** Each site has 2 simplex bonds (1 internal + 1 external for NNs; 2 internal for the center). The internal-bond mode unfolding gives the truncated bond direction; the external bond stays at its original dimension and doesn't participate in this gate's truncation.

**Step 5: External λ extraction (un-absorbing).**

For each NN site, divide out `√λ_external` on its 1 external simplex bond. The center has no external simplex bond. The 2 simplex tensors have no external legs.

**Step 6: Writeback.**

Write back the 5 new site tensors, 2 new simplex tensors, and 6 new λ spectra. Update `state.bond_inds` for affected internal bonds (fresh `Index` objects of the new truncated dimensions). Lambda dictionaries update both ends of each affected bond, sharing the same vector by reference.

**Step 7: Diagnostics.**

```julia
SimpleUpdateDiagnostics(
    discarded_weight = sum over 6 internal bonds of (relative truncated-singular-value-mass),
    affected_bonds = list of 6 internal bond keys,
    output_bond_dims = list of 6 post-truncation dims,
)
```

### Unique-decomposition caveat

HOSVD of a 7-tensor network is gauge-non-unique: the simplex tensors have a residual gauge group acting on each bond. We fix the gauge by demanding that each new λ be sorted descending with `‖λ‖ = √dim(λ)` (matching the existing convention) and that singular vectors be chosen with a fixed sign convention (e.g., first nonzero entry positive). This makes the algorithm deterministic.

### Translational invariance with the 9-site UC

The 9-site UC was chosen specifically so that 5-star positions are 5 distinct reps. So the writeback is direct: each new tensor goes to its rep. **No sibling merging required.** This was the structural fix for the triangular failure.

For consistency, when a star is applied at a different center within the same color layer, the writeback updates a different set of reps. After all 3 color layers (one Trotter sweep), every rep has been touched exactly the same number of times by translational invariance.

## Validation Tests

Three layers, mirroring the triangular plan structure.

### Layer 1: Kernel ED tests on the 5-site cluster (fast, run every commit)

In `test/test_kagome_simple_update.jl`. Uses helper `cluster_vector_from_pess(state, center)` that returns the 32-dim vector by contracting the 5 site tensors and 2 simplex tensors of the star, with external simplex bonds traced against all-ones environments. Defined in `test/util_kagome_finite_ed.jl`.

Tests:

- **Identity gate, D=2 random PESS input**: tensors unchanged, `discarded_weight == 0`, bond dims unchanged.
- **Random unitary gate, D=2 input, no truncation (`maxdim = 64`)**: `cluster_vector_after ≈ G * cluster_vector_before` up to global phase, within `1e-10`.
- **Random unitary gate, D=2 input, with truncation (`maxdim = 2`)**: `discarded_weight > 0`, bond dims `<= 2`, finite norms.
- **Imaginary-time projected PXP gate, D=2 random input**: tensors stay finite, hermiticity-preserving (real input → real output within tolerance).
- **Site-product 5-site gate via the general path**: agrees with a direct site-by-site application within `1e-10`.

### Layer 2: Finite-torus integration test (indicator, not gate)

**Status: indicator-only.** Demoted from acceptance gate per the SU-prototype framing in Purpose. We run it every CI cycle, record the result, and use the result to interpret whether SU dynamics is acceptable for scientific use. Pass within tolerance → continue with first scarfinder runs; fail by orders of magnitude → halt scientific claims and queue the NTU PR before publishing any results.

In `test/test_kagome_evolution.jl`. Uses helpers `build_kagome_3x3_torus_pxp_hamiltonian()` and `kagome_torus_local_z_per_sublattice(vec)` in `test/util_kagome_finite_ed.jl`.

Setup: a small kagome torus matched to the natural 3-site UC (12-site = 2×2 supercell, Hilbert space `2^12 = 4096`, trivial ED). The iPEPS uses the 9-site enlarged UC, which represents a state with translation invariance on a coarser period than the torus's; this is acceptable for a short-time reference because the iPEPS is *more* translation-invariant than the torus, not less. The torus is the comparison oracle, not a faithful representation of the iPEPS state.

If the 12-site torus turns out not to wrap consistently with the kagome lattice geometry chosen, the implementation plan falls back to a 9-site or 18-site torus. The validation criterion is geometry-agnostic: per-sublattice `<Z>` from the iPEPS matches the per-sublattice `<Z>` from the torus ED to within `1e-3` absolute over the first 3 Trotter steps at `dt = 0.01, maxdim = 4`.

Test: `dt = 0.01`, `maxdim = 4`, second-order, 3 steps. Compare per-sublattice `<Z>` to torus full-ED. Tolerance `1e-3` absolute.

### Layer 3: Kagome stabilizer benchmark (per CI run, true 2D)

Closed-form solvable: `K_v = X_v ∏_{u ∈ NN(v)} Z_u`. All `K_v` commute → no Trotter error. Initial state `:up`. Closed form `<Z_c>(t)` can be derived analogously to the triangular cluster Hamiltonian.

Run real-time evolution, `dt = 0.05`, `maxdim = 4`, 4 steps. Compare to analytic. Tolerance `1e-6`.

## API Surface

Public:

```julia
KagomeCoord, KagomeTriangleCoord, KagomeBondKey
NineSiteKagomeUC                         # the only supported UC for now
unit_cell_representatives(uc)
neighbor(c, d)                           # d ∈ 1:4
star_sites(c)                            # returns 5 KagomeCoord
star_color(c)                            # returns 1, 2, or 3 (sublattice)
disjoint_stars(a, b)                     # bool

KagomePESS, product_pess, random_pess
site_tensor, simplex_tensor, phys_index, simplex_bond_index, simplex_bond_lambda

pxp_kagome_star_hamiltonian(projector = projector_down(), flip = pauli_x())  # 32×32
kagome_blockade_projector()              # 32×32
kagome_cluster_star_hamiltonian()        # 32×32

dense_kagome_gate, projected_kagome_gate
first_order_kagome_colors, second_order_kagome_colors, kagome_schedule_layers
kagome_color_canonical_center(color)

apply_star_gate_simple_update_pess!(state, gate, center; cutoff, maxdim) -> SimpleUpdateDiagnostics
projected_pxp_step_kagome!(state, dt; order, maxdim, cutoff, evolution) -> ProjectedPXPStepDiagnostics
imaginary_projected_pxp_step_kagome!(state, dτ; order, maxdim, cutoff)
run_projected_pxp_kagome!(state, dt, nsteps; order, maxdim, cutoff, evolution)

local_expectation_kagome(state, c, op)   # single-site, exact for D=1, proxy for D>1
mean_blockade_violation_kagome(state, centers)
dense_star_blockade_violation_kagome(state, center)

KagomeScarFinderConfig, KagomeScarCandidate
scar_search_kagome(config; seed)
rank_candidates_kagome(candidates)
```

Private kernels (same file-private pattern as the triangular code):

```julia
_absorb_lambda_into_kagome_star_tensors(state, center)
_build_kagome_cluster_with_gate(absorbed_sites, absorbed_simplices, gate)
_hosvd_decompose_kagome_cluster(cluster, ...; cutoff, maxdim)
_extract_and_writeback_kagome!(state, center, new_sites, new_simplices, new_lambdas, ...)
```

The `apply_star_gate_simple_update_pess!` function dispatches:

1. Identity gate (numerical detection) → no-op.
2. Site-product gate (via the existing `_try_factorize_product_gate` from triangular code, generalized to 5 sites) → per-rep single-site apply.
3. D=1 dense product oracle (analogous to triangular path, on 5 sites instead of 7) → tested exact path.
4. General path → cluster contraction + HOSVD decomposition.

## Files To Create

In `src/`:

- `src/KagomeGeometry.jl`
- `src/KagomePESS.jl`
- `src/KagomeModels.jl`
- `src/KagomeGates.jl`
- `src/KagomeSchedules.jl`
- `src/KagomeSimpleUpdate.jl`
- `src/KagomeObservables.jl`
- `src/KagomeEvolution.jl`
- `src/KagomeScarFinder.jl`

Modify:

- `src/TriangularPEPSDynamics.jl` — `include` and re-export the new modules. Rename or repurpose the package later (out of scope).
- `src/SolvableModels.jl` — add `kagome_cluster_center_z_expectation_exact(t)`.

In `test/`:

- `test/test_kagome_geometry.jl`
- `test/test_kagome_pess.jl`
- `test/test_kagome_models.jl`
- `test/test_kagome_gates.jl`
- `test/test_kagome_schedules.jl`
- `test/test_kagome_simple_update.jl`
- `test/test_kagome_observables.jl`
- `test/test_kagome_evolution.jl`
- `test/test_kagome_scar_finder.jl`
- `test/util_kagome_finite_ed.jl`
- `test/runtests.jl` — include the new test files and the helper.

In `Notes/`:

- README updates and a short `2026-05-09-kagome-pivot-status.md` once the work lands.

## Risk Register

- **Cluster scaling at large D.** `2^5 × D^4` is comfortable through D = 16 in raw cluster size, but per-iteration HOSVD costs grow as `D^k` for some k > 4 depending on contraction order. Realistically D = 8 is comfortable on a workstation; D = 12 needs care; D ≥ 16 is out of scope for the first PR.
- **HOSVD gauge non-uniqueness.** Could introduce numerical drift across many gate applications. Mitigation: enforce sign conventions and norm normalization on every λ extraction; spot-check with the kernel ED tests after long runs.
- **No published kagome PXP dynamics reference.** Validation depends entirely on Layer 1 (kernel ED), Layer 2 (small-torus ED), and Layer 3 (analytic stabilizer). If physics intuition disagrees with results, no peer-validated numbers to cross-reference. Mitigation: be conservative with claims; lean on the analytic stabilizer benchmark as a hard accuracy floor.
- **9-site UC correctness.** The claim "9-site enlargement makes 5-star positions 5 distinct reps" must be verified by enumeration in tests, not just asserted. Layer-1 setup includes this enumeration test.
- **Cluster reconstruction when site has external dangling bonds.** The "trace external bonds with all-ones" recipe used by the cluster-vector helper must be consistent before and after the gate. Same recipe, same external bond identities (which don't change during this gate's update). Mitigation: tested explicitly with identity-gate round-trip in Layer 1.
- **Triangulation of the gate.** The 5-site PXP gate has 6 internal blockade edges (4 center-NN + 2 NN-NN within each triangle the star spans). Getting the edge set wrong means the projected gate is wrong. Mitigation: test the projector against hand-computed reference matrix elements for a few specific 5-site configurations.
- **Reusable triangular helpers may not generalize cleanly.** `_try_factorize_product_gate` was written for `n=7`; generalizing to `n=5` is straightforward but needs its own test. The D=1 oracle similarly.

## Acceptance Criteria

Hard gates (must pass for the PR to land):

- All existing triangular tests still pass (no regression in the triangular code path).
- Kagome geometry tests pass: 9-site UC has 5 distinct reps per star; 3-coloring partitions are vertex-disjoint.
- Layer 1 kernel ED tests pass at D = 2 and D = 4.
- Layer 3 stabilizer benchmark passes within `1e-6` at D = 4.
- `scar_search_kagome` on `NineSiteKagomeUC` with `maxdim = 4` produces a ranked candidate list with finite scores, finite diagnostics, nonzero discarded weight on at least one layer (proving real truncation activates), and bond dims that grow during evolution (proving real bond growth).
- README has a new "Kagome status" section honestly documenting that kagome PESS PXP dynamics is implemented as an SU prototype, with NTU as the planned next step, and Layer-2 dynamics fidelity recorded but not used for scientific claims.
- No silent fallbacks, no placeholder updates, no test that only checks `isfinite`.

Indicator (recorded, not gating):

- Layer 2 torus integration test result. Target: per-sublattice `<Z>` within `1e-3` absolute. Outcome interpretation:
  - Pass within tolerance → SU dynamics is acceptable for first scientific scarfinder runs; the NTU follow-up PR can be queued at lower priority.
  - Fail by a factor of 2-10× tolerance → SU dynamics is borderline; explicitly disclaim scientific accuracy in any paper draft and prioritize the NTU PR.
  - Fail by orders of magnitude → SU dynamics is unfit for this problem; halt scientific use of the SU implementation and dispatch the NTU PR immediately.

## Out-Of-Scope Follow-Ups

- **NTU on kagome PESS** — the planned immediate follow-on PR. Sketch in `Notes/2026-05-09-kagome-pess-ntu-followup.md`. Replaces the Simple Update kernel with a neighborhood-tensor-update kernel; reuses all infrastructure built in this PR. Priority depends on Layer-2 indicator outcome.
- BP-gauge maintenance (Tindall-Fishman 2023). Vanilla λ environment first; orthogonal layer that can stack with either SU or NTU.
- Loop-corrected / cluster-BP corrections for triangle 3-cycles.
- 3-site or other aliased UCs, with sibling-merge semantics.
- Full Update / CTMRG.
- PEPS expectation values via boundary MPS or environment contraction.
- Imaginary-time energy correction in ScarFinder.
- D > 8 performance work.
- Triangular code path revival.
- Package rename from `TriangularPEPSDynamics` to something kagome-aware.
