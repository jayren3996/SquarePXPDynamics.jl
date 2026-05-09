# D>1 General Dense-Star Simple Update Design

## Purpose

Replace the existing rank-1 mean-field placeholder in `_apply_general_star_gate_simple_update!` with a real Simple Update path: proper SVD-based bond truncation, lambda spectrum updates, and reported discarded weight. The current placeholder runs (so D>1 evolution does not throw), but it approximates every star gate by a per-rep rank-1 single-site map computed from each rep's dominant local profile; it never updates bond dimensions or lambda spectra and is not a Simple Update. This is the blocker that prevents the existing ScarFinder loop, blockade diagnostics, and step-driver scaffolding from being exercised on a state whose entanglement structure can actually grow.

The bar for "done" is **internally validated against finite-cluster ED for short times plus invariants** — not a research-grade reproduction of published triangular-PEPS Simple Update benchmarks. Long-time accuracy and external benchmark reproduction are NTU's job.

## Scope

In scope:

- General (non-product, non-identity) dense 7-site star gates at `D >= 2` on `ThreeSiteUnitCell` triangular iPEPS.
- Lambda absorption / extraction in the Vidal-form symmetric square-root convention already shaped into `States.jl`.
- Sequential center-anchored peel-split decomposition with deterministic spoke order.
- `cutoff` and `maxdim` truncation, with `discarded_weight` reporting.
- Three layers of validation tests: kernel ED, torus integration, stabilizer benchmark.
- Relaxation of the `ScarFinder.scar_search` `D == 1` guard to `unitcell isa ThreeSiteUnitCell || maxdim == 1`.
- README status section update.

Out of scope:

- `OneSiteUnitCell` for non-product D>1 gates. The 7 star positions alias to one rep and the cluster is degenerate; this is a separate semantic question and stays as an explicit `ArgumentError` from `_apply_general_star_gate_simple_update!`.
- Larger unit cells (4-site, 7-site, etc.).
- Performance work beyond `D <= 4`. The cluster scaling is `2^7 * D^18`; D=2 is comfortable, D=4 is at the edge, D>=6 will be slow.
- PEPS-level energy expectation. Requires CTMRG or boundary MPS and is a separate item.
- Imaginary-time energy correction inside ScarFinder (depends on the energy expectation above).
- NTU. This design is the bridge to NTU, not its replacement.

## Architecture

The new code replaces the body of `_apply_general_star_gate_simple_update!` in `src/SimpleUpdate.jl` and introduces 3-4 file-private helpers. The dead helpers used by the rank-1 placeholder (`_dominant_site_vector`, `_product_projection_targets`, `_representative_target`, `_regularized_physical_map`, `_apply_physical_map!`, `_normalize_site_tensor!`, `_relative_residual`) are deleted along the way. No new files in `src/`. No new module. The public API (`apply_star_gate_simple_update!`, `projected_pxp_step!`, `run_projected_pxp!`, `scar_search`) does not change shape. The standalone `truncate_state!` in `States.jl` is orthogonal and stays — it is used by `scar_search` to project from `dynamics_maxdim` down to `scar_maxdim` between iterations and is unrelated to the in-step gate-application truncation we are adding here.

The dispatch order in `apply_star_gate_simple_update!` already does the right thing:

1. Identity gate detected by Frobenius distance to `I`: no-op (existing).
2. Site-product factorization via `_try_factorize_product_gate`: per-rep single-site apply (existing, works at any D).
3. D=1 dense product oracle when `_is_d1_state(state)` is true (existing).
4. Else: `_apply_general_star_gate_simple_update!` (this is what gets filled in).

Inside the new function, the flow is:

```text
absorb_lambda_into_star_tensors
  -> build_cluster_network
  -> contract_with_gate
  -> sequentially peel-split spokes 1..6 with truncation
  -> extract new lambdas, divide-out external sqrt(lambda) on each new spoke tensor
  -> write back tensors and lambdas, respecting opposite-bond shared-reference invariants
  -> assemble SimpleUpdateDiagnostics
```

## Data Flow

For `ThreeSiteUnitCell` and a star centered at `c` (rep `r0`):

- Spokes at directions `{1, 3, 5}` land on rep `r1`.
- Spokes at directions `{2, 4, 6}` land on rep `r2`.
- Each spoke has 1 internal bond (to center) + 2 internal bonds (to adjacent spokes around the hexagonal ring) + 3 external bonds (to non-star sites).
- The center has 6 internal bonds (one to each spoke), 0 external bonds within the star.
- Internal bonds in the cluster: 6 center-spoke + 6 spoke-spoke ring = 12.
- External bonds dangling out of the cluster: 6 spokes * 3 external = 18.

After contracting the gate `Geff` (an `ITensor` with 7 in-phys + 7 out-phys indices) into the cluster, the worst-case materialized tensor has 7 out-phys indices (each dim 2) and 18 external virtual indices (each dim D), total `2^7 * D^18` complex entries:

- D=2: ~33M entries, ~530 MB. Comfortable.
- D=3: ~50G entries. Not feasible to materialize in one block.
- D=4: ~9G entries. ~140 GB. Not feasible.

ITensors' chosen contraction order generally avoids materializing this worst case in one block, but it is the upper bound and we cap CI tests at `D <= 4`. If memory bites earlier than expected, the fallback is an interleaved peel-as-we-contract mode (added later only if needed).

## Algorithm

### Step 1: Lambda absorption

For each of the 7 star sites, multiply its tensor on every bond leg by `sqrt(lambda)` of that bond. This produces "absorbed" copies of the site tensors. The original `state.tensors` and `state.lambdas` are not yet mutated.

The symmetric square-root convention is already what `States.jl` shapes the `lambdas` dictionary for: each bond carries a single shared lambda vector between two sites, with `norm(lambda) == sqrt(length(lambda))` and entries nonnegative. Absorbing `sqrt(lambda)` on each end means the full `lambda` shows up exactly once when the bond is contracted.

### Step 2: Cluster contraction

Form a network of:

- 7 absorbed star site tensors connected by their 12 internal bonds.
- The gate `Geff` represented as an `ITensor` of shape `[2,2,2,2,2,2,2] x [2,2,2,2,2,2,2]` with the `in` indices contracted into the 7 phys legs and 7 fresh `out` indices left as the new physicals.

Contract over: all 12 internal virtual bonds and all 7 in-phys indices. The result is a single tensor with 7 out-phys legs and 18 external virtual legs. Use ITensors' default contraction order; do not force a single materialized intermediate.

### Step 3: Sequential peel-split

Choose a fixed deterministic spoke order: `1, 2, 3, 4, 5, 6` (radial direction order from `TRIANGULAR_DIRECTIONS`). For `d` in that order:

- Identify the indices on the spoke side of the cut: spoke `d`'s out-phys index plus spoke `d`'s 3 external virtual indices.
- All other indices are on the rest side.
- SVD the cluster tensor across that bipartition with the requested `cutoff` and `maxdim`. Use `ITensors.svd(cluster, spoke_inds; cutoff, maxdim)`.
- The left singular vectors `U` become the new (still-absorbed) site tensor for spoke `d`. They carry the spoke `d`'s out-phys index, its 3 external bonds (now still carrying `sqrt(lambda)` from absorption), and a new fresh bond index toward the center.
- The new center-spoke `d` bond carries the singular values `S` directly; this is the new `lambda` for that bond.
- The right singular vectors times `S` become the new "rest" cluster tensor for the next iteration. The next iteration's cluster has one fewer spoke phys index and one more "internal" bond index where `d`'s split was performed.
- Track per-step `discarded_weight_d = sum(sigma^2 for sigma cut by cutoff or maxdim) / sum(sigma^2 over all)`.

After 6 peels, the residual rest tensor has exactly 1 phys index (the center's new out-phys), 6 internal bond indices (one per spoke, with the freshly truncated dimensions), and 0 external indices. This is the new center site tensor (still in absorbed form).

### Step 4: Lambda extraction and writeback

For each new spoke site tensor:

- Divide out `sqrt(lambda_external)` on each of the 3 external bonds (un-absorb). The external bond lambdas were not modified by this update, so this restores the Vidal-form representation on those bonds.
- Do NOT divide out anything on the new center-spoke bond: that bond's new `lambda` is the singular vector `S` from the SVD.

For the new center site tensor: it has 6 internal bond legs (one per spoke). No external bonds, no further division.

Update `state.tensors[r0]`, `state.tensors[r1]`, `state.tensors[r2]` with the new tensors (not three separate writes per spoke; collect contributions and write per rep). For the 3-site UC where three spokes share the same rep (e.g., directions `{1, 3, 5}` all map to `r1`), the rep tensor must absorb all three spoke contributions consistently. The three spokes-of-rep-`r1` were treated as three independent legs in the cluster, but in the iPEPS they are all the same physical site. After the peel-split, the three new "spoke tensors at rep `r1`" must agree as single-site operators (translational invariance).

For projected PXP on the 3-site UC the three sibling tensors are equal up to numerical noise because the gate respects sublattice symmetry; the policy below is a guard, not the common path:

- Bring the three sibling tensors into a common index labeling (rename their fresh center-spoke bond indices to a canonical name `b_d` per direction, and reorder so the canonical leg order is `(phys, b_external_1, b_external_2, b_external_3, b_center)`).
- Compute the average tensor.
- Compute the residual: `max_i ||T_i - T_avg||_F / max(||T_avg||_F, 1)`.
- If the residual exceeds `1e-8`, raise an `ArgumentError` naming the rep, the residual magnitude, and a hint that the gate may not respect the unit cell's sublattice symmetry.
- Otherwise use the average. Same procedure for rep-`r2`.

This errors loudly on gates that don't fit the assumed symmetry, instead of silently averaging incompatible tensors.

For the 6 new center-spoke `lambda` vectors: normalize each to nonnegative, descending order, with `norm(lambda) == sqrt(length(lambda))`. Update `state.lambdas` for both `(r0, d)` and the opposite-bond entry on the spoke rep, sharing the same vector by reference (per the `States.jl` invariant verified in `test_states.jl`).

### Step 5: Diagnostics

Build and return:

```julia
SimpleUpdateDiagnostics(
    discarded_weight = sum of per-spoke discarded_weight_d,
    affected_bonds = [(r0, d) for d in 1:6],
    output_bond_dims = [dim(state.bond_inds[(r0, d)]) for d in 1:6],
)
```

`projected_pxp_step!` already wraps this into `ProjectedPXPStepDiagnostics` with per-layer breakdowns.

## Translational-Invariance Handling on `ThreeSiteUnitCell`

The 6 spokes of a single star straddle two reps (`r1` and `r2`), with 3 spokes per rep. After the peel-split, each spoke produces an independent new tensor. The three new tensors at rep `r1` must represent the same physical site in the translationally invariant iPEPS, so they must agree up to gauge.

Resolution policy is fully specified in Step 4 above. Summary: rename siblings into a canonical leg labeling, average them, and raise `ArgumentError` if the per-sibling residual exceeds `1e-8`. For projected PXP at the 3-site UC the symmetry holds exactly, so this branch should never trip in normal use.

## Validation Tests

Three layers, all in the existing `test/` tree.

### Layer 1: Kernel ED tests (fast, run every commit)

Located in `test/test_simple_update.jl`. New `@testset "general star simple update kernel"` block.

Setup helper: `_finite_cluster_vector(state, center)` reads the 7 star tensors, contracts them as if the cluster were finite (no wrap-around), and returns the 128-dim flattened state vector. This is well-defined for `ThreeSiteUnitCell` because the 7 star sites span 3 distinct reps and no internal bond aliases.

Tests:

- **Identity gate, D=2 random input.** Construct `random_ipeps(ThreeSiteUnitCell(), 2; seed=...)`. Capture per-tensor copies. Apply `apply_star_gate_simple_update!(state, I_128, c; maxdim=2)`. Assert tensors unchanged (within numeric tolerance), `discarded_weight == 0`, bond dims unchanged.
- **Random unitary gate, D=2 input, no truncation.** Build a Haar-random 128x128 unitary `U`. Apply via the general path with `maxdim` large enough that no truncation is forced (e.g., `maxdim = 64`). Compare `_finite_cluster_vector` before/after to direct `U * vector_before`. Assert match within `1e-10` of vector norm.
- **Random unitary gate, D=2 input, with truncation.** Same as above with `maxdim = 2`. Assert `discarded_weight > 0`, bond dims `<= 2`, and the truncated cluster vector still has norm close to 1.
- **Site-product gate via the general path.** Construct a known site-product gate but skip the factorization branch (call the private `_apply_general_star_gate_simple_update!` directly). Assert the result agrees with the existing site-product path within `1e-10`.
- **Hermiticity preservation under imaginary-time gate.** Build `Geff = exp(-tau * H_pxp_star)` with small `tau`. Apply to a real-valued D=2 input. Assert resulting tensors are real within tolerance and norms are finite.
- **Order-of-spokes invariance check.** Apply the same general gate twice from the same starting state, with peel order 1..6 and with peel order 6..1 (parametrized). Compare local Z expectations. Assert agreement within a documented tolerance (~1e-6) — this is a regression check on order bias, not a strict requirement.

### Layer 2: Torus integration test (runs once per CI)

Located in `test/test_evolution.jl`. New `@testset "projected PXP iPEPS matches 3x3 torus ED at short times"`.

Setup helper: `test/util_finite_ed.jl` (new file) provides:

```julia
build_3x3_torus_pxp_hamiltonian()  -> 512x512 dense Hermitian
build_3x3_torus_initial_all_down() -> 512-dim vector
torus_local_z_per_sublattice(vec)  -> NamedTuple{(:r0, :r1, :r2)}
```

The 3x3 torus has 9 sites = 3 unit cells along each lattice vector. Sublattices match `ThreeSiteUnitCell` exactly: 3 sites per sublattice. The torus PXP Hamiltonian is the sum over all 9 stars with periodic wrap, restricted to the blockade-allowed subspace via `P_blockade` factors per term.

Test:

- Initial state: all-down product on both sides.
- Trotter parameters: `dt = 0.01`, `maxdim = 4`, `order = :second`, 3 steps.
- iPEPS side: `run_projected_pxp!(state, 0.01, 3; order=:second, maxdim=4, cutoff=1e-12)`.
- Torus side: full ED `exp(-i * 3 * dt * 7 * H_torus / N_layers)` is wrong because Trotter splitting is what's being tested. Instead apply the same color-by-color projected gate sequence to the 512-dim vector. (`util_finite_ed.jl` provides a helper that mirrors the iPEPS schedule layer-by-layer.)
- Compare per-sublattice `<Z>` after each step. Tolerance: `1e-3` absolute for the integrated test.

### Layer 3: Stabilizer benchmark (runs once per CI, true 2D)

Located in `test/test_evolution.jl`. New `@testset "cluster stabilizer benchmark at D=4"`.

Initial state: all-`Z+` product, which in the repo convention is `:up` because `|up> = |0>` and `pauli_z() * |0> = +|0>`. Hamiltonian: `cluster_star_hamiltonian()` — all `K_i` commute so Trotter error is zero.

Run real-time evolution with `dt = 0.05`, `maxdim = 4`, second-order, 4 steps. Measure center `<Z_c>(t)` after each step. Compare to `cluster_center_z_expectation_exact(t)`. Tolerance: `1e-6`. Any drift IS truncation error, and the stabilizer evolution stays in a low-D manifold so D=4 should be more than enough for short times.

## API Changes

- `apply_star_gate_simple_update!`: no signature change. The fall-through general-gate branch now does proper Simple Update instead of the rank-1 placeholder. The `maxdim === nothing` check inside `_apply_general_star_gate_simple_update!` stays — `maxdim` remains required for general updates.
- `_apply_general_star_gate_simple_update!`: rejects `OneSiteUnitCell` with an explicit `ArgumentError("general star updates for OneSiteUnitCell are not yet supported")`. Removes the existing `current_dim <= maxdim` precondition because the new path can now grow bond dimensions up to `maxdim` and truncate back.
- `ScarFinder.scar_search`: no change required. The current code already has no D=1 guard. `dynamics_maxdim`/`scar_maxdim` separation stays as-is. After this change the workflow will actually exercise bond growth and truncation per-step rather than running rank-1 placeholders.
- `ScarFinder._seed_state`: no change. `random_ipeps(uc, D)` already supports D>1.
- `TriangularPEPSDynamics.jl`: no new exports.
- The dead helpers listed in the Architecture section are removed in the same PR; their callsites disappear when the rank-1 path body is replaced.

## Files To Modify

- `src/SimpleUpdate.jl` — replace the body of `_apply_general_star_gate_simple_update!`, add 3-4 file-private helpers (`_absorb_lambda_into_star_tensors`, `_build_cluster_network`, `_peel_split_spoke!`, `_extract_and_writeback!`), delete the rank-1 placeholder helpers (`_dominant_site_vector`, `_product_projection_targets`, `_representative_target`, `_regularized_physical_map`, `_apply_physical_map!`, `_normalize_site_tensor!`, `_relative_residual`).
- `src/ScarFinder.jl` — no change.
- `src/TriangularPEPSDynamics.jl` — no change.
- `test/test_simple_update.jl` — add Layer 1 tests under a new `@testset`. Update or remove any existing tests that pin the rank-1 placeholder behavior (residual ≈ ?, no bond growth, etc.) — those tests have to be reframed as "old behavior, replaced".
- `test/test_evolution.jl` — add Layer 2 and Layer 3 tests under new `@testset`s.
- `test/util_finite_ed.jl` — new file. Finite-torus ED helpers; not exported from the package.
- `test/runtests.jl` — `include("util_finite_ed.jl")` near the top so test files can use it.
- `README.md` — update `## Simple Update Status` and `## ScarFinder Status` to reflect that ThreeSiteUnitCell at D>1 now does proper Simple Update with bond truncation. Document the OneSiteUnitCell limitation explicitly.

## Risk Register

- **Cluster contraction memory at D=4.** Worst-case `2^7 * 4^18 ≈ 140 GB`. Mitigation: ITensors' contraction order generally avoids the worst case; cap CI at `D <= 4`; document the scaling honestly. If it bites earlier than expected, add an interleaved peel-as-we-contract path.
- **Order-dependence of peel sequence.** The fixed `1..6` order introduces a bias smaller than truncation error in practice. Mitigation: Layer 1 includes a regression test comparing two peel orders to keep an eye on the bias magnitude. If the bias ever exceeds the truncation error, add an option to symmetrize over orders.
- **Translational-invariance enforcement on rep-`r1`/`r2`.** The "average three new spoke tensors per rep" step assumes the gate respects sublattice symmetry. For projected PXP at the 3-site UC, it does. For arbitrary user-supplied gates it may not, in which case the explicit `ArgumentError` fires. Mitigation: documented in the error message.
- **Lambda normalization drift.** Repeated absorb/extract cycles can compound floating-point drift in `norm(lambda) == sqrt(length(lambda))`. Mitigation: re-normalize after each writeback; existing tests already check this invariant.
- **Lazy contraction surprises.** ITensors may pick a contraction order that materializes a large intermediate even when smaller orders exist. Mitigation: if observed, force a specific contraction sequence via explicit `contract` calls rather than relying on the network-level `*` operator.

## Acceptance Criteria

This PR is acceptable when:

- All existing tests still pass, after explicitly retiring or rewriting any test that was pinning the rank-1 placeholder's behavior (zero discarded weight, no bond growth, residual reporting).
- The rank-1 placeholder helpers listed in the Architecture section are deleted, not just unused.
- Layer 1 kernel ED tests pass at D=2 and D=4.
- Layer 2 torus integration test passes within `1e-3` absolute on per-sublattice `<Z>` at `dt=0.01`, 3 steps, `maxdim=4`.
- Layer 3 stabilizer benchmark passes within `1e-6` at `dt=0.05`, 4 steps, `maxdim=4`.
- `scar_search` on `ThreeSiteUnitCell` with `maxdim = 2` and a non-product seed produces ScarCandidate diagnostics whose `discarded_weight` is nonzero on at least one layer (proving the new path actually truncates) and whose bond dimensions evolve up to `maxdim` (proving the new path actually grows bonds).
- README accurately reflects the new boundary: ThreeSiteUnitCell at D>1 supported; OneSiteUnitCell at D>1 with non-product gates explicitly errors.
- No silent fallbacks, no placeholder updates, no test that only checks `isfinite`.

## Out-Of-Scope Follow-Ups

These are noted here so they don't drift away. None are blockers for this PR.

- `OneSiteUnitCell` aliasing semantics for non-product gates.
- PEPS energy expectation at D>1 (CTMRG or boundary MPS).
- Imaginary-time energy correction in ScarFinder.
- NTU backend with the same `apply_gate!` interface.
- 4-site, 7-site, larger unit cells.
- Symmetrization of the peel order as a built-in option.
- Performance work for D>=4 (peel-as-we-contract or structure-aware decomposition).
