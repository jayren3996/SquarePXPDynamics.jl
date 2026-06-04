# Native (PEPSKit) backend: D>1 validation status

Date: 2026-06-04. Closes the "native engine is unfalsifiable at D>1" gap flagged
in the migration review. Companion plan:
`docs/superpowers/plans/2026-06-04-native-d-gt-1-validation-oracle.md`.

## Why D=1 was insufficient

At `maxdim = 1` every star truncates back to a product state, so the native and
legacy engines are bit-for-bit identical **by construction** and the entire
truncation / SVD / leg-routing / `sqrt(λ)`-gauge machinery is bypassed. A
coordinate transpose/reflection, a center↔leaf leg swap (invisible because all
four neighbour projectors are the identical `P↓`, so the gate is permutation-
invariant while the bond routing is not), or a `√λ`-vs-`λ` double-count would all
have passed every pre-existing native test. The collapse-and-revival physics is a
D>1 (entanglement) phenomenon, so D=1 validates plumbing, not dynamics.

## What is NOW trusted (gated by `test/test_pepskit_native_dgt1.jl`)

- **Native exact finite oracle.** `exact_density_finite(::PXPIPEPSState)` (and
  `dense_state_finite`, `exact_one_site_expectation_finite`) contract ONLY
  `measurement_peps(state)` (convention B) on the finite torus, reusing the
  ED-validated legacy statevector readout. Validated **== the legacy exact oracle
  at D=2 to 1e-9** (T1, keystone). The legacy oracle is itself ED-validated
  short-time (`test_d_convergence.jl`), so this is an independent reference.
- **Native star-update routing/convention at D>1.** A single *lossless* native
  `project_star_pepskit!` (maxdim large, cutoff 0, rel_floor 0) on an **anisotropic
  D=2** cell (Lx≠Ly) matches the legacy `project_star!` exact density to **1e-9**
  (T3). Catches center↔leaf leg swaps invisible on C4-symmetric / D=1 states.
- **Convention B — `measurement_peps` is the source of truth.** Corrupting the
  `SUWeight` ledger leaves the exact density **unchanged to 1e-12**, while the
  mean-field `measure_simple` density **does** change (T2). Proves the exact
  contraction neither double-counts nor misses `λ` and does not depend on the
  ledger.
- **Coordinate map (no reflection/transpose).** Per-site native-vs-legacy density
  on a **vertically-asymmetric** checkerboard (even Ly) agrees site-by-site to
  **1e-9** (T4, all sites). A `row=Ly−y+1` ↔ `row=y` bug would mirror the profile.
- **Native evolve! pipeline tracks legacy.** Over a short D=2 trajectory the first
  step agrees to **<1e-5** (scheduling + single-step parity) and the whole
  trajectory stays within **2e-3** (the engines keep marginally different 2-D
  subspaces under truncation; measured drift ~2.5e-7 → ~5e-4 over t∈[0,0.4]) (T5,
  extended).
- **D=1 absolute correctness** (T0): all-down→0, x_plus→1/2, checkerboard→pattern
  mean, exact==simple.

## What is STILL NOT trusted (out of scope for this PR)

- **Native D>1 vs ED over the full revival trajectory.** T1/T3 compare native vs
  the *legacy exact oracle* (ED-validated short-time), not native vs ED across the
  whole `n(t)` revival [0,2.8] at D=2,3,4. The oracle now EXISTS to do this, but
  the headline D-ladder gate (`test_d_convergence.jl`) still builds legacy states.
  Wiring native into that trajectory-RMS gate is the next validation step.
- **D≥3 native dynamics.** Tests use D=2; the dense `2^N` native oracle is a
  small-cell/small-D tool, so D≥3 is not yet gated.
- **Native CTM observables at D>1** remain D=1-only, and the CTM convergence
  residual is still mislabeled (`truncation_error`, not `η`) — separate issue.
- **The `:boundary` oracle at D≥5** (native or legacy) is still unverified; the
  native oracle is `:dense` only (≤16 sites).
- **CTM-aware evolution (Stage C/D)** — not attempted here.

## Implementation note

The oracle lives in `src/PEPSKitObservables.jl` (extends the
`FiniteIPEPSObservables` generics for `PXPIPEPSState`); it reuses
`_local_expectation_from_state`/`_dense_index`/`_local_positions` so only the
native statevector *construction* is new code (an `ncon` over `measurement_peps`
with explicit E↔W / N↔S bond routing and `cell.reps`-ordered physical legs).
