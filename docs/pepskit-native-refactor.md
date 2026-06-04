# PEPSKit-native backend refactor

Status: **in progress** (incremental). The legacy ITensors Γ–λ backend remains
the production path; the PEPSKit-native backend is being introduced in parallel
and proven on small cases before any public API is migrated.

## Why PEPSKit owns the tensor-network backend

The repo historically carried a hand-rolled ITensors Γ–λ iPEPS state
(`SquareIPEPSState`), custom link weights, a custom five-site square-star simple
update, and used PEPSKit only as a CTMRG measurement adapter
(`to_pepskit_infinitepeps`). Maintaining a second tensor engine is expensive and
duplicates functionality PEPSKit already provides robustly:

- `InfinitePEPS` — the iPEPS state container.
- `SUWeight` — Schmidt bond weights for simple/cluster update (the λ ledger).
- `CTMRGEnv` + `leading_boundary` — boundary environments.
- `LocalOperator` + `expectation_value` — observables, including multi-site
  cross terms.
- `correlator` / `correlation_length`.
- `absorb_weight` — Vidal-gauge weight absorption/removal.

The target is: PEPSKit owns the backend; this repo owns only PXP-specific
physics (the star Hamiltonian/gate, the five-colour schedule, PXP observables,
ED/CTM validation, ScarFinder).

## Why the PXP star update stays custom

The PXP star term is a genuine **five-site cross**:

```
    P_right P_up P_left P_down X_center
```

PEPSKit's registered `trotterize`/`time_evolve` simple-update path supports
1-site, nearest-neighbour 2-site, and (via `simpleupdate3site`) L-shaped 3-site
clusters only — `is_nearest_neighbour(H)` and the 3-site MPO machinery cannot
represent a cross without an arbitrary snake-path embedding that biases the
truncation. We therefore keep a custom five-site star simple update, but rewrite
it to operate on PEPSKit `InfinitePEPS + SUWeight` and to reuse PEPSKit's
`absorb_weight` and gauge conventions, rather than on the bespoke ITensors
state. A snake-path `LocalCircuit`/MPO backend is deferred to an explicitly
experimental milestone (M7) and must not replace the star update unless
benchmarks show it is accurate and not path-biased.

## Milestone 0 audit (this environment)

- `Project.toml` pins `PEPSKit = "0.7.0"`, `TensorKit = "0.15.3"`, `julia = "1.12"`.
  These resolve together in the committed `Manifest.toml` (MPSKit `0.13.8`).
- Public/exported PEPSKit API used by the native backend:
  `InfinitePEPS`, `SUWeight`, `CTMRGEnv`, `leading_boundary`, `LocalOperator`,
  `expectation_value`, `correlator`, `correlation_length`, `absorb_weight`,
  `SimpleUpdate`.
- Internal (unexported) but stable helpers referenced via `PEPSKit.`:
  `truncrank` (already used by the legacy adapter). We avoid copying PEPSKit
  internals; where an unexported name is needed it is accessed through the
  `PEPSKit.` namespace and documented here.
- `InfiniteWeightPEPS` does **not** exist in 0.7.0, so this repo defines its own
  lightweight `PXPIPEPSState` bundling `InfinitePEPS + SUWeight` plus
  bookkeeping refs.

### PEPSKit structural facts (verified by probe)

- `InfinitePEPS` has one field `A::Matrix{TensorMap}`; `size(peps) == (Nr, Nc) ==
  (Ly, Lx)`; element `peps.A[row, col]` is a `(1,4)` `TensorMap`
  `P ← N ⊗ E ⊗ S ⊗ W`.
- `SUWeight` has one field `data::Array{DiagonalTensorMap,3}` of shape
  `(2, Nr, Nc)`. `weights[1, r, c]` is the **east** (x) bond of `T[r,c]`;
  `weights[2, r, c]` is the **north** (y) bond of `T[r,c]`.
- `absorb_weight(t, weights, row, col, ax; inv)` with `ax ∈ 1:4 = (N,E,S,W)`
  absorbs `sqrt(weight)` (or its inverse) onto that axis, reading the correct
  SUWeight entry for the bond (N→`[2,r,c]`, E→`[1,r,c]`, S→`[2,next(r),c]`,
  W→`[1,r,prev(c)]`).

## Known cost: 5-site CTM expectation compile time

PEPSKit's generic `expectation_value` for a multi-site `LocalOperator`
contracts the operator with the CTM environment via `@tensor`, whose
contraction-order search (`optimaltree`) is expensive *at compile time* for the
five-site cross network (~a dozen tensors). The first such call costs minutes;
subsequent calls reuse the compiled method. 1-site (density) and 2-site
(blockade) expectations are cheap. Consequence: PXP energy density via the
per-center 5-site sum is correct but compile-heavy, so the routine test suite
exercises the star update through 1-/2-site CTM observables and keeps the
5-site energy/star parity in `scripts/dev/verify_*.jl`. A consolidation
follow-up could compute energy from the simple-update reduced density matrices
(no CTM) for the fast path, reserving the 5-site CTM term for validation. The
legacy adapter shares this cost (same PEPSKit code path), so it is not a
regression introduced here.

## Coordinate convention

`SquareCoord(x, y)` (custom geometry, `x` = column, `y` = row increasing
upward) maps to PEPSKit matrix coordinates `CartesianIndex(row, col)` by

```
    row = Ly - y + 1,   col = x
```

so the custom `:up` neighbour is PEPSKit **north** (`row - 1`), `:right` is
PEPSKit **east** (`col + 1`), `:down` is south (`row + 1`), `:left` is west
(`col - 1`). This is identical to the legacy `to_pepskit_infinitepeps` adapter,
so the two backends share one convention and parity tests are meaningful.

## Physical basis convention

Preserved from the legacy backend:

```
    basis 1 = Rydberg / up / excited
    basis 2 = down / vacuum
```

The density operator is `n = |excited⟩⟨excited| = [1 0; 0 0]`. An all-down
product state has density 0; all-up has density 1.

## Virtual-space convention

Physical space `ComplexSpace(2)`; virtual space `ComplexSpace(D)`. PEPS tensor
domain order is `N ⊗ E ⊗ S ⊗ W` with **N, E non-dual** and **S, W dual**
(`ComplexSpace(D)'`), matching the legacy adapter and PEPSKit's neighbour
contraction convention.

## Simple-update gauge convention

The backend uses PEPSKit's native simple-update gauge ("convention B"): every
stored `InfinitePEPS` Γ tensor carries `sqrt(λ)` on each bond, so two
neighbouring tensors reproduce the full bond weight λ and the `InfinitePEPS` is
the measurement state directly (the `SUWeight` ledger is used only by the
update, not for measurement — `measurement_peps(state) === state.peps`). One
star update:

1. absorbs `sqrt(λ)` of each leaf's three *environment* bonds and normalises;
2. QR-reduces each leaf (3 environment legs → isometry X, reduced core →
   physical + center bond + internal q);
3. contracts the center with the four reduced cores and applies the five-site
   star gate on the physical legs;
4. sequentially `tsvd`-splits into center + four leaf cores, truncating each new
   center-leaf bond (`truncdim(maxdim) & truncbelow(cutoff)`, plus a relative
   condition floor `rel_floor`), and absorbs `sqrt(s)` symmetrically
   (`PEPSKit.absorb_s`, with `PEPSKit.flip_svd` to keep the stored weight
   non-dual as `SUWeight` requires);
5. reattaches the QR isometries, removes the environment weights, and writes the
   new `InfinitePEPS` tensors and `SUWeight` bonds.

A fully blockade-forbidden projected star annihilates the patch; for D=1 this is
the same physics as the legacy `project_star!`. Parity is verified: one D=1
short-time projected real-time star step reproduces the legacy ITensors
backend's CTMRG density to < 1e-6, and matches the analytic `sin²(dt)` single-
site rotation.

## What has been migrated

- `src/PEPSKitBackend.jl` — `PXPIPEPSState` wrapper + product/checkerboard
  constructors, `copy_state`, `state_version`, `log_norm`, `mark_mutated!`,
  `add_log_norm!`, coordinate maps, `measurement_peps`, and the shared
  dense-operator→`TensorMap` helpers (`_dense_index`,
  `_dense_operator_tensor_map`) imported by the native gate and observables. (M1)
- `src/PEPSKitStarUpdate.jl` — PEPSKit-native five-site star gate
  (`pxp_star_hamiltonian_matrix`, `pxp_star_gate_tensormap`) and the
  `project_star_pepskit!` simple update (the simple-update algorithm, formerly
  the separate `pepskit_star_update_impl.jl`, is now inlined here). (M3, M4)
- `src/PEPSKitEvolution.jl` — `evolve_pepskit!` Trotter driver preserving the
  five-colour palindromic schedule (`trotter_sequence_pk`); it calls
  `project_star_pepskit!` per star and deliberately does NOT use PEPSKit's
  generic `time_evolve`. (M5)
- `src/PEPSKitObservables.jl` — both native observable layers in one module:
  PXP `LocalOperator`s with CTMRG-backed measurements (M2), and the simple/local
  observables with the backend-agnostic `measure_simple`/`evolve!`/
  `reverse_evolve!` adapters + the `pxpipeps_from_square_ipeps` converter, so
  ScarFinder runs on either engine by dispatch (M6). (Merged from the former
  `PEPSKitPXPObservables.jl` + `PEPSKitSimpleObservables.jl`.)
- `src/PEPSKitStarSnake.jl` — experimental exact MPO factorization of the
  five-site star gate (`StarGateMPO`, `star_gate_mpo`, `reconstruct_star_gate`,
  `star_mpo_bond_dims`). Not wired into production. (M7)
- Tests: `test/test_pepskit_native.jl` (constructors, observables/parity, gate,
  star-update parity/reversibility/maxdim, schedule) and
  `test/test_pepskit_m6.jl` (native measure_simple parity at D = 1 and D > 1,
  evolve! parity, end-to-end ScarFinder parity, audit tensor-engine metadata).

## What remains legacy / next steps

- `SquareIPEPSState`, `StarSimpleUpdate`, `IPEPSEvolution`, and the
  `to_pepskit_infinitepeps` measurement adapter remain in place and are the
  production path. They are NOT removed in this PR (per the incremental
  constraint).
- **M6 (in progress — increment 1 landed):** a tensor-engine selector
  (`:legacy_itensors` | `:pepskit_simple`) threaded through ScarFinder and the
  audit.
  - `src/PEPSKitObservables.jl` (simple-observable layer) bridges the
    PEPSKit-native state onto the
    shared generics: `measure_simple(::PXPIPEPSState)` (local density / blockade
    / five-site energy / bond entropy, computed on the convention-"B" patches
    with the squared-Schmidt cut-leg environment, no CTMRG) and
    `evolve!(::PXPIPEPSState; params, protocol)::EvolutionLog` (wraps
    `evolve_pepskit!`). A `pxpipeps_from_square_ipeps` converter reuses the
    measurement adapter so native/legacy simple observables can be compared in
    the same Vidal gauge for `D > 1`.
  - `scarfinder!` and its state-changing helpers now dispatch on
    `Union{SquareIPEPSState, PXPIPEPSState}`, so the *same* orchestration runs on
    either engine. The PEPSKit modules are now `include`d before `ScarFinder.jl`.
  - `ScarFinderAuditConfig` gains a `tensor_engine` field (recorded in every CSV
    row and the JSON `config`), and `_build_state` routes to the native
    constructors. Parity is proven end-to-end: `scarfinder!` on the legacy and
    native engines (D = 1, `SimpleBackend`) agree per-iteration on observables
    and score (`test/test_pepskit_m6.jl`).
  - *Staged for increment 2:* trusted-CTM rows, per-iteration compression, the
    imaginary-time energy correction, and JLD2 snapshots are not yet supported on
    `:pepskit_simple` (the audit constructor and the relevant methods reject
    these loudly). These need native CTM-trust, native compression, and native
    snapshot support respectively.
- **M7 (experimental — exact MPO landed):** `src/PEPSKitStarSnake.jl` provides
  `StarGateMPO` / `star_gate_mpo(dt; projected, evolution)` — an *exact* matrix-
  product-operator factorization of the five-site square-star gate along the
  star-site ordering `(center, right, up, left, down)`, built by successive SVD.
  Only exactly-zero singular values are dropped, so `reconstruct_star_gate`
  reproduces `pxp_star_gate_dense` to machine precision while
  `star_mpo_bond_dims` exposes the genuine operator Schmidt rank: ≤ 2 across
  every cut (since `H = X_c⊗P^4` obeys `(X_c⊗P^4)² = I_c⊗P^4`, the gate is
  `I + (cosθ-1)P^4 - i sinθ X_cP^4` — two operator terms per star-site cut).
  Tests in `test/test_pepskit_star_snake.jl`. This module is *not* wired into
  production evolution: the genuine five-site cross gate is still applied by
  `project_star_pepskit!`. Applying the MPO onto an `InfinitePEPS` along a
  spatial snake (the path-biased research question) is deliberately deferred and
  must only be attempted — and only kept — once the native star update is fully
  production-validated and benchmarks show it is accurate and not path-biased.
- Before migrating public APIs onto the native backend, it must additionally
  pass: ED short-time comparison, D-convergence smoke test, CTM-χ sensitivity,
  multi-step reversibility audit, and snapshot/load — matching the legacy
  validation suite. Progress:
  - **ED short-time comparison — done.** `PXPValidationConfig` gained a
    `tensor_engine` selector (`:legacy_itensors` | `:pepskit_simple`), so
    `validate_pxp_ed_ipeps` runs the same finite-ED-versus-iPEPS trajectory on
    the native backend through the shared `evolve!`/`measure_simple` generics.
    At `D = 1` the native trajectory matches the ED-validated legacy engine
    bit-for-bit and is exact against ED at `t = 0`
    (`test/test_pxp_validation.jl`, "runs on the :pepskit_simple engine"). The
    native engine rejects `exact_finite_observables` and CTM trust sweeps (still
    legacy-only).
  - **Multi-step reversibility audit — done (small D).** `evolve_pepskit!`
    gained a `reverse` flag (apply each forward step's substeps in inverse order
    with negated layer times) and `reverse_evolve!(::PXPIPEPSState, …)` extends
    the shared generic, so `validate_pxp_reversibility` runs on the native state.
    From all-down at `maxdim = 1` the single-layer round trip returns to machine
    precision (drift `~1e-11`); multi-step at `maxdim = 2` with no truncation
    drifts only `~1e-6`, the simple-update integrator's intrinsic non-exact
    reversibility (`test/test_pxp_validation.jl`). NOTE: the simple-observable
    star-energy patch is a *dense* contraction (`~4^12` external legs at
    `D = 4`), so `measure_simple` — and hence this audit — is only feasible at
    `D ≲ 2–3`; large-`D` reversibility will need a CTM/MPO-environment energy.
  - **Native exact finite oracle + D>1 traps — done (small cell/D).**
    `exact_density_finite(::PXPIPEPSState)`, `dense_state_finite(::PXPIPEPSState)`,
    and `exact_one_site_expectation_finite(::PXPIPEPSState, …)` extend the
    `FiniteIPEPSObservables` generics for the native state (`src/PEPSKitObservables.jl`),
    contracting ONLY `measurement_peps` (convention B) — the result is independent
    of the `SUWeight` ledger. This makes the native backend falsifiable at D>1.
    `test/test_pepskit_native_dgt1.jl` gates it: native exact density == the
    ED-validated legacy exact oracle at D=2 (1e-9); a lossless native star update
    matches legacy on an anisotropic D=2 cell (1e-9); the exact density is
    invariant under ledger corruption while the mean-field `measure_simple`
    density is not (convention-B source-of-truth); a per-site coordinate-reflection
    trap on a vertically-asymmetric checkerboard; and a native-vs-legacy short-time
    trajectory. See `docs/superpowers/plans/2026-06-04-native-d-gt-1-validation-oracle.md`
    and `notes/stage1-dynamics/native-d-gt-1-validation.md`.
  - *Remaining:* native D>1 vs **ED over the full revival trajectory** (the oracle
    now exists but the D-ladder gate still builds legacy states); native at **D≥3**;
    CTM-χ sensitivity and native CTM observables at D>1; snapshot/load on the native
    engine; a native `:boundary` oracle (current native oracle is `:dense`, ≤16 sites).

## Consolidation roadmap (the payoff)

The file count grew *temporarily* because the new backend is introduced in
parallel with the old one (so parity can be proven before anything is deleted).
The end state is substantially **smaller**, with a single tensor engine and no
ITensors dependency. Six files (~4,200 lines) currently carry the ITensors
backend; the largest single file, `PEPSKitMeasurements.jl` (1,267 lines), is
mostly the ITensors→PEPSKit *conversion adapter* that becomes dead code once
states are already PEPSKit-native.

### Delete once parity is proven

| File | Lines | Why it goes |
|------|------:|-------------|
| `SquareIPEPS.jl` | 591 | legacy Γ-λ state → `PEPSKitBackend.PXPIPEPSState` |
| `StarSimpleUpdate.jl` | 609 | legacy star update → `PEPSKitStarUpdate` |
| `IPEPSEvolution.jl` | 538 | legacy Trotter (538 lines) → `PEPSKitEvolution` (184) |
| `StarModels` (in `PXPModel.jl`) | ~80 | ITensors gate wrappers; the dense gate moves to `pxp_star_gate_tensormap` |
| conversion half of `PEPSKitMeasurements.jl` | ~870 | `to_pepskit_infinitepeps`, weight absorption, dense→ITensor site builders — no ITensors states left to convert |
| simple-on-legacy half of `Observables.jl` | ~450 | mean-field observables on `SquareIPEPSState`; reimplement small native versions from the SUWeight |

Then drop `ITensors` (and likely `Strided`) from `Project.toml` — removing an
entire second tensor engine, its precompile cost, and the dual-maintenance
burden. (`EDKit` stays for the ED reference; `TensorKit`/`PEPSKit` become the
only tensor stack.)

### Merge / dedup

- **One observables module — DONE.** `PEPSKitPXPObservables.jl` (CTM
  observables) and `PEPSKitSimpleObservables.jl` (simple observables + the
  `measure_simple`/`evolve!`/`reverse_evolve!` adapters + converter) are merged
  into a single `src/PEPSKitObservables.jl`. *Remaining:* fold the surviving CTM
  logic of `PEPSKitMeasurements.jl` (params, context, `leading_boundary`,
  `expectation_value`, diagnostics, `correlator`/`correlation_length`, χ-sweep,
  threading, trust integration) into it once the legacy/ITensors conversion half
  of that file is deleted.
- **One star-update file — DONE.** `pepskit_star_update_impl.jl` is inlined into
  `PEPSKitStarUpdate.jl`.
- **Shared TensorKit helpers.** `_dense_operator_tensor_map`, `_dense_index`,
  `squarecoord_to_cartesianindex`, `_star_sites_cartesian` are currently
  duplicated across the measurement, observable, and star-update modules; factor
  them into one place (`PEPSKitBackend`) and import. (Low-risk; can be done even
  before legacy removal.)
- **Retarget, don't rewrite.** `PXPValidation`, `ScarFinder*`,
  `ScarFinderSupport` (snapshots/compression), and `CTMTrust` keep their logic
  but operate on `PXPIPEPSState`; the dual-path/backend-selector scaffolding
  (M6) is *itself temporary* and is removed once the migration completes — we do
  not keep `:legacy_itensors` forever.

### Projected end state

~20 files / ~11 k lines  →  ~14 files / ~6.5–7 k lines (≈40 % smaller), a single
state type (`PXPIPEPSState`), a single star update, a single evolution driver, a
single CTM path, and **no ITensors dependency**. The net simplification is in
code volume and dependency surface, not merely file count.
