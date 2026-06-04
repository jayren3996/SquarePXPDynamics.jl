# Plan: native D>1 validation oracle

> Status: proposed (2026-06-04). Prereq for the PEPSKit-native migration; design in
> `docs/superpowers/specs/2026-06-04-pepskit-native-ctm-aware-migration-design.md`.
> REQUIRED SUB-SKILL: test-first (red → green), one focused PR, no new physics.

**Goal:** Make the existing PEPSKit-native backend (`PXPIPEPSState`) validatable at D>1
against the trusted exact oracle, and add the three D>1 traps that the D=1-only suite
cannot catch. This closes the single biggest validation gap (the native engine is
currently unvalidated at D>1; every native end-to-end test is gauge-trivial D=1, and
`verify_m4.jl` — cited as the D>1 CTM-parity proof — is itself D=1).

**Architecture:** Give `exact_density_finite` / the `exact_*_finite` family a
`PXPIPEPSState` method by converting the native state to per-site dense arrays (the simple
layer already reconstructs these in `PEPSKitObservables.jl:323-333`) and reusing the
existing `:dense`/`:boundary` contractor in `Observables.jl` — no new contraction code.
Then add native D>1 trap tests. Compare **observables only, never raw tensors** (the
engines use different λ conventions: legacy full-λ, native √λ).

**Tech Stack:** Julia 1.12, TensorKit 0.15.3, PEPSKit 0.7.0, EDKit; existing
`FiniteIPEPSObservables` contractor; `test/runtests.jl` (+ `SQUAREPXP_EXTENDED_TESTS`).

**Non-goals:** No CTM-aware evolution (Stage C/D). No public-API migration. No ITensors
removal. No change to `project_star_pepskit!` physics.

---

## Tasks

1. **Add failing `test/test_pepskit_native_dgt1.jl` coverage** (red first):
   - *Anisotropic D=2 single-star round-trip vs dense.* Build a D=2 native state with a
     **non-symmetric** seed (distinct √λ on the right vs up bond, not all-down/checkerboard),
     apply one `project_star_pepskit!`, and assert the 1-site density and the 1-/2-site
     reduced density matrices match an independent dense `ncon` of the *same* native tensors
     to `atol 1e-9`; FAIL if any 1-site density differs >1e-8. (Catches transpose/reflection,
     leg-routing swap, √λ-vs-λ mismatch — all invisible at D=1.)
   - *Convention-B equivalence.* For a native D=2 state, compute density/blockade via
     `measurement_peps` directly, then independently reabsorb the `SUWeight` λ into the bare
     Γ (undo √λ on each leg) and contract the full-weight PEPS; assert agreement to `atol
     1e-10`. (Proves "stored Γ carries exactly √λ".)
   - *Coordinate-reflection trap.* Build a state excited only on the bottom row (NOT
     row-reflection-symmetric), assert its row-resolved density profile matches an explicit
     `row = Ly − y + 1` reference contraction; FAIL if vertically mirrored. Add the
     east/west analogue.
   - *Native vs ED, D=2,3,4 trajectory* (behind `SQUAREPXP_EXTENDED_TESTS`): 4×4 Néel,
     n(t) over [0,2.8] sampled every 0.2, measured by `exact_density_finite` (NO CTM);
     assert `rms[D] < 3.5e-2` and monotone `rms3 ≤ rms2+2e-3`, `rms4 ≤ rms3+2e-3` —
     mirroring `test_d_convergence.jl`, but on the native engine.

2. **Type the exact oracle for native states.** Add `exact_density_finite(::PXPIPEPSState;
   …)` (and the rest of the `exact_*_finite` family as needed) that reconstructs per-site
   dense arrays from `measurement_peps(state)` + the `SUWeight` ledger (reuse `_site_array`
   logic) and dispatches into the existing `:dense`/`:boundary`/`:auto` contractor. No new
   contraction algorithm — only a state→array adapter.

3. **Wire native into `validate_pxp_ed_ipeps` with the exact oracle.** Lift the
   `:pepskit_simple` rejection of `exact_finite_observables` (`PXPValidation.jl:176-180`)
   now that the oracle accepts native states; keep the CTM-trust rejection (still
   legacy-only — out of scope). Add a native D=2 `density_error_exact_finite < 1e-3` at
   t=0.1 gate, mirroring the legacy trap test.

4. **Extend `boundary == dense` to D=5** (and a small-cell D=6 if memory permits) in
   `test_finite_ipeps_observables.jl` so the sole D≥5 oracle is numerically verified at the
   regime where it becomes load-bearing.

5. **Make the tests pass (green).** Fix any convention/leg-routing/√λ bug the traps surface.
   If a trap fails, that is a real native bug — fix `project_star_pepskit!` /
   `PEPSKitObservables`, do not weaken the test.

6. **Docs:** export any new public symbol with a docstring (`test_public_docs.jl` enforces
   it); note the new native oracle path in
   `docs/superpowers/specs/2026-06-04-pepskit-native-ctm-aware-migration-design.md` and
   `docs/pepskit-native-refactor.md` ("What has been migrated").

Verify:

```
julia --project=. -e 'using Pkg; Pkg.test()'
SQUAREPXP_EXTENDED_TESTS=1 julia --project=. test/runtests.jl test_pepskit_native_dgt1.jl test_finite_ipeps_observables.jl
```

Commit as `feat: native D>1 exact-finite oracle + D>1 validation traps`.
