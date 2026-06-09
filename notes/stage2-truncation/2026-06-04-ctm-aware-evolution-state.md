# CTM-aware (loop-carrying) evolution — full state snapshot (2026-06-04)

**Purpose of this note:** a complete handoff of *where we are* in implementing a
CTM-aware / full-environment real-time evolution step. Read this together with
`improvement-roadmap.md` (priority #3) and `../methodology/revival-validation.md`.
No new calculations were run to write this; it records the state as of 2026-06-04.

---

## 1. The goal (unchanged)

Implement a **CTM-aware** (full, loop-carrying environment) real-time evolution of
the 5-site square-PXP iPEPS so that the genuine 5-site cross gate is
applied/truncated against the *real* environment instead of the crude single-site
λ² Bethe/tree mean-field. Target: **remove the revival-rise lag** (and, as a bonus,
fix the D=2 blowup, see §4).

"CTM-aware" here means **a loop-carrying environment**, as opposed to the tree
mean-field that severs the closed plaquette loops which rebuild Néel order. CTMRG
is one way to get such an environment; it is not the only one (see §6 engine fork).

## 2. Why (the structural lag — established earlier)

- 4×4 Néel quench: collapse (t<1.4, ED-near-exact) then revival rise (t≈1.4–2.6)
  with a **structural** mean-field lag. n(t) trajectory density RMS over [0,2.8]:
  **D2 2.62e-2, D3 1.82e-2, D4 9.79e-3** (reproduced this session, §3).
- The lag is **not truncation-limited**: across D the discarded weight and the lag
  *decouple*. Legacy probes (`probe_d2/d4_trajectory.out`): at t=2.2, going D2→D4
  the discarded weight *increases* ~4.2× yet the lag *drops* ~5×. At t=2.0, |err|
  vs truncerr is 346× (D2) / 95× (D4). ⇒ the bottleneck is the **environment**
  (severed loops), not the bond-dimension budget.

## 3. What is VALIDATED / reproduced (with exact numbers + file pointers)

All data files live in `notes/stage2-truncation/data/`.

1. **Plan figures reproduced** (legacy ITensors backend probes):
   `probe_d2_trajectory.out`, `probe_d4_trajectory.out`. D4 trajectory RMS
   = 9.788e-3 (matches the plan's 9.79e-3). The truncerr-vs-lag decoupling above.

2. **The CTM-aware bond-truncation PRIMITIVE works** (single horizontal bond) —
   `ctm_bond_proto.out`, script `scripts/g0/ctm_bond_proto.jl`:
   - Frozen native D=2 state at t=1.0, lossless 5-site gate on one star, CTMRG env
     (χ=24), truncate the center→right bond back to D=2.
   - `n_lossless = 0.15183976`
   - `n_CTM-aware = 0.15183974` (|dn| = **1.52e-8**)
   - `n_bare-SVD  = 0.15183936` (|dn| = **3.97e-7**) → **CTM-aware is 26× better.**
   - `bond_truncate` fid = 1.0; `cond(Z'Z) = 3.4e20` (env is near-singular on the
     rise → `fixgauge_benv` will be needed for stability, see §5/§7).

3. **Full-star geometry: horizontal bonds PASS, vertical bonds blocked** —
   `ctm_star_truncate.out`, script `scripts/g0/ctm_star_truncate.jl`. No-op
   self-checks (maxdim=BIG=64 must reproduce n_lossless exactly):
   - right (H): dev = 2.15e-10 ✅   left (H): dev = 2.79e-10 ✅
   - up/down (V, via `rotl90`): **FAILS** with
     `SpaceMismatch("((ℂ^8)' ⊗ (ℂ^6)') ≠ ((ℂ^8)' ⊗ (ℂ^8)')")` inside `bondenv_fu`.
   - **Diagnosis (confirmed by the space dump in the run):** the ℂ^6 is NOT in any
     peps/X/Y/CTM-edge space (all are 2,8,16,24). It comes from a CTM **corner**
     whose χ-bond CTMRG truncated to 6 at that position. The vertical path built a
     **fresh** CTMRG (random init) on the rotated peps, which converges to a
     *different* fixed point with a per-position truncation pattern that is
     internally inconsistent for `bondenv_fu`'s 2×3 window.
   - **Fix identified but NOT YET TESTED:** PEPSKit's own `su_iter`
     (`simpleupdate3site.jl:530`) never rebuilds — it does
     `state2, env2 = rotl90(state2), rotl90(env2)`, rotating the *same* converged
     fixed point so χ stays consistent. Lines ~104 and ~132 of
     `ctm_star_truncate.jl` were edited to use `rotl90(env)` instead of a fresh
     CTMRG on the rotated peps. **The run was killed by the user before it produced
     output, so the rotl90(env) fix is unverified.** (`rotl90(::CTMRGEnv)` and
     `rotl90(::InfinitePEPS)` both exist.)

## 4. Native D=2 blowup — root-caused, and why it matters here

- `native_relfloor_test.out`: native D=2 evolution blows up on the rise
  (t=1.6→1.7: 0.143→**0.396**, ED says ~0.10→0.10) **regardless of rel_floor**
  (1e-3 vs 0 give the same blowup; min bond dim stays 2). ⇒ NOT a rel_floor issue.
- `native_d34_trajectory.out`: native **D=3 and D=4 track ED to ~3 digits**
  (t=1.8: native D4 0.2336, legacy D4 0.2324, ED 0.1556 — i.e. native D4 ≈ legacy
  D4, both lag ED as expected; no blowup). ⇒ the blowup is **D=2-specific**, a
  truncation-subspace defect, NOT a fundamental native bug. Native D≥3 is sound, so
  the CTM-aware machinery can be built and validated at D=4.
- **Expectation:** CTM-aware truncation is gauge-independent (minimizes the true
  norm error), so it should *also* fix the D=2 bad-subspace blowup as a bonus.

- **LOAD-BEARING INSIGHT for how to demonstrate the fix** (`g0_frozen_step.out`,
  script `scripts/g0/frozen_step_diagnostic.jl`): at t=1.8, D=2, the **frozen**
  state is already n=0.434 (ED 0.156) while a **single** step's truncation error is
  only max|dn_trunc| = 5.4e-4 and gate error 9.2e-4. ⇒ **the blowup is CUMULATIVE.**
  Freezing an already-corrupted state and truncating it well cannot un-corrupt it.
  **Therefore a frozen-step test CANNOT demonstrate the fix.** The only convincing
  demonstration is a **full trajectory** run with the better-environment truncation
  applied at *every* step from t=0. This is decisive for choosing the engine (§6):
  the engine must be trajectory-feasible.

## 5. The PEPSKit primitive idiom (so it isn't rediscovered)

Per-BOND CTM-aware truncation with PEPSKit 0.7.0 (validated, §3.2):

```
X, a, b, Y = PEPSKit._qr_bond(A, B)                  # X,Y = env-facing isometries
a_bt = permute(a, ((1,2),(3,)))                      # bond_truncate wants {2,1}/{1,2}
b_bt = permute(b, ((1,),(2,3)))
benv = PEPSKit.bondenv_fu(row, col, X, Y, env)       # CTM bond env, {2,2} on q-bonds
# (optionally) Z = PEPSKit.positive_approx(benv); fixgauge_benv for cond(Z'Z)~1e20
a1, s1, b1, info = PEPSKit.bond_truncate(a_bt, b_bt, benv, ALSTruncation(trunc,maxiter,tol))
a1f, b1f = PEPSKit.absorb_s(a1, s1, b1)
A2, B2 = PEPSKit._qr_bond_undo(X, permute(a1f,((1,2),(3,))), permute(b1f,((1,),(2,3))), Y)
```

Key facts learned:
- `bondenv_fu` is **HORIZONTAL-only** (bond (r,c)-(r,c+1)); vertical bonds go via
  `rotl90` of *both* peps and env (§3.3).
- `BondEnv = AbstractTensorMap{T,S,2,2}` on the two bond-end virtual legs. **`bond_truncate`/`ALSTruncation`/`FullEnvTruncation` accept ANY such {2,2} env** —
  CTMRG (`bondenv_fu`) is just one producer. This is the door to a non-CTMRG
  loop-carrying environment (§6).
- PEPSKit 0.7.0 ships these building blocks **with tests but NO full-update
  driver**, and **no NTU/cluster/patch env builder** (searched: only `bondenv_fu`
  exists under `src/algorithms/contractions/bondenv/`).
- The norm-fidelity environment is **PSD in real time** (the "indefinite metric"
  worry is just regularization), so env-weighted truncation is well-posed.

## 6. THE OPEN ARCHITECTURAL FORK — which engine? (decision pending)

Per-bond CTM truncation is easy and validated. The hard part is the **per-step
environment for a genuine 5-site gate in a trajectory** (this is the plan's "P4").
The 5-site cross gate enlarges **all 4 center-leaf bonds simultaneously**, so you
cannot reduce it to 4 independent 2-site internal-bond enlargements with fixed
X/Y env-isometries (the standard 2-site full-update trick). Candidate engines:

| Engine | Loop-carrying? | Trajectory-feasible? | Status / blockers |
|---|---|---|---|
| **Per-star full CTMRG** | yes (large χ) | **NO** — ~56k–112k CTMRG solves; + the vertical/rotation/χ-consistency pain (§3.3); + 5-site dim-mismatch | validation/accuracy ceiling only. Plan demoted it. |
| **NTU finite-patch** | yes (local plaquettes) | yes (cheap, no CTMRG) | plan's PRIMARY engine; must hand-contract a {2,2} patch env (no PEPSKit helper) and feed `bond_truncate`. Not yet built. |
| **Exact finite-cluster env** | yes (ALL loops on the torus) | yes on 4×4 (the oracle already contracts this every measurement) | NEW idea this session: build the {2,2} bond env by the SAME finite contraction as `_native_dense_state`/`exact_density_finite` but leave the two bond-end virtual legs open. Gold-standard env on the finite torus; reuses validated code; sidesteps CTMRG + NTU + rotation pain. Best first target to *prove the thesis*; NTU/CTM are cheaper approximations to chase later for larger cells. Not yet built. |
| **Warm-started per-step CTMRG** | yes | maybe (1–2 iters/step if env persisted across dt=0.02) | under-explored; still hits the 5-site simultaneous-enlargement dim-mismatch. |

**Recommendation on the table (not yet agreed with author):** build the
**exact finite-cluster bond environment** first — it is the most accurate
loop-carrying env, reuses the already-validated oracle contraction, is feasible on
the 4×4 trajectory, and directly answers "does a full loop-carrying environment
remove the lag / fix the D=2 blowup". If yes, it sets the target and NTU becomes
the cheaper production engine for larger cells.

Oracle contraction to adapt: `_native_dense_state` in
`src/PEPSKitObservables.jl:474-510` (single-layer `ncon` over `measurement_peps`,
PEPSKit leg order [P,N,E,S,W]; bonds via `Hbond`/`Vbond`). The bond env is a
double-layer (bra⊗ket) version of this with two virtual legs left open and the
physical legs closed → a PSD {2,2} `BondEnv`.

## 7. Standing constraints (still in force)

- Genuine **5-site cross gate** (center,right,up,left,down; basis 1=up/excited,
  2=down). **NO** 2-site gate substitution. **NO** PEPSKit generic `time_evolve`.
- **Convention B**: `measurement_peps(state) === state.peps`; every Γ carries √λ;
  the `SUWeight` ledger is update-aid + diagnostic only; the exact oracle is
  ledger-independent.
- Pin **PEPSKit 0.7.0**. No PEPSKit private internals without a shim + contract
  test (deferred for prototyping — the `_qr_bond`/`bondenv_fu`/`bond_truncate`
  calls above are private and currently used directly in throwaway scripts).
- Compare **observables (density), never raw tensors**. Judge by the **n(t)
  trajectory RMS over [0,2.8]**, never a single time. **D=1 is never validation.**
- Do not delete the legacy ITensors backend. Durable knowledge → `notes/`.
- `fixgauge_benv` will be needed (env `cond(Z'Z) ~ 3.4e20` on the rise, §3.2).

## 8. Concrete next steps (in order)

1. **Decide the engine** (§6 fork) with the author. Leaning: exact finite-cluster
   env first (prove thesis), NTU second (production).
2. (If continuing CTM-per-star for validation only) verify the **untested
   `rotl90(env)` fix** in `ctm_star_truncate.jl` closes the vertical-bond ℂ^6
   mismatch; complete the 4-bond no-op self-check; then the D=2 single-shot
   CTM-vs-bare comparison. NB: this validates the *primitive on the full star* but
   is trajectory-infeasible, so it cannot by itself show the lag/blowup fix (§4).
3. Build the chosen env as a {2,2} `BondEnv`, wire it into the 4-bond split of
   `project_star_pepskit!` (replacing the bare-SVD `_split_bond`), behind an
   **opt-in flag** in `project_star_pepskit!` / `evolve_pepskit!`. Add `fixgauge_benv`.
4. **Correctness test**: no-op at maxdim=BIG reproduces the lossless density; then
   a **full D=2 and D=4 trajectory** vs ED/legacy — must (a) fix the D=2 blowup,
   (b) reduce the D=4 lag (trajectory RMS over [0,2.8]).
5. Convention-B writeback + weight-ledger bookkeeping for the new split.

## 9. File map

**Scripts** (`scripts/g0/`): `ctm_bond_proto.jl` (✅ single-bond CTM primitive),
`ctm_star_truncate.jl` (4-bond star; H ✅, V blocked, rotl90(env) fix untested),
`frozen_step_diagnostic.jl` (cumulative-blowup insight), `native_relfloor_test.jl`
(rel_floor ruled out), `native_d34_trajectory.jl` (D=2-specific), legacy
`probe_d2_trajectory.jl` / `probe_d4_trajectory.jl`. `native_evolve_anomaly.jl` was
redundant (killed/superseded by native_d34).

**Data** (`notes/stage2-truncation/data/`): the `*.out` companions of the above.

**Source touch-points**: `src/PEPSKitStarUpdate.jl` (`project_star_pepskit!`,
`_split_bond`), `src/PEPSKitEvolution.jl` (`evolve_pepskit!`),
`src/PEPSKitObservables.jl` (`exact_density_finite`, `_native_dense_state`,
`pepskit_native_context`), `src/Lattice.jl` (`square_star_sites`,
`PeriodicSquareUnitCell`).

**PEPSKit (0.7.0, `~/.julia/packages/PEPSKit/zz9BS/`)**: bond primitives in
`src/algorithms/time_evolution/evoltools.jl` (`_qr_bond`, `_qr_bond_undo`),
`src/algorithms/contractions/bondenv/{benv_ctm.jl,als_solve.jl,gaugefix.jl}`
(`bondenv_fu`, `bond_truncate`, `positive_approx`, `fixgauge_benv`), rotation idiom
in `src/algorithms/time_evolution/simpleupdate3site.jl:530`.
