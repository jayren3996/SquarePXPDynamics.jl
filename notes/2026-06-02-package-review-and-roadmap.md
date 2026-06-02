# SquarePXPDynamics.jl — consolidated review + 3-stage roadmap

Date: 2026-06-02. Author: review campaign (4-agent audit + independent read +
empirical ground-truth). Supersedes scattered status in `memory/short_term/`.

## Executive summary

The package is **heavier than it needs to be (~3–4k of 12k src lines are
removable)** and its *physics core is actually sound where it counts* — but the
reliability is gated by three concrete, fixable defects, one per user stage.
Crucially, several "known facts" in `memory/` are **stale**:

- **The test suite was completely broken** (`Pkg.test()` could not even load:
  `runtests.jl` uses `Distributed`+`Printf` but they were undeclared in
  `test/Project.toml`). Fixed 2026-06-02; full suite now passes **32/32 files**.
  The memory's "65356 passing" predated the parallel-`pmap` harness switch.
- **Néel/checkerboard does NOT crash** (memory said it did). It evolves fine on
  even-tileable cells with `schedule=:serial`. Only `:z_up` (all-up, genuinely
  blockade-forbidden) and odd-cell checkerboard *seams* crash.
- **The iPEPS evolution is correct.** Exact-finite contraction matches ED to
  ~1e-7 for D=1–4 at short time (truncerr ~6e-29). The D>1 "error" users see is
  a **measurement bug** in `measure_simple`, not a dynamics bug.

## Ground-truth measurements (this campaign)

| Check | Result |
|---|---|
| `Pkg.test()` full suite | **32/32 PASS** after `Distributed`+`Printf` fix (~10 min, 32 workers) |
| 3×3 t=0.02 D=1 | exact_finite & simple both match ED to 7e-7 ✓ |
| 3×3 t=0.02 D=2,3,4 | **exact_finite matches ED to ~1e-7**, but `simple` off by **1.9e-4** |
| `:z_up` evolve | crash: "zero or non-finite singular spectrum" (blockade-forbidden) |
| `:checkerboard` evolve (even cell, serial) | **OK** |

Slowest tests: `test_ipeps_compression.jl` 582s, `test_pepskit_measurements.jl`
511s, `test_ctm_gauge_readiness.jl` 217s — too slow for an iterate gate.

## Stage 1 — reliable iPEPS dynamics matching ED

**HARD RULE (owner, 2026-06-02): D=1 is NEVER validation.** iPEPS reliability is
judged ONLY by the trusted observable (`exact_density_finite`/CTM) converging
toward ED across a D-ladder (D≥2,3,4) in an entangled regime; error vs ED must be
non-increasing in D within tolerance and shrink toward ED where D matters; a
larger-D-worse-than-smaller-D result is a HARD REGRESSION. The simple/local
observable is mean-field (wrong for D>1) and is never the convergence metric. See
`memory/stage1_d_convergence_rule.md`; enforced by `test/test_d_convergence.jl`.

D-ladder evidence (3x3, dt=0.02, cutoff 1e-12, rel_floor 1e-4; |exact_finite−ED|):
t=0.1 → D1 1.6e-3, D2/3/4 2.8e-6 (D=1 useless, D≥2 converged). t=0.3 → D1 6.0e-2,
D2/3/4 1.6e-4. t=0.5 (entangled) → D1 1.5e-1, D2 1.21e-3, **D3 1.24e-3 (slightly
worse — residual non-monotonicity)**, **D4 7.3e-4 (better than D2 — D helps)**.
The D=3/t=0.5 wrinkle is the next Stage-2 target (tune rel_floor / Vidal √λ).

**Status: largely achieved. Evolution is correct; the catastrophic Stage-2
conditioning defect is fixed (rel_floor); the D>1 "simple error" is an inherent
diagnostic limit, not a bug; one small residual D-non-monotonicity remains at long time.**

1. **D>1 `measure_simple` density is NOT a bug — UPDATED 2026-06-02 (CASE B,
   rigorously settled).** `local_density_simple` computes the canonical Vidal
   single-site mean-field RDM *exactly* (a hand-built λ² RDM matches it to machine
   zero). The ~2e-4 gap vs ED is the inherent single-site Bethe/mean-field failure
   on a loopy lattice: 1-site off 3.5e-4, but 2-site matches to 3.8e-27, 5-site
   star to 2.7e-7. **Do not edit `Observables.jl`.** D>1 density already routes
   through `exact_density_finite`/CTM in the validation pipeline. The
   `density_error_simple > 1e-4` lower-bound tests correctly *document* this
   offset — keep them; they are not traps. See `memory/stage1_simple_density_inherent.md`.
2. **Long-time conditioning drift FIXED** by the `rel_floor` condition floor
   (default 1e-4): at t=0.2, D=4, cutoff 1e-12 the exact_finite error drops
   9.7e-3 → 8e-12 (D=2/3/4 collapse to the ED-converged value). No need to
   loosen cutoff. The principled Vidal √λ successor remains optional Stage-2 work.
3. **Genuinely missing:** a real-`measure_ctm`-vs-ED D=2 test (all existing
   CTM-value tests use fake closures). This is the one Stage-1 test gap worth adding.

## Stage 2 — bond truncation + proper regauging  ← the crux

**Status: bare/mean-field SVD truncation with NO regauging. The named
environment-aware capability (`fix_bond_gauge!`/CTMGaugeReadiness, 662 lines) is
DEAD CODE — never called by `evolve!`/`project_star!`/`compress_to_target_maxdim!`/
ScarFinder.**

Root cause of the larger-D-worsens instability (critical): `project_star!`
absorbs the **full** λ then `deabsorb_link_weight` divides the full λ back out
with only a 1e-14 floor (`StarSimpleUpdate.jl:231/344`, `SquareIPEPS.jl:403-408`).
At cutoff 1e-12, retained λ near ~1e-6–1e-12 get inverted to ~1e6–1e12 on the
next neighboring star update → poisons the state. The diagnostics that measure
exactly this (`min_lambda`, `touched_min_lambda`, `diagonal_condition`) are
recorded but never act as a control.

Fixes (in order):
- **(S, stopgap, high payoff)** Relative singular-value condition floor in
  `_split_reduced_theta`: drop retained λ_i below `rel_floor·λ_max`, caps the
  bond condition number, kills the inversion blowup *without* loosening cutoff.
- **(L, principled)** Vidal √λ absorb/deabsorb so the next update only ever
  divides √λ of *other* bonds, never the freshly-truncated tiny λ.
- **Decision needed**: wire `fix_bond_gauge!` into a post-compression regauge
  (reconcile its √λ-vs-full-λ convention mismatch first) **or delete** the
  662-line orphan.

## Stage 3 — ScarFinder, better-than-Néel initial state

**Status: premise narrowly unblocked, but the search is NOT expressible today.**

- `scarfinder!` is a **single-seed** evolve-(compress)-measure-rank loop. It
  ranks ONE state's trajectory snapshots. There is **no search over initial
  states**, no optimizer, and **no Néel baseline** — "find a better initial
  state than Néel" has no API surface.
- The "**revival**" objective is just instantaneous staggered magnetization
  `|n_even − n_odd|` at one time — not a revival. **No overlap / return-amplitude
  / Loschmidt observable exists anywhere** in src. Rewarding instantaneous
  imbalance is maximized at t=0 by Néel itself.
- The `:z_up`/seam crash throws opaquely deep in the star split. Needs a legible
  error or graceful skip.
- Bug: the audit silently drops reversibility for **all checkerboard rows**
  (`ScarFinderAudit.jl:338` calls `product_square_ipeps` instead of `_build_state`).
- Bug: `scripts/run_scarfinder_audit.jl` defaults `initial_state=:z_up` → driver
  broken out of the box.

Path to real Stage 3: (1) legible crash handling; (2) a real return-amplitude/
revival observable over the trajectory time series; (3) an `scarfinder_search`
outer layer over a candidate-state family ranked vs the Néel baseline;
(4) optionally an energy-variance proxy (the standard scar selector).

## Heaviness — concrete slim list (~3–4k lines)

| Target | Lines | Action |
|---|---|---|
| `CTMGaugeReadiness.jl` (`fix_bond_gauge!`) | 662 | wire into Stage 2 **or** delete (orphan) |
| TFIM axis: `Benchmarks.jl`+`FiniteTFIM*`+TFIM model/obs | ~950 src +700 test | move to `examples/` or delete (orthogonal to PXP) |
| Duplicate PXP campaign drivers + struct families + serialization | ~550 | merge into one grid runner |
| Placeholder ScarFinder objectives (TargetEnergy, LowVariance) + JSON store | ~150 | delete until real |
| `*_simple` individual exports (D>1-divergent) | ~15 exports | demote from public API |
| `SquarePEPS.jl` (lends 2 method names), `GaugeDiagnostics.jl` (unused) | ~275 | delete/inline |
| Dead exports (census) + 245-symbol surface | ~30-40 exports | strip; relax `test_public_docs` to core API |
| 30 design docs incl. kagome/PESS (never built) | — | archive |

## Test integrity — fix the traps before fixing the physics

- **Golden tests PIN the D>1 simple breakdown as "expected"** (`>1e-4` lower
  bounds, `@test_broken` on the correct values) → these actively trap any fix.
  Convert to `@test_broken` + add `exact_finite < 1e-6` upper-bound gates.
- **No real-CTM-vs-ED D>1 test** — every CTM-value test uses fake closures.
- Move Aqua + schema/serialization smoke to a nightly tier.

## Iterate loop design

- **Fast gate (every pass, seconds in a warm process):**
  `test_pxp_validation.jl` + `test_pxp_d2_localization.jl` +
  `test_ipeps_compression.jl` + `test_star_simple_update.jl`. Holds the D=1 ED
  match, the D=2 exact-finite-vs-simple separation, and compression bookkeeping.
  Kill the ~75s compile tax with a warm Julia daemon / sysimage.
- **Nightly:** full suite incl. Aqua, real-CTM-vs-ED, long-time D=4 conditioning.
- **Improvement focus order:** (1) Stage-2 relative floor [crux] → (2) Stage-3
  legible crash + audit/script fixes → (3) de-trap tests → (4) Stage-1 simple
  measurement bug → (5) Stage-3 revival observable + search layer → (6) heaviness
  slimming (gated on owner decisions) → (7) Vidal √λ / CTM regauge.

## Owner decisions that fork the work

1. **"Better than Néel" metric**: longer-lived staggered-mag revival vs lower
   energy variance vs larger overlap-return amplitude? Search space: per-sublattice
   product-angle family vs low-D seeds?
2. **Heaviness disposition**: aggressive delete (TFIM, fix_bond_gauge!, dup
   drivers) vs conservative demote-to-experimental?
3. **Stage-3 geometry**: fixed even-Lx/Ly + `:serial` (Néel works today) vs
   `:five_color` (forces 5-divisible odd-tiling cells that can't host Néel)?
