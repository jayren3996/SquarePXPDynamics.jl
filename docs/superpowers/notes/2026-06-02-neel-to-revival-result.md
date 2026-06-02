# Néel-to-revival reliability result (2026-06-02)

The decisive Stage-1/2 benchmark: 4×4 Néel, PXP, simulate to the FIRST n(t)
revival and check D-convergence vs ED. Initial state checkerboard, schedule
:serial, dt=0.02, order 2, cutoff 1e-12, rel_floor 1e-4 (default), CTM chi=8.
ED: n collapses to 0.101 at t=1.4, first revival peak 0.483 at t=2.6.
Artifact: `artifacts/neel_to_revival_4x4.json`.

## Result — reaches the revival; D-monotone in buildup; BREAKS at the peak

|iPEPS_CTM − ED| through the trajectory:

| time | ED | D=2 | D=3 | D=4 |
|---|---|---|---|---|
| collapse t=1.4 | 0.101 | 0.0138 | 0.0017 | **0.0002** |
| rise t=2.0 | 0.337 | 0.0579 | 0.0363 | **0.0129** |
| **revival peak t=2.6** | 0.483 | 0.0052 | **0.0006** | **0.0221 ← worse** |
| end t=3.0 | 0.408 | 0.0537 | **0.0039** | 0.0120 |

max_truncerr: D2 1.3e-4, D3 9.1e-5, **D4 5.3e-4 (highest)**.

- **Positive:** the iPEPS reproduces the qualitative scar revival (collapse +
  return) at all D, and through the collapse/rise the error decreases
  monotonically with D (D=4 best) — the Stage-1 D-convergence rule holds in the
  buildup. So the revival window (t≈2.6) IS reachable.
- **Problem (the Stage-2 failure, exactly where physics matters):** at the
  revival peak and beyond, D=4 error (0.022) is WORSE than D=2/D=3 (<0.005) — a
  D-non-monotonicity at the most-entangled point. This is the sharp form of the
  t=0.5 D=3 wrinkle.

## Two candidate causes (opposite fixes) — chi-sweep disentangles
1. **CTM chi=8 under-measures** the more-entangled D=4 state at the revival
   (measurement artifact). Fix: larger chi.
2. **rel_floor/simple-update over-truncates or mis-regauges** D=4 at high
   entanglement (evolution artifact; supported by D4's highest max_truncerr).
   Fix: Vidal √λ / environment-aware regauge, or an adaptive/looser floor.

**VERDICT (chi-sweep, settled):** chi=8 → chi=16 did NOT reduce the D=4 revival
error (0.0221 → 0.0241, unchanged; identical max_truncerr 5.3e-4). The evolution
is identical across chi; only measurement resolution changed and the error did
not move. **So the D=4 revival error is an EVOLUTION/regauge limit, NOT a CTM
measurement limit.** Larger bond dimension gives WORSE dynamics at the revival —
the simple-update mean-field environment failing at high entanglement, with
rel_floor=1e-4 cutting directions D=4 needs.

**Consequence:** rel_floor reaches the revival and is D-monotone through the
buildup, but is NOT sufficient for D-monotone reliability AT the revival. The
remaining Stage-2 fix must be environment-aware truncation/regauging (Vidal √λ
canonicalization and/or full-update with the CTM environment), not a measurement
knob. This is now the empirically-justified TOP Stage-2 priority — above slimming.
Sub-question being probed: is rel_floor=1e-4 over-truncating D=4 (fix = looser/
adaptive floor, easy) or is it the simple-update environment itself (fix =
full-update, hard)? See `artifacts/` rel_floor sweep at the revival.

## rel_floor sweep at the revival (D=4, t=2.6, ED peak 0.4825) — regime-dependent

| rel_floor | revival error | max_truncerr |
|---|---|---|
| 1e-2 | 0.0740 (over-truncated) | 9.6e-5 |
| 1e-3 | **0.0057 (good)** | 1.8e-3 |
| **1e-4 (default)** | **0.0215 (bad pocket)** | 5.3e-4 |
| 1e-5 | 0.0225 (bad) | 9.5e-4 |
| 1e-6 | **0.0056 (good)** | 6.7e-4 |
| 0 | **0.0035 (best)** | 8.4e-4 |

Non-monotonic and regime-dependent: the rel_floor that is OPTIMAL at short time
(1e-4: collapses D=2/3/4 to 8e-12 at t=0.2/cutoff 1e-12) is in a BAD POCKET at
the revival; rel_floor=0 (which FAILS at short time with the conditioning blowup)
is BEST at the revival. **No single rel_floor is optimal across regimes.**

## DECISIVE revival D-ladder across floors (2026-06-02, current code)

Re-ran 4x4 Neel to the revival t=2.6 (CTM chi=8 vs ED 0.4825) for the FULL
D-ladder at several floors (scripts/dev_relfloor_revival_confirm.jl +
dev_relfloor_revival_lowfloor.jl):

| rel_floor | D=2 | D=3 | D=4 | worst-D | worst D |
|---|---|---|---|---|---|
| 1e-4 (old default) | 4.7e-3 | 4.3e-6 | 2.2e-2 | 2.2e-2 | D=4 catastrophe |
| 1e-3 | 3.4e-3 | 2.6e-3 | 5.7e-3 | 5.7e-3 | D=4 mild |
| 1e-6 | 4.7e-3 | 1.3e-2 | 5.6e-3 | 1.3e-2 | D=3 |
| 0    | 4.7e-3 | 1.3e-2 | 3.5e-3 | 1.3e-2 | D=3 |

**No floor is D-monotone at the revival, and WHICH D is worst shifts with the
floor (D=4 at 1e-4, D=3 at 0/1e-6).** This erratic shift confirms the revival
non-monotonicity is the mean-field-ENVIRONMENT limit, not a conditioning effect:
no rel_floor (constant OR adaptive) fixes it. FULL UPDATE is required for true
D-monotone reliability at the revival.

Floor consequences (a floor is a robustness knob, not an accuracy fix):
- **1e-3 is the most ROBUST floor** — tightest D-band (2.6e-3..5.7e-3), kills the
  D=4 catastrophe (2.2e-2 -> 5.7e-3). Cost: it over-truncates short time, capping
  all D at ~8e-5 (vs 1e-4's 2.5e-5) and flattening short-time D-resolution.
- **rel_floor=0 is NOT robustly best** (D=3 -> 1.3e-2). So adaptive-toward-0 is
  WRONG; the only defensible ramp is 1e-4 (low entanglement) -> 1e-3 (high), which
  is mostly captured by just using the constant 1e-3.
- Regauging (canonicalize_simple!) was separately ruled out (gauge transform;
  cannot change truncation quality). See 2026-06-02-stage2-regauge-map.md.

## Takeaway for the roadmap
The true acceptance test (D-converged Néel dynamics to ≥ first revival) is
PARTIALLY met: reachable + D-monotone in buildup, but not D-monotone at the
revival peak with the current default. Key conclusions:
1. `rel_floor` is a useful STOPGAP for the short-time conditioning instability,
   but it is a crude, regime-dependent knob — NOT a principled reliability fix.
   The default 1e-4 should be revisited (it is bad at the revival); a robust
   default likely needs the floor to be adaptive or replaced.
2. The principled Stage-2 fix is environment-aware truncation/regauging — Vidal
   √λ canonicalization and/or full-update with the CTM bond environment — which
   would give D-monotone reliability across regimes WITHOUT a hand-tuned floor.
   This is the empirically-justified TOP Stage-2 priority (above slimming).
3. Immediate cheap step: re-verify a robust rel_floor (or an adaptive rule) on
   BOTH the short-time D-ladder AND the revival D-ladder before changing the
   default; never tune to one regime alone (that is how 1e-4 got its bad pocket).

## EXACT-ORACLE DECONTAMINATION (2026-06-02) — the benchmark is now CTM-free

The chi-sweep above argued "env-limit, not measurement" from chi 8→16 not moving
the error. The exact 16-site contraction (`exact_density_finite(max_sites=16)`,
zero environment approximation) now settles it directly and REFINES that argument.
`scripts/dev_exact_oracle_decontaminate.jl` (refactor bit-identical to old per-site
formula, |diff|=0). 4×4 Néel to t=2.6, CTM chi=8 vs EXACT-16, vs ED 0.4825:

| rel_floor | D | CTM chi=8 err | EXACT-16 err | CTM contamination |
|---|---|---|---|---|
| 1e-4 | 2 | 4.7e-3 | 1.01e-2 | 5.4e-3 |
| 1e-4 | 3 | **4.3e-6** | 2.89e-3 | 2.89e-3 |
| 1e-4 | 4 | 2.15e-2 | **3.45e-2** | 1.30e-2 |
| **1e-3** | 2 | 3.4e-3 | **5.91e-3** | 2.50e-3 |
| **1e-3** | 3 | 2.6e-3 | **5.48e-3** | 2.91e-3 |
| **1e-3** | 4 | 5.7e-3 | **9.65e-3** | 3.92e-3 |

1. **CTM chi=8 contaminated this benchmark by 2.5e-3–1.3e-2 — the SAME ORDER as the
   signal.** It FLATTERED every result (exact error is always larger). The
   near-perfect chi=8 D=3/1e-4 (4e-6) was a measurement artifact; the true error is
   2.9e-3. chi 8→16 (0.0221→0.0241) was creeping toward the exact 0.0345 — CTM was
   under-converged, not converged. So the old "the error did not move" sub-claim was
   imprecise; the env-limit CONCLUSION stands and is now exact-confirmed.
2. **The revival D-non-monotonicity is REAL evolution error.** D=4 is worst at the
   peak under exact measurement too, and BIGGER than CTM showed. The mean-field
   environment ceiling is not a CTM artifact. (See
   `memory/stage2_meanfield_environment_ceiling.md`.)
3. **1e-3 > 1e-4 holds under exact measurement** (D=4 9.65e-3 vs 3.45e-2
   catastrophe) — the rel_floor default decision was not a CTM illusion.

The exact-16 oracle is now the trusted benchmark measurement (committed as the
`test_d_convergence.jl` revival D-ladder, `@test_broken err[4]<=err[3]`). Score
every future env fix on it, NOT on CTM. The implementation is cheap: the existing
`dense_state_finite` already contracts the 16-site cell exactly (~4s, ~const in D);
`exact_density_finite` was just refactored to build that state once (not 16×).
