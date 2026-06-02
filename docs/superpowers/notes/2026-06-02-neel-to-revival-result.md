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
