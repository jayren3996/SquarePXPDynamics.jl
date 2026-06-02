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

Running D=4 at chi=16,24 to settle it (`artifacts/neel_chi{16,24}_D4.json`). If
the revival error shrinks with chi → measurement; if it persists → evolution and
Stage-2 needs the principled regauge (rel_floor alone is insufficient at the
revival). Either way: **rel_floor reaches the revival but is not sufficient for
full D-monotone reliability through it** — this reprioritizes Stage-2
(Vidal √λ / full-update or chi scaling) above slimming.

## Takeaway for the roadmap
The true acceptance test (D-converged Néel dynamics to ≥ first revival) is
PARTIALLY met: reachable + D-monotone in buildup, but not D-monotone at the
revival peak. Closing that gap is the central remaining Stage-2 problem.
