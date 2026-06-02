# Validating PXP revival dynamics

The load-bearing methodology for this project. Detailed record:
`../stage2-truncation/2026-06-02-neel-to-revival-result.md`.

## 1. Judge by the n(t) TRAJECTORY, never a single time

The benchmark is the 4×4 Néel quench: n(t) collapses (~0.10 at t≈1.2–1.4) and
revives (peak 0.483 at t≈2.6). **Always compare the full n(t) trajectory error
(max/RMS of |n_iPEPS(t) − n_ED(t)| over [0, ~2.8]), not the value at one time —
least of all the revival peak t=2.6.**

t=2.6 is a CROSSING POINT: every D-curve passes near the ED peak there, so per-D
errors converge AND scramble. Evaluating only at t=2.6 produced a chain of FALSE
conclusions this project had to retract ("D=4 bad pocket", "D=5 best", "mean-field
ceiling"). By the trajectory metric the ladder is instead CLEAN MONOTONE:

| D | rms\|err\|_traj | err@2.6 (misleading) |
|---|---|---|
| 2 | 2.6e-2 | 6.5e-3 |
| 3 | 1.8e-2 | 6.0e-3 |
| 4 | **9.8e-3 (best)** | 1.0e-2 (looked worst) |

The error lives in the RISE back to the revival (t≈1.4–2.6) — the iPEPS lags ED
slightly; higher D reduces the lag. The collapse (t<1.4) is tracked near-perfectly
at all D. Scripts: `scripts/dev_revival_trajectory.jl`, `plot_revival_trajectory.py`.

## 2. The D-convergence HARD RULE

Judge iPEPS reliability ONLY by a TRUSTED observable converging toward ED across a
D-ladder (D≥2,3,4) in an entangled regime. **D=1 is a product state and is NEVER
validation.** The trusted observable is the exact finite contraction
(`exact_density_finite`) or CTM — NOT the simple/mean-field observable (wrong for
D>1, see `../stage1-dynamics/simple-density-inherent.md`). Error vs ED must be
non-increasing in D within tolerance; a larger-D-worse-than-smaller-D result (at
fixed dt/cutoff/time, by the TRAJECTORY metric) is a hard regression. Enforced by
`test/test_d_convergence.jl`.

## 3. Use the EXACT oracle, not CTM, for the headline benchmark

CTM chi=8 was found to *flatter* the 4×4 revival error by ~3–13e-3 (it
under-measures). The trusted measurement is the exact 16-site contraction
`exact_density_finite(max_sites = 16)` (~4s/call, no environment approximation).
Caveat: the dense contractor OOMs at the D=5 trajectory and D=6 — reaching higher
D needs a memory-efficient boundary-MPS contractor (bond ≤ D⁴), the current top
infrastructure priority.

## 4. Benchmark must reach ≥ the first revival

Short benchmarks (t≤0.5) prove nothing about scar dynamics — they only see the
collapse. Always simulate to at least the first n(t) revival (t≈2.6 for 4×4 Néel;
6×6 collapse t≈1.3, revival t≈2.6).
