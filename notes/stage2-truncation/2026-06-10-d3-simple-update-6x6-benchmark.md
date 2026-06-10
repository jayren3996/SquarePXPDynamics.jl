# Plain simple update solves the PXP quench — 4×4 closed, 6×6 benchmarked (2026-06-10)

**TL;DR.** The 2D PXP Néel-quench goal ("reliable n(t) through collapse and revival,
benchmarked against exact ED") is met by **plain native simple update with `rel_floor=0`
and no environment machinery at all** — on the original 4×4 testbed *and* on a 6×6 torus
against a 36-site exact Krylov ED oracle. The entire CTM/exact-cluster/NTU environment
program of 2026-06-04..09 was rescuing a D=2-specific pathology; at D≥3 there is nothing
to rescue. Figure: `artifacts/neel_6x6_d3_vs_ed.png`.

## 1. 4×4 baseline at rel_floor=0 (`native_relfloor0_baseline.out`)

D=3 / D=4, dt=0.02 (+ dt=0.01 control), t=0→3.0, vs `artifacts/neel_to_revival_4x4.json`:

| regime | D=3 | D=4 |
|---|---|---|
| collapse+onset (t≤1.7) | ≤5e-5 | ≤4.4e-4 |
| rise (t=2.0–2.4) | −0.005…−0.010 | −0.005…−0.006 |
| peak t=2.6 | −0.016 (−3.4%) | **−0.008 (−1.7%)** |
| t=2.8 | −0.036 | −0.018 |
| RMS (0.2–2.8 grid) | 1.13e-2 | **5.9e-3** (best to date) |

- **dt is converged**: dt=0.01 reproduces dt=0.02 to ≲1e-3 everywhere → the late-time
  deficit is bond-dimension truncation (entanglement growth after the onset), not Trotter.
- **D-ladder is monotone** at rel_floor=0 (D3→D4 halves the late error); no conditioning
  crash (the old "rel_floor=0 crashes at D=4" was the legacy backend).
- The old flat-ladder / "D=3 lags ED by ~1e-2 on the rise" numbers were the
  `rel_floor=1e-3` floor artifact (see addendum in `rel-floor.md`).

## 2. The D=2 instability is size-independent (`native_6x6_d2_control.out`)

6×6 D=2 control: tracks ED to t=1.6 (−0.009), then the same catastrophic jump as 4×4
(t=1.7: n=0.440 vs ED 0.198), revival destroyed (trajectory RMS 0.18). Constraint
leakage stays ≤2e-6 throughout — the blowup is an **in-sector truncation-subspace
pathology of D=2**, not constraint violation and not a 4×4 finite-size artifact.

## 3. Observable for >16 sites: the sector contraction (the one new build)

The dense single-layer oracle caps at 16 sites (2^N). Two measured dead ends, then the fix:

- **Legacy `:boundary` ring contraction — dead for this regime.** The native→legacy
  adapter itself is bit-exact (`legacy_boundary_gate.out` phase A: diffs ~1e-16), but the
  exact (cutoff 1e-14) ring rank-explodes on entangled states (measured: ~48,600-dim ring
  bonds, 74 GB, >1 h for ONE 4×4 density at t=0.5; the 6×6 t=0.2 call never finished in
  11 h). Root cause: the torus wrap legs ride the ring uncompressed, inflating the "true
  Schmidt rank" regardless of physical entanglement. A `boundary_cutoff` knob was added to
  `exact_density_finite` (src/Observables.jl) and calibrated
  (`boundary_cutoff_calibration.out`): accurate (−4e-5 at cutoff 1e-5) but still 354 s …
  2.6 h per 4×4 density — structurally too slow. Kept as documented reference.
- **`sector_density` (scripts/g0/sector_density.jl) — the production observable.**
  Meet-in-the-middle contraction over the hard-square-constrained sector: single layer
  (bond D, not D²), physical configs fixed per constrained row-stack, row-transfer
  matrices cached per (row, ring-config), halves tree-shared, joined by one BLAS GEMM
  over both torus seams. Exact in-sector density + seam-leakage diagnostic; no SVD, no
  cutoff, no gauge. 6×6 D=3: **~80 s / ~41 GB per evaluation** (2,406,862 valid configs).
  Gated at machine precision against an independent dense-sector reference on 4×4
  (|Δn| ≤ 9.4e-16 at t=0.5/1.6/2.6; per-config amplitudes ≤2.4e-17) and to **1.8e-9**
  against the 6×6 Krylov oracle at t=0.2 (`sector_density_gate.out`); adversarially
  reviewed (zero blocking findings; combinatorics independently recomputed).
  Limits: even Nr; memory ~(stacks)·D^(2Nc) → D=4 6×6 needs GEMM slicing (~800 GB dense,
  feasible out-of-core in ~10 min/point if ever needed).

## 4. The 6×6 result (`native_6x6_d3_{even,odd}.out`, merged 0.1 grid, 30 points)

D=3, dt=0.02, rel_floor=0, two staggered lanes vs `artifacts/neel_ed_krylov_6.json`:

| t | 0.1–1.3 | 1.7 | 2.0 | 2.3 | **2.6 (peak)** | 2.8 | 3.0 |
|---|---|---|---|---|---|---|---|
| Δn | ≤2e-5 | −4e-4 | −0.0047 | −0.0084 | **−0.0102 (−2.1%)** | −0.0155 | −0.0273 |

- **RMS over [0.1, 2.8] = 5.5e-3** — equal to the best 4×4 number (D=4), at 36 sites
  with only D=3. Max leakage 3.2e-6. Lanes evolved independently and interleave smoothly
  (no grid/lane artifact).
- The revival is **size-converged**: 6×6 ED peak 0.4825 vs 4×4 ED 0.4831 at the same
  t=2.6, and the PEPS reproduces it — the t≈2.6 revival is not a 4×4 artifact (answers
  external-review question #2 of `2026-06-09-external-review-status.md`).
- Cost of the whole benchmark: ~40 min/lane wall (evolution ~25 s per 0.1; observable
  ~80 s/point), i.e. the 36-site oracle-validated quench is now a routine run.

## 5. What this supersedes

- The **CTM/exact-cluster/NTU environment program** (2026-06-04..09) is superseded *for
  the dynamics goal*: its effect-A result (a loop env stabilizes D=2 to ±2%) stands as
  the autopsy of the D=2 pathology, but the production answer is "use D≥3", not "fix
  D=2's environment". Status note `2026-06-09-external-review-status.md` remains the
  record of that program.
- `rel-floor.md`'s "keep 1e-3" verdict was legacy-backend-scoped; for the **native**
  backend the floor is an accuracy artifact and rel_floor=0 is correct (addendum added).
- The 4×4 16-site dense-oracle ceiling: `sector_density` is the trusted finite-torus
  observable beyond it (gated; constraint-sector-restricted, leakage reported).

## 6. Open levers (not goals)

- D=4 at 6×6 (sector observable needs sliced GEMM; or revisit a cached/looser ring).
- Longer time (second revival t≈5.2): extend the Krylov oracle (`run_neel_ed_krylov.jl`,
  ~minutes) and run the same lanes longer; D-deficit will grow with entanglement.
- 8×8: Krylov oracle infeasible (~1.5^64); sector observable scales as stacks·D^(2·8)
  (3^16≈43M seam — needs slicing); validation would be 6×6-anchored + D/dt convergence.

## 7. Data index

`notes/stage2-truncation/data/`: `native_relfloor0_baseline.out` (4×4 D3/D3-dt01/D4),
`native_6x6_d3_even.out` + `native_6x6_d3_odd.out` (the benchmark),
`native_6x6_d2_control.out`, `sector_density_gate.out` (observable gate),
`legacy_boundary_gate.out` + `boundary_cutoff_calibration.out` (ring-boundary autopsy).
Scripts: `scripts/g0/{native_relfloor0_baseline,sector_density,sector_density_gate,native_6x6_d3_trajectory,native_to_legacy,boundary_cutoff_calibration}.jl`.
Figure: `artifacts/neel_6x6_d3_vs_ed.png`.
