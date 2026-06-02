# PXP iPEPS+CTM Larger-D Debug — Root Cause

Date: 2026-05-19

## Problem

For `3x3`, `dt = 0.02`, `total_time = 0.2`, the CTM-measured density error
appeared to grow with iPEPS bond dimension `D` at the original cutoff
`1e-12`, contradicting the naive expectation of monotone convergence.
Final-time site-averaged density from
`artifacts/m3-systematic/threaded-ctm-trajectory-3x3-t020-dt002-D2D3D4-chi2chi4.{csv,json}`:

- ED reference (symmetric PBC, global site average): `0.0374514509`
- D=2 CTM density: `0.0374430017` (error `≈ -8.45e-6`)
- D=3 CTM density: `0.0368134433` (error `≈ -6.38e-4`)
- D=4 CTM density: `0.0277543865` (error `≈ -9.70e-3`)

`ctm_trust_status` was `unaccepted_diagnostics` for these low-`chi` runs,
which is consistent with finite-`chi` CTM not being trusted but does not
itself explain the D-monotone worsening.

## Root cause

Numerical / algorithmic — **not** a CTM density-operator code bug. The
larger-D failure localizes to **simple-update truncation cutoff stability**:
at a very tight `cutoff = 1e-12`, the QR/SVD split inside `project_star!`
retains ill-conditioned near-zero singular values, and those directions
poison subsequent steps via gauge conditioning. The CTM measurement faithfully
tracks the evolved state.

## Evidence

Exact finite iPEPS contraction (`exact_density_finite`) on the evolved 3×3
state tracks CTM density, so the discrepancy is in evolution, not measurement.
Sweeping `cutoff ∈ {1e-12, 1e-10, 1e-9, 1e-8}` at `t=0.2`, all from the same
ED reference `0.0374514509`
(source: `artifacts/m3-systematic/cutoff-exact-3x3-t020-dt002-D2D3D4.csv`):

| D | cutoff | ipeps_exact_finite_density | error vs ED |
|---|--------|----------------------------|-------------|
| 2 | 1e-12 | `0.0374264910` | `-2.50e-5` |
| 2 | 1e-9  | `0.0374264910` | `-2.50e-5` |
| 3 | 1e-12 | `0.0368066348` | `-6.45e-4` |
| 3 | 1e-10 | `0.0373725084` | `-7.89e-5` |
| 3 | 1e-9  | `0.0374264910` | `-2.50e-5` |
| 3 | 1e-8  | `0.0374264910` | `-2.50e-5` |
| 4 | 1e-12 | `0.0277497509` | `-9.70e-3` |
| 4 | 1e-10 | `0.0374012639` | `-5.02e-5` |
| 4 | 1e-9  | `0.0374264910` | `-2.50e-5` |
| 4 | 1e-8  | `0.0374264910` | `-2.50e-5` |

At cutoff ≥ `1e-9`, D=2/3/4 all collapse to the same final density and match
ED within Trotter/finite-size tolerance. At cutoff `1e-12`, D=4 is dominated
by ill-conditioned directions; cf. simple-gauge diagnostics, where the
maximum diagonal condition number drops from `~2.6e4` at `1e-12` to `~14`
at `1e-9` (decision log 2026-05-17).

## Fix / next experiment

No code change required: the dense five-site PXP gate, simple-update star
split, CTM density operator, and exact-finite contraction are all consistent.

Operational recommendation, already captured in
`memory/mid_term/decision_log.md` (2026-05-17 entry "Treat Larger-D PXP
Simple-Update Cutoff As A Stability Parameter"):

- For larger-D PXP probes, sweep `cutoff ∈ {1e-12, 1e-10, 1e-9, 1e-8}` and
  treat cutoff as a stability knob, not a fixed accuracy target.
- For tiny cells (`3x3`), report `exact_density_finite` alongside CTM
  density to separate evolution error from finite-`chi` measurement error.

Reproducer command pattern (same one used to produce the cutoff sweep
artifacts under `artifacts/m3-systematic/`):

```bash
SQUAREPXP_LARGERD_N=3 \
SQUAREPXP_LARGERD_DT=0.02 \
SQUAREPXP_LARGERD_D=2,3,4 \
SQUAREPXP_LARGERD_CUTOFF=1e-12,1e-10,1e-9,1e-8 \
SQUAREPXP_LARGERD_TOTAL_TIME=0.20 \
SQUAREPXP_LARGERD_EXACT_FINITE=true \
SQUAREPXP_LARGERD_EXACT_FINITE_MAX_SITES=9 \
SQUAREPXP_LARGERD_JSON=artifacts/<stem>.json \
SQUAREPXP_LARGERD_CSV=artifacts/<stem>.csv \
julia --project=. scripts/pxp_larger_d_ed_benchmark.jl
```

## What this does *not* claim

- It does not validate finite-`chi` CTM trust. Low-`chi` runs (e.g. `chi=2,4`)
  are still expected to be rejected by `ctm_trust_status` until the trust
  policy has enough sweep points to converge — that is a separate axis from
  the cutoff stability issue resolved here.
- It does not extend to longer times or larger cells. The reproducer is
  scoped to `3x3`, `t = 0.20`; behavior past Trotter saturation or at larger
  cells must be re-validated with the same cutoff grid.

## Sources

- `prompts/pxp-d-debug-prompt.txt`
- `memory/mid_term/decision_log.md` (2026-05-17 entries)
- `artifacts/m3-systematic/threaded-ctm-trajectory-3x3-t020-dt002-D2D3D4-chi2chi4.{json,csv}`
- `artifacts/m3-systematic/cutoff-exact-3x3-t020-dt002-D2D3D4.{json,csv}`
- `src/IPEPSEvolution.jl`, `src/StarSimpleUpdate.jl`,
  `src/FiniteIPEPSObservables.jl`, `src/PEPSKitMeasurements.jl`
