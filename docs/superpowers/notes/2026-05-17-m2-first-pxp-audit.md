# M2 First PXP Audit Artifact

Source artifacts:

- `artifacts/pxp_audit_noctm.csv`
- `artifacts/pxp_audit_noctm.json`

Requested Stage 2 CTM artifacts were intentionally not generated because the
Stage 1 no-CTM gate was not clean. The no-CTM audit already exposed a
short-time D=2/simple-update pathology, so running the CTM grid with
`D_values = [1, 2]` would mix finite-chi diagnostics with a pre-CTM
simple/local diagnostic boundary.

Update: follow-up localization in
`docs/superpowers/notes/2026-05-17-d2-measurement-localization.md` showed that
the D=2 star update itself matches the exact dense serial-star reference under
full finite contraction. The no-CTM density anomaly comes from treating
`density_simple` as an exact finite observable on a loopy D>1 PEPS.

## Command

```bash
SQUAREPXP_AUDIT_N=3 \
SQUAREPXP_AUDIT_DT=0.02,0.01,0.005 \
SQUAREPXP_AUDIT_D=1,2 \
SQUAREPXP_AUDIT_CUTOFF=1e-10,1e-12 \
SQUAREPXP_AUDIT_TOTAL_TIME=0.02 \
SQUAREPXP_AUDIT_JSON=artifacts/pxp_audit_noctm.json \
SQUAREPXP_AUDIT_CSV=artifacts/pxp_audit_noctm.csv \
julia --project=. scripts/pxp_audit_campaign.jl
```

The first attempt failed before running the audit because the fresh worktree had
not been instantiated and `EDKit` was unavailable. Running
`julia --project=. -e 'using Pkg; Pkg.instantiate()'` fixed the environment;
the same audit command then completed.

## Summary

The `3 x 3`, all-down, `total_time = 0.02` audit produced 12 no-CTM rows:
three `dt` values, two bond dimensions, and two cutoffs. The two cutoff values
produced identical rows at this scale.

Key values:

| D | dt | max density error | max truncerr | log-norm abs delta | reversibility density drift |
|---|---:|---:|---:|---:|---:|
| 1 | 0.02 | 7.436e-7 | 1.599e-7 | 1.077e-5 | 1.015e-9 |
| 1 | 0.01 | 1.299e-6 | 1.596e-7 | 2.201e-5 | 4.060e-9 |
| 1 | 0.005 | 2.544e-6 | 1.587e-7 | 6.122e-5 | 1.599e-8 |
| 2 | 0.02 | 1.888e-4 | 1.058e-28 | 7.278e0 | 1.337e-31 |
| 2 | 0.01 | 2.443e-4 | 1.727e-15 | 2.599e1 | 3.433e-20 |
| 2 | 0.005 | 2.442e-4 | 1.077e-15 | 6.342e1 | 2.570e-20 |

## Bottleneck Classification

Superseded classification: this note originally called the dominant bottleneck
D=2 update/normalization behavior in the simple-update path. Follow-up exact
finite contractions localized the issue more narrowly to the D>1 no-CTM
simple/local measurement interpretation, not the QR/SVD star-update path, CTM
finite-chi, or reversibility.

Evidence:

- D=1 tracks ED density closely through `t = 0.02`: final-density errors remain
  between about `7.4e-7` and `2.5e-6`.
- D=2 is much worse: final-density error is about `1.9e-4` to `2.4e-4`, roughly
  two orders of magnitude larger than D=1.
- D=2 truncation error is tiny (`~1e-28` to `~1e-15`), so the problem is not
  ordinary truncation pressure.
- D=2 log-norm deltas are large (`7.28`, `25.99`, `63.42`) and grow as `dt`
  is refined because more update intervals are applied, but a prior scalar
  normalization trial did not change the density anomaly.
- Reversibility drift is tiny for D=2, so the forward/reverse path is internally
  round-tripping despite producing a poor ED density match.

## Follow-Up

Before using no-CTM D>1 audit density as a finite validation signal, use exact
finite contractions for tiny cells or trusted CTM observables for production
diagnostics. The current `max_abs_density_error_simple` field should be read as
a simple/local diagnostic offset, not as an exact finite PEPS density error.

Do not start CTM-aware/full-update design yet. The first bottleneck is a
short-time D=2 audit interpretation issue in the existing no-CTM simple
measurement stack.
