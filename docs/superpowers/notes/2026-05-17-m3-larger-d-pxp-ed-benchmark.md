# M3 Larger-D PXP ED Benchmark

Date: 2026-05-17

## Contract

The M3 benchmark compares iPEPS dynamics against finite PBC ED through
observable modes that are explicitly labeled in JSON and CSV artifacts.

- `exact_finite`: tiny-cell exact finite contraction of the current iPEPS state.
- `simple_diagnostic`: simple/local environment diagnostic, not exact for D>1.
- `ctm_trusted`: CTM-backed value only when finite-chi trust is attached.
- `symmetric_pbc_ed_global`: finite PBC ED global site-averaged density and ED
  return probability in the selected symmetric sector.

The symmetric PBC ED path does not provide central 3x3 observables. PBC has no
physical center, and the current ED basis is fully symmetry reduced by default.

## Default Fast Campaign

```bash
SQUAREPXP_LARGERD_N=3 \
SQUAREPXP_LARGERD_DT=0.02 \
SQUAREPXP_LARGERD_D=1,2,3,4 \
SQUAREPXP_LARGERD_CUTOFF=1e-12 \
SQUAREPXP_LARGERD_TOTAL_TIME=0.02 \
SQUAREPXP_LARGERD_EXACT_FINITE=true \
SQUAREPXP_LARGERD_EXACT_FINITE_MAX_SITES=9 \
julia --project=. scripts/pxp_larger_d_ed_benchmark.jl
```

## Larger Odd PBC Campaign

Use `n = 5` for the first larger odd benchmark. Use `n = 7` only as a manual
capacity/runtime boundary probe.

## Postponed

Unreduced PBC 5x5 ED can support local operators, but a central region is still
not physically privileged under PBC. Open-boundary ED is the cleaner route for
literal central 3x3 observables and should be planned as a separate milestone.
