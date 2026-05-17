# D2 Measurement-Path Localization

Date: 2026-05-17

## Question

Does the `3 x 3`, all-down, `t = 0.02`, `dt = 0.02`, serial PXP D=2 anomaly
come from the D=2 star update path or from the D>1 simple/local measurement
path?

## Harness

Added `test/test_pxp_d2_localization.jl` with test-only helpers:

- exact dense `2^9` serial-star PXP reference in `PeriodicSquareUnitCell(3, 3)`
  representative order `[SquareCoord(x, y) for y = 1:3 for x = 1:3]`;
- exact finite contraction of a `3 x 3` `SquareIPEPSState` to a dense vector by
  absorbing every periodic link weight once and contracting all site tensors;
- per-star D=1/D=2 trace comparing dense reference density, dense-contracted
  iPEPS density, `measure_simple` density, `log_norm`, `norm_factors`,
  `keptdims`, `min_lambda`, and truncation diagnostics.

The exact dense serial-star final density is
`0.00039962698926202146`, matching the independent diagnostic reference.

## Red Regression

The initial regression required both D=2 dense-contracted density and
`measure_simple` density to match the dense reference after every serial star.
It failed at the simple-measurement assertion:

- first dense-contracted D=2 divergence: `nothing`
- first D=2 `measure_simple` versus dense-contracted divergence: star `3`

## First Divergent Star

At star `3`, center `SquareCoord(3, 1)`:

- dense reference density: `0.00013326224449912612`
- dense-contracted D=2 iPEPS density: `0.00013326224449912625`
- D=2 `measure_simple` density: `0.000111054099003352`
- D=2 `log_norm`: `1.7328679513998648`
- D=2 `keptdims`: `Dict(:left => 2, :right => 2, :up => 1, :down => 1)`
- D=2 `norm_factors`:
  `Dict(:left => 1.4142135623730956, :right => 1.0000000000000004, :up => 1.4142135623730956, :down => 2.000000000000001)`
- D=2 `min_lambda`:
  `Dict(:left => 0.00999883273437467, :right => 0.0003997867841006291, :up => 1.0, :down => 1.0)`
- D=2 max truncation error: `5.9362248149573574e-30`

## Final Star

At star `9`, center `SquareCoord(3, 3)`:

- dense reference density: `0.0003996269892620213`
- dense-contracted D=2 iPEPS density: `0.0003996269892620211`
- D=2 `measure_simple` density: `0.00021094978264193137`
- D=2 `log_norm`: `7.278045395879476`
- D=2 max truncation error: `3.651022102560827e-29`

## Conclusion

Confirmed: the D=2 QR/SVD star-update path is not the source of this tiny-run
density anomaly as measured by exact finite contraction. The D=2 state produced
by grow-on-demand serial `project_star!` matches the exact dense serial-star
circuit through all nine stars to tight tolerance.

Confirmed: the no-CTM audit anomaly is in the D>1 simple/local measurement
path or, more precisely, in treating `density_simple` as an exact finite `3 x 3`
observable after overlapping D=2 serial updates. The current `Observables.jl`
documentation already describes simple/local observables as cheap
simple-update local-environment diagnostics, not CTMRG-quality measurements.

No production observable was changed in this slice. The regression keeps the
exact `density_simple == dense contraction` expectation as `@test_broken`, so
future work can either make the intended exactness contract true or keep audit
interpretation conservative.

Postponed:

- CTM Stage 2;
- CTM-aware/full-update design;
- new CTM observables;
- tensor persistence.
