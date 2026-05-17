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

## Exact Observable Boundary

Added test-only exact finite observable helpers, later promoted as public
`_finite` APIs:

- `exact_one_site_expectation_finite(psi, c, O)`;
- `exact_nearest_neighbor_expectation_finite(psi, c, dir, O)`;
- `exact_star_expectation_finite(psi, center, O)`.

Each helper contracts the full `3 x 3` network through the dense-vector path
after absorbing each periodic link weight exactly once.

For the first divergent D=2 state after star `3`, the first per-site one-site
density mismatch is:

- site: `SquareCoord(2, 1)`
- exact finite `<n>`: `0.0003997867121725804`
- `local_density_simple`: `0.0001999133387617944`
- simple minus exact: `-0.000199873373410786`

The average exact finite one-site density remains
`0.00013326224449912612`, matching the dense serial-star reference, while the
average `density_simple` is `0.000111054099003352`.

Nearest-neighbor blockade density `<n_i n_j>` stays zero within tolerance on
canonical bonds for this state, so it does not expose the boundary. A generic
two-site ZZ observable does:

- bond: `SquareCoord(1, 1)` `:right`
- exact finite `<ZZ>`: `0.9984005332366328`
- simple local `<ZZ>`: `0.9988001200421318`
- simple minus exact: `0.00039958680549900816`

The first star-patch mismatch found with a center-density operator embedded in
star order `(center, right, up, left, down)` is:

- center: `SquareCoord(1, 1)`
- exact finite star-center density: `0.00039994666951102603`
- `star_expectation_simple` value: `0.00040006586165461806`
- simple minus exact: `1.1919214359202906e-7`

The square-PXP star Hamiltonian expectation itself still agrees within the
diagnostic tolerance on this state; it is less sensitive than the local density
and ZZ probes here.

## Conclusion

Confirmed: the D=2 QR/SVD star-update path is not the source of this tiny-run
density anomaly as measured by exact finite contraction. The D=2 state produced
by grow-on-demand serial `project_star!` matches the exact dense serial-star
circuit through all nine stars to tight tolerance.

Confirmed: the no-CTM audit anomaly is in treating D>1 simple/local
measurements as exact finite `3 x 3` observables after overlapping D=2 serial
updates. The exact observable helpers show mismatches in one-site, generic
two-site, and star-patch probes. This is expected for simple/local environments
on a loopy finite PEPS and is not evidence of a concrete lambda-counting bug in
`Observables.jl`.

No production observable formula was changed in this slice. The regression
keeps exact simple-environment equality expectations as `@test_broken`, and
the README / validation docstrings now state that D>1 no-CTM simple-density
errors are local diagnostic offsets, not exact finite-PEPS validation errors.

Postponed:

- CTM Stage 2;
- CTM-aware/full-update design;
- new CTM observables;
- tensor persistence.

## Follow-Up Implementation Plan

Plan `docs/superpowers/plans/2026-05-17-d2-exact-finite-observables.md`
promotes the test-only exact finite helpers into a size-limited module and
wires exact finite density into PXP validation/audit as an opt-in field. The
plan keeps simple/local observable formulas unchanged.
