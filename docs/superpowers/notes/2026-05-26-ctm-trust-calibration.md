# CTM trust calibration framework

Purpose: turn `CTMTrustPolicy` from a placeholder defaults object into a
data-derived gate. Without calibrated thresholds, `assess_ctm_trust` only
catches gross failures (`min_points`, `nonmonotonic_sweep`,
`unaccepted_diagnostics`) and the numerical thresholds
(`max_density_delta = 1e-3`, `max_blockade_delta = 1e-4`,
`max_energy_delta = 1e-3`) are arbitrary.

## Available trust policies

| Constructor                              | Use when                                                                                          |
|------------------------------------------|---------------------------------------------------------------------------------------------------|
| `CTMTrustPolicy()`                       | Smoke testing, regression coverage. Conservative defaults — passes most stable iPEPS states.      |
| `tight_ctm_trust_policy()`               | `D >= 4`, where iPEPS truncation is `< 1e-5` and CTM error must not dominate.                     |
| `calibrated_ctm_trust_policy(window; …)` | Data-derived. Reads observed adjacent χ-drift from a sensitivity sweep and applies a safety margin. |

`tight_ctm_trust_policy()` returns `CTMTrustPolicy(3, true, 1e-5, 1e-6, 1e-5, 1e-6)`.
The `min_points=3` raise means at least three increasing-χ measurements are
compared, which gives a more stable drift estimate than the default 2.

`calibrated_ctm_trust_policy(window; safety=3.0, floor=1e-8)` reads adjacent
drift in `density`, `blockade_violation`, and `pxp_energy_density` from
`window` (a vector of `CTMValidationPoint`s or `CTMObservableSummary`s) and
sets each threshold to `max(safety * observed_max_adjacent_delta, floor)`. The
default `safety=3` gives roughly an order-of-magnitude margin above measured
drift; the `floor` prevents collapse to zero when the window is too short or
too easy.

## Calibration procedure

1. Pick a representative state. For the ScarFinder context this means: a
   PXP-evolved iPEPS at the `D` and approximate evolution time you intend to
   sample (e.g. `D=2, evolve_time=0.5` from `:z_up`). A near-product state
   gives a floor-dominated calibration that is *not* representative.
2. Run `scripts/run_ctm_chi_sensitivity.jl` with χ values spanning at least
   one octave (e.g. `8,16,32,64`).
3. Read `calibrated_policy.max_*_delta` from the JSON output.
4. Compare against the iPEPS truncation error you expect at that `D`. The
   CTM trust thresholds should be **tighter** than truncation, otherwise CTM
   error dominates the audit budget.
5. Either adopt the calibrated thresholds for the audit campaign, or pick a
   compromise constant policy that is consistent across all states the audit
   will visit.

## Smoke validation (2026-05-26)

`scripts/run_ctm_chi_sensitivity.jl` was smoke-run at:
- cell `3x3`, `D=2`, `evolve_time=0.02`, χ ∈ `{4, 8, 12, 16}`
- output: `artifacts/ctm_chi_sensitivity_smoke.csv`,
  `artifacts/ctm_chi_sensitivity_smoke.json`

Result: density `≈ 3.996e-4` was χ-constant to ~1e-16; blockade and energy
were at machine precision. All three policies (default, tight, calibrated)
trusted the sweep. **These numbers are not a calibration** — at
`evolve_time = 0.02` the state is still essentially `:down`, so χ-drift is
unmeasurable. They confirm the harness runs end-to-end. A real calibration
must be repeated on a longer-evolved state on the full server.

## Recommended production calibration runs

After Phase 1 throughput tuning lands, run the sensitivity sweep on the
states the Phase 5 audit will actually rank:

```bash
JULIA_NUM_THREADS=42 \
SQUAREPXP_CTM_BLAS_THREADS=1 \
SQUAREPXP_CTM_STRIDED_THREADS=42 \
SQUAREPXP_CTM_STRIDED_THREADED_MUL=true \
SQUAREPXP_CTM_PEPSKIT_SCHEDULER=dynamic \
SQUAREPXP_SENS_CELL_LX=3 SQUAREPXP_SENS_CELL_LY=3 \
SQUAREPXP_SENS_D=2 SQUAREPXP_SENS_EVOLVE_TIME=0.5 \
SQUAREPXP_SENS_CHI_VALUES=8,16,32,64 \
SQUAREPXP_SENS_LABEL=d2-t05 \
SQUAREPXP_SENS_OUTPUT_CSV=artifacts/ctm-sensitivity-d2-t05.csv \
SQUAREPXP_SENS_OUTPUT_JSON=artifacts/ctm-sensitivity-d2-t05.json \
julia --project=. scripts/run_ctm_chi_sensitivity.jl
```

Repeat for `D=3` and `D=4`. Record observed drifts and adopted thresholds
under "Selected thresholds" below before moving to Phase 5.

## Selected thresholds

### D = 2, evolve_time = 0.5 (cell 3x3, `:down` start)

Production calibration run on 2026-05-26 (`artifacts/v1/ctm_sensitivity_d2_t05.json`):

```
chi=8   density=1.715412e-01  blockade=1.019e-11  energy=-1.48e-15
chi=16  density=1.715412e-01  blockade=1.019e-11  energy=-1.51e-15
chi=32  density=1.715412e-01  blockade=1.019e-11  energy=-1.57e-15
chi=64  density=1.715412e-01  blockade=1.019e-11  energy=-1.43e-15
```

Observed adjacent χ-drift was below the calibration floor (1e-8) on all three
observables. `default`, `tight`, and `calibrated` policies all return
`trusted = true`. **At `D=2, t=0.5` on a 3x3 cell, CTMRG converges by χ=8.**

Conclusion: for `D = 2` audit runs in this regime, use
`tight_ctm_trust_policy()` — its thresholds (1e-5, 1e-6, 1e-5) give physically
meaningful gating that will catch drift if a different state in the audit grid
*isn't* fully converged, while staying loose enough to pass converged points
without false rejections.

The `calibrated_ctm_trust_policy(window; safety=3, floor=1e-8)` output is too
tight (all thresholds at floor) — adopting it would risk false rejections in
audit points where chi-drift is non-negligible. Use the calibrated policy as a
**lower bound check** only.

### D = 1 (any state)

Trivial — at D=1 the iPEPS is a product state and CTMRG is exact. Use
`CTMTrustPolicy()`.

### D >= 3 (pending)

Re-run `scripts/run_ctm_chi_sensitivity.jl` with `SQUAREPXP_SENS_D=3` and a
representative evolve_time before launching D=3 audits. If observed drift is
similar to D=2, `tight_ctm_trust_policy()` remains the recommendation.
