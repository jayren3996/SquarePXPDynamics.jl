# S7 CTM Trust And Gauge Readiness Design

Date: 2026-05-16

## Goal

State S7 should make CTMRG-backed measurements trustworthy enough to support
future gauge-conditioned updates. This design covers the first implementation
slice, S7a: CTMRG trust hardening, finite-chi validation policy, D>1
gauge-invariant regressions, and read-only gauge diagnostics. Mutating
environment-based gauge fixing is deliberately deferred until S7a exit criteria
are met. S7a trust is a measurement-validation signal; it is not by itself
permission to mutate gauges.

## Current Context

The repository already has:

- A custom ITensors `SquareIPEPSState` in Gamma-lambda simple-update gauge.
- QR-reduced five-site star updates through `project_star!`.
- PEPSKit conversion and CTMRG-backed density, blockade, and PXP energy
  measurements in `src/PEPSKitMeasurements.jl`.
- `CTMRGDiagnostics`, `CTMValidationPoint`, `validate_ctm_sweep`, and CSV
  export for CTM validation records.
- Fast and extended tests split by `SQUAREPXP_EXTENDED_TESTS`.

The missing S7a layer is not raw CTMRG access. It is a policy and regression
layer that decides when CTM data is acceptable for downstream ranking or for
later gauge-fixing algorithms.

## Non-Goals

- Do not implement repo-native CTMRG.
- Do not replace PEPSKit as the CTMRG environment backend.
- Do not implement mutating `fix_bond_gauge!` in S7a.
- Do not claim CTMRG values are final physics-quality results from one chi.
- Do not compare raw tensors after gauge-changing or split-order-changing
  operations.
- Do not change ScarFinder ranking or logging semantics in S7a. Trust
  assessment is an external validation workflow until a later ScarFinder
  integration design chooses how to consume it.
- Do not export a gauge-update readiness predicate in S7a. S7b should introduce
  that only after it can validate state freshness, CTM assessment provenance,
  and local environment norm-matrix quality.

## Architecture

S7a adds two focused pieces:

1. `CTMTrust`: a small policy layer over existing CTM sweep records.
2. `GaugeDiagnostics`: read-only local diagnostics for the current
   Gamma-lambda simple-update gauge.

`PEPSKitMeasurements` remains the CTMRG adapter. It should keep doing
conversion, PEPSKit context creation, operator construction, and CTM
measurements. The new trust layer should consume its public types rather than
moving more policy logic into the adapter.

The top-level module should export only stable user-facing policy and
diagnostic APIs. Internal helpers for extracting fields, comparing sweep
points, or contracting local norm matrices should remain unexported.

## File And Module Structure

Create `src/CTMTrust.jl` as a leaf module included after
`src/PEPSKitMeasurements.jl`. It should depend on the public
`PEPSKitMeasurements` records (`PEPSKitCTMRGParams`, `CTMRGDiagnostics`,
`CTMObservableSummary`, and `CTMValidationPoint`) and standard Julia libraries
only. It should not import `PEPSKit` or `TensorKit` directly.

Create `src/GaugeDiagnostics.jl` as a leaf module included after
`src/SquareIPEPS.jl` and `src/SquareUnitCells.jl`. It should depend on the
custom ITensors iPEPS representation and local linear algebra only. It should
not depend on `PEPSKitMeasurements` or `CTMTrust`.

Update `src/SquarePXPDynamics.jl` to include both files, import stable public
symbols with `using .CTMTrust: ...` and `using .GaugeDiagnostics: ...`, and
export only those stable symbols. Do not export the submodules themselves. Add
docstrings for every exported symbol because `test/test_public_docs.jl`
requires public documentation.

Add focused test files to `test/runtests.jl`:

- `test_ctm_trust.jl`
- `test_gauge_diagnostics.jl`

## CTM Trust Policy

Add a `CTMTrustPolicy` value that defines software-level acceptance thresholds:

```julia
struct CTMTrustPolicy
    min_points::Int
    require_accepted_diagnostics::Bool
    max_density_delta::Float64
    max_blockade_delta::Float64
    max_energy_delta::Float64
    max_residual::Union{Float64,Nothing}
end
```

The validating constructor should require:

- `min_points >= 2`,
- each maximum observable delta finite and nonnegative,
- `max_residual === nothing` or a finite nonnegative residual threshold.

The default constructor should use:

```julia
CTMTrustPolicy(
    2,
    true,
    1e-3,
    1e-4,
    1e-3,
    nothing,
)
```

These defaults are regression-trust thresholds, not universal physics
thresholds. Callers doing publication-quality comparisons should pass stricter
project-specific values and still inspect finite-chi trends.

Add a `CTMTrustAssessment` record:

```julia
struct CTMTrustAssessment
    trusted::Bool
    reason::Symbol
    message::String
    compared_points::Int
    finite_chi_density_delta::Union{Float64,Nothing}
    finite_chi_blockade_delta::Union{Float64,Nothing}
    finite_chi_energy_delta::Union{Float64,Nothing}
    observed_max_residual::Union{Float64,Nothing}
end
```

Allowed `reason` values are stable public API:

- `:trusted`
- `:too_few_points`
- `:nonmonotonic_sweep`
- `:missing_diagnostics`
- `:unaccepted_diagnostics`
- `:missing_residual`
- `:residual_too_large`
- `:density_delta_too_large`
- `:blockade_delta_too_large`
- `:energy_delta_too_large`

The `CTMTrustAssessment` constructor should validate that `reason` is in this
set, `compared_points >= 0`, all present delta/residual fields are finite and
nonnegative, and `trusted == true` only appears with `reason === :trusted`.

The public function:

```julia
assess_ctm_trust(points; policy = CTMTrustPolicy())::CTMTrustAssessment
```

must preserve input order and assess the final `policy.min_points` sweep
points. Within that window, `params.chi` must be strictly increasing and
`params.tol` must be nonincreasing. Observable trust deltas are the maximum
absolute adjacent CTM-to-CTM drift in that final window:

```julia
finite_chi_density_delta =
    maximum(abs(points[i].measurement.density - points[i - 1].measurement.density)
            for i in window_indices[2:end])
finite_chi_blockade_delta =
    maximum(abs(points[i].measurement.blockade_violation -
                points[i - 1].measurement.blockade_violation)
            for i in window_indices[2:end])
finite_chi_energy_delta =
    maximum(abs(points[i].measurement.pxp_energy_density -
                points[i - 1].measurement.pxp_energy_density)
            for i in window_indices[2:end])
```

These fields must not use `CTMValidationPoint.delta_*`, because those existing
fields mean CTM-minus-simple-reference.

`assess_ctm_trust` should reject with `trusted = false` for ordinary trust
failures:

- fewer than `policy.min_points` points,
- nonmonotonic final-window `chi` or `tol`,
- missing or unaccepted diagnostics when required,
- missing residuals when `policy.max_residual` is set,
- residuals above `policy.max_residual` when that field is set,
- density, blockade, or energy deltas above policy thresholds.

The function should return a structured rejection reason instead of throwing
for ordinary non-trust outcomes. It should throw `ArgumentError` only for
malformed policy values, wrong point element types, or malformed records such
as nonfinite measurement summaries.

## Finite-Chi Sweep Behavior

The existing `validate_ctm_sweep` compares CTM outputs to a simple/local
reference. S7a should keep that behavior and add a trust assessment step over
the produced CTM measurements.

Recommended usage should become:

```julia
points = validate_ctm_sweep(psi; params = [
    PEPSKitCTMRGParams(4, 1e-6, 50, 0),
    PEPSKitCTMRGParams(8, 1e-8, 100, 0),
])
assessment = assess_ctm_trust(points)
```

CSV output should include enough information to audit trust decisions without
changing the existing validation CSV format. Add:

```julia
write_ctm_trust_csv(points, path; policy = CTMTrustPolicy())
```

The writer should compute the assessment internally from `points` and `policy`
so the CSV cannot pair points with an unrelated assessment. It should write one
row per sweep point and repeat the sweep-level assessment fields on every row
so the file remains self-contained when opened as a flat table.

Required trust CSV columns:

```text
chi,tol,maxiter,verbosity,
ctm_density,ctm_blockade_violation,ctm_pxp_energy_density,
ctm_iterations,ctm_residual,ctm_converged,ctm_accepted,
trust_policy_min_points,trust_policy_require_accepted_diagnostics,
trust_policy_max_density_delta,trust_policy_max_blockade_delta,
trust_policy_max_energy_delta,trust_policy_max_residual,
trust_trusted,trust_reason,trust_compared_points,
trust_finite_chi_density_delta,trust_finite_chi_blockade_delta,
trust_finite_chi_energy_delta,trust_observed_max_residual
```

The CSV should omit free-form assessment messages and use stable `reason`
symbols for machine-readable audit trails. If future CSV fields contain
strings, they must use CSV quoting/escaping rather than simple comma joining.

## Gauge Diagnostics

Add read-only diagnostics for the existing Gamma-lambda simple-update gauge:

```julia
struct SimpleGaugeDiagnostic
    bond::BondKey
    deviation::Float64
    frobenius_norm::Float64
    diagonal_min::Float64
    diagonal_max::Float64
    diagonal_condition::Float64
end

gauge_diagnostic_simple(psi, c, dir)::SimpleGaugeDiagnostic
gauge_deviation_simple(psi, c, dir)::Float64
all_gauge_deviations_simple(psi)::Dict{BondKey,Float64}
```

`gauge_diagnostic_simple` should require `psi.gauge === :simple`. It should
accept all four directions by canonicalizing `(c, dir)` through
`bondkey(psi.unitcell, c, dir)`. The returned `bond` field should be the
canonical `BondKey`, so `:left`/`:down` aliases report the same canonical bond
as their `:right`/`:up` counterpart.

For a canonical bond `(a, dir)` with endpoint `b = neighbor(cell, a, dir)`,
build a local D x D bond norm matrix in the current lambda basis:

1. Copy the endpoint tensors `A = psi.tensors[a]` and `B = psi.tensors[b]`.
2. Do not absorb the target bond lambda.
3. Absorb the three external simple-update link weights into each endpoint
   tensor using the same full-`lambda` convention as `absorb_link_weight`.
4. For each endpoint, contract the tensor with its conjugate over the physical
   leg and the three external virtual legs, leaving only the target bond ket
   and bra indices. This gives endpoint Gram matrices `GA` and `GB`.
5. Form the local bond norm matrix as the elementwise product
   `N = GA .* GB` after aligning both matrices to the canonical target bond
   basis.
6. Symmetrize numerically as `(N + N') / 2` before extracting diagnostics.

Let `normN = norm(N)`. If `normN == 0` or is nonfinite, throw
`ArgumentError`. The diagnostic deviation is:

```julia
norm(N - Diagonal(diag(N))) / normN
```

`diagonal_min` and `diagonal_max` are computed from `real.(diag(N))`.
`diagonal_condition` is `diagonal_max / diagonal_min` when
`diagonal_min > 0`, and `Inf` otherwise. `gauge_deviation_simple` returns
`gauge_diagnostic_simple(...).deviation`.

For D=1 product states the deviation must be exactly zero within numerical
tolerance. For deterministic D=2 seeded states with active off-diagonal virtual
sectors, at least one chosen bond must have deviation greater than `1e-10`.
For a deliberately diagonal D=2 fixture, all deviations should be zero within
tolerance.

This diagnostic is not a canonical-gauge certificate. It is an early warning
that the current lambda basis is a poor local truncation basis and that future
environment-aware gauge conditioning may matter.

## Environment Contract For Future Gauge Fixing

Document the S7b contract now, but do not implement `fix_bond_gauge!` yet.

Future mutating gauge-fixing code must:

- accept a fresh `PEPSKitMeasurementContext` or explicitly build one,
- require measurement trust plus separate local CTM norm-matrix checks before
  mutating tensors,
- mutate only after all local factorization steps have succeeded,
- increment the state version and invalidate old CTM contexts,
- compare observables and CTM summaries, not raw tensors, in tests,
- preserve or clearly update the Gamma-lambda representation invariants.

S7a should explicitly not export `ctm_ready_for_gauge_updates`. S7b should add
that predicate only after it can check all required inputs together: the current
`psi`, fresh `PEPSKitMeasurementContext`, CTM trust assessment provenance,
local CTM bond norm matrices, Hermiticity residuals, positive-semidefinite
eigenvalue floors, condition numbers/rcond, finite entries, and bond-direction
coverage.

## Error Handling

Malformed user inputs should throw `ArgumentError`:

- invalid trust thresholds,
- wrong CTM point element types or missing required fields,
- nonfinite CTM measurement summaries,
- invalid bond directions,
- unsupported state/index layouts in gauge diagnostics.

Non-converged CTMRG, empty or too-short sweep point collections, nonmonotonic
final sweep windows, excessive finite-chi deltas, and unaccepted diagnostics
are expected trust outcomes and should return
`CTMTrustAssessment(trusted = false, ...)` rather than throw.

## Testing Strategy

Fast CI tests should cover:

- `CTMTrustPolicy` validation.
- `assess_ctm_trust` accepts synthetic stable sweep records.
- `assess_ctm_trust` rejects too few points, nonmonotonic sweep windows,
  unaccepted diagnostics, missing residuals when residual policy is enabled,
  excessive residual, and excessive observable deltas.
- Synthetic cases where CTM is stable across chi but far from the simple
  reference are trusted, proving trust does not use `CTMValidationPoint.delta_*`.
- Synthetic cases where CTM is close to the simple reference but drifts across
  chi are rejected.
- CSV trust metadata round trips on synthetic records.
- `gauge_deviation_simple` is zero for D=1 product states.
- `gauge_diagnostic_simple` reports finite nonnegative fields for seeded D=2
  states, and at least one chosen off-diagonal D=2 fixture has deviation above
  `1e-10`.
- A deliberately diagonal D=2 fixture has zero gauge deviation within
  tolerance.
- Existing stale-context behavior remains unchanged, including a stale
  aggregate `pxp_energy_density_ctm` call after state mutation.

Extended tests behind `SQUAREPXP_EXTENDED_TESTS=1` should cover:

- Actual PEPSKit CTMRG sweeps on D=1 product/checkerboard states.
- At least one short-evolved D=1 state with two chi values.

D=2 PEPSKit CTMRG is not required for S7a because it can dominate local and CI
runtime. S7a D>1 coverage comes from fast gauge diagnostics and
gauge-invariant simple-observable regressions. Extended CTMRG tests should not
assert exact iteration counts or residual values because PEPSKit CTMRG starts
from a random initial environment; they should assert finite summaries,
well-formed diagnostics, and trust outcomes only for cases shown empirically
stable.

## Documentation

Update README CTM guidance so users see the trust workflow:

1. Run `validate_ctm_sweep`.
2. Run `assess_ctm_trust`.
3. Optionally write `write_ctm_trust_csv`.
4. Inspect diagnostics and finite-chi deltas before ranking or claiming energy
   comparisons.

Documentation should state clearly that S7a is gauge-readiness work and that
mutating full-update gauge fixing is S7b.

## Exit Criteria

S7a is complete when:

- Trust-policy APIs and gauge diagnostics are exported and documented.
- Fast tests pass without requiring extended CTMRG solves.
- Extended CTM tests pass manually with `SQUAREPXP_EXTENDED_TESTS=1`.
- README shows the new trust workflow.
- The repository has
  `docs/superpowers/notes/2026-05-16-s7b-gauge-fixing-handoff.md` with these
  sections: prerequisites from S7a, local CTM norm-matrix requirements,
  mutation contract, required gauge-invariant tests, and known blockers for
  mutating `fix_bond_gauge!`.
