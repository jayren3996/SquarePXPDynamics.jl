# S7 CTM Trust And Gauge Readiness Design

Date: 2026-05-16

## Goal

State S7 should make CTMRG-backed measurements trustworthy enough to support
future gauge-conditioned updates. This design covers the first implementation
slice, S7a: CTMRG trust hardening, finite-chi validation policy, D>1
gauge-invariant regressions, and read-only gauge diagnostics. Mutating
environment-based gauge fixing is deliberately deferred until S7a exit criteria
are met.

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
    density_delta::Union{Float64,Nothing}
    blockade_delta::Union{Float64,Nothing}
    energy_delta::Union{Float64,Nothing}
    max_residual::Union{Float64,Nothing}
end
```

The public function:

```julia
assess_ctm_trust(points; policy = CTMTrustPolicy())::CTMTrustAssessment
```

must compare the final two sweep points after preserving input order. It should
reject:

- fewer than `policy.min_points` points,
- nonfinite CTM summaries,
- missing or unaccepted diagnostics when required,
- residuals above `policy.max_residual` when that field is set,
- density, blockade, or energy deltas above policy thresholds.

The function should return a structured rejection reason instead of throwing
for ordinary non-trust outcomes. It should throw only for malformed policy
values or impossible records such as nonfinite summary fields.

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
write_ctm_trust_csv(points, assessment, path)
```

The trust CSV should write one row per sweep point and repeat the sweep-level
assessment fields on every row so the file remains self-contained when opened
as a flat table.

## Gauge Diagnostics

Add read-only diagnostics for the existing Gamma-lambda simple-update gauge:

```julia
gauge_deviation_simple(psi, c, dir)::Float64
all_gauge_deviations_simple(psi)::Dict{BondKey,Float64}
```

`gauge_deviation_simple` should build a local two-site bond norm matrix in the
current lambda basis, normalize it by its Frobenius norm, and report the
relative Frobenius norm of its off-diagonal part. For D=1 product states the
deviation must be exactly zero within numerical tolerance. For D>1 seeded
states it must be finite and nonnegative.

This diagnostic is not a canonical-gauge certificate. It is an early warning
that the current lambda basis is a poor local truncation basis and that future
environment-aware gauge conditioning may matter.

## Environment Contract For Future Gauge Fixing

Document the S7b contract now, but do not implement `fix_bond_gauge!` yet.

Future mutating gauge-fixing code must:

- accept a fresh `PEPSKitMeasurementContext` or explicitly build one,
- require a trusted CTM assessment before mutating tensors,
- mutate only after all local factorization steps have succeeded,
- increment the state version and invalidate old CTM contexts,
- compare observables and CTM summaries, not raw tensors, in tests,
- preserve or clearly update the Gamma-lambda representation invariants.

Add an exported readiness helper:

```julia
ctm_ready_for_gauge_updates(assessment::CTMTrustAssessment)::Bool
```

This function should be equivalent to `assessment.trusted` and exist mainly to
make downstream S7b call sites readable.

## Error Handling

Malformed user inputs should throw `ArgumentError`:

- invalid trust thresholds,
- empty or malformed CTM point collections where a trust assessment cannot be
  constructed,
- invalid bond directions,
- unsupported state/index layouts in gauge diagnostics.

Non-converged CTMRG, too few sweep points, excessive finite-chi deltas, and
unaccepted diagnostics are expected outcomes and should return
`CTMTrustAssessment(trusted = false, ...)` rather than throw.

## Testing Strategy

Fast CI tests should cover:

- `CTMTrustPolicy` validation.
- `assess_ctm_trust` accepts synthetic stable sweep records.
- `assess_ctm_trust` rejects too few points, unaccepted diagnostics, excessive
  residual, and excessive observable deltas.
- CSV trust metadata round trips on synthetic records.
- `gauge_deviation_simple` is zero for D=1 product states.
- `gauge_deviation_simple` is finite and nonnegative for seeded D=2 states.
- Existing stale-context behavior remains unchanged.

Extended tests behind `SQUAREPXP_EXTENDED_TESTS=1` should cover:

- Actual PEPSKit CTMRG sweeps on D=1 product/checkerboard states.
- At least one short-evolved D=1 state with two chi values.

D=2 PEPSKit CTMRG is not required for S7a because it can dominate local and CI
runtime. S7a D>1 coverage comes from fast gauge diagnostics and
gauge-invariant simple-observable regressions.

## Documentation

Update README CTM guidance so users see the trust workflow:

1. Run `validate_ctm_sweep`.
2. Run `assess_ctm_trust`.
3. Inspect diagnostics and finite-chi deltas before ranking or claiming energy
   comparisons.

Documentation should state clearly that S7a is gauge-readiness work and that
mutating full-update gauge fixing is S7b.

## Exit Criteria

S7a is complete when:

- Trust-policy APIs and gauge diagnostics are exported and documented.
- Fast tests pass without requiring extended CTMRG solves.
- Extended CTM tests pass manually with `SQUAREPXP_EXTENDED_TESTS=1`.
- README shows the new trust workflow.
- The repository has a written S7b follow-up note or design pointer for
  mutating `fix_bond_gauge!`.
