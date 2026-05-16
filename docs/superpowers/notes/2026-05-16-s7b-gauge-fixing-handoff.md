# S7b Gauge-Fixing Handoff

Date: 2026-05-16

## Prerequisites From S7a

- `assess_ctm_trust(points)` must be available for finite-chi measurement
  validation.
- `gauge_diagnostic_simple(psi, c, dir)` must be available for read-only local
  simple-gauge diagnostics.
- CTM contexts must reject stale state usage through state id/version checks.
- S7a trust remains a measurement-validation signal. It is not sufficient by
  itself to authorize mutating gauge updates.

## Local CTM Norm-Matrix Requirements

Future `fix_bond_gauge!` work must define and test environment-backed local
bond norm matrices before mutating tensors. Required diagnostics:

- finite entries,
- Hermiticity residual,
- positive-semidefinite eigenvalue floor,
- condition number or reciprocal condition number,
- bond-direction coverage over canonical `:right` and `:up` links,
- clear failure behavior when the local environment is singular or indefinite.

## Mutation Contract

Future gauge-fixing code must:

- accept a fresh `PEPSKitMeasurementContext` or build one explicitly,
- verify measurement trust and local norm-matrix quality before mutation,
- perform all local factorization steps before writing into `psi`,
- increment `state_version(psi)` after mutation,
- invalidate old CTM contexts by relying on the existing state-version guard,
- preserve Gamma-lambda invariants or document any representation change.

## Required Gauge-Invariant Tests

S7b tests must compare observables and CTM summaries, not raw tensor entries.
Required checks:

- D=1 product states remain unchanged under any no-op gauge path.
- D=2 seeded states keep finite simple observables after gauge conditioning.
- Fresh CTM context measurements before/after a pure gauge change agree within
  documented tolerance.
- Old CTM contexts throw after gauge mutation.
- Singular or ill-conditioned norm matrices fail without partially mutating
  the state.

## Known Blockers For Mutating `fix_bond_gauge!`

- S7a does not compute CTM local bond norm matrices.
- S7a does not define whitening or ALS/full-update truncation factors.
- PEPSKit environment internals needed for local norm matrices still need a
  focused feasibility check.
- The project has not chosen whether `fix_bond_gauge!` should update only the
  current Gamma tensors or introduce a richer gauge state.
