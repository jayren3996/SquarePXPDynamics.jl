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

Slice 4 defines and tests environment-backed local bond norm diagnostics before
mutating tensors. Required diagnostics:

- finite entries,
- Hermiticity residual,
- positive-semidefinite eigenvalue floor,
- condition number or reciprocal condition number,
- bond-direction coverage over canonical `:right` and `:up` links,
- clear failure behavior when the local environment is singular or indefinite.

## Mutation Contract

The `fix_bond_gauge!` entry point enforces freshness, trust, and norm
diagnostics before mutation. The D=1 product path is a no-op; D>1 uses PEPSKit
bond-environment gauge factorization, then writes the two conditioned tensors
back into the custom ITensors Gamma-lambda state after deabsorbing stored link
weights.

- It accepts a fresh `PEPSKitMeasurementContext`.
- It verifies measurement trust and local norm-matrix quality before mutation.
- It performs all local factorization and Gamma conversion steps before writing
  into `psi`.
- It increments `state_version(psi)` after D>1 mutation.
- It invalidates old CTM contexts by relying on the existing state-version
  guard.
- It preserves the Gamma-lambda representation with the existing link weights.

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

## Remaining Follow-Ups After S7b

- The current D>1 path conditions the existing two tensors on a bond; it does
  not add a full ALS truncation/update solver.
- Gauge conditioning currently requires positive link weights on all legs of
  the two tensors being converted back from PEPSKit.
- Production ScarFinder validation still needs physics-facing CTM workflows on
  top of the completed S0-S7 infrastructure.
