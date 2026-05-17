# GPT PXP Roadmap Completion

## Summary

Confirmed: The GPT review follow-up roadmap was implemented, merged to
`main`, verified, pushed to GitHub, and branch-cleaned on 2026-05-17.

Source: `git log`; `README.md`;
`docs/superpowers/plans/2026-05-16-complete-gpt-pxp-roadmap.md`

## Completed

- First-class ScarFinder measurement backends:
  `MeasurementBackend`, `SimpleBackend`, `TrustedCTMBackend`, and
  `measure_scarfinder`.
- Trusted ScarFinder ranking path:
  `scarfinder!` accepts `measurement = TrustedCTMBackend(...)` and
  `require_trusted_ctm = true`.
- Physics-objective API:
  `RevivalObjective`, `TargetEnergyObjective`, `LowVarianceObjective`, and
  `CompositeObjective`.
- Scar-oriented observables:
  sublattice imbalance and checkerboard structure factor in simple and CTM
  summaries.
- Reproducible CTMRG initialization:
  `PEPSKitCTMRGParams(...; seed = ...)` with seed metadata in validation
  outputs.
- PXP ED/iPEPS validation and convergence reports:
  `PXPValidationConfig`, `PXPConvergenceConfig`,
  `validate_pxp_ed_ipeps`, `validate_pxp_convergence`, and JSON writers.
- Real-time reversibility diagnostics:
  `reverse_evolve!`, `PXPReversibilityReport`, and
  `validate_pxp_reversibility`.
- Projection semantics clarification:
  `projected_square_pxp_gate(...; projection = :left | :sandwich)` with tests
  pinning default `P * U`, explicit `P * U * P`, and equivalence for the
  current square-star PXP Hamiltonian.
- ScarFinder candidate metadata persistence:
  `CandidateStore`, `NoCandidateStore`, and `JSONCandidateStore`.
- PEPSKit helper compatibility wrapper:
  `pepskit_private_full_update_available()` and guarded use of PEPSKit
  full-update helper names.
- README and memory decision-log updates for the trusted ScarFinder workflow.

## Verification

- Feature-branch verification:
  `julia --project=. -e 'using Pkg; Pkg.test()'` passed `65356/65356` in
  `7m28.3s`.
- Post-merge verification on `main`:
  `julia --project=. -e 'using Pkg; Pkg.test()'` passed `65356/65356` in
  `4m57.8s`.
- `git diff --check HEAD~1 HEAD` passed.
- Pushed head:
  `98e1ad7 Merge GPT PXP roadmap completion` on `origin/main`.

## Still Not Shipped

- CTM-aware/full-update evolution.
- Full tensor snapshot persistence for ScarFinder candidates.
- Publication-grade physics audit sweeps across `dt`, `D`, `chi`, cutoff, unit
  cell, and update scheme.
- True CTM-compatible many-body return/fidelity proxy.
- Richer CTM observables such as two-point correlations, transfer-matrix
  correlation length, and energy-variance-quality diagnostics.
- A large production ScarFinder campaign.

## Recommended Next Choice

Start with a small CTM-trusted ScarFinder audit campaign using the current
simple-update plus trusted-measurement stack. This tests the new control flow
before investing in full tensor snapshot persistence or CTM-aware/full-update
mutation.
