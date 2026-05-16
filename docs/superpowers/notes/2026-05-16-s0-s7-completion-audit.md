# S0-S7 Completion Audit

Date: 2026-05-16

Branch: `codex/s0-s7-completion`

## Result

S0-S7 is complete for the current project architecture. The original backend
facade and broad measurement facade items are explicitly superseded by the
custom ITensors iPEPS update path plus PEPSKit CTM measurement/readiness
adapters. Production ScarFinder validation remains outside this completion
claim.

## Evidence Map

- S0 baseline conventions are locked by `src/SquarePXP.jl`,
  `src/SquareIPEPS.jl`, `src/SquarePEPS.jl`, and regression tests under
  `test/test_square_pxp.jl`, `test/test_square_ipeps.jl`, and
  `test/test_square_peps.jl`.
- S0.5 is complete by decision: `AbstractProjectionBackend` and
  `projection_backend` are superseded until a second update backend exists.
  See `docs/superpowers/specs/2026-05-16-s0-s7-completion-design.md` and
  `memory/mid_term/decision_log.md`.
- S1 periodic unit-cell and iPEPS state helpers are shipped in
  `src/SquareIPEPS.jl` and exported from `src/SquarePXPDynamics.jl`.
  Coverage: `test/test_square_ipeps.jl`.
- S2 gate/link-weight invariants are shipped through ITensor gate wrappers,
  `normalize_link_weights!`, entropy helpers, and `:x_plus` product-state
  support. Coverage: `test/test_square_ipeps_s2.jl`,
  `test/test_observables_evolved.jl`, and related square-state tests.
- S3 QR-reduced five-site star simple update is complete with
  pre-update touched-link minima and post-split lambda diagnostics in
  `src/StarSimpleUpdate.jl`. Coverage: `test/test_star_simple_update.jl`.
- S4 deterministic iPEPS evolution is complete with model/protocol metadata
  and log-normalization diagnostics in `src/IPEPSEvolution.jl`.
  Coverage: `test/test_ipeps_evolution.jl`.
- S5 observables are complete through explicit simple and CTM APIs:
  `measure_simple`, `measure_ctm`, CTM validation sweeps, and trust
  assessment. The broad `measure`/environment facade is intentionally not
  shipped. Coverage: `test/test_observables*.jl`,
  `test/test_pepskit_measurements.jl`, and `test/test_ctm_trust.jl`.
- S6 ScarFinder is complete for the non-production pipeline with opt-in
  guarded simple-energy correction, transactional rejection of non-improving
  correction attempts, and CSV/JSON diagnostics in `src/ScarFinder.jl`.
  Coverage: `test/test_scarfinder.jl`.
- S7 is complete through S7a CTM trust plus S7b CTM local bond norm matrices,
  readiness checks, D=1 transactional no-op behavior, and D>1 PEPSKit
  bond-environment gauge conditioning in `src/CTMTrust.jl`,
  `src/GaugeDiagnostics.jl`, and `src/CTMGaugeReadiness.jl`.
  Coverage: `test/test_ctm_trust.jl`, `test/test_pepskit_measurements.jl`,
  and `test/test_ctm_gauge_readiness.jl`.

## Verification

- `julia --project=. test/runtests.jl test_ctm_gauge_readiness.jl`
  passed: 98/98.
- `SQUAREPXP_EXTENDED_TESTS=1 julia --project=. test/runtests.jl test_ctm_gauge_readiness.jl`
  passed: 102/102.
- `julia --project=. -e 'using Pkg; Pkg.test()'`
  passed: 65149/65149.
- `git diff --check` passed before the Slice 5 code commit and after the
  final audit documentation changes.
