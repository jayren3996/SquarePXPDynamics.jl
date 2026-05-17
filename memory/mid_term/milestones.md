# Milestones

## Completed On Current Main

- Confirmed: S0-S7 prototype pipeline exists on local `main`: dense model
  definitions, finite and periodic PEPS/iPEPS state containers, QR-reduced
  five-site star updates, deterministic Trotter evolution, simple/local and
  CTM-backed diagnostics, benchmark/reference paths, ScarFinder orchestration,
  CTM trust, and S7b gauge-conditioning readiness.
- Source: `README.md`
- Source: `src/SquarePXPDynamics.jl`
- Source: `test/runtests.jl`

- Confirmed: GPT PXP roadmap completion was merged and pushed to `main` on
  2026-05-17. It adds first-class ScarFinder measurement backends, trusted CTM
  ranking gates, explicit physics objectives, scar-oriented observables,
  candidate metadata persistence, reproducible CTMRG seeding, ED/iPEPS
  validation and convergence reports, reverse-evolution validation, projection
  semantics clarification, and PEPSKit helper compatibility guards.
- Source: `docs/superpowers/plans/2026-05-16-complete-gpt-pxp-roadmap.md`
- Source: `docs/superpowers/notes/2026-05-17-gpt-roadmap-completion.md`
- Source: `README.md`
- Source: `src/ScarFinder.jl`
- Source: `src/PXPValidation.jl`

- Confirmed: CTM diagnostics and context safety were hardened after review,
  including diagnostics metadata, stale-context guards, CTM validation sweeps,
  and README warnings.
- Source: `README.md`
- Source: `notes/2026-05-15-gpt-pro-ctm-scarfinder-revision-notes.md`
- Source: `src/PEPSKitMeasurements.jl`

- Confirmed: S0-S7 completion has been reconciled against the current
  architecture. Slice 1 adds public iPEPS helper APIs, link-weight
  normalization, and PXP `:x_plus` observable regression coverage while
  preserving the custom ITensors update plus PEPSKit measurement boundary.
- Source: `docs/superpowers/specs/2026-05-16-s0-s7-completion-design.md`
- Source: `src/SquareIPEPS.jl`
- Source: `test/test_square_ipeps.jl`
- Source: `test/test_square_ipeps_s2.jl`
- Source: `test/test_observables_evolved.jl`

- Confirmed: S0-S7 Slice 2 strengthens S3/S4 auditability by recording
  pre-update touched-link minima in `StarUpdateInfo` and model/protocol
  metadata in `EvolutionLog`.
- Source: `docs/superpowers/plans/2026-05-16-s0-s7-slice2-diagnostics.md`
- Source: `src/StarSimpleUpdate.jl`
- Source: `src/IPEPSEvolution.jl`
- Source: `test/test_star_simple_update.jl`
- Source: `test/test_ipeps_evolution.jl`

- Confirmed: S0-S7 Slice 3 completes the non-CTM ScarFinder plan by adding
  opt-in guarded simple-energy correction, correction outcome fields, and
  CSV/JSON logging of correction diagnostics while preserving default S6-lite
  behavior.
- Source: `docs/superpowers/plans/2026-05-16-s0-s7-slice3-scarfinder-energy-correction.md`
- Source: `src/ScarFinder.jl`
- Source: `test/test_scarfinder.jl`
- Source: `README.md`

- Confirmed: S0-S7 Slice 4 adds S7b CTM local bond norm diagnostics,
  `ctm_ready_for_gauge_updates`, and a transactional D=1 product/no-op
  `fix_bond_gauge!` path.
- Source: `docs/superpowers/plans/2026-05-16-s0-s7-slice4-ctm-gauge-readiness.md`
- Source: `src/CTMGaugeReadiness.jl`
- Source: `test/test_ctm_gauge_readiness.jl`
- Source: `docs/superpowers/notes/2026-05-16-s7b-gauge-fixing-handoff.md`

- Confirmed: S0-S7 Slice 5 completes S7b by adding transactional D>1
  PEPSKit bond-environment gauge conditioning with Gamma-lambda writeback,
  state-version invalidation, D=2 finite-simple tests, and extended fresh-CTM
  before/after regression coverage.
- Source: `docs/superpowers/plans/2026-05-16-s0-s7-slice5-d2-gauge-conditioning.md`
- Source: `src/CTMGaugeReadiness.jl`
- Source: `test/test_ctm_gauge_readiness.jl`

## Completed And Merged Locally

- Confirmed: `codex/infinite-tfim-benchmark` implemented the v1 TFIM benchmark
  framework with star-model abstraction, protocol-aware evolution, TFIM simple
  observables, finite schedule checks, benchmark runner, JSON/CSV writers, docs,
  and full-suite verification. The branch was merged into local `main` on
  2026-05-15.
- Source: `docs/superpowers/notes/2026-05-15-current-work-infinite-tfim-benchmark.md`
- Source: `git log` in
  `/Users/ren/.config/superpowers/worktrees/iPEPS/codex-infinite-tfim-benchmark`

- Confirmed: `codex/s0-s7-completion` was fast-forward merged into local
  `main` on 2026-05-16 after a focused S7b test rerun. The branch adds S0-S7
  reconciliation, guarded ScarFinder simple-energy correction, CTM local bond
  norm diagnostics, readiness checks, and transactional D>1 gauge conditioning.
- Source: `docs/superpowers/notes/2026-05-16-s0-s7-completion-audit.md`
- Source: `git log`
- Source: `test/test_ctm_gauge_readiness.jl`

## Future Milestones

- Confirmed: Local `main` is aligned with `origin/main` at
  `98e1ad7 Merge GPT PXP roadmap completion`.
- Open question: Add the broader Tier 2 TFIM smoke matrix across multiple
  `h/J` values, initial states, `D`, and `dt`.
- Open question: Run the first CTM-trusted ScarFinder audit campaign using the
  current simple-update plus trusted-measurement stack.
- Open question: Add full tensor snapshot persistence for ScarFinder candidates.
- Open question: Expand CTM-backed observables and design CTM-aware/full-update
  evolution.
- Source: `docs/superpowers/specs/2026-05-15-infinite-tfim-benchmark-design.md`
- Source: `docs/superpowers/notes/2026-05-15-current-work-infinite-tfim-benchmark.md`
- Source: `docs/superpowers/notes/2026-05-16-s0-s7-completion-audit.md`
- Source: `docs/superpowers/notes/2026-05-17-gpt-roadmap-completion.md`
- Source: `README.md`
