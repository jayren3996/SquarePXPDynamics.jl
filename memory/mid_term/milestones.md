# Milestones

## Completed On Current Main

- Confirmed: S0-S6 prototype pipeline exists: dense model definitions, finite
  and periodic PEPS/iPEPS state containers, QR-reduced five-site star updates,
  deterministic Trotter evolution, simple/local observables, experimental CTM
  measurement hooks, and ScarFinder-lite orchestration.
- Source: `README.md`
- Source: `src/SquarePXPDynamics.jl`
- Source: `test/runtests.jl`

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

## Completed And Merged Locally

- Confirmed: `codex/infinite-tfim-benchmark` implemented the v1 TFIM benchmark
  framework with star-model abstraction, protocol-aware evolution, TFIM simple
  observables, finite schedule checks, benchmark runner, JSON/CSV writers, docs,
  and full-suite verification. The branch was merged into local `main` on
  2026-05-15.
- Source: `docs/superpowers/notes/2026-05-15-current-work-infinite-tfim-benchmark.md`
- Source: `git log` in
  `/Users/ren/.config/superpowers/worktrees/iPEPS/codex-infinite-tfim-benchmark`

## Future Milestones

- Confirmed: The local integration path chosen on 2026-05-15 was direct merge
  into `main`, followed by verification and push to GitHub.
- Open question: Add a full finite Hilbert-space TFIM schedule reference if it
  is required beyond the current coefficient/non-overlap/mapping tests.
- Open question: Add the broader Tier 2 TFIM smoke matrix across multiple
  `h/J` values, initial states, `D`, and `dt`.
- Open question: Make CTMRG validation production-like enough for physics
  claims and energy-targeted ScarFinder ranking.
- Source: `docs/superpowers/specs/2026-05-15-infinite-tfim-benchmark-design.md`
- Source: `docs/superpowers/notes/2026-05-15-current-work-infinite-tfim-benchmark.md`
- Source: `README.md`
