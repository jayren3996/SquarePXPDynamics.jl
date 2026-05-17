# Experiments And Results

## Code Quality Audit

- Confirmed: A 2026-05-15 audit ran the baseline test suite, formatting,
  Aqua checks, exploratory JET probes, and performance smoke checks.
- Confirmed: The audit fixed link-weight validation, nonfinite tolerance
  rejection, star-update diagnostic validation, nested evolution diagnostic
  validation, invalid-input coverage, D=2 repeated-update smoke coverage, and
  README/module status updates.
- Source: `notes/2026-05-15-code-quality-audit.md`

## PEPSKit Feasibility Probe

- Confirmed: PEPSKit `0.7.0` could instantiate minimal `InfinitePEPS` states,
  build CTMRG environments, and measure 1-site, 2-site, and custom 5-site local
  operators.
- Confirmed: PEPSKit's registered simple-update path was not suitable as a
  ready backend for the custom five-site PXP update.
- Source: `notes/2026-05-15-pepskit-backend-feasibility.md`
- Source: `scripts/dev/pepskit_feasibility_probe.jl`

## CTM And ScarFinder Hardening

- Confirmed: Later CTM work added diagnostics, stale-context protection,
  validation sweeps, CTM metadata in summaries/logs, and README warnings about
  CTM trust and simple-score limitations.
- Source: `README.md`
- Source: `src/PEPSKitMeasurements.jl`
- Source: `src/ScarFinder.jl`
- Source: `notes/2026-05-15-gpt-pro-ctm-scarfinder-revision-notes.md`

## M3 ED And iPEPS/CTM Benchmark Cleanup

- Confirmed: The M3 systematic ED/iPEPS campaign used all-down/all-zero,
  translationally symmetric PBC dynamics and symmetry-reduced ED global
  excitation density/return probability only.
- Confirmed: ED trajectories/artifacts exist through `3x3..6x6`; `7x7` did
  not produce a useful artifact in the stopped capacity probe and was
  deprioritized by user direction.
- Confirmed: For `3x3`, exact finite iPEPS density shows strong improvement
  from `D=1` to `D=2`, while `D=3/4` are not monotone improvements over `D=2`.
- Confirmed: D>1 simple/local density is not an exact finite or thermodynamic
  observable. A direct CTM probe at `3x3`, `t = 0.02`, `chi = 2` showed D=2
  CTM density matches exact finite density while simple density does not.
- Confirmed: Current iPEPS+CTM performance investigation found that PEPSKit
  CTMRG needs Julia threads for its `dtmap` parallel regions; BLAS thread
  settings alone do not maximize CTM CPU utilization.
- Source: `docs/superpowers/notes/2026-05-17-m3-systematic-larger-d-results.md`
- Source: `artifacts/m3-systematic/ctm-direct-3x3-t002.json`
- Source: `src/PEPSKitMeasurements.jl`

## TFIM Benchmark Branch Verification

- Confirmed: Branch `codex/infinite-tfim-benchmark` was verified with focused
  TFIM/benchmark tests passing `609/609`, a benchmark/public-doc/Aqua subset
  passing `48/48`, and the full suite passing `64704/64704` in `4m43.0s`.
- Confirmed: A manual `J = 0` smoke benchmark produced three samples and final
  `<Z> = 0.9992001066609787`.
- Source: `docs/superpowers/notes/2026-05-15-current-work-infinite-tfim-benchmark.md`

## S0-S7 Completion Verification

- Confirmed: Before merging `codex/s0-s7-completion` into `main`, the focused
  S7b gauge-readiness test command
  `julia --project=. test/runtests.jl test_ctm_gauge_readiness.jl` passed
  `98/98` in the feature worktree on 2026-05-16.
- Confirmed: After fast-forwarding `main` and refreshing memory files,
  `julia --project=. -e 'using Pkg; Pkg.test()'` passed `65149/65149` in
  `4m35.1s` on 2026-05-16.
- Confirmed: The S0-S7 audit records prior full-suite verification for the
  branch: default package tests passed `65149/65149`, focused S7b tests passed
  `98/98`, and extended S7b tests passed `102/102`.
- Source: `docs/superpowers/notes/2026-05-16-s0-s7-completion-audit.md`
- Source: `test/test_ctm_gauge_readiness.jl`

## GPT PXP Roadmap Completion Verification

- Confirmed: The GPT PXP roadmap completion branch was verified before merge
  with `julia --project=. -e 'using Pkg; Pkg.test()'` passing `65356/65356`
  in `7m28.3s`.
- Confirmed: After merging to `main`,
  `julia --project=. -e 'using Pkg; Pkg.test()'` passed `65356/65356` in
  `4m57.8s`.
- Confirmed: `git diff --check HEAD~1 HEAD` passed.
- Confirmed: `main` was pushed to GitHub, and `origin/main` points to
  `98e1ad7 Merge GPT PXP roadmap completion`.
- Source: `docs/superpowers/notes/2026-05-17-gpt-roadmap-completion.md`
- Source: `git log`
- Source: `README.md`
