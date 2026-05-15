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

## TFIM Benchmark Branch Verification

- Confirmed: Branch `codex/infinite-tfim-benchmark` was verified with focused
  TFIM/benchmark tests passing `609/609`, a benchmark/public-doc/Aqua subset
  passing `48/48`, and the full suite passing `64704/64704` in `4m43.0s`.
- Confirmed: A manual `J = 0` smoke benchmark produced three samples and final
  `<Z> = 0.9992001066609787`.
- Source: `docs/superpowers/notes/2026-05-15-current-work-infinite-tfim-benchmark.md`
