# Current Work Note: Infinite TFIM Benchmark

Date: 2026-05-15

## Summary

We are working on the infinite-system TFIM benchmark framework for the iPEPS
codebase, based on the design in
`docs/superpowers/specs/2026-05-15-infinite-tfim-benchmark-design.md`.

The implementation was executed with a subagent-driven workflow in an isolated
worktree:

- Worktree: `/Users/ren/.config/superpowers/worktrees/iPEPS/codex-infinite-tfim-benchmark`
- Branch: `codex/infinite-tfim-benchmark`
- Current HEAD: `3a73d88 fix: record benchmark split order metadata`
- Status: implementation complete, committed, and worktree clean

The original checkout remains at `/Users/ren/Codex/iPEPS` on `main`, currently
ahead of `origin/main` by three documentation/design commits.

## What Landed

- Added JSON3 dependency for deterministic benchmark serialization.
- Added a square-star model abstraction with `PXPStarModel`, `TFIMStarModel`,
  `StaticModel`, and explicit star gate helpers.
- Threaded model protocols through `project_star!`, `evolve!`, and Trotter
  scheduling while preserving legacy PXP call shapes.
- Added product-state aliases for benchmark initialization, including
  `:z_up`, `:z_down`, and `:x_plus`.
- Added simple TFIM observables, exact-limit checks, and energy-density
  diagnostics comparing star-patch and decomposed estimates.
- Added finite schedule reference coverage for Trotter coefficients, color
  disjointness, and center/site mapping.
- Added benchmark runner APIs:
  `BenchmarkSpec`, `run_benchmark`, `write_benchmark_json`, and
  `write_benchmark_csv`.
- Added deterministic benchmark metadata, including package version,
  protocol/model type, TFIM parameters, cell size, Trotter settings,
  `split_order`, and measurement cadence.
- Documented a TFIM smoke benchmark in `README.md`, including the caveat that
  v1 results are simple-update regression records rather than CTMRG-quality
  physics claims.

## Verification

The implementation branch was verified with:

- Focused TFIM/benchmark tests: `609/609` passed.
- Benchmark/public-doc/Aqua subset after the final metadata fix: `48/48`
  passed.
- Full package test suite: `64704/64704` passed in `4m43.0s`.
- Manual `J = 0` smoke benchmark produced three samples and final
  `<Z> = 0.9992001066609787`.
- `git diff --check` passed.
- Implementation worktree is clean.

## Residual Caveats

- The finite schedule reference follows the executable implementation plan:
  coefficient, non-overlap, and mapping checks. It is not yet a full finite
  Hilbert-space TFIM simulator.
- The manual smoke run covers the planned `J = 0` case. The broader Tier 2
  benchmark matrix remains future work.

## Likely Next Step

Decide how to integrate `codex/infinite-tfim-benchmark`: merge locally, push and
open a PR, or add the future finite Hilbert-space reference / broader smoke
matrix before review.
