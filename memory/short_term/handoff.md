# Agent Handoff

## Current Objective

Active implementation branch: `codex/m3-larger-d-pxp-ed-benchmark`.

Current objective: improve iPEPS+CTM observable efficiency and CPU utilization.
Do not spend more time optimizing ED unless the user explicitly redirects.

## What Was Just Done

M3 ED/iPEPS benchmark work established that ED through `3x3..6x6` is enough for
the current campaign and that `7x7` should be stopped/deprioritized. The user
then redirected the active work to iPEPS+CTM observables and CPU utilization.

Current code now exposes CTM threading controls through
`configure_ctm_threading!`, `configure_ctm_threading_from_env!`, and
`SQUAREPXP_CTM_*` environment variables consumed by
`scripts/pxp_larger_d_ed_benchmark.jl`.

## Important Files Touched

- `Project.toml`
- `src/PEPSKitMeasurements.jl`
- `src/SquarePXPDynamics.jl`
- `scripts/pxp_larger_d_ed_benchmark.jl`
- `test/test_pepskit_measurements.jl`
- `memory/short_term/*.md`

## Commands Run

- `julia --project=test test/runtests.jl test_pepskit_measurements.jl`
- `julia --project=test test/runtests.jl test_public_docs.jl`
- Script smoke with `SQUAREPXP_CTM_*` variables and `total_time=0`.
- Threading API smoke:
  `JULIA_NUM_THREADS=8 ... julia --project=. -e 'using SquarePXPDynamics; println(configure_ctm_threading_from_env!())'`

## Tests/Results

- `test_pepskit_measurements.jl` passed `93/93` in `6m16.1s`.
- `test_public_docs.jl` passed `8/8` in `0.6s`.
- Script smoke printed:
  `(julia_threads = 1, blas_threads = 1, strided_threads = 1, strided_threaded_mul = false, pepskit_scheduler = :default)`.
- Environment smoke with `JULIA_NUM_THREADS=8` printed:
  `(julia_threads = 8, blas_threads = 1, strided_threads = 8, strided_threaded_mul = true, pepskit_scheduler = :dynamic)`.

## Known Problems

- CTM throughput still needs a clean warmed timing matrix; earlier simultaneous
  CTM probes were polluted by first-use compilation in separate Julia sessions.
- PEPSKit CTMRG parallelism is Julia-thread based; `OPENBLAS_NUM_THREADS` alone
  is insufficient.
- Strided threaded matrix splitting is now configurable, but its best setting
  for this CTM workload has not yet been measured.
- CTM-aware/full-update evolution is not implemented.
- Full tensor snapshot persistence for ScarFinder candidates is not implemented.
- Publication-grade physics audit sweeps across `dt`, `D`, `chi`, cutoff, unit
  cell, and update scheme remain future work.
- Return/fidelity proxy, richer CTM observables, transfer-matrix correlation
  length, and energy-variance-quality observables remain future work.

## Next Recommended Actions

Run the warmed CTM threading/timing matrix and record the best command pattern
before launching larger iPEPS+CTM observable sweeps.

## Things Not To Do

- Do not treat simple/local observables as CTMRG-quality physics measurements.
- Do not claim publication-grade ScarFinder validation from trusted CTM plumbing
  alone; convergence/audit campaigns are still required.
- Do not restart `7x7` ED by default; the user explicitly said to stop at
  `6x6` and focus on iPEPS+CTM.
- Do not rerun `project-memory-curator` unless explicitly requested.
- Do not delete older notes just because they are summarized here.
