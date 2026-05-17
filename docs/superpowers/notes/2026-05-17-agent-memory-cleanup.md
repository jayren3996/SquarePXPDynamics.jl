# Agent Memory Cleanup

Date: 2026-05-17

## Why This Note Exists

The project memory had several stale short-term entries from the previous
merged `main` state. They incorrectly described the active workspace as
`/Users/ren/Codex/iPEPS`, the checkout as clean local `main`, and the next task
as broad milestone selection. Current work is instead happening in
`/data/djxg096/SquarePXPDynamics.jl` on
`codex/m3-larger-d-pxp-ed-benchmark`.

## Current Priority

The active priority is iPEPS+CTM observable throughput and CPU utilization.
ED work is complete enough through `3x3..6x6`; the user explicitly redirected
away from `7x7` ED and toward improving iPEPS CTM/environment observables.

## Updated Problematic Entries

- Short-term memory now records the active branch, dirty worktree, current CTM
  threading changes, and the iPEPS+CTM performance priority.
- The old "stop before CTM" memory is marked superseded. The D=2 issue is now
  understood as a simple/local observable limitation for D>1 loopy iPEPS, not a
  reason to avoid CTM.
- The old `7x7` ED target memory is marked superseded for the current campaign.
  The current ceiling is `6x6` unless the user explicitly reopens `7x7`.
- Mid-term architecture and experiment memory now mention the direct `Strided`
  dependency and CTM threading controls.

## Current Technical State

New CTM threading controls:

- `configure_ctm_threading!`
- `configure_ctm_threading_from_env!`
- `SQUAREPXP_CTM_BLAS_THREADS`
- `SQUAREPXP_CTM_STRIDED_THREADS`
- `SQUAREPXP_CTM_STRIDED_THREADED_MUL`
- `SQUAREPXP_CTM_PEPSKIT_SCHEDULER`

The benchmark script `scripts/pxp_larger_d_ed_benchmark.jl` applies the
`SQUAREPXP_CTM_*` variables and prints the applied threading tuple.

## Current Evidence

- `test_pepskit_measurements.jl` passed `93/93`.
- `test_public_docs.jl` passed `8/8`.
- A direct CTM probe at `3x3`, `t = 0.02`, `chi = 2` showed D=2 CTM density
  matching exact finite density while simple density did not.
- Earlier simultaneous CTM timing probes were polluted by first-use compilation
  in separate Julia sessions. Future timing should use warmed sessions.

## Recommended Immediate Next Step

Run a warmed CTM timing matrix before launching larger CTM sweeps. First
candidate full-server command shape:

```bash
JULIA_NUM_THREADS=42 \
SQUAREPXP_CTM_BLAS_THREADS=1 \
SQUAREPXP_CTM_STRIDED_THREADS=42 \
SQUAREPXP_CTM_STRIDED_THREADED_MUL=true \
SQUAREPXP_CTM_PEPSKIT_SCHEDULER=dynamic \
julia --project=. scripts/pxp_larger_d_ed_benchmark.jl
```
