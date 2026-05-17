# Next Steps

1. Finish the CTM throughput cleanup: review the new threading API, make sure
   it is documented clearly enough for benchmark commands, and decide whether
   README updates are needed.
2. Run a small CTM-only timing matrix in warmed sessions, not simultaneous
   first-compile probes. Compare at least:
   `JULIA_NUM_THREADS=1/8/42`, `SQUAREPXP_CTM_BLAS_THREADS=1`, and
   `SQUAREPXP_CTM_STRIDED_THREADED_MUL=false/true`.
3. Prefer `JULIA_NUM_THREADS=42`,
   `SQUAREPXP_CTM_BLAS_THREADS=1`,
   `SQUAREPXP_CTM_STRIDED_THREADS=42`,
   `SQUAREPXP_CTM_STRIDED_THREADED_MUL=true`,
   `SQUAREPXP_CTM_PEPSKIT_SCHEDULER=dynamic` as the first full-server CTM
   utilization experiment.
4. After threading is understood, rerun a small iPEPS+CTM observable campaign
   using environment/CTM density fields and reuse cached ED/iPEPS artifacts
   where possible.
5. Only after the CTM performance path is stable, return to broader milestones:
   CTM-trusted ScarFinder audit, tensor snapshot persistence, expanded CTM
   observables, or CTM-aware/full-update design.

Source: `memory/short_term/current_state.md`
Source: `memory/mid_term/open_questions.md`
Source: `src/PEPSKitMeasurements.jl`
Source: `scripts/pxp_larger_d_ed_benchmark.jl`
