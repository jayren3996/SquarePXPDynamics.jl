# Active Tasks

- Active priority: improve iPEPS+CTM observable throughput and CPU utilization.
- Current implementation branch: `codex/m3-larger-d-pxp-ed-benchmark`.
- Current code change: expose CTM tensor-threading controls for PEPSKit,
  Strided/TensorOperations, and BLAS through `configure_ctm_threading!`,
  `configure_ctm_threading_from_env!`, and `SQUAREPXP_CTM_*` script variables.
- Current benchmark guidance: start CTM-heavy runs with Julia threads, for
  example `JULIA_NUM_THREADS=42`, and avoid assuming `OPENBLAS_NUM_THREADS`
  alone will make PEPSKit CTMRG use all cores.
- Current scientific guidance: use CTM/environment observables for D>1 local
  iPEPS density comparisons; simple/local observables are diagnostics only.
- Current ED guidance: ED work is complete enough through `6x6`; do not restart
  `7x7` ED unless the user explicitly asks.
- Source: user request on 2026-05-17
- Source: `src/PEPSKitMeasurements.jl`
- Source: `scripts/pxp_larger_d_ed_benchmark.jl`
- Source: `git status`
