# Current State

- Confirmed: Current workspace is `/data/djxg096/SquarePXPDynamics.jl`.
- Confirmed: Current checkout is branch `codex/m3-larger-d-pxp-ed-benchmark`,
  not local `main`.
- Confirmed: The worktree is intentionally dirty with current iPEPS+CTM
  throughput changes and untracked benchmark artifacts/logs. Do not clean or
  revert unrelated artifacts without explicit user direction.
- Confirmed: The active priority is iPEPS+CTM observable performance, not ED.
  ED runs through `3x3..6x6` are sufficient for the current campaign; `7x7` ED
  dynamics were explicitly deprioritized by the user.
- Confirmed: Current CTM work added direct `Strided` dependency and exported
  CTM threading controls:
  `configure_ctm_threading!` and `configure_ctm_threading_from_env!`.
- Confirmed: `scripts/pxp_larger_d_ed_benchmark.jl` now applies
  `SQUAREPXP_CTM_*` threading environment variables and prints the active CTM
  threading tuple.
- Confirmed: Tests run after the CTM threading change:
  `julia --project=test test/runtests.jl test_pepskit_measurements.jl` passed
  `93/93`; `julia --project=test test/runtests.jl test_public_docs.jl` passed
  `8/8`.
- Confirmed: Direct CTM probe at `3x3`, `t = 0.02`, `chi = 2` showed
  environment/CTM density fixes the D=2 local-observable mismatch:
  D=2 `density_simple ≈ 0.0002109498`, exact finite density
  `≈ 0.0003996270`, and CTM density `≈ 0.0003996270`.
- Source: `git status`
- Source: `src/PEPSKitMeasurements.jl`
- Source: `scripts/pxp_larger_d_ed_benchmark.jl`
- Source: `artifacts/m3-systematic/ctm-direct-3x3-t002.json`
- Source: current 2026-05-17 session
