# Current State

- Confirmed: Current workspace is `/data/djxg096/SquarePXPDynamics.jl`.
- Confirmed: Current checkout is local branch `main` at commit `8c1eed7` as of
  2026-05-19, ahead of the earlier `98e1ad7` reference point by several
  refactor/docs commits.
- Confirmed: Recent main-branch commits since `98e1ad7` (most recent first):
  - `8c1eed7 docs: root-cause note for pxp-d-debug — cutoff stability at larger D`
  - `c371dcd refactor: inline single-use _ctm_* unwrap helpers in PXPValidation`
  - `f1f51eb refactor: drop unreachable peps-nothing checks in CTM measurement helpers`
  - `7761280 refactor: drop redundant private wrappers in PEPSKitMeasurements`
  - `4432af2 refactor: drop unused legacy outer constructors in PXPValidation`
  - `43e86c2 refactor: collapse redundant CTMObservableSummary constructors and dead throw`
  - `1cde3bf fix: include offending value in env-var ArgumentError messages`
- Confirmed: The active priority remains iPEPS+CTM observable performance, not
  ED. ED runs through `3x3..6x6` are sufficient for the current campaign;
  `7x7` ED dynamics were explicitly deprioritized by the user.
- Confirmed: A multi-agent code-quality + correctness review was run on
  2026-05-19. Source: this session.
- Confirmed: CTM threading controls are exported (`configure_ctm_threading!`,
  `configure_ctm_threading_from_env!`) and use a direct `Strided` dependency.
- Confirmed: Direct CTM probe at `3x3`, `t = 0.02`, `chi = 2` previously showed
  CTM density matching exact finite density to ~1e-15 against simple-update
  D=2 diverging by ~3.8e-4. Recorded baseline:
  D=2 `density_simple ≈ 0.0002109498`, exact finite density
  `≈ 0.0003996270`, CTM density `≈ 0.0003996270`.
- Source: `git log`
- Source: `git rev-parse --abbrev-ref HEAD`
- Source: `src/PEPSKitMeasurements.jl`
- Source: `scripts/pxp_larger_d_ed_benchmark.jl`
- Source: `artifacts/m3-systematic/ctm-direct-3x3-t002.json`
