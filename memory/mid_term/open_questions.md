# Open Questions

## TFIM Schedule Reference Scope

- Resolved on 2026-05-16: A dense finite-Hilbert-space TFIM reference was added
  for small periodic cells such as `3 x 3`; larger sparse/EDKit-backed
  references remain future scope if needed.
- Former open question: The broad TFIM design mentions a tiny dense/sparse finite
  Hilbert-space schedule reference. The executed implementation plan delivered
  coefficient, non-overlap, and mapping schedule checks. Decide whether a full
  finite Hilbert-space TFIM simulator is required after the initial merge.
- Source: `docs/superpowers/specs/2026-05-15-infinite-tfim-benchmark-design.md`
- Source: `docs/superpowers/notes/2026-05-15-current-work-infinite-tfim-benchmark.md`
- Source: `src/FiniteTFIMReference.jl`

## TFIM Smoke Matrix

- Open question: The broad TFIM design lists a Tier 2 smoke matrix across
  `J = 0`, `h = 0`, small field, near-critical field, and large field. The
  executed plan's manual smoke covered the planned `J = 0` case. Decide whether
  to add the broader matrix as follow-up work.
- Source: `docs/superpowers/specs/2026-05-15-infinite-tfim-benchmark-design.md`
- Source: `docs/superpowers/notes/2026-05-15-current-work-infinite-tfim-benchmark.md`

## CTMRG Trust Policy

- Partially resolved: S7a provides a software trust policy over finite-chi CTM
  sweep records, S7b requires that trust before gauge-changing updates, and
  ScarFinder can now require trusted CTM measurements for candidate ranking.
- Open question: What physics-facing convergence thresholds, finite-chi
  sensitivity policy, and benchmark evidence are sufficient before CTM values
  can drive energy-oriented ScarFinder ranking or external physics claims?
- Source: `README.md`
- Source: `notes/2026-05-15-gpt-pro-ctm-scarfinder-revision-notes.md`
- Source: `src/CTMTrust.jl`
- Source: `src/CTMGaugeReadiness.jl`
- Source: `src/ScarFinder.jl`

## CTM Throughput Tuning

- Open question: What CTM threading layout gives the best throughput on the
  current multi-core server for representative iPEPS+CTM observable runs?
- Current hypothesis: use Julia threads for PEPSKit `dtmap` parallel regions
  and Strided/TensorOperations work, keep BLAS threads small initially, and
  test Strided threaded matrix splitting explicitly.
- Required evidence: a warmed timing matrix that avoids first-compile noise and
  records wall time, CPU percentage, `JULIA_NUM_THREADS`,
  `SQUAREPXP_CTM_BLAS_THREADS`, `SQUAREPXP_CTM_STRIDED_THREADS`,
  `SQUAREPXP_CTM_STRIDED_THREADED_MUL`, and PEPSKit scheduler.
- Source: `src/PEPSKitMeasurements.jl`
- Source: current 2026-05-17 session

## Production ScarFinder Validation

- Partially resolved: ScarFinder now has explicit objective objects,
  `TrustedCTMBackend`, `require_trusted_ctm`, scar-oriented observables,
  candidate metadata persistence, and convergence-report infrastructure.
- Open question: What first production/audit campaign should be run, and what
  acceptance thresholds should define a publishable candidate trajectory?
- Open question: Should the next implementation milestone prioritize full
  tensor snapshot persistence, expanded CTM observables, or CTM-aware/full
  updates?
- Source: `README.md`
- Source: `docs/superpowers/notes/2026-05-16-s0-s7-completion-audit.md`
- Source: `docs/superpowers/notes/2026-05-16-s7b-gauge-fixing-handoff.md`
- Source: `docs/superpowers/notes/2026-05-17-gpt-roadmap-completion.md`

## ScarFinder Candidate Persistence

- Partially resolved: `JSONCandidateStore` writes candidate metadata, scores,
  trust fields, and rejection reasons.
- Open question: What tensor snapshot format should be used for exact candidate
  reruns: JLD2/HDF5, ITensors-native serialization, or a custom JSON metadata
  plus binary tensor payload layout?
- Source: `src/ScarFinder.jl`
- Source: `README.md`

## CTM Observable Roadmap

- Open question: Which CTM-backed observables should be implemented next:
  return/fidelity proxy, two-point correlations, transfer-matrix correlation
  length, structure-factor variants, or energy-variance-quality diagnostics?
- Source: `README.md`
- Source: `docs/superpowers/notes/2026-05-17-gpt-roadmap-completion.md`

## PEPSKit Public API Boundary

- Partially resolved: Private/helper PEPSKit full-update names are now guarded
  by `pepskit_private_full_update_available()`.
- Open question: Should PEPSKit/TensorKit-facing code remain core public API,
  become experimental-but-exported API, or move behind a package extension once
  project boundaries settle?
- Source: `notes/2026-05-15-code-quality-audit.md`
- Source: `src/CTMGaugeReadiness.jl`

## Superseded Context To Watch

- Superseded: Older notes describe the repo as having only a few modules and
  not yet having production Simple Update, evolution, or ScarFinder scaffolding.
  Current `README.md` and source show the S0-S7 prototype now exists, including
  S7b CTM gauge-readiness and conditioning APIs.
- Superseded/problematic: Older short-term memory claimed the active workspace
  was `/Users/ren/Codex/iPEPS`, local `main`, clean, and ready for milestone
  selection. Current work is `/data/djxg096/SquarePXPDynamics.jl` on
  `codex/m3-larger-d-pxp-ed-benchmark`, with active iPEPS+CTM performance
  changes.
- Superseded/problematic: Older guidance to stop before CTM because of the D=2
  simple-density anomaly is no longer the active interpretation. The anomaly is
  a simple/local measurement limitation for D>1; CTM/environment observables
  are now the correct path for local-density comparisons.
- Superseded/problematic: Older 7x7 ED target language should not drive current
  work. The user explicitly redirected to stop at `6x6` and improve iPEPS+CTM
  observables/performance.
- Source: `notes/2026-05-15-ipeps-literature-code-algorithm-notes.md`
- Source: `README.md`
- Source: `src/SquarePXPDynamics.jl`
- Source: `memory/short_term/current_state.md`
- Source: `docs/superpowers/notes/2026-05-17-m3-systematic-larger-d-results.md`
