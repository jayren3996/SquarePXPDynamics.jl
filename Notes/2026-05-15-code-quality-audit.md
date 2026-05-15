# Code Quality Audit - 2026-05-15

Pre-S5b cleanup and hardening audit for the S0-S6 square-lattice PXP prototype.

## Tools Run

- Baseline: `julia --project=. -e 'using Pkg; Pkg.test()'`
  - Julia 1.12.6
  - Result: 10856/10856 passing
- Formatting: temporary `JuliaFormatter` environment at `/tmp/squarepxp-format`, `format(".")`
- Post-format tests: `julia --project=. -e 'using Pkg; Pkg.test()'`
  - Result: 12399/12399 passing
- Aqua: temporary `/tmp/squarepxp-aqua` environment
  - `Aqua.test_all(SquarePXPDynamics)` and ambiguity check passed
- JET: temporary `/tmp/squarepxp-jet` environment
  - `JET.@report_call` probes for `density_simple`, `measure_simple`, and `evolve!`
  - No report output from the requested calls
- Performance smoke: warm elapsed timings on a 10x10, D=2 state
  - `project_star!`: about 0.24 s
  - `evolve!`: about 0.18 s for one warmed step
  - `measure_simple`: about 0.49 s

Manifest.toml is present. CI is present at `.github/workflows/ci.yml` and runs Julia 1.12 with `Pkg.test()`.

## Subagent Findings

- `repo-structure-agent`: README/module status was stale relative to the current public surface; `SquareIPEPS.jl` mixes state, link-weight, and ITensor gate-wrapper responsibilities; optional module split should be deferred unless kept mechanical.
- `julia-quality-agent`: Aqua checks passed; no undefined exports or missing public docstrings found; JET was exploratory; PEPSKit context fields remain untyped and TensorKit compat is pinned.
- `numerics-agent`: link weights could be read from internally corrupted state without revalidation; `deabsorb_link_weight` accepted non-finite `atol`; star-update diagnostics needed explicit finite/normalization checks.
- `tests-agent`: source-inspection tests are fragile but retained; invalid-input coverage and D=2 repeated-update coverage could be improved safely; CTMRG tests may become slow for normal cleanup iteration.
- `karpathy-guidelines-agent`: keep changes boring and local; avoid backend abstractions, CTMRG work, broad module splits, or new physics algorithms in this cleanup.

## Fixed

- Revalidated stored link weights at `link_weight`/`link_weight_tensor` boundaries, so directly corrupted internal spectra fail before use.
- Rejected non-finite `deabsorb_link_weight` tolerances.
- Added explicit star-update singular-spectrum, normalized-link, truncation-error, and norm-factor validation before committing updates.
- Added nested `StarUpdateInfo` diagnostic validation when building `EvolutionLog`.
- Added invalid-input tests for corrupted link weights, non-finite `atol`, non-finite entropy inputs, Trotter parameters, and evolution total time.
- Added repeated D=2 star-update smoke coverage with finite observable and normalized-link checks.
- Updated the module docstring and README status/example/warnings.

## Deferred

- Do not split `SquareIPEPS.jl` into state/link-weight/gate-wrapper modules in this pass; the dependency surface is broad.
- Decide separately whether the existing PEPSKit/TensorKit-facing code is intended to be core public API, experimental S5b surface, or moved behind an optional boundary.
- Clarify or rename `pxp_energy_density_ctm`, because it currently uses a simple/local fallback for energy.
- Consider an export allowlist test once the intended public surface is settled.
- Consider gating slow CTMRG tests if CI timing becomes a problem.
- Consider future log-normalization tracking if repeated updates show norm-factor underflow/overflow.

## Recommendations For S5b

- Settle the PEPSKit/TensorKit boundary before adding more CTMRG measurements.
- Keep simple/local diagnostics clearly separated from future `*_ctm` diagnostics.
- Add CTMRG adapter tests around context/state compatibility and periodic boundary coordinates before relying on CTM summaries.
- If performance becomes relevant, first investigate caching dense/ITensor gates, update centers, physical-index tuples, and dense local operators.
