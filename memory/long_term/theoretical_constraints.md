# Theoretical Constraints

## Measurement Quality

- Confirmed: Simple/local observables are regression diagnostics and should not
  be used for final physics claims.
- Confirmed: CTMRG-backed measurements require convergence diagnostics and
  finite-chi sensitivity checks before being used for quantitative ranking or
  literature comparison.
- Confirmed: S7b gauge conditioning is a readiness/conditioning layer over
  trusted CTM contexts and local bond norm diagnostics; it is not by itself a
  full-update solver or production ScarFinder validation.
- Source: `README.md`
- Source: `notes/2026-05-15-gpt-pro-ctm-scarfinder-revision-notes.md`
- Source: `src/PEPSKitMeasurements.jl`
- Source: `src/CTMGaugeReadiness.jl`

## Blockade And Projection

- Confirmed: Projected local gates do not guarantee the truncated PEPS/iPEPS
  state remains globally blockade-constrained, so blockade leakage must be
  monitored.
- Source: `README.md`
- Source: `notes/README.md`

## Scope Boundaries

- Confirmed: The repository should not grow into a general tensor-network
  package unless a concrete ScarFinder need justifies it.
- Confirmed: Avoid arbitrary graph support, broad Hamiltonian packaging, GPU
  backends, symmetry machinery, and broad module splitting unless they solve a
  specific near-term project problem.
- Source: `notes/README.md`
- Source: `notes/2026-05-15-code-quality-audit.md`
- Source: `docs/superpowers/specs/2026-05-15-infinite-tfim-benchmark-design.md`

## PEPSKit Boundary

- Confirmed: PEPSKit `0.7.0` can support `InfinitePEPS`, CTMRG environments,
  and custom local measurements, including five-site star observables.
- Confirmed: PEPSKit's registered simple-update/time-evolution API did not
  expose a reliable ready path for the custom five-site square-star update, so
  the project keeps custom local update logic outside PEPSKit.
- Source: `notes/2026-05-15-pepskit-backend-feasibility.md`

## Testing Constraints

- Confirmed: D=1 dense/product-state references are useful for signs,
  conventions, and exact limits.
- Confirmed: D=2 update tests should compare gauge-invariant observables, not
  raw tensor equality, because gauge-changing factorizations can alter tensor
  representatives.
- Confirmed: S7b tests should compare observables, CTM summaries, freshness
  guards, and transactional mutation behavior rather than raw tensor entries.
- Source: `notes/2026-05-15-gpt-pro-ctm-scarfinder-revision-notes.md`
- Source: `notes/2026-05-15-ipeps-literature-code-algorithm-notes.md`
- Source: `docs/superpowers/notes/2026-05-16-s7b-gauge-fixing-handoff.md`
