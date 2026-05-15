# Open Questions

## Integration State

- Open question: Should `codex/infinite-tfim-benchmark` be merged locally,
  pushed as a PR, or extended before review?
- Source: `docs/superpowers/notes/2026-05-15-current-work-infinite-tfim-benchmark.md`

## TFIM Schedule Reference Scope

- Open question: The broad TFIM design mentions a tiny dense/sparse finite
  Hilbert-space schedule reference. The executed implementation plan delivered
  coefficient, non-overlap, and mapping schedule checks. Decide whether a full
  finite Hilbert-space TFIM simulator is required before integration.
- Source: `docs/superpowers/specs/2026-05-15-infinite-tfim-benchmark-design.md`
- Source: `docs/superpowers/notes/2026-05-15-current-work-infinite-tfim-benchmark.md`

## TFIM Smoke Matrix

- Open question: The broad TFIM design lists a Tier 2 smoke matrix across
  `J = 0`, `h = 0`, small field, near-critical field, and large field. The
  executed plan's manual smoke covered the planned `J = 0` case. Decide whether
  to add the broader matrix before publishing the branch.
- Source: `docs/superpowers/specs/2026-05-15-infinite-tfim-benchmark-design.md`
- Source: `docs/superpowers/notes/2026-05-15-current-work-infinite-tfim-benchmark.md`

## CTMRG Trust Policy

- Open question: What convergence thresholds and finite-chi sensitivity policy
  are sufficient before CTM values can drive energy-oriented ScarFinder ranking
  or physics claims?
- Source: `README.md`
- Source: `notes/2026-05-15-gpt-pro-ctm-scarfinder-revision-notes.md`

## PEPSKit Public API Boundary

- Open question: Should PEPSKit/TensorKit-facing code remain core public API,
  become experimental-but-exported API, or move behind a package extension once
  project boundaries settle?
- Source: `notes/2026-05-15-code-quality-audit.md`

## Superseded Context To Watch

- Superseded: Older notes describe the repo as having only a few modules and
  not yet having production Simple Update, evolution, or ScarFinder scaffolding.
  Current `README.md` and source show the S0-S6 prototype now exists.
- Source: `notes/2026-05-15-ipeps-literature-code-algorithm-notes.md`
- Source: `README.md`
- Source: `src/SquarePXPDynamics.jl`
