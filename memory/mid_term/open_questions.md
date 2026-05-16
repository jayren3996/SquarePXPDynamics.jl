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
  sweep records, and S7b requires that trust before gauge-changing updates.
- Open question: What physics-facing convergence thresholds, finite-chi
  sensitivity policy, and benchmark evidence are sufficient before CTM values
  can drive energy-oriented ScarFinder ranking or external physics claims?
- Source: `README.md`
- Source: `notes/2026-05-15-gpt-pro-ctm-scarfinder-revision-notes.md`
- Source: `src/CTMTrust.jl`
- Source: `src/CTMGaugeReadiness.jl`

## Production ScarFinder Validation

- Open question: How should the completed S0-S7 infrastructure be assembled
  into production ScarFinder validation runs with CTM-trusted energy/ranking,
  finite-chi sweeps, benchmark comparisons, and acceptance criteria?
- Source: `README.md`
- Source: `docs/superpowers/notes/2026-05-16-s0-s7-completion-audit.md`
- Source: `docs/superpowers/notes/2026-05-16-s7b-gauge-fixing-handoff.md`

## PEPSKit Public API Boundary

- Open question: Should PEPSKit/TensorKit-facing code remain core public API,
  become experimental-but-exported API, or move behind a package extension once
  project boundaries settle?
- Source: `notes/2026-05-15-code-quality-audit.md`

## Superseded Context To Watch

- Superseded: Older notes describe the repo as having only a few modules and
  not yet having production Simple Update, evolution, or ScarFinder scaffolding.
  Current `README.md` and source show the S0-S7 prototype now exists, including
  S7b CTM gauge-readiness and conditioning APIs.
- Source: `notes/2026-05-15-ipeps-literature-code-algorithm-notes.md`
- Source: `README.md`
- Source: `src/SquarePXPDynamics.jl`
