# Decision Log

## 2026-05-15 - Keep The Project Square-Lattice Focused

Decision:

Keep `SquarePXPDynamics` focused on square-lattice PXP dynamics and directly
related PEPS/iPEPS tooling for ScarFinder.

Reason:

The near-term goal is a clean, testable square-lattice baseline before harder
geometries or broader tensor-network abstractions.

Alternatives considered:

General tensor-network package, arbitrary graph/lattice support, broad
Hamiltonian packaging.

Consequences:

New features should justify themselves through evolve/project/diagnose/rank
workflows for square-lattice PXP or explicitly scoped benchmarks.

Source:

`notes/README.md`; `notes/2026-05-15-ipeps-literature-code-algorithm-notes.md`

Status: active

## 2026-05-15 - Dense Five-Site PXP Is The Source Of Truth

Decision:

Use the dense 32x32 five-site square-star PXP Hamiltonian and gates as the
source of truth for PXP local physics.

Reason:

The local star operator is small enough to test exhaustively, and it provides a
stable convention anchor for ITensor gate conversion and PEPS/iPEPS updates.

Alternatives considered:

Building directly around a PEPS backend's local operator abstraction or a broad
Hamiltonian framework.

Consequences:

Basis order, star order, blockade projection, and dense-to-ITensor conversion
must remain directly tested.

Source:

`README.md`; `notes/README.md`; `src/SquarePXP.jl`; `test/test_square_pxp.jl`

Status: active

## 2026-05-15 - Use A Custom Five-Site Update And PEPSKit For CTM Measurements

Decision:

Keep the custom five-site star update in this package and use PEPSKit/TensorKit
as an experimental CTMRG measurement adapter rather than the update backend.

Reason:

PEPSKit `0.7.0` supported `InfinitePEPS`, CTMRG, and custom five-site
measurements, but its registered simple-update/time-evolution API did not
provide a reliable ready path for arbitrary connected five-site gates.

Alternatives considered:

Build the whole backend directly on PEPSKit, or avoid PEPSKit in the main
package and leave CTMRG validation external.

Consequences:

The package keeps PEPSKit as a measurement/validation dependency while
maintaining custom QR-reduced star-update code.

Source:

`notes/2026-05-15-pepskit-backend-feasibility.md`; `src/PEPSKitMeasurements.jl`;
`src/StarSimpleUpdate.jl`

Status: active

## 2026-05-15 - Simple Observables Are Diagnostics Only

Decision:

Treat simple/local observables and ScarFinder simple scores as diagnostics and
regression signals, not CTMRG-quality physics claims.

Reason:

Simple update lacks the full infinite environment needed for reliable
quantitative observables, especially near criticality or for energy ranking.

Alternatives considered:

Use simple/local energy-like diagnostics directly as physics-quality ranking
metrics.

Consequences:

README and benchmark docs must warn users to inspect CTM diagnostics and finite
chi sensitivity before trusting energy comparisons or physics claims.

Source:

`README.md`; `notes/2026-05-15-gpt-pro-ctm-scarfinder-revision-notes.md`;
`src/ScarFinder.jl`

Status: active

## 2026-05-15 - TFIM Benchmark V1 Uses Static Star Models And Simple-Update Records

Decision:

For the v1 TFIM benchmark branch, add a narrow star-model/protocol layer and
serialize simple-update benchmark trajectories to JSON/CSV.

Reason:

The goal is to validate the existing square-star iPEPS machinery on a model
with stronger external references while preserving PXP behavior and avoiding
premature CTMRG physics claims.

Alternatives considered:

General Hamiltonian framework, arbitrary lattice support, local-quench
spectral-function tooling, or CTMRG-quality TFIM comparisons in v1.

Consequences:

The branch records reproducible simple-update diagnostics, includes exact-limit
tests, and defers CTM-backed TFIM validation and broader smoke matrices.

Source:

`docs/superpowers/specs/2026-05-15-infinite-tfim-benchmark-design.md`;
`docs/superpowers/notes/2026-05-15-current-work-infinite-tfim-benchmark.md`

Status: active

## 2026-05-15 - Maintain Repository-Local Project Memory

Decision:

Create a repository-local `memory/` system with long-term, mid-term, and
short-term layers plus a decision log and handoff file.

Reason:

Future agent sessions need a concise recovery path for scientific context,
project state, decisions, and active work without rereading every note and
source file.

Alternatives considered:

Continue relying only on scattered dated notes, README text, and chat context.

Consequences:

Future substantial work should begin by reading `memory/README.md` and the
listed mid/short-term files. The `project-memory-curator` skill should only run
when explicitly requested.

Source:

User request on 2026-05-15; `/Users/ren/.codex/skills/project-memory-curator/SKILL.md`

Status: active

## 2026-05-16 - S7a Uses Separate CTM Trust And Gauge Diagnostics Modules

Decision:

Keep CTM finite-chi trust policy in `src/CTMTrust.jl` and read-only
simple-gauge diagnostics in `src/GaugeDiagnostics.jl`, with no mutating
gauge-fixing API in S7a.

Reason:

CTM measurement trust and gauge-update readiness have different correctness
requirements. S7a validates measurement stability and local simple-gauge
diagnostics, while S7b must separately validate environment norm matrices
before mutating tensors.

Consequences:

S7a exports `assess_ctm_trust`, `write_ctm_trust_csv`, and read-only
`gauge_diagnostic_simple` helpers. Mutating `fix_bond_gauge!` or
full-update-style conditioning remains deferred until S7b.

Source:

`docs/superpowers/specs/2026-05-16-s7-ctm-trust-gauge-readiness-design.md`;
`src/CTMTrust.jl`; `src/GaugeDiagnostics.jl`

Status: active
