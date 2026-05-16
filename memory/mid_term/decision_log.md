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

## 2026-05-16 - Separate S7b Readiness From D>1 Gauge Mutation

Decision:

Add CTM local bond norm diagnostics and a public readiness predicate before
attempting D>1 mutating gauge conditioning. Slice 4 initially made
`fix_bond_gauge!` transactional and no-op for D=1 product bonds, with D>1
mutation added later in Slice 5.

Reason:

The PEPSKit CTM environment supplies a usable bond-environment contraction for
readiness diagnostics. Writing gauge-conditioned D>1 factors back into the
custom ITensors Gamma-lambda state needed a separate representation decision
and gauge-invariant regression suite, which Slice 5 added.

Consequences:

ScarFinder and future full-update work can now gate gauge-changing updates on
fresh CTM contexts, finite-chi trust, bond coverage, Hermiticity, PSD floor,
and reciprocal condition number before D=1 no-op or D>1 gauge conditioning.

Source:

`docs/superpowers/plans/2026-05-16-s0-s7-slice4-ctm-gauge-readiness.md`;
`src/CTMGaugeReadiness.jl`; `test/test_ctm_gauge_readiness.jl`

Status: active

## 2026-05-16 - Use PEPSKit Bond-Environments For D>1 Gauge Conditioning

Decision:

Implement D>1 `fix_bond_gauge!` by using PEPSKit's bond-environment
factorization on the converted absorbed PEPS tensors, then converting the two
conditioned tensors back into the custom ITensors Gamma-lambda representation
by deabsorbing the stored link weights.

Reason:

This reuses PEPSKit's tested CTM bond-environment gauge machinery instead of
inventing an independent full-update factorization. Writing back only after
successful factorization and conversion preserves the project's transactional
mutation contract and invalidates stale CTM contexts through `state_version`.

Consequences:

The S7b gauge-fixing API now mutates D>1 bonds when readiness passes and link
weights on the affected tensor legs are positive. The implementation remains a
gauge-conditioning layer, not a full ALS/full-update truncation solver.

Source:

`docs/superpowers/plans/2026-05-16-s0-s7-slice5-d2-gauge-conditioning.md`;
`src/CTMGaugeReadiness.jl`; `test/test_ctm_gauge_readiness.jl`

Status: active

## 2026-05-16 - Add Serial Star-Sweep Evolution Schedule

Decision:

Support a serial square-star Trotter schedule in addition to the existing
five-color batched schedule. The default remains `:five_color`; users can set
`schedule = :serial` in `TrotterParams` or the `evolve!` convenience keyword.

Reason:

The five-color schedule is efficient and layer-parallel but imposes rectangular
unit-cell dimensions compatible with the coloring. A serial sweep applies one
star gate per layer, so overlapping stars are ordered rather than batched and
non-five-color-compatible cells such as `4 x 4` can be evolved.

Alternatives considered:

Replace the five-color scheduler entirely, or relax the five-color unit-cell
assertions while still batching overlapping stars.

Consequences:

Serial sweeps make the center order part of the algorithm. First-order serial
evolution visits unit-cell representatives in stored order; second-order serial
evolution uses a forward half sweep, a full final-center step, and the reversed
half sweep. Existing five-color behavior is preserved as the default.

Source:

Current user discussion; `src/IPEPSEvolution.jl`;
`test/test_ipeps_evolution.jl`

Status: active

## 2026-05-16 - TFIM Benchmark Uses Serial 3x3 Smoke Scheme

Decision:

Switch the TFIM smoke benchmark path from a five-color-compatible unit cell to
a `3 x 3` unit cell using `schedule = :serial`, and record the Trotter schedule
in benchmark metadata and flattened CSV output.

Reason:

The current serial star-sweep scheme is the relevant comparison target after
the square-star scheduling discussion. A `3 x 3` cell is intentionally not
five-color-compatible, so the benchmark now exercises the serial path directly
instead of silently relying on the older five-color batching constraint.

Consequences:

Benchmark JSON/CSV records now distinguish `:five_color` and `:serial`
evolution schedules. Existing benchmark calls that omit `schedule` still keep
the default `:five_color`, but the documented TFIM smoke run and realistic TFIM
dynamics script use the serial 3x3 configuration.

Source:

Current user request; `src/Benchmarks.jl`; `test/test_benchmarks.jl`;
`scripts/realistic_tfim_dynamics.jl`; `README.md`

Status: active

## 2026-05-16 - Add Dense Finite TFIM Reference For Small Benchmarks

Decision:

Add a dense finite-Hilbert-space TFIM reference path for small periodic cells,
starting with the `3 x 3` benchmark cell. The implementation is repository
local and dependency-light; it does not add `EDKit.jl` as a runtime dependency.

Reason:

For nine spins the dense Hamiltonian is only `512 x 512`, so a local exact
reference is simpler and pins exactly the same basis, bond, and observable
conventions as the package. `EDKit.jl` remains a possible future option if the
reference scope grows beyond small dense cells, but adding a GitHub-only ED
dependency now would make this benchmark less reproducible without solving a
current size problem.

Consequences:

The TFIM dynamics script now compares simple-update observables against exact
finite `3 x 3` TFIM time evolution. The reference is a finite-size baseline,
not an infinite-system physics target, and should be used mainly for convention
checks, short-time drift, and regression diagnostics.

Source:

Current user request; `src/FiniteTFIMReference.jl`;
`test/test_tfim_finite_reference.jl`; `scripts/realistic_tfim_dynamics.jl`

Status: active

## 2026-05-16 - Add Open-Boundary 6x6 TFIM MPS Reference

Decision:

Add an ITensorMPS-backed finite-MPS TFIM reference path for open-boundary square
lattices, with a `6 x 6` smoke script using a snake mapping and TDVP.

Reason:

Dense exact diagonalization is not feasible for `6 x 6` (`2^36` states), while
an MPS reference provides a practical finite-size comparison tier above the
exact `3 x 3` dense reference. Open boundaries avoid adding long periodic wrap
bonds to the first MPS benchmark and make the 2D-to-1D mapping explicit.

Consequences:

The `6 x 6` MPS trajectory is an approximate finite-size reference. It should
be compared through trends, short-time observables, energy drift, norm drift,
and bond-dimension convergence, not treated as an exact or infinite-system
target.

Source:

Current user request; `src/FiniteMPSTFIMReference.jl`;
`test/test_finite_mps_tfim_reference.jl`; `scripts/finite_mps_tfim_6x6.jl`

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

## 2026-05-16 - Compare 6x6 TDVP Reference Against 3x3 Serial iPEPS Through t=0.30

Decision:

Add a standalone comparison script that aligns the open-boundary 6x6
ITensorMPS TDVP trajectory with the infinite periodic 3x3 serial-sweep iPEPS
trajectory on the shared time grid `0:0.05:0.30`.

Reason:

The 6x6 TDVP reference reaches `maxdim = 64` at `t = 0.30`, so comparing only
through this time avoids extrapolating beyond the first saturated TDVP sample.
The iPEPS run keeps the finer `dt = 0.01` Trotter step and reports every five
steps to match the TDVP sampling grid.

Consequences:

The comparison is a short-time local-observable benchmark, not an equality test:
TDVP is finite and open-boundary, while the iPEPS result is infinite and
periodic. In the first run, `<Z>` remained near zero in both methods and energy
density stayed within about `1e-2`, but `<X>` separated by about `0.235` at
`t = 0.30`; this should motivate bond-dimension and update-scheme convergence
checks before interpreting the discrepancy physically.

Source:

`scripts/compare_tfim_tdvp_ipeps_t03.jl`

Status: active

## 2026-05-16 - Target 7x7 As The Practical Symmetry-Reduced PXP ED Ceiling

Decision:

Use `7 x 7` as the realistic maximum square-lattice PXP ED benchmark target
when a custom constrained EDKit-compatible basis can enumerate translation and
square-point-group orbits directly. Treat `5 x 5` as the out-of-the-box EDKit
ceiling when relying on `basis(...; f, symmetries=...)`, because the current
generic constrained-basis path scans the full `2^N` product space.

Reason:

Exact independent-set counts on periodic square lattices give about `1.19e6`
states in the fully symmetric `7 x 7` space-group orbit basis, but about
`4.17e8` states for `8 x 8`. EDKit's adaptive Krylov implementation stores a
Lanczos basis with many full state vectors, so the `8 x 8` sector would require
hundreds of GB of Krylov-vector storage even before Hamiltonian/application
overheads.

Consequences:

The first ED benchmark should use EDKit's `Operator`/`timeevolve` APIs but add
a PXP-specific constrained orbit basis rather than asking EDKit's generic
`ProjectedBasis`/`AbelianBasis` constructors to discover the constraint by
full Hilbert-space scanning. The fully symmetric sector is suitable for
short-time dynamics from a symmetrized initial state. Broader momentum or
point-group sector coverage remains a later extension.

Source:

Current EDKit feasibility investigation; EDKit `0.5.1` source/docs;
periodic hard-square independent-set and orbit-count calculations in this
thread.

Status: active

## 2026-05-16 - Implement PXP ED Benchmark Through EDKit Krylov

Decision:

Add a repository-local PXP ED benchmark path using EDKit as the operator and
adaptive-Krylov backend. The implementation introduces a PXP-specific
space-group basis that enumerates periodic hard-square states directly and then
uses EDKit `Operator`, sparse conversion, and `KrylovEvolutionCache` for
short-time dynamics.

Reason:

EDKit's generic `basis(...; f, symmetries=...)` path is correct for small
systems but discovers the PXP constraint by scanning the full product space.
That is not viable for the `7 x 7` target. A narrow constrained-orbit basis
keeps the benchmark square-PXP-specific while still relying on EDKit for the
ED-facing APIs requested for the benchmark.

Consequences:

`scripts/pxp_ed_7x7_benchmark.jl` is the large benchmark entry point. Normal
tests validate the same code path only on `3 x 3` and `4 x 4` sectors. The
default initial state is all-down because odd periodic square lattices do not
support a perfect checkerboard product state.

Source:

Current user request; `src/FinitePXPEEDBenchmark.jl`;
`test/test_pxp_ed_benchmark.jl`; `scripts/pxp_ed_7x7_benchmark.jl`

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

## 2026-05-16 - Reconcile S0-S7 Completion With Current Architecture

Decision:

Fulfill the original S0-S7 trajectory by preserving the verified custom
ITensors iPEPS update stack, completing remaining helper, diagnostic,
ScarFinder, and S7b readiness gaps, and explicitly superseding stale
backend-abstraction requirements that no longer match the chosen architecture.

Reason:

The current implementation has passed default and extended verification and
already embodies an active decision to use custom five-site updates with
PEPSKit as an experimental CTMRG measurement adapter. Reintroducing speculative
single-backend abstractions would add complexity without improving ScarFinder
or S7 gauge-readiness workflows.

Consequences:

S0.5/S1 `AbstractProjectionBackend` and `projection_backend` requirements are
treated as superseded unless a second update backend is introduced. Immediate
completion work should focus on small state helper APIs, link-weight
normalization, stronger update/evolution diagnostics, guarded ScarFinder
energy correction, and S7b CTM norm-matrix/gauge-conditioning infrastructure.

Source:

`docs/superpowers/specs/2026-05-16-s0-s7-completion-design.md`;
current S0-S7 audit in this thread; `notes/2026-05-15-ipeps-superpowers-multistage-plan.md`

Status: active

## 2026-05-15 - Include S7 Gauge-Fixed CTMRG Work In Next Steps

Decision:

Include State S7, CTMRG and gauge-fixed full-update infrastructure, in the next
planning/implementation steps.

Reason:

The current iPEPS update path intentionally uses the simplest Gamma-lambda
simple-update gauge and does not yet select an environment-informed gauge. Prior
plans already identified S7 as the stage where CTMRG-backed measurements,
environment stability, and gauge-conditioned truncation should be introduced.

Alternatives considered:

Continue with only simple-update benchmark/scarfinder diagnostics, or add local
QR regauging without first hardening the CTMRG/environment layer.

Consequences:

Next-step planning should explicitly cover CTMRG convergence diagnostics,
finite-chi sensitivity, gauge-invariant D>1 regression tests, and only then
`fix_bond_gauge!`/full-update-style conditioning of local norm matrices.
Simple-update results remain diagnostic records until that validation layer is
trustworthy.

Source:

User confirmation in current thread; `notes/2026-05-15-chatgpt-pro-ipeps-review-plan.md`;
`notes/2026-05-15-ipeps-superpowers-multistage-plan.md`;
`memory/long_term/literature_context.md`

Status: active

## 2026-05-16 - Promote CTM-Trusted PXP Validation Reports

Decision:

Add a focused `PXPValidation` layer that composes existing CTM sweep, CTM
trust, iPEPS evolution, and finite PXP ED APIs into a machine-readable
validation report.

Reason:

The next reliability step is to make trusted measurement and ED comparison a
normal workflow before changing ScarFinder ranking. A narrow report layer
preserves the existing explicit `measure_simple` / `measure_ctm` boundary and
avoids a broad measurement facade before a second production backend exists.

Consequences:

Short-time PXP runs can now produce reproducible JSON artifacts with ED
density references, iPEPS diagnostics, optional finite-`chi` CTM trust, and
run metadata. ScarFinder can consume this report shape in a follow-on slice.

Source:

`src/PXPValidation.jl`; `test/test_pxp_validation.jl`;
`scripts/validate_pxp_ed_ipeps.jl`

Status: active
