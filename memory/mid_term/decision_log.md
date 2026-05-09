# Decision Log

## 2026-05-08 - Keep PEPS as internal root-package tooling

Decision:

Keep the PEPS/iPEPS layer inside the root `TriangularPEPSDynamics` Julia project. Do not create a nested package, nested `Project.toml`, or separate PEPS environment.

Reason:

The repository exists for 2D triangular-lattice PXP ScarFinder work, and PEPS features should serve that workflow rather than becoming a standalone tensor-network library.

Alternatives considered:

Nested Julia package or independent PEPS environment.

Consequences:

New PEPS code belongs under `src/`; tests belong under `test/`; run with `julia --project=. -e 'using Pkg; Pkg.test()'`.

Source:

`AGENTS.md`; `README.md`; `Project.toml`

Status: active

## 2026-05-08 - Preserve center-first 7-site star and basis conventions

Decision:

Use `|up> = |0>`, `|down> = |1>`, and dense star ordering `center, neighbor_1, ..., neighbor_6` in triangular direction order.

Reason:

These conventions pin dense gate construction, blockade projectors, tests, and iPEPS writeback semantics.

Alternatives considered:

Changing basis labels or neighbor ordering to match other physics/tensor-network conventions.

Consequences:

All future gates, diagnostics, and tests must respect the existing order.

Source:

`AGENTS.md`; `Notes/implemented_peps_algorithm_detail.md`; `src/Models.jl`; `src/States.jl`

Status: active

## 2026-05-08 - Use local projected gates for constrained PXP evolution

Decision:

Build constrained real- and imaginary-time PXP gates as `P_blockade * U` and `P_blockade * G` on the dense 7-site star.

Reason:

This keeps blockade enforcement explicit and local while the PEPS projection backend is still developing.

Alternatives considered:

Penalty evolution, post-hoc cleanup only, or relying on candidate rejection.

Consequences:

Projection/truncation may still leak outside the constrained manifold, so blockade diagnostics remain mandatory.

Source:

`Notes/implementation_roadmap.md`; `Notes/current_peps_evolution_solution.md`; `src/Gates.jl`; `src/Observables.jl`

Status: active

## 2026-05-08 - Make Simple Update the first projection backend

Decision:

Use Simple Update as the first fixed-bond-dimension projection backend, with dense-star correctness paths and local diagnostics before stronger update schemes.

Reason:

Simple Update is the cheapest stable path to an evolve-project loop and avoids CTMRG/full-environment infrastructure while gate ordering and diagnostics are still being validated.

Alternatives considered:

Starting directly with NTU, CTMRG/full update, or a PESS ansatz refactor.

Consequences:

Current `D>1` diagnostics are local approximations; NTU or better local factorization remains future work.

Source:

`Notes/literature_review.md`; `Notes/implementation_roadmap.md`; `src/SimpleUpdate.jl`

Status: active

## 2026-05-08 - Separate dynamics and ScarFinder truncation dimensions

Decision:

Use `dynamics_maxdim` during Trotter evolution and optionally hard-truncate to `scar_maxdim` after each projection interval.

Reason:

ScarFinder benefits from evolving in a slightly larger working manifold before returning to the lower-dimensional search manifold.

Alternatives considered:

Single fixed dimension for both evolution and candidate manifold.

Consequences:

The split does not remove per-layer iPEPS projection; it only changes the target dimension used during those local updates versus later hard truncation.

Source:

`Notes/current_peps_evolution_solution.md`; `src/ScarFinder.jl`; `test/test_scar_finder.jl`

Status: active
