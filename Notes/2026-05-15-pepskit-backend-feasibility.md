# PEPSKit Backend Feasibility

Date: 2026-05-15

## Scope

This note checks whether `SquarePXPDynamics` should use PEPSKit for square-lattice iPEPS/CTMRG work around the 5-site square-star PXP gate. I did not modify package source code. The only exploratory code added is `scripts/dev/pepskit_feasibility_probe.jl`, which activates a temporary Julia environment.

Probe environment:

- Julia `1.12.6`
- PEPSKit `0.7.0`
- TensorKit `0.15.3`
- MPSKit `0.13.8`

## Findings

1. Minimal `InfinitePEPS` works.

   The probe successfully instantiated spin-1/2 square-lattice iPEPS states with `ComplexSpace(2)` physical legs and `ComplexSpace(1)` virtual legs for both `(1, 1)` and `(2, 2)` unit cells. It also built `SUWeight`, `CTMRGEnv`, and ran a minimal CTMRG contraction with `leading_boundary(...; alg = :simultaneous)`.

2. PEPSKit does not currently expose a ready 5-site update path in the registered API.

   A 5-site star `LocalOperator` can be constructed, but `time_evolve(..., ::SimpleUpdate, ...)` failed because the simple-update machinery expects 2-site bond tuples in `is_equivalent_bond`. PEPSKit dev docs also state that `trotterize` currently supports only 1-site terms, nearest-neighbor 2-site terms, and next-nearest-neighbor 2-site terms. The dev docs mention `LocalCircuit`, but `PEPSKit.LocalCircuit` and `PEPSKit.trotterize` were not defined in the registered `0.7.0` API used by the probe, so they are not a reliable dependency target without pinning/validating a newer PEPSKit revision.

3. CTMRG observables are useful, including custom 5-site terms.

   `expectation_value(peps, O::LocalOperator, env::CTMRGEnv)` worked for:

   - 1-site local operator on a `(1, 1)` unit cell.
   - nearest-neighbor 2-site local operator on a periodic `(1, 1)` unit cell.
   - custom 5-site star operator on a `(3, 3)` unit cell.

   Important detail: use tuple site keys, not vector site keys. Vector keys construct a `LocalOperator`, but PEPSKit `0.7.0` CTMRG observable dispatch expects tuple keys such as `(CartesianIndex(2, 2), CartesianIndex(2, 3), ...)`.

## Recommendation

Choose **C: implement a thin adapter between dense PXP gates and TensorKit/PEPSKit**.

Do not build the whole backend directly on PEPSKit yet. PEPSKit is strong enough for `InfinitePEPS`, CTMRG environments, and 1-/2-/5-site local measurements, but its registered simple-update/time-evolution path is not generic enough for the square-star 5-site PXP gate.

The practical architecture should be:

- Keep the existing dense 32x32 PXP gate as the source of truth.
- Add a small TensorKit adapter that maps dense gate/operator arrays into `TensorMap`s with the intended site ordering.
- Use PEPSKit `InfinitePEPS` and CTMRG `expectation_value` as the measurement/reference backend.
- Keep custom star-gate update logic outside PEPSKit's built-in `SimpleUpdate` unless/until PEPSKit exposes a generic local-circuit update that supports arbitrary connected multi-site gates.

Option B is a reasonable fallback if the project wants to avoid a PEPSKit dependency in the main package for now, but it leaves CTMRG validation as external tooling. Option A is not recommended because the 5-site update is the central PXP operation and is not supported by PEPSKit's current registered simple-update interface.

## References

- PEPSKit `InfinitePEPS` constructors and square unit-cell convention: https://quantumkithub.github.io/PEPSKit.jl/dev/lib/lib/
- PEPSKit `LocalOperator` and `LocalCircuit` docs: https://quantumkithub.github.io/PEPSKit.jl/dev/lib/lib/
- PEPSKit CTMRG `expectation_value` docs: https://quantumkithub.github.io/PEPSKit.jl/dev/lib/lib/
- PEPSKit `trotterize` support limitation in dev docs: https://quantumkithub.github.io/PEPSKit.jl/dev/lib/lib/
