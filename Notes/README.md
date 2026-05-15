# Project Purpose

This repository is internal Julia tooling for studying square-lattice PXP dynamics with PEPS/iPEPS methods. The long-term goal is to support a ScarFinder workflow: generate low-entanglement candidate states, evolve them under constrained PXP dynamics, diagnose blockade leakage and truncation error, and rank promising scar-like trajectories.

The project has intentionally moved back to the square lattice. The immediate goal is not to solve every PEPS/iPEPS algorithmic problem at once; it is to build a clean, testable square-lattice baseline before revisiting harder geometries.

## Active Scope

- Square-lattice PXP dynamics.
- PEPS/iPEPS tooling only when it directly supports ScarFinder.
- Root Julia project only; no nested PEPS package or separate Julia environment.
- Local dense gates and local diagnostics first.
- Fixed-bond-dimension evolve-project loops before full-environment methods.

This code should not grow into a general tensor-network package. Avoid arbitrary graph support, broad Hamiltonian packaging, CTMRG infrastructure, GPU backends, or symmetry machinery unless a concrete ScarFinder need justifies them.

## Conventions

- Physical basis: `|up> = |0>`, `|down> = |1>`.
- Square coordinate: `(x, y)`.
- Direction order: right, up, left, down.
- Dense square star ordering: center first, then right, up, left, down.
- Local square PXP term:

```text
h_c = X_c P_down(right) P_down(up) P_down(left) P_down(down)
```

- Projected constrained gates should stay explicit:

```text
U_eff = P_blockade * exp(-i dt h_c)
G_eff = P_blockade * exp(-dt h_c)
```

`P_blockade` removes local output support with an up center adjacent to an up neighbor. Approximate PEPS projection can still leak outside the constrained manifold, so blockade diagnostics remain mandatory.

## Current Code Foundation

- `src/SpinOps.jl`: dense spin-1/2 operators and small Kronecker helpers.
- `src/SquareGeometry.jl`: square coordinates, nearest neighbors, 5-site stars, and a 5-color disjoint star schedule.
- `src/SquarePXP.jl`: dense 32x32 square-star Hamiltonian, blockade projector, real/imaginary gates, and projected gates.
- `src/SquarePEPS.jl`: minimal ITensors-backed square PEPS product-state container.
- `src/SquarePXPDynamics.jl`: public module and exports.

The current PEPS state container is only a clean starting point. It is not yet a production Simple Update, NTU, or ScarFinder implementation.

## Candidate Algorithm Path

The likely first working algorithm should be a local Simple Update style evolve-project loop on square PEPS.

1. Build a finite or unit-cell square PEPS with one physical index and four virtual bond indices per site.
2. Sweep over disjoint square-star color classes. The current 5-color rule is:

```text
color(x, y) = mod(x + 2y, 5) + 1
```

Same-color radius-1 stars are disjoint, so a first-order Trotter step can apply colors `1:5`. A second-order step can sweep `1:5` and then `5:1` with half-step layers.

3. For each center, form the local 5-site star cluster: center, right, up, left, down.
4. Absorb nearby bond weights or local gauges into the cluster.
5. Apply the dense projected PXP gate on the five physical legs.
6. Refactor the updated cluster back into five PEPS site tensors with fixed virtual dimension `D`.
7. Record truncation residuals, discarded weights, norm changes, and blockade diagnostics.

The first refactorization backend can be deliberately modest:

- Start with a product-state or `D=1` dense oracle path to verify gate signs and blockade projection.
- Add a local SVD/HOSVD-style split for small `D`.
- Use Simple Update bond spectra as the first approximation to the environment.
- Move to NTU only after local Simple Update diagnostics are trustworthy and a specific failure mode appears.

CTMRG/full update should be treated as later accuracy infrastructure, not as a prerequisite for the first ScarFinder-facing loop.

## Validation Ladder

Use small, concrete checks before running ScarFinder searches.

1. Dense 5-site gate tests:
   - Hamiltonian has size `32 x 32`.
   - `P_blockade` is Hermitian and idempotent.
   - projected gates remove forbidden local output support.
   - real-time gates compose correctly for one local term.

2. Product-state PEPS checks:
   - all-up/all-down product states have expected local amplitudes;
   - neighboring square sites share the same ITensor link index;
   - boundary legs are dimension one.

3. Local update checks:
   - identity gates leave tensors unchanged up to gauge;
   - `D=1` update agrees with direct dense product-state evolution;
   - truncation residuals are finite and deterministic.

4. Short-time dynamics checks:
   - compare small square clusters against exact diagonalization where feasible;
   - monitor local `Z`, blockade violation, norm drift, and discarded weight.

5. ScarFinder smoke path:
   - deterministic seed set;
   - repeated real-time evolve-project;
   - optional imaginary-time or projection cleanup;
   - candidate ranking by blockade violation, truncation residual, and low-entanglement proxy.

## Near-Term Milestones

1. Add square nearest-neighbor blockade diagnostics.
2. Add local observable helpers for product states and small PEPS states.
3. Implement the first `D=1` dense/product oracle update.
4. Implement a small-`D` square-star Simple Update or HOSVD refactorization path.
5. Add real- and imaginary-time projected PXP step drivers.
6. Add a minimal ScarFinder candidate loop once diagnostics and evolution are stable.

## Working Rule

Every new feature should answer one ScarFinder-facing question: does it help evolve, project, diagnose, or rank square-lattice PXP candidate states? If not, leave it out for now.
