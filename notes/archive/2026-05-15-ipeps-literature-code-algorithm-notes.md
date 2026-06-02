# Square iPEPS Literature, Codebase, And Algorithm Notes

Date: 2026-05-15

These notes combine external web/literature search, open-source repository survey, and local codebase inspection for the square-lattice PXP/iPEPS implementation plan.

## Local Repo Constraints

The current repo is a small Julia package named `SquarePXPDynamics` with one root module and four implementation modules:

- `src/SpinOps.jl`: dense spin-1/2 operators and Kronecker helpers.
- `src/SquareGeometry.jl`: square coordinates, nearest neighbors, 5-site stars, five-color scheduling, and existing simple unit-cell wrappers.
- `src/SquarePXP.jl`: dense 32x32 square-star PXP Hamiltonian, blockade projector, and projected real/imaginary gates.
- `src/SquarePEPS.jl`: finite PEPS product-state container backed by ITensors.

The entrypoint `src/SquarePXPDynamics.jl` re-exports the public API. The public docstring test in `test/test_public_docs.jl` means every new exported symbol needs a docstring.

Conventions that must not drift:

- Physical basis: `|up> = |0>` and `|down> = |1>`.
- Rydberg/excited state: `:up`, basis index `1`.
- Vacancy/unexcited state: `:down`, basis index `2`.
- Direction order: `:right`, `:up`, `:left`, `:down`.
- Dense square-star order: `(center, right, up, left, down)`.
- Local PXP term: `X_center * P_down(right) * P_down(up) * P_down(left) * P_down(down)`.
- Existing finite PEPS tensor index order: physical, left, right, up, down.

Implementation implication: add a separate iPEPS stack. Do not mutate the finite `SquarePEPSState` into a periodic/iPEPS state.

## Core Literature

### ScarFinder

Source: https://arxiv.org/abs/2504.12383

ScarFinder is a variational evolve-project loop. The paper introduces a framework that iteratively evolves and projects states within a low-entanglement variational manifold, and applies it to PXP dynamics. The arXiv page currently lists PRX Quantum 6, 040332 (2025).

Implementation implications:

- ScarFinder should not own tensor-index logic.
- The projection step in this repo should be the fixed-bond-dimension iPEPS truncation.
- Energy drift after projection is expected; the driver needs energy diagnostics and optional imaginary-time correction.
- Early scoring should rank by energy stability, blockade leakage, truncation error, and entanglement proxies before attempting expensive revival-quality checks.

### iPEPS Gauge Fixing And Full Update

Source: https://arxiv.org/abs/1503.05345

Phien, Bengua, Tuan, Corboz, and Orus discuss fast full update and local gauge fixing for iPEPS. The paper emphasizes that gauge fixing improves stability and accelerates ALS convergence in full-update-style algorithms.

Implementation implications:

- Do not call the simple-update Gamma-lambda form canonical.
- Use `:simple` or `:simple_update` as the first gauge label.
- Full-update gauge fixing requires an environment backend and should wait until CTMRG exists.
- Local gauge diagnostics should be added early so future environment-based gauge improvements have a baseline.

### CTMRG For iPEPS Contraction

Source: https://arxiv.org/abs/0905.3225

Orus and Vidal revisit infinite-lattice PEPS contraction using CTMRG. The abstract states that CTMRG provides the environment for iPEPS contraction and improves estimates near criticality relative to older transfer-matrix approaches.

Implementation implications:

- Simple-update observables are useful for development but not the final measurement layer.
- CTMRG should become the environment source for reliable energy, correlations, fidelity density, and full-update/gauge-fixed truncations.
- CTMRG is large enough to be its own milestone.

### Neighborhood Tensor Update

Source: https://arxiv.org/abs/2107.06635

Dziarmaga compares simple update, full update, and neighborhood tensor update (NTU). The abstract frames SU and FU as the two standard paradigms and NTU as an intermediate local-environment method.

Implementation implications:

- A QR-reduced local star update is a sensible first step toward environment-aware updates.
- If simple update fails for the PXP star gate, NTU-like nearest-neighbor environments may be the next algorithmic improvement before a full CTMRG/FU stack.

### PXP Scar Context

Useful sources:

- Bernien et al. Rydberg scar experiment: https://arxiv.org/abs/1707.04344
- Turner et al. quantum many-body scars: https://www.nature.com/articles/s41567-018-0137-5
- Emergent SU(2) / deformed PXP dynamics: https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.122.220603

Implementation implications:

- The local square-lattice Hamiltonian term naturally has star geometry, `h_i = X_i prod_j P_j`.
- The local gate is small enough to keep dense 32x32 matrices as test truth.
- Blockade leakage is a required diagnostic because projected local gates do not guarantee a truncated PEPS remains globally constrained.

## Open-Source Code References

### ITensors.jl

Source: https://github.com/ITensor/ITensors.jl

ITensors.jl provides named-index tensor algebra, contraction, and SVD/factorization primitives. The current repo already depends on ITensors 0.9.

Implementation implications:

- Use ITensors as the low-level tensor engine, not as a full PEPS algorithm library.
- Use named indices and explicit replacement of SVD-generated links to avoid index-ordering bugs.
- Keep dense gate-to-ITensor conversion tested on every computational basis vector.

### ITensorNetworks.jl

Source: https://github.com/ITensor/ITensorNetworks.jl

ITensorNetworks.jl is a general higher-dimensional tensor-network package in the ITensor ecosystem.

Implementation implications:

- It may provide design ideas for graph/network organization, but this repo should stay focused on square-lattice PXP rather than importing broad graph abstractions prematurely.

### PEPSKit.jl

Sources:

- Repo: https://github.com/QuantumKitHub/PEPSKit.jl
- Simple-update example: https://quantumkithub.github.io/PEPSKit.jl/stable/examples/heisenberg_su/

PEPSKit.jl is the closest Julia reference for iPEPS workflows. Its simple-update Heisenberg example initializes an `InfinitePEPS`, creates `SUWeight(peps)`, evolves with `SimpleUpdate`, monitors `|Delta lambda|`, then computes observables with CTMRG.

Implementation implications:

- The proposed Gamma-lambda design matches the reference pattern: tensors plus separate simple-update weights.
- Track link-weight changes and truncation errors as first-class diagnostics.
- Use decreasing time steps and convergence on bond spectra for imaginary-time workflows.
- For serious measurements, normalize tensors and converge a CTMRG environment after simple update.

### peps-torch

Source: https://github.com/jurajHasik/peps-torch

peps-torch is a Python iPEPS/CTM/AD optimization codebase with square and Kagome examples, CTMRG routines, complex tensors, and symmetry-aware variants through YASTN.

Implementation implications:

- Good reference for CTM layout, unit-cell tilings, and measurement APIs.
- Less directly useful for ITensor syntax because it is PyTorch-first.
- Do not import its broad optimization architecture into this repo before the simple-update engine works.

### YASTN fPEPS CTM

Source: https://yastn.github.io/yastn/fpeps/environment_ctm.html

YASTN documents CTM environments with local corner/edge tensors, update methods, measurements, and bond metrics.

Implementation implications:

- Useful API reference for `CTMEnvironment` layout and method naming.
- Its environment operations suggest future functions: `measure_1site`, `measure_2site`, `measure_nsite`, `transfer_matrix_spectrum`, and `bond_metric`.

### tensors.net Julia iPEPS

Source: https://www.tensors.net/j-peps

tensors.net hosts Julia iPEPS TEBD examples and index-ordering conventions.

Implementation implications:

- Useful as a compact TEBD/simple-update benchmark reference.
- Keep index-ordering diagrams and local tensor leg conventions explicit in tests.

## Algorithm Notes

### Unit Cells

Use a true periodic rectangular unit cell:

```julia
struct PeriodicSquareUnitCell <: SquareUnitCell
    Lx::Int
    Ly::Int
    reps::Vector{SquareCoord}
end
```

Use one-based representatives:

```julia
[SquareCoord(x, y) for y in 1:Ly for x in 1:Lx]
```

Reasons:

- Matches the existing finite PEPS coordinate convention.
- Avoids conflating the five-color schedule with the physical ansatz unit cell.
- Makes checkerboard parity and five-color compatibility explicit.

Initial compatibility rules:

- Five-color disjoint star sweeps: require `Lx % 5 == 0` and `Ly % 5 == 0`.
- Checkerboard ansatz in same rectangular cell: require even `Lx` and `Ly`.
- Default robust cell: `10 x 10`.

### Bond Keys And Link Storage

Store each undirected periodic bond once:

```julia
struct BondKey
    site::SquareCoord
    dir::Symbol
end
```

Canonical directions:

- `:right`
- `:up`

Map noncanonical directions:

- `(c, :left)` maps to `(neighbor(cell, c, :left), :right)`.
- `(c, :down)` maps to `(neighbor(cell, c, :down), :up)`.

Keep endpoint access with `link_indices[(c, dir)]` so tensor construction and local algorithms can ask for all four directional legs naturally.

### Gamma-Lambda Simple Gauge

Represent the network as:

```text
... Gamma -- lambda -- Gamma ...
```

State fields:

```julia
tensors::Dict{SquareCoord,ITensor}
physical_indices::Dict{SquareCoord,Index}
link_indices::Dict{Tuple{SquareCoord,Symbol},Index}
link_weights::Dict{BondKey,Vector{Float64}}
gauge::Symbol # :simple
```

Rules:

- Link weights are nonnegative and normalized.
- Absorb lambdas only for local computations.
- Deabsorb external lambdas after local updates using a safe inverse.
- Never silently invert a tiny lambda.

### ITensor PXP Gate

Dense source of truth:

```julia
projected_square_pxp_gate(step; evolution = :real)
```

ITensor wrapper:

```julia
square_pxp_gate_itensor(step, phys; evolution = :real, projected = true)
```

Index order:

```text
output: p_center', p_right', p_up', p_left', p_down'
input:  p_center,  p_right,  p_up,  p_left,  p_down
```

Tests should compare every computational basis vector against the dense matrix.

### QR-Reduced Star Update

Avoid the full raw patch:

```text
2^5 D^12
```

Use QR/factorization on the four leaves:

```text
leaf tensor -> Q_external * R_active
```

Then apply the gate only to:

```text
center Gamma + four R_active tensors
```

Sequentially split leaves back out with SVD, storing each new singular spectrum as the corresponding center-leaf lambda.

Required diagnostics:

```julia
StarUpdateInfo(
    center,
    max_truncerr,
    truncerrs,
    keptdims,
    min_lambda,
    norm_factors,
)
```

### Observables

Two layers:

- Simple-update observables for early diagnostics and tests.
- CTMRG observables for production measurements.

Minimum simple observables:

```julia
density_simple
sublattice_densities
blockade_violation_simple
energy_density_simple
mean_bond_entropy
max_bond_entropy
```

Important testing principle: computational-basis product states have zero PXP energy because `X` has zero diagonal expectation.

### ScarFinder

ScarFinder should be a thin orchestration module:

```text
evolve! -> normalize -> optional energy_correct! -> measure -> rank/log
```

It should not:

- construct low-level ITensor gate indices,
- absorb/deabsorb lambda weights,
- split tensors,
- own CTMRG implementation details.

Early ranking should use:

```text
energy error
blockade violation
mean/max bond entropy
max truncation error
```

Revival quality can be a separate validation trajectory after the engine produces stable candidate states.

## Practical Milestone Order

1. Periodic unit cells and `SquareIPEPSState`.
2. Link weights and ITensor gate conversion.
3. QR-reduced five-site star update.
4. Trotter evolution driver.
5. Simple observables and diagnostics.
6. ScarFinder orchestration.
7. CTMRG and gauge-fixed full update.

## High-Risk Bugs To Test Against

- Mixed physical convention, especially assuming `|1>` is Rydberg when this repo uses `|0>`.
- Reusing `FiveSiteSquareUC` as a physical iPEPS unit cell.
- Double-absorbing an internal lambda in the star update.
- Contracting external virtual legs into the update core and exploding scaling.
- Losing canonical link indices after ITensor SVD.
- Comparing raw tensors after gauge-changing operations.
- Returning observables for unsupported gauge/state cases instead of throwing a clear error.
- Treating simple-update observables as final quantitative results.

## Source Index

- ScarFinder: https://arxiv.org/abs/2504.12383
- iPEPS fast full update/gauge fixing: https://arxiv.org/abs/1503.05345
- CTMRG for iPEPS: https://arxiv.org/abs/0905.3225
- NTU: https://arxiv.org/abs/2107.06635
- Bernien et al. Rydberg scar experiment: https://arxiv.org/abs/1707.04344
- Turner et al. scars: https://www.nature.com/articles/s41567-018-0137-5
- Deformed PXP / SU(2): https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.122.220603
- ITensors.jl: https://github.com/ITensor/ITensors.jl
- ITensorNetworks.jl: https://github.com/ITensor/ITensorNetworks.jl
- PEPSKit.jl: https://github.com/QuantumKitHub/PEPSKit.jl
- PEPSKit simple-update example: https://quantumkithub.github.io/PEPSKit.jl/stable/examples/heisenberg_su/
- peps-torch: https://github.com/jurajHasik/peps-torch
- YASTN CTM docs: https://yastn.github.io/yastn/fpeps/environment_ctm.html
- tensors.net Julia iPEPS: https://www.tensors.net/j-peps
