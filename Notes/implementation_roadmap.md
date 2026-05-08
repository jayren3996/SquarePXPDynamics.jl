# Implementation Roadmap From Literature

This note translates the literature review into project-local decisions.

## Algorithm Contract

The ScarFinder-facing map should be:

```julia
candidate = scarfind(
    model = TriangularPXP(...),
    manifold = FixedDPEPS(unit_cell, D),
    DeltaT,
    dt;
    update = :simple,
    projected = true,
    target_energy = ...,
    seeds = ...
)
```

Conceptual loop:

```text
for seed in seeds
  psi = initialize_low_entanglement_peps(seed, unit_cell, D)
  for k in 1:Nproject
    psi = evolve_projected_real_time(psi, H, DeltaT; dt, update)
    diagnostics = measure_energy_constraints_truncation(psi)
    psi = maybe_energy_correct_imag_time(psi, diagnostics, target_energy)
    reject early if constraints or energy correction fail
  end
  score candidate by revival, entanglement proxy, constraints, energy drift
end
```

## Local Gate Semantics

Preserve repo conventions:

- physical basis: `|up> = |0>`, `|down> = |1>`;
- dense star ordering: center first, then six triangular neighbors in direction order;
- central star Hamiltonian accepts explicit neighbor projector and center flip operators;
- projected gate is local and explicit.

Use two projector helpers:

```text
central_flip_projector: checks the six neighbors of the center.
full_star_blockade_projector: checks all 12 triangular-star local edges.
```

For projected dense gates:

```text
U = exp(-im * dt * H_star)
Ueff = P_full_star_blockade * U
```

For imaginary time:

```text
G = exp(-dtau * H_star)
Geff = P_full_star_blockade * G
```

Tests should cover:

- all valid local basis states pass `P_full_star_blockade`;
- all 12 invalid edge pairs are rejected;
- `P_blockade * U` has no forbidden output support;
- `P_blockade * U * P_blockade` agrees with `P_blockade * U` on constrained inputs;
- dense matrix exponential for a 7-site star matches direct construction.

## Scheduling

Use a 7-color triangular-lattice schedule for radius-1 star centers. Each active layer must satisfy:

```text
no two centers are equal, adjacent, or share a nearest neighbor
```

Second-order evolution should build half-step gates per layer:

```text
S2(dt) = colors(1:7, dt/2) followed by colors(7:1, dt/2)
```

Do not reuse a prebuilt full-step gate for second-order layers.

## Simple Update Target

Simple Update is the first fixed-D projection backend because it is the minimum viable PEPS projection for ScarFinder.

Required diagnostics:

- discarded weight per decomposed bond;
- max discarded weight across the star;
- updated bond spectra;
- norm change;
- local projection residual;
- optional condition/gauge warnings.

Implementation phases:

1. exact `D=1` dense-star product update;
2. dense 7-site apply on a small star cluster for general `D`;
3. star-to-site refactor by SVD/HOSVD path;
4. lambda absorption and update;
5. stable truncation tests using solvable 7-site gates.

## NTU Target

NTU should share the same high-level API as Simple Update:

```julia
apply_star_gate!(state, gate, center; update = :ntu, maxdim, cutoff)
```

NTU adds a local neighborhood metric around the star replacement. It should be introduced only after Simple Update produces reliable diagnostics and the local energy/blockade tests are in place.

Do not build CTMRG or full update before NTU unless a specific ScarFinder failure requires it.

## ScarFinder Candidate Ranking

Return a table-like record per candidate:

```text
seed
unit_cell
D
DeltaT
dt
Nproject
energy_density
target_energy_error
blockade_violation
max_discarded_weight
mean_discarded_weight
norm_drift
revival_score
observable_period
three_sublattice_density_contrast
entanglement_proxy
status
```

Early rejection conditions:

- blockade violation exceeds tolerance;
- energy correction repeatedly increases target-energy error;
- norm or truncation residual diverges;
- candidate collapses to a trivial low-energy state when target energy is finite;
- frequency filter rejects an attractor from a different unit-cell sector.

## Energy Correction

Implement the ScarFinder-style correction as an optional policy:

```text
DeltaE = E(psi) - E_target
dtau = DeltaE / n
for i in 1:n
  psi_i = projected_imaginary_time_step(psi_{i-1}, dtau)
  keep the state with energy closest to E_target
end
```

For early Simple Update, label energy estimates as local/star-energy estimates. Later NTU or environment contraction can replace the estimator.

## Benchmarks

Use the repo design report's benchmark hierarchy:

1. `H = 0` no-op gate;
2. one-site transverse-field product rotations;
3. dense 7-site PXP star against `128 x 128` matrix exponential;
4. full-star projected PXP mask;
5. triangular cluster/stabilizer star Hamiltonian;
6. small-cluster ED for short-time PXP;
7. ScarFinder pilot on small `D`, small unit cells.

## Non-Goals For Now

- nested Julia package or separate PEPS environment;
- arbitrary graph PEPS;
- CTMRG/full update before NTU;
- PESS ansatz refactor;
- broad Rydberg Hamiltonian package with long-range interactions;
- symmetry-aware or GPU tensor backend before dense correctness tests pass.
