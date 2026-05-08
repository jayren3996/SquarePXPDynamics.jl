# Triangular PEPS Dynamics Design

## Purpose

Build high-performance Julia PEPS-based tensor-network tooling for translationally invariant triangular lattices. This tooling is an internal subpackage-style layer for the 2D triangular-lattice PXP ScarFinder project, not a standalone general-purpose PEPS package. It should use the repository root `Project.toml` and root Julia environment.

The first target is dense spin-1/2 real-time and imaginary-time evolution for a two-dimensional PXP-type model on the triangular lattice, where each Hamiltonian term acts on a 7-site star: one central site and its six nearest neighbors.

The implementation should be research-grade from the start: reproducible, benchmarkable, and structured so that higher-accuracy truncation schemes can be added without rewriting the lattice or gate layers. However, scope should remain ScarFinder-driven: do not add broad PEPS-library abstractions unless they directly support triangular PXP dynamics or ScarFinder workflows.

The package should also be suitable as a backend for ScarFinder-style searches for low-entanglement scar trajectories in the 2D triangular PXP model. This means the evolution and truncation pipeline must be callable as an explicit "evolve then project back to PEPS" map, not only as a forward simulator.

## Scope

The first build targets:

- Julia package scaffold using `ITensors.jl` as the tensor backend.
- Dense spin-1/2 physical sites with no conserved quantum numbers.
- Native triangular iPEPS tensors with six virtual legs and one physical leg.
- Translationally invariant unit cells, starting with 1-site and 3-site cells.
- Real-time evolution for quench and dynamical simulations.
- Imaginary-time evolution for ground-state search, state preparation, and validation.
- 7-site star Hamiltonian terms suitable for triangular-lattice PXP dynamics.
- Simple Update for the first working truncation path.
- Neighborhood Tensor Update as the main production truncation scheme.
- ScarFinder support: fixed-bond-dimension projection loops, target-energy correction, blockade-constraint diagnostics, and random low-entanglement initial-state ensembles.

Out of scope for the first milestone:

- Fermionic PEPS.
- Non-Abelian or U(1) symmetric tensors.
- Production-grade Full Update with CTMRG.
- Arbitrary lattice graphs.
- Finite PEPS boundary contractions.

## Scientific Model Target

The eventual model is a triangular-lattice two-dimensional PXP system. In a spin-1/2 basis, each local term has the form of a spin flip on a central site dressed by projectors on its six neighboring sites. A representative term is

```text
h_c = P_1 P_2 P_3 P_4 P_5 P_6 X_c
```

with model-specific choices for whether `P` projects onto the ground, unexcited, or allowed Rydberg state. The code should not hard-code one convention. Instead, it should expose a star-term builder that accepts the local flip and neighbor projector operators.

The 7-site term structure is central to the design. The gate, scheduling, and truncation interfaces should all understand a star update instead of reducing everything to nearest-neighbor gates.

## Architecture

The PEPS code should remain inside this repository as internal modules under the existing root Julia project. Do not create a nested package or independent PEPS-specific environment.

Core modules:

- `Geometry`: triangular lattice directions, unit-cell coordinates, bond labels, and 7-site star neighborhoods.
- `States`: iPEPS tensor containers, bond spectrum containers, tensor initialization, gauge normalization, and checkpoint serialization.
- `Models`: local spin operators and Hamiltonian-term builders, including triangular PXP star terms.
- `Gates`: construction and compression of real-time and imaginary-time 7-site gates.
- `Schedules`: first-order and second-order Trotter schedules using a coloring of star centers.
- `Updates`: Simple Update and Neighborhood Tensor Update implementations.
- `Observables`: one-site, bond, star, norm, and diagnostic measurements.
- `ScarFinder`: iterative evolve-project search drivers, target-energy correction, candidate ranking, and ensemble orchestration.
- `Benchmarks`: microbenchmarks for contractions, SVDs, update kernels, and scaling with bond dimension.
- `SolvableModels`: analytically tractable 2D benchmark Hamiltonians with exact real-time and imaginary-time checks.

## Data Model

Each PEPS tensor represents one triangular-lattice site and has indices:

```text
(phys, n0, n1, n2, n3, n4, n5)
```

where `phys` has dimension 2 and the six virtual legs follow a fixed counterclockwise direction convention.

The state container stores:

- unit-cell shape and tensor map;
- virtual bond indices by direction and cell coordinate;
- optional diagonal bond spectra for Simple Update;
- metadata for physical dimension, bond dimension, time, time step, and backend device;
- stable labels for reconstructing neighboring tensors and star neighborhoods.

The initial implementation uses dense `ITensor` objects. The design should avoid APIs that would prevent later block-sparse or GPU tensor storage.

## Time Evolution

Real-time evolution applies

```text
U_c(dt) = exp(-im * dt * h_c)
```

to each 7-site star. Imaginary-time evolution applies

```text
G_c(dτ) = exp(-dτ * h_c)
```

followed by normalization and convergence checks. Same-color star centers must be disjoint so all updates in one color layer can be applied independently.

Use second-order Trotter by default:

```text
S2(dt) = S1(dt / 2, colors = 1:K) * S1(dt / 2, colors = K:-1:1)
```

where `K` is the number of star-color layers. The implementation should include first-order schedules for debugging and performance comparisons.

The triangular-lattice 7-site star coloring must be generated and tested explicitly. The scheduler should verify that no two stars in a layer overlap for the configured unit cell.

Imaginary-time runs should support step-size schedules, for example a coarse-to-fine sequence of `dτ` values. The state should be normalized or gauge-conditioned after each layer or full step so tensor norms do not dominate convergence behavior.

## Gate Compression

For spin-1/2, a dense 7-site gate has dimension `2^7 x 2^7`, which is small enough to form for early validation. However, applying it naively to a PEPS cluster creates large intermediate tensors.

The gate layer should support two representations:

- `DenseStarGate`: direct dense 7-site gate, used for tests and small bond dimensions.
- `CompressedStarGate`: tree or MPO-like factorization over the central site and six neighbors, used for production runs.
- `ProjectedStarGate`: constrained PXP gate that applies the local blockade projector after the dense or compressed gate.

The PXP star gate is structured because it is a projector-controlled central flip. The implementation should exploit this structure where possible, avoiding unnecessary generic exponentiation for every time step or imaginary-time step. For constrained PXP runs, projected gates are the default.

## ScarFinder Compatibility

ScarFinder iterates a map of the form:

```text
|psi> -> Project_M[exp(-i H DeltaT) |psi>]
```

where `M` is a low-entanglement variational manifold. In this package, `M` is the fixed-bond-dimension translationally invariant triangular PEPS manifold. The projection is implemented by the PEPS truncation backend, initially Simple Update and later NTU.

The ScarFinder driver should support:

- projection interval `DeltaT`, distinct from the microscopic Trotter step `dt`;
- repeated evolve-project iterations;
- random dense PEPS initialization at fixed unit cell and bond dimension;
- product-state and low-bond-dimension seeds such as triangular charge-density-wave patterns;
- candidate ranking by entanglement growth, revival fidelity, local-observable periodicity, and constraint violation;
- target energy density or target energy per site;
- optional imaginary-time energy correction after projection, using projected imaginary-time gates for constrained PXP runs;
- stopping rules for convergence, divergence, or excessive constraint violation.

Energy correction follows the ScarFinder idea: after projection, measure `DeltaE = <H> - E_target`; apply a short imaginary-time correction within the same fixed-`D` PEPS manifold, and keep the correction step that best returns the state toward `E_target`. The implementation should expose this as a policy because imaginary-time correction can be approximate in PEPS and may fail for some initial states.

The PXP blockade constraint is especially important. The default enforcement strategy should be local projected evolution: replace each unconstrained local update by a projected gate,

```text
U_eff(dt) = P_blockade U(dt)
G_eff(dτ) = P_blockade G(dτ)
```

where `P_blockade` projects the affected local neighborhood back into the allowed subspace after applying the real-time or imaginary-time gate. This is the simplest first-line constraint enforcement for PXP dynamics and should be integrated directly into the gate layer.

Bond truncation can still generate weight on forbidden nearest-neighbor `|up up>` configurations outside the updated star or through approximate projection. The package should therefore include:

- nearest-neighbor blockade violation observables;
- `ProjectedStarGate` and projected imaginary-time gate variants;
- optional post-projection cleanup using additional local projected gates;
- optional non-Hermitian penalty evolution with `-im * mu * sum_<ij> |up up><up up|_ij` only as a diagnostic or fallback mode for ScarFinder searches;
- diagnostics that reject candidates whose constraint violation exceeds a configured tolerance.

For the triangular lattice, ScarFinder should not assume bipartite order or known scar states. The design must keep unit cells flexible and support ensembles over different unit-cell choices, including 1-site, 3-site, 4-site, and larger cells when computationally feasible.

## Solvable 2D Benchmarks

The main nontrivial true-2D benchmark should be a triangular-lattice cluster/stabilizer star Hamiltonian,

```text
K_i = X_i prod_{j in NN(i)} Z_j
H_cluster = -sum_i K_i
```

Each `K_i` has the same 7-site star support as the triangular PXP term: one central `X` operator and six neighbor operators. Unlike PXP, all stabilizers commute. The exact real-time evolution therefore factorizes as

```text
U(t) = prod_i exp(i t K_i)
```

and Trotter ordering should introduce no physical Trotter error. Any discrepancy comes from tensor truncation, gate compression, or implementation mistakes. This model should be used as the primary nontrivial 2D regression benchmark.

The benchmark layer should also include:

- `H = 0` and identity gates for no-op update tests.
- Single-site transverse-field dynamics, `H = sum_i X_i`, for exact product-state rotations.
- Commuting triangular Ising dynamics, `H = sum_<ij> Z_i Z_j`, for diagonal entangling gates.
- Commuting 7-site diagonal star dynamics, `H = sum_i Z_i prod_{j in NN(i)} Z_j`, for star-gate imaginary-time checks.
- Exact dense 7-site PXP star tests against a `128 x 128` matrix exponential.
- Projected PXP star-gate tests verifying that `P_blockade U` removes forbidden local configurations.
- Small-cluster exact diagonalization tests for short-time triangular PXP dynamics.
- Quasi-1D triangular-cylinder ScarFinder reference tests where useful, while keeping true 2D benchmarks separate.

For the cluster/stabilizer benchmark, the package should expose exact observables such as stabilizer expectations `<K_i>`, simple Pauli-string expectations, and norm checks. Since the stabilizer evolution is exactly solvable, it should be the first test used to distinguish algorithmic truncation error from physics.

## Truncation Strategy

### Simple Update

Simple Update is the first implementation target. It should:

- absorb local bond spectra around the 7-site star;
- apply the compressed star gate;
- decompose the updated cluster back into seven site tensors;
- truncate each affected bond to maximum bond dimension `D`;
- update local bond spectra;
- report discarded weight and local truncation diagnostics.

This path is the fast baseline and the main vehicle for early testing.

### Neighborhood Tensor Update

Neighborhood Tensor Update is the production target. It should:

- construct a small exactly contractible neighborhood around the updated star;
- define a Hermitian positive local error measure for the truncated tensors;
- optimize the replacement tensors using alternating least squares or local linear solves;
- reuse Simple Update output as the initial guess;
- parallelize across disjoint stars in the same Trotter color layer.

NTU should be written as a separate update backend sharing the same `apply_gate!` interface as Simple Update.

### Future Full Update

Full Update with CTMRG is a later accuracy layer. The first design should leave room for:

- triangular CTMRG environments;
- environment reuse across time steps;
- variational truncation with a full infinite environment.

No Full Update implementation is required in the first milestone.

## Observables And Diagnostics

The package should measure:

- local magnetizations and spin-flip expectation values;
- nearest-neighbor correlators;
- 7-site PXP star energy terms;
- blockade-constraint violation on all nearest-neighbor bonds;
- return fidelity or local proxy fidelities for scar revivals;
- entanglement-spectrum diagnostics used to rank ScarFinder candidates;
- norm drift or local normalization diagnostics;
- discarded weights per update;
- entanglement spectra from Simple Update bonds;
- wall-clock timing per Trotter layer and per truncation kernel.

For real-time dynamics, norm and local observable stability are mandatory regression diagnostics. For imaginary-time dynamics, energy density, step-to-step state change, truncation error, and convergence with decreasing `dτ` are mandatory diagnostics.

## Performance Plan

Initial high-performance choices:

- use `ITensors.jl` named indices to prevent contraction mistakes;
- cache contraction orders and intermediate index groupings;
- minimize index reconstruction inside hot update loops;
- use in-place state updates where practical;
- keep update kernels type-stable through function barriers;
- benchmark dense CPU first before adding GPU paths;
- parallelize independent star updates within a Trotter color.

GPU support should be added after the CPU implementation is correct and benchmarked. The design should keep tensor movement explicit so the user can control CPU/GPU placement.

## Testing Plan

Unit tests:

- triangular direction and opposite-direction consistency;
- star-neighborhood construction;
- Trotter color layers contain non-overlapping stars;
- PXP local operator and star Hamiltonian construction;
- blockade projector construction, projected-gate construction, blockade violation observable, and optional penalty-term construction;
- triangular cluster/stabilizer star construction and commutation checks;
- dense 7-site gate agrees with exact matrix exponentiation;
- gate compression reconstructs the dense gate within tolerance;
- PEPS tensor index order and bond matching;
- Simple Update preserves expected dimensions and truncates to `D`.

Integration tests:

- `D=1` product-state evolution against exact local-star behavior where possible;
- one or two Trotter steps at small `D` with deterministic seeds;
- real-time norm/diagnostic regression tests;
- imaginary-time energy monotonicity or convergence smoke tests where the Hamiltonian and update approximation make that expectation meaningful;
- coarse-to-fine imaginary-time schedule regression tests;
- nontrivial true-2D triangular cluster/stabilizer dynamics with zero expected Trotter error;
- ScarFinder smoke test on a tiny fixed-`D` manifold: evolve, project, energy-correct, and rank several seeded candidates deterministically.

Benchmarks:

- one star update as a function of `D`;
- one full Trotter step as a function of unit-cell size;
- Simple Update versus NTU wall-clock and memory use.

## Milestones

1. Package scaffold, geometry, state container, and PXP star-term builder.
2. Dense 7-site gate construction, solvable benchmark models, and exact small-gate tests.
3. First-order and second-order Trotter schedules with star coloring checks.
4. Simple Update implementation for dense triangular iPEPS.
5. Observables, real-time diagnostics, and imaginary-time convergence diagnostics.
6. Compressed star-gate backend.
7. ScarFinder driver with ensemble initialization, candidate ranking, target-energy correction, and blockade diagnostics.
8. NTU backend.
9. CPU performance tuning and benchmarks.
10. Optional GPU backend experiments.

## Open Technical Decisions

- Exact unit-cell defaults: 1-site is minimal, but 3-site may be more natural for triangular-lattice ordering and PXP dynamics.
- Exact PXP projector convention and sign convention.
- Exact stabilizer benchmark conventions, including Pauli normalization and sign of `H_cluster`.
- Default imaginary-time step-size schedule and stopping criteria.
- Default ScarFinder ranking objective for 2D triangular PXP: entanglement growth, fidelity revivals, local-observable periodicity, or a weighted combination.
- Whether projected gates are sufficient at target bond dimensions or whether a structurally constrained PEPS ansatz is needed later.
- Whether the first NTU neighborhood should include only the updated star shell or an additional ring of neighboring tensors.
- Whether to use `ITensor` throughout hot loops or lower selected kernels to `NDTensors.Tensor` after correctness is established.

## Acceptance Criteria

The first complete version is acceptable when it can:

- initialize a dense spin-1/2 triangular iPEPS at configurable bond dimension `D`;
- construct a triangular-lattice PXP 7-site star gate;
- run real-time second-order Trotter evolution for multiple time steps;
- run imaginary-time second-order Trotter evolution with normalization and convergence diagnostics;
- truncate bonds through Simple Update with reported discarded weights;
- compute local and star observables;
- run a ScarFinder-style fixed-`D` evolve-project loop with target-energy correction and candidate diagnostics;
- measure and report nearest-neighbor blockade violations;
- pass exact 7-site gate tests and geometry tests;
- run a true-2D triangular cluster/stabilizer benchmark with zero expected Trotter error and reported truncation error;
- provide benchmark output for one-step runtime and memory scaling.
