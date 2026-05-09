# Current PEPS Evolution Solution

This note records the current repository solution for PEPS evolution in the
triangular-lattice PXP ScarFinder work. The PEPS layer remains internal project
tooling, not a standalone tensor-network package.

## Algorithm In Use

The evolution algorithm is a local star-gate evolve-project scheme on a
translationally invariant triangular iPEPS:

```text
state
  -> apply scheduled 7-site triangular PXP star gates
  -> project the updated local state back to the current PEPS manifold
  -> collect truncation, norm, local-observable, and blockade diagnostics
```

For the PXP model, each local gate acts on the dense 7-site star ordered as:

```text
center first, then the six triangular nearest neighbors in direction order
```

Real-time and imaginary-time gates are built as dense `128 x 128` matrices:

```text
U(dt) = exp(-im * dt * H_star)
G(dτ) = exp(-dτ * H_star)
```

Constrained PXP evolution uses local projected gates:

```text
U_eff(dt) = P_blockade * U(dt)
G_eff(dτ) = P_blockade * G(dτ)
```

This keeps blockade enforcement explicit at the local gate level. Truncation can
still leak outside the constrained manifold, so the evolution path also reports
nearest-neighbor and dense-star blockade screening diagnostics.

## Trotter Scheduling

The triangular lattice uses a 7-color schedule for radius-1 star centers. A
color layer contains disjoint star updates, so layer updates can be treated as
independent in the translational schedule.

The supported schedules are:

```text
first order:
  colors 1, 2, 3, 4, 5, 6, 7 with scale 1.0

second order:
  colors 1, 2, 3, 4, 5, 6, 7 with scale 0.5
  colors 7, 6, 5, 4, 3, 2, 1 with scale 0.5
```

The Hamiltonian-based evolution API constructs a fresh gate for each layer
using `dt * layer.scale`. The older prebuilt-gate API converts a second-order
full-step gate to a half-step gate before applying the symmetric schedule.

## Implemented Update Backends

### Simple Update Path

`Simple Update` is the baseline projection strategy and the main correctness
path right now. The current implementation supports:

- identity gates as no-ops;
- factorized product gates `u_1 ⊗ ... ⊗ u_7`, applied directly to physical
  indices;
- dense non-product 7-site gates only for `D = 1` product iPEPS, where the star
  update is applied exactly to the dense 7-site product vector and then projected
  back to one-site product factors.

The `D = 1` dense-star path is the current oracle for projected PXP evolution.
It is intentionally narrow: it verifies gate construction, ordering, schedules,
unit-cell wrapping, and diagnostics before the general fixed-`D` truncation path
is trusted.

For `D > 1`, the current Simple Update path uses a local product projection:
it extracts representative physical profiles, applies the dense star gate to
their product vector, projects the result back to representative one-site
profiles, and updates the tensors by local physical maps. This preserves the
existing virtual bond dimensions up to `dynamics_maxdim`. It is a local
approximation; the later production target is still an SVD/HOSVD, ring SU, or
NTU-style update with stronger truncation diagnostics.

Important boundary: this does not allow bond dimension to grow without
truncation throughout a full projection interval. iPEPS gate application still
needs a projection/truncation step at each scheduled gate layer. The current
two-tier scheme only separates the larger working PEPS dimension used during
Trotter dynamics from the smaller ScarFinder manifold dimension used after a
projection interval.

## Hard Truncation

The explicit hard truncation API is:

```julia
truncate_state!(state, target_maxdim) -> StateTruncationDiagnostics
```

It operates on the stored virtual bond spectra and tensor indices:

```text
for each unique logical bond
  deduplicate opposite-direction aliases
  sort the shared lambda vector by magnitude
  keep the top target_maxdim entries
  allocate a fresh virtual Index of dimension target_maxdim
  slice the two adjacent site tensors to the kept rows
  update the bond Index entries
  replace the lambda vector with the kept values
  record discarded lambda weight
end
```

For one-site unit cells, opposite directions can share lambda vectors even when
the tensor uses distinct virtual `Index` objects. The truncation code therefore
uses object identity (`===`) before writing the opposite-direction lambda entry,
so it preserves true lambda aliases without assuming every opposite direction
also shares the same tensor index.

## ScarFinder Driver Status

The current ScarFinder loop samples seed states and repeatedly calls projected
PXP evolution:

```text
for each seed
  initialize product or random triangular iPEPS
  repeat niterations
    run projected PXP for projection_interval Trotter steps at dynamics_maxdim
    if scar_maxdim < dynamics_maxdim
      hard-truncate the state to scar_maxdim
    end
    record diagnostics
  rank by discarded weight, blockade violation, and lambda entropy proxy
```

`ScarFinderConfig` keeps separate dimensions:

- `dynamics_maxdim`: the working maximum bond dimension used during Trotter
  evolution;
- `scar_maxdim`: the hard truncation dimension for the ScarFinder manifold.

For backward compatibility, the positional `maxdim` argument sets both values
unless `scar_maxdim` is supplied as a keyword. This two-tier split does not
delay all truncation until the end of an interval; iPEPS gate application still
projects at each gate layer, but it can project to the larger
`dynamics_maxdim` before the ScarFinder hard truncation step.

The active update backend is `apply_star_gate_simple_update!`. Full Update,
Fast Full Update, and NTU are not shipped active backends in this solution;
they remain future alternatives if the Simple Update/ring-update route is not
accurate enough.

Current candidate diagnostics include:

- per-layer discarded weights;
- maximum and mean bond dimension;
- lambda spectrum summaries;
- tensor norms;
- local `Z`, `X`, and `|up><up|` expectations;
- local blockade screening.

Energy targeting and imaginary-time energy correction are still planned
ScarFinder features, not a complete production workflow.

## Implementation Map

The main code paths are:

- `src/Gates.jl`: dense real-time, dense imaginary-time, and projected star
  gate construction;
- `src/Schedules.jl`: first- and second-order 7-color layer schedules;
- `src/Evolution.jl`: `evolve_step!`, `projected_pxp_step!`,
  `imaginary_projected_pxp_step!`, and `run_projected_pxp!`;
- `src/SimpleUpdate.jl`: identity/product-gate updates and the `D = 1`
  dense-star product-state oracle, plus the local product-projection `D > 1`
  Simple Update path;
- `src/States.jl`: iPEPS containers, unit-cell wrapping, and hard state
  truncation;
- `src/Observables.jl`: local and dense-star blockade diagnostics;
- `src/ScarFinder.jl`: seed loop, repeated evolve-project calls, and candidate
  ranking with `dynamics_maxdim` / `scar_maxdim` separation.

## Current Boundary

What is implemented now:

- dense projected 7-site PXP gates;
- first- and second-order 7-color star schedules;
- real- and imaginary-time projected PXP step helpers;
- exact `D = 1` projected PXP product-state evolution;
- local product-projection `D > 1` Simple Update evolution;
- hard state truncation from `dynamics_maxdim` to `scar_maxdim`;
- local diagnostics and ScarFinder candidate ranking;

What remains future work:

- production SVD/HOSVD or ring Simple Update for non-product star gates;
- Neighborhood Tensor Update;
- CTMRG/full-environment observables;
- target-energy correction for ScarFinder;
- robust fixed-bond-dimension ScarFinder searches beyond prototype scale.
