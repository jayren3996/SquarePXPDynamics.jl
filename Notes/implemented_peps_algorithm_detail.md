# Implemented PEPS Evolution Algorithm Detail

This note records the algorithm that this repository has chosen and currently
implements for the triangular-lattice PXP ScarFinder PEPS workflow. It is based
on the Notes references and the actual code paths in `src/`.

The short version is:

```text
triangular iPEPS seed
  -> scheduled dense 7-site projected PXP star gates
  -> per-layer Simple Update-style local projection
  -> optional hard truncation from dynamics_maxdim to scar_maxdim
  -> local blockade, norm, lambda, and discarded-residual diagnostics
  -> ScarFinder candidate ranking
```

This is deliberately not a general PEPS package. The implementation is a
ScarFinder-facing internal tool for constrained triangular PXP dynamics.

## Reference Decision

The Notes references point to the standard PEPS time-evolution problem:
applying a local Trotter gate increases the effective tensor-network bond
dimension, so the state must be projected back to a fixed maximum `D` after
local gates.

The relevant algorithm families are:

- Simple Update: local, cheap, stable, and based on local bond weights or
  lambda spectra. This is the chosen baseline because it gives an immediate
  fixed-`D` evolve-project loop without requiring a full environment.
- Neighborhood Tensor Update: a later target. NTU uses an exactly contractible
  nearest-neighbor environment metric and sits between Simple Update and Full
  Update in cost and accuracy.
- Fast Full Update and CTMRG-based Full Update: future accuracy infrastructure.
  They require environment contraction, gauge fixing, regularization, and ALS
  stability work that is not needed for the first ScarFinder path.
- PESS/HOSVD-style triangular updates: useful reference material for
  simplex-aware triangular lattices, but this repository keeps the current PEPS
  ansatz with one physical leg and six virtual legs per triangular site.

The current code therefore chooses a project-local star-gate Simple Update
path first, with explicit blockade projection and detailed diagnostics. The
implemented `D > 1` dense-star path is a local product-projection approximation,
not a production SVD/HOSVD, ring-SU, NTU, or Full Update backend.

## Core Conventions

The implementation preserves these conventions throughout:

- physical basis: `|up> = |0>`, `|down> = |1>`;
- local physical dimension: `d = 2`;
- triangular coordinate: axial `Coord(q, r)`;
- triangular directions:

```text
1: ( 1,  0)
2: ( 0,  1)
3: (-1,  1)
4: (-1,  0)
5: ( 0, -1)
6: ( 1, -1)
```

- dense star ordering: center first, then neighbors in directions `1:6`;
- dense star Hilbert-space size: `2^7 = 128`;
- a site tensor has indices `(phys, n1, n2, n3, n4, n5, n6)`;
- opposite virtual directions are `1 <-> 4`, `2 <-> 5`, and `3 <-> 6`.

The triangular PXP blockade projector treats adjacent `|up>|up>` pairs as
forbidden. On the dense star, it checks the 12 local star edges:

```text
center-neighbor edges:
  (1,2), (1,3), (1,4), (1,5), (1,6), (1,7)

neighbor-ring edges:
  (2,3), (3,4), (4,5), (5,6), (6,7), (7,2)
```

This is intentionally stricter than only checking the center flip constraint:
the dense star projector screens the full local triangular star.

## State Representation

`TriangularIPEPS` stores a translational triangular iPEPS:

```julia
TriangularIPEPS(
    unitcell,
    phys_inds,
    bond_inds,
    tensors,
    lambdas,
)
```

The important implementation details are:

- `phys_inds[rep]` is the physical ITensors `Index` for a unit-cell
  representative.
- `bond_inds[(rep, d)]` is the virtual `Index` leaving representative `rep` in
  triangular direction `d`.
- Distinct representatives connected by the same logical bond share the same
  `Index` object.
- Self-loop bonds, such as the one-site unit cell wrapping back to itself, use
  distinct tensor `Index` objects so a tensor never repeats the same `Index`.
- `lambdas[(rep, d)]` stores the local diagonal bond spectrum. Opposite
  directions on the same logical bond share the same `Vector{Float64}` object.

The supported unit cells are:

- `OneSiteUnitCell`: every coordinate wraps to `Coord(0, 0)`;
- `ThreeSiteUnitCell`: wraps by `(q - r) mod 3`, useful for three-sublattice
  triangular structure;
- `SevenSiteUnitCell`: wraps by the seven-color star schedule,
  `star_color(c) = (q + 3r) mod 7 + 1`.

Seeds come from:

- `product_ipeps(uc, :up; D)`;
- `product_ipeps(uc, :down; D)`;
- `random_ipeps(uc, D; seed)`.

Product states put the requested local vector at the all-ones virtual-bond
corner. Random states fill every tensor entry with complex Gaussian data.

## Gate Construction

The PXP star Hamiltonian is built in `src/Models.jl` as:

```text
H_star = X_center kron P_down_neighbor_1 kron ... kron P_down_neighbor_6
```

with the star order `center, neighbor_1, ..., neighbor_6`. In code this is:

```julia
H = pxp_star_hamiltonian(projector_down(), pauli_x())
```

Dense gates are then constructed in `src/Gates.jl`:

```text
real time:      U(dt)  = exp(-im * dt * H_star)
imaginary time: G(dtau) = exp(-dtau * H_star)
```

Projected gates are left-projected by the full dense star blockade projector:

```text
real time:      U_eff(dt) = P_blockade * U(dt)
imaginary time: G_eff(dtau) = P_blockade * G(dtau)
```

This means the local gate removes forbidden output support. For inputs already
inside the local blockade subspace, `P_blockade * U` agrees with
`P_blockade * U * P_blockade`.

`src/Evolution.jl` also caches constructed dense gates by:

```text
(objectid(H), Float64(step), evolution, projected, objectid(projector or 0))
```

This avoids rebuilding identical `128 x 128` matrix exponentials during a
Trotter run.

## Trotter Scheduling

A radius-1 triangular star centered at `c` is colored by:

```text
star_color(c) = (c.q + 3c.r) mod 7 + 1
```

Centers with the same color have disjoint radius-1 stars. For translational
iPEPS, each color layer is represented by one canonical center:

```text
color_canonical_center(color) = Coord(color - 1, 0)
```

The implemented schedules are:

```text
first order:
  (1, scale 1.0), (2, scale 1.0), ..., (7, scale 1.0)

second order:
  (1, scale 0.5), (2, scale 0.5), ..., (7, scale 0.5),
  (7, scale 0.5), (6, scale 0.5), ..., (1, scale 0.5)
```

The Hamiltonian API constructs a fresh cached gate for each layer using:

```text
layer_step = dt * layer.scale
```

The prebuilt-gate API cannot know the Hamiltonian step directly, so for
second-order evolution it uses the matrix square root of the full-step gate as
the layer gate.

## Evolution Entry Points

The main evolution functions are:

```julia
evolve_step!(state, gate; order, update, maxdim, cutoff)
evolve_step!(state, H, dt; order, evolution, projected, maxdim, cutoff)
projected_pxp_step!(state, dt; order, maxdim, cutoff, evolution)
imaginary_projected_pxp_step!(state, dtau; order, maxdim, cutoff)
run_projected_pxp!(state, dt, nsteps; order, maxdim, cutoff, evolution)
```

Only `update = :simple` is accepted. Some update-shape keywords
(`chi`, `maxiter`, `tol`, `regularization`) are already threaded through the
public API, but `_update_config` currently rejects every backend except
`:simple` and returns no backend-specific configuration.

One projected PXP step is:

```text
function projected_pxp_step!(state, dt)
  H = pxp_star_hamiltonian(projector_down(), pauli_x())
  layer_diagnostics = []

  for layer in schedule_layers(order)
    gate = projected_gate(H, dt * layer.scale; evolution)
    center = color_canonical_center(layer.color)
    diag = apply_star_gate_simple_update!(state, gate, center; maxdim, cutoff)
    push!(layer_diagnostics, diag)
  end

  return _step_diagnostics(state, layer_diagnostics)
end
```

The important boundary is that projection happens at every scheduled gate
layer. The two-tier `dynamics_maxdim` / `scar_maxdim` split in ScarFinder does
not postpone all projection until the end of an interval; it only allows the
per-layer projection to keep a larger working dimension before a later hard
truncation.

## Simple Update Dispatch

`apply_star_gate_simple_update!` takes a dense `128 x 128` star gate and a
star center. It reports `SimpleUpdateDiagnostics`:

```julia
SimpleUpdateDiagnostics(
    discarded_weight,
    affected_bonds,
    output_bond_dims,
)
```

`affected_bonds` is currently the six outgoing bonds of the wrapped center
representative. `output_bond_dims` reports the corresponding post-update bond
dimensions.

The function chooses one of four paths:

1. identity or near-identity gate;
2. factorized product gate `u_1 kron ... kron u_7`;
3. dense non-product gate on a `D = 1` product iPEPS;
4. dense non-product gate on a `D > 1` iPEPS using local product projection.

### Identity Path

The gate is treated as identity when:

```text
norm(G - I_128) <= 1e-12 * max(norm(I_128), 1.0)
```

The state is unchanged. The diagnostic discarded weight is zero.

### Product-Gate Path

The code tries to factorize the dense star gate into seven `2 x 2` matrices.
It peels factors from right to left using SVD reshapes that match Julia's
`kron` storage order, then reverses the list so the final factors satisfy:

```text
G ~= kron(u_1, kron(u_2, ... kron(u_6, u_7)...))
```

If factorization succeeds, the factors are grouped by wrapped unit-cell
representative. When multiple star positions wrap to the same representative,
their factors must be scalar multiples of a common one-site operator; otherwise
the translational update would be inconsistent and an error is raised.

For each representative, the chosen factor is normalized to a canonical
Frobenius scale and applied only to the physical index:

```text
T_rep <- noprime(ITensor(u, prime(phys), phys) * T_rep)
```

No virtual bond index or lambda spectrum changes in this path.

### Dense `D = 1` Product-State Path

When the gate is not factorized and every virtual index has dimension one, the
code uses the dense 7-site product-state oracle:

```text
star = star_sites(center)
v_i = local two-component vector at wrapped star site i
psi = kron(v_1, v_2, ..., v_7)
phi = G * psi
```

It then tries to factorize the output vector `phi` as a product state:

```text
phi ~= kron(w_1, w_2, ..., w_7)
```

If exact product-state factorization succeeds, each wrapped representative is
assigned the corresponding one-site vector. If several star positions wrap to
the same representative, `_common_factor_for_positions` checks whether their
vectors are scalar multiples and averages the aligned factors; if they are not,
it falls back to the first factor.

If `phi` is entangled and cannot be factorized as a product state, the code
projects it back to one-site factors by reduced density matrices:

```text
rho_i = Tr_{all sites except i} |phi><phi| / <phi|phi>
w_i   = dominant eigenvector of rho_i
```

For a representative that appears multiple times in the same wrapped star, the
current implementation uses the first position in that representative group in
this non-factorized branch. That is a useful implementation detail: the `D = 1`
path is an oracle for dense gate ordering and projection behavior, but it is
still a product-manifold projection when the exact output leaves the product
manifold.

The final writeback rebuilds a normalized `D = 1` site tensor with the selected
local vector and all virtual indices fixed to index value `1`.

### Dense `D > 1` Local Product-Projection Path

For non-product dense gates at `D > 1`, the current implementation does not
perform true Simple Update truncation by SVD/HOSVD. Instead it uses a local
product projection in physical space while preserving all existing virtual
bond dimensions.

The guardrails are:

```text
maxdim must be supplied
maxdim >= 1
cutoff >= 0
every affected current bond dimension must be <= maxdim
```

If an affected bond already exceeds `maxdim`, the function errors because
dimension-changing star writeback is not implemented in this path.

The algorithm is:

```text
star = star_sites(center)
reps = wrapped representatives for the seven star positions

for each star position i:
  T = tensor for reps[i]
  reshape T as a 2 x (all virtual entries) matrix M
  rho = M * M'
  v_i = dominant eigenvector of rho / tr(rho)

psi = kron(v_1, ..., v_7)
phi = G * psi

for each star position i:
  rho_i = one-site reduced density matrix of phi
  target_i = dominant eigenvector of rho_i

for each unique representative rep:
  collect all target_i whose star position wraps to rep
  phase-align them to the first target
  average and normalize
  build a regularized physical map old_rep -> target_rep
  apply the map to the physical leg of T_rep
  normalize T_rep by its Frobenius norm
```

The regularized physical map is:

```text
oldn    = old / ||old||
targetn = target / ||target||
projector = oldn * oldn'
op = targetn * oldn' + eps_reg * (I - projector)
```

with `eps_reg = 1e-10` in the helper default. This maps the current dominant
physical profile toward the projected target while keeping a small component
on the orthogonal subspace.

The diagnostic `discarded_weight` in this path is a relative product-projection
residual, not an SVD discarded weight:

```text
projected = kron(target_1, ..., target_7)
scale = <projected|phi> / <projected|projected>
residual = ||phi - scale * projected||^2 / ||phi||^2
```

Current limitations of this path:

- no bond dimension growth;
- no SVD/HOSVD factorization of the dense 7-site output;
- no lambda update from a true truncation spectrum;
- no ring environment, neighborhood metric, CTMRG, or ALS solve;
- local physical profiles are extracted from single-site tensor density
  matrices, not from a contracted PEPS environment.

This path is useful for prototype ScarFinder ranking and diagnostics at
`D > 1`, but the note should not be read as claiming production PEPS accuracy.

## Hard Truncation

`truncate_state!(state, target_maxdim)` is a separate hard truncation step used
by ScarFinder when `scar_maxdim < dynamics_maxdim`.

It does not run after every gate layer by default. Instead, ScarFinder calls it
after a `projection_interval` block of projected PXP steps.

For each unique logical bond, the truncation algorithm:

```text
key = (rep, d)
opp_key = (wrap_coord(neighbor(rep, d)), opposite_direction(d))
lambda_old = state.lambdas[key]

sort lambda entries by descending abs(lambda)
keep top target_maxdim entries
drop the rest
discarded_weight = sum(abs2(drop)) / sum(abs2(lambda_old))

create a fresh virtual Index of dimension target_maxdim
slice the tensor at key to the kept rows
slice the opposite tensor consistently
replace lambda vectors while preserving true aliasing
record output dimension and discarded weight
```

The alias handling matters. In multi-representative bonds, the two sides may
share the same `Index` object. In one-site self-loop cases, the lambda vectors
can be aliased while the tensor indices are deliberately distinct. The code
therefore checks object identity for both bond indices and lambda vectors before
deciding whether to reuse or separately replace the opposite-side data.

This truncation is a mechanical bond-spectrum truncation. It is not a
variational reoptimization of the state after truncation.

## Step Diagnostics

`ProjectedPXPStepDiagnostics` reports:

```julia
layer_diagnostics
discarded_weights
max_bond_dim
mean_bond_dim
lambda_summaries
blockade_violation
tensor_norms
local_z
local_x
local_projector_up
```

The diagnostics are assembled after each scheduled projected PXP step.

`lambda_summaries` stores, for each stored bond key:

```text
min(lambda), max(lambda), norm(lambda)
```

Local observables are computed as single-site tensor contractions with a
trivial environment:

```text
<op>_local = (Tdag * op * T) / (Tdag * T)
```

This is exact for `D = 1` product states. For `D > 1`, it is a bounded local
screening diagnostic, not a full PEPS expectation value.

The nearest-neighbor blockade diagnostic is:

```text
local_blockade_violation(c, d) =
  clamp(<P_up>_c * <P_up>_neighbor(c,d), 0, 1)
```

`mean_blockade_violation` averages this over unit-cell representatives and six
directions. `dense_star_blockade_violation` separately builds a dense 7-site
product vector from dominant local physical profiles and evaluates the exact
dense blockade projector on that vector.

## ScarFinder Loop

`ScarFinderConfig` stores:

```julia
dt
projection_interval
niterations
maxdim
dynamics_maxdim
scar_maxdim
cutoff
unitcell
seed_count
blockade_tolerance
update
```

The positional `maxdim` is the dynamics dimension. For compatibility it is also
stored as `maxdim`, and `scar_maxdim` defaults to the same value unless supplied
as a keyword. The config enforces:

```text
dt > 0
projection_interval >= 1
niterations >= 0
dynamics_maxdim >= 1
scar_maxdim >= 1
scar_maxdim <= dynamics_maxdim
cutoff >= 0
seed_count >= 1
blockade_tolerance >= 0
update == :simple
```

The seed policy is deterministic:

```text
seed_index 1: product_down
seed_index 2: product_up
seed_index >= 3: random(seed + seed_index)
```

The candidate loop is:

```text
for seed_index in 1:seed_count
  kind, state = seed_state(seed_index)
  diagnostics = []
  failed_projection = false

  for iteration in 1:niterations
    try
      append!(diagnostics,
        run_projected_pxp!(
          state, dt, projection_interval;
          order = :second,
          maxdim = dynamics_maxdim,
          cutoff,
          evolution = :real,
          update = :simple,
        )
      )

      if scar_maxdim < dynamics_maxdim
        truncate_state!(state, scar_maxdim)
      end
    catch ArgumentError
      failed_projection = true
      break
    end
  end

  discarded = sum(all per-layer discarded weights)
  blockade = mean_blockade_violation(state, unit_cell_representatives(unitcell))
  entropy = lambda_entropy_proxy(state)
  score = discarded + blockade + entropy + failed_projection_penalty
  accepted = !failed_projection && blockade <= blockade_tolerance
end

return candidates sorted by:
  (score, blockade_violation, entanglement_proxy, seed_index)
```

The lambda entropy proxy deduplicates aliased lambda vectors by `objectid`,
normalizes `abs2(lambda)`, and averages the Shannon entropy over unique logical
bond spectra.

## What Is Implemented Now

The active shipped algorithm currently includes:

- dense `128 x 128` real- and imaginary-time star gates;
- explicit dense full-star blockade projector;
- first- and second-order 7-color triangular star schedules;
- dense gate caching for repeated Hamiltonian layers;
- one-site, three-site, and seven-site unit cells;
- product and random iPEPS seeds;
- identity and product-gate Simple Update writebacks;
- `D = 1` dense-star product-state oracle;
- `D > 1` local product-projection writeback for dense non-product gates;
- hard truncation from `dynamics_maxdim` to `scar_maxdim`;
- local tensor, lambda, observable, and blockade diagnostics;
- deterministic ScarFinder seed loop and ranking.

## What Is Not Implemented Yet

The current code does not yet implement:

- a true SVD/HOSVD 7-site star refactor for general `D`;
- ring Simple Update with meaningful lambda absorption and re-emission for
  non-product star gates;
- Neighborhood Tensor Update;
- CTMRG contraction;
- Full Update or Fast Full Update;
- environment-based observables;
- target-energy correction or imaginary-time energy targeting in ScarFinder;
- a production-quality finite-`D` ScarFinder search.

## Practical Reading Of The Algorithm

The safest way to describe the current solution is:

```text
We use dense projected PXP star gates and a local Simple Update-style
evolve-project shell. The `D = 1` path is the dense correctness oracle. The
`D > 1` path is a local product-profile projection that preserves virtual
dimensions and provides prototype diagnostics, but it is not yet a true PEPS
truncation algorithm. ScarFinder can run with separate dynamics and scar
dimensions, but projection still occurs at every scheduled star layer.
```

This wording matches both the Notes literature direction and the actual code.
It also keeps the next algorithmic upgrade clear: replace the `D > 1` local
product-projection writeback with an SVD/HOSVD/ring-SU update or NTU behind the
same `apply_star_gate_simple_update!`-style interface, then promote diagnostics
from local screening metrics to environment-aware quantities only when CTMRG or
another contraction backend is available.
