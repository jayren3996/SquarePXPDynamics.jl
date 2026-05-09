# TriangularPEPSDynamics.jl

Internal Julia tooling for PEPS-based dynamics on translationally invariant triangular lattices, built for the 2D triangular-lattice PXP ScarFinder project.

This is not intended to become an independent general-purpose PEPS package. Keep it as a subpackage-style module inside this repository, using the root `Project.toml` and root Julia environment.

The implemented layer provides:

- triangular axial geometry and 7-site star neighborhoods;
- 7-color non-overlapping star schedules;
- dense spin-1/2 operators;
- dense PXP, projected PXP, and cluster/stabilizer star gates;
- analytically solvable true-2D benchmark helpers, including a narrow stabilizer benchmark helper;
- translational iPEPS containers for one-site and three-site unit cells;
- correctness-checked Simple Update paths for identity gates, site-product gates, `D=1` dense product-state projected PXP updates, and local product-projection `D>1` updates;
- projected real- and imaginary-time PXP step helpers with diagnostics;
- local and dense-star blockade screening diagnostics;
- an internal `ScarFinder` loop with separate dynamics and ScarFinder truncation dimensions, deterministic seed generation, repeated evolve-project iteration, tolerance flagging, and candidate ranking.

The next implementation layers should stay focused on ScarFinder needs: reliable constrained PXP evolution, fixed-bond-dimension projection, blockade diagnostics, and candidate-ranking workflows.

## Current Supported Workflows

The code supports dense 7-site projected PXP gates and PEPS evolution using the root Julia environment.

```julia
using TriangularPEPSDynamics

state = product_ipeps(OneSiteUnitCell(), :down; D = 1)
diag = projected_pxp_step!(state, 0.01; order = :second, maxdim = 1, cutoff = 1e-12)

history = run_projected_pxp!(state, 0.01, 4; order = :first, maxdim = 1)
```

For a minimal ScarFinder-facing search:

```julia
config = ScarFinderConfig(0.005, 1, 2, 1, 1e-12, OneSiteUnitCell(), 1, 1e-6)
candidates = scar_search(config; seed = 123)
ranked = rank_candidates(candidates)
```

## Simple Update Status

The `D=1` product-state path applies dense non-product projected star gates through a dense 7-site oracle and remains the exact regression path for product iPEPS.

For `D>1`, the current Simple Update path applies dense non-product star gates by projecting the updated local star back onto representative physical profiles while preserving the existing virtual bond dimensions up to `dynamics_maxdim`. This is a local approximation, not a replacement for a later SVD/HOSVD or ring/NTU-style update.

Lambda spectra are kept nonnegative and normalized with `norm(lambda) == sqrt(length(lambda))`.

## Diagnostics

Projected PXP step diagnostics report per-layer discarded weights, max/mean bond dimension, lambda summaries, tensor norms, local `Z`, local `X`, local `projector_up`, and blockade screening values.

Blockade diagnostics include:

- nearest-neighbor local screening via one-site expectations;
- aggregate unit-cell bond screening;
- dense-star blockade violation from dominant local physical profiles.

These PEPS diagnostics are exact for product states and bounded screening metrics for `D>1`; they are not substitutes for a contracted environment.

## ScarFinder Status

`ScarFinder` is currently internal to this root module. It provides deterministic seed handling, repeated projected-PXP evolve-project iterations, candidate diagnostics, blockade tolerance flagging, and deterministic ranking by discarded weight, blockade violation, and a lambda entropy proxy.

Target-energy correction, full search orchestration, NTU, and environment-based observables remain future work.

## Test

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```
