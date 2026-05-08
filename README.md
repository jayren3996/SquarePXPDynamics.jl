# TriangularPEPSDynamics.jl

Internal Julia tooling for PEPS-based dynamics on translationally invariant triangular lattices, built for the 2D triangular-lattice PXP ScarFinder project.

This is not intended to become an independent general-purpose PEPS package. Keep it as a subpackage-style module inside this repository, using the root `Project.toml` and root Julia environment.

The first implemented layer provides:

- triangular axial geometry and 7-site star neighborhoods;
- 7-color non-overlapping star schedules;
- dense spin-1/2 operators;
- dense PXP, projected PXP, and cluster/stabilizer star gates;
- analytically solvable true-2D benchmark helpers, including a narrow stabilizer benchmark helper.

The next implementation layers should stay focused on ScarFinder needs: reliable constrained PXP evolution, fixed-bond-dimension projection, blockade diagnostics, and candidate-ranking workflows.

## Current ScarFinder-Facing Boundary

The code supports dense 7-site projected PXP gates and an initial PEPS evolution path. The Simple Update backend is being built in stages: exact `D=1` non-product star updates first, then fixed-`D` SVD truncation with lambda updates. PEPS blockade diagnostics are currently local-environment diagnostics and should be treated as screening metrics until a stronger environment contraction is added.

## Test

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```
