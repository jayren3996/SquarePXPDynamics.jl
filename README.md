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

## Test

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```
