# TriangularPEPSDynamics.jl

High-performance Julia foundations for PEPS-based dynamics on translationally invariant triangular lattices.

The first implemented layer provides:

- triangular axial geometry and 7-site star neighborhoods;
- 7-color non-overlapping star schedules;
- dense spin-1/2 operators;
- dense PXP, projected PXP, and cluster/stabilizer star gates;
- analytically solvable true-2D benchmark helpers, including a narrow stabilizer benchmark helper.

The next implementation layer will add native triangular iPEPS state containers, Simple Update truncation, observables, and ScarFinder evolve-project drivers.

## Test

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```
