# KagomePXPDynamics.jl

Internal Julia tooling for PEPS-based dynamics on the 2D kagome lattice, built for the 2D kagome-lattice PXP ScarFinder project.

## Status

In active development. The first implementation pass is the kagome PESS Simple Update prototype — see [`docs/superpowers/plans/2026-05-09-kagome-pess-pxp-dynamics.md`](docs/superpowers/plans/2026-05-09-kagome-pess-pxp-dynamics.md).
NTU follow-up is sketched in [`Notes/2026-05-09-kagome-pess-ntu-followup.md`](Notes/2026-05-09-kagome-pess-ntu-followup.md).

## Currently shipped

- Generic spin-1/2 operators (`src/SpinOps.jl`).
- Solvable benchmark helpers (`src/SolvableModels.jl`).

## Planned (per the active spec/plan)

- KagomeGeometry, KagomePESS state container, kagome PXP gate, 3-color schedule
- Simple Update on PESS via cluster-and-HOSVD decomposition
- Three-layer validation: kernel ED, finite torus integration (indicator), kagome stabilizer benchmark
- ScarFinder driver

## History

An earlier attempt targeted the triangular lattice. It was abandoned due to two structural blockers in cluster-and-split Simple Update at D>1 (sublattice aliasing in 3-site UC writeback and 2^7*D^18 cluster scaling). See git log before commit `eef0ead` for the triangular source tree.

## Test

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```
