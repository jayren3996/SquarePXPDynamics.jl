# Literature Review: Tensor-Network Algorithms For PEPS Evolution

Prepared for the repository at `/Users/ren/Codex/PEPs`.

## Scope

This note keeps only tensor-network algorithm literature that is directly useful
for fixed-bond-dimension PEPS evolution and projection in this repository. It
excludes scar, Rydberg-platform, and broad model-background papers.

The implementation target remains project-local triangular PEPS dynamics:

- dense local gates followed by fixed-`D` projection;
- real- and imaginary-time evolution;
- local truncation and residual diagnostics;
- a path from Simple Update to stronger local and environment-aware updates.

## Immediate Algorithm References

- Jiang, Weng, Xiang, "Accurate Determination of Tensor Network State of Quantum Lattice Models in Two Dimensions", PRL 101, 090603, 2008. DOI: https://doi.org/10.1103/PhysRevLett.101.090603, arXiv: https://arxiv.org/abs/0806.3719. Local PDF: `literature/simple_update_jiang_0806.3719.pdf`.
- Dziarmaga, "Time evolution of an infinite projected entangled pair state: a neighborhood tensor update", PRB 104, 094411, 2021. DOI: https://doi.org/10.1103/PhysRevB.104.094411, arXiv: https://arxiv.org/abs/2107.06635. Local PDF: `literature/ntu_ipeps_2107.06635.pdf`.
- Phien, Bengua, Tuan, Corboz, Orus, "The iPEPS algorithm, improved: fast full update and gauge fixing", PRB 92, 035142, 2015. DOI: https://doi.org/10.1103/PhysRevB.92.035142, arXiv: https://arxiv.org/abs/1503.05345. Local PDF: `literature/fast_full_update_1503.05345.pdf`.
- Czarnik, Dziarmaga, "Time Evolution of an Infinite Projected Entangled Pair State: an Algorithm from First Principles", PRB 98, 045110, 2018. DOI: https://doi.org/10.1103/PhysRevB.98.045110, arXiv: https://arxiv.org/abs/1804.03872. Local PDF: `literature/first_principles_ipeps_1804.03872.pdf`.

## Main Takeaways

All PEPS evolution algorithms face the same local bottleneck: applying a gate
increases the effective bond dimension, then the network must be projected back
to a chosen maximum `D`.

Simple Update is the right baseline for this repository:

- it is cheap and stable;
- it tracks diagonal bond spectra on virtual bonds;
- it gives useful local discarded-weight diagnostics;
- it avoids introducing CTMRG and full-environment conditioning too early.

The limitation is also clear: Simple Update treats the environment as local
bond weights, so it misses loop and long-range environment effects. That is
acceptable for the first correctness path, but it should not be treated as a
final accuracy benchmark.

Neighborhood Tensor Update is the best next algorithmic target. It optimizes
the replacement tensors with an exactly contractible local-neighborhood metric,
remaining local enough to share infrastructure with star-gate diagnostics while
capturing more of the surrounding tensor network than Simple Update.

Fast Full Update and first-principles iPEPS time evolution are important
reference points for a later variational projection backend. They bring CTMRG,
gauge fixing, regularization, and ALS stability concerns, so they should come
after the local PEPS update path is correct and after diagnostics show that
Simple Update or NTU is insufficient.

## Environment And Variational Infrastructure

- Orus, Vidal, "Simulation of two dimensional quantum systems on an infinite lattice revisited: corner transfer matrix for tensor contraction", PRB 80, 094403, 2009. DOI: https://doi.org/10.1103/PhysRevB.80.094403, arXiv: https://arxiv.org/abs/0905.3225. Local PDF: `literature/ctm_ipeps_orus_vidal_0905.3225.pdf`.
- Vanderstraeten et al., "Fast convergence of imaginary time evolution tensor network algorithms by recycling the environment", PRB 91, 125136, 2015. DOI: https://doi.org/10.1103/PhysRevB.91.125136, arXiv: https://arxiv.org/abs/1411.0391. Local PDF: `literature/environment_recycling_1411.0391.pdf`.
- Poilblanc, Mambrini, Alet, "Tensor network variational optimizations for real-time dynamics: application to the time-evolution of spin liquids", SciPost Phys. 15, 158, 2023. DOI: https://doi.org/10.21468/SciPostPhys.15.4.158, arXiv: https://arxiv.org/abs/2304.13184. Local PDF: `literature/variational_realtime_peps_2304.13184.pdf`.

These are not first-pass implementation targets, but they are directly relevant
to the current algorithm boundary:

- CTMRG is the natural contraction backend for full-environment observables,
  full update, and variational overlap objectives.
- Environment recycling is relevant if the project introduces repeated
  environment-aware imaginary-time or real-time updates.
- Variational real-time PEPS comparisons help define what a stronger projection
  backend should improve over Simple Update: overlap fidelity, longer stable
  real-time windows, and diagnostics that separate truncation error from
  entanglement growth.

## Update-Algorithm References

- Verstraete, Cirac, "Renormalization algorithms for Quantum-Many Body Systems in two and higher dimensions", arXiv: https://arxiv.org/abs/cond-mat/0407066. Local PDF: `literature/peps_origin_0407066.pdf`.
- Jordan, Orus, Vidal, Verstraete, Cirac, "Classical Simulation of Infinite-Size Quantum Lattice Systems in Two Spatial Dimensions", PRL 101, 250602, 2008. DOI: https://doi.org/10.1103/PhysRevLett.101.250602, arXiv: https://arxiv.org/abs/cond-mat/0703788. Local PDF: `literature/ipeps_jordan_orus_vidal_0703788.pdf`.
- Wang, Verstraete, "Cluster update for tensor network states", arXiv: https://arxiv.org/abs/1110.4362. Local PDF: `literature/cluster_update_1110.4362.pdf`.
- Lubasch, Cirac, Banuls, "Algorithms for finite Projected Entangled Pair States", PRB 90, 064425, 2014. DOI: https://doi.org/10.1103/PhysRevB.90.064425, arXiv: https://arxiv.org/abs/1405.3259. Local PDF: `literature/finite_peps_algorithms_1405.3259.pdf`.
- Czarnik, Dziarmaga, "Projected entangled pair states at finite temperature: Iterative self-consistent bond renormalization for exact imaginary time evolution", PRB 92, 035120, 2015. DOI: https://doi.org/10.1103/PhysRevB.92.035120, arXiv: https://arxiv.org/abs/1411.6778. Local PDF: `literature/self_consistent_bond_renorm_1411.6778.pdf`.
- Czarnik, Dziarmaga, Corboz, "Time Evolution of an Infinite Projected Entangled Pair State: an Efficient Algorithm", PRB 99, 035115, 2019. DOI: https://doi.org/10.1103/PhysRevB.99.035115, arXiv: https://arxiv.org/abs/1811.05497. Local PDF: `literature/efficient_ipeps_evolution_1811.05497.pdf`.
- Hubig, Cirac, "Time-dependent study of disordered models with infinite projected entangled pair states", SciPost Phys. 6, 031, 2019. DOI: https://doi.org/10.21468/SciPostPhys.6.3.031, arXiv: https://arxiv.org/abs/1812.03801. Local PDF: `literature/realtime_disordered_ipeps_1812.03801.pdf`.
- Dziarmaga, "Time evolution of an infinite projected entangled pair state: a gradient tensor update in the tangent space", PRB 106, 014304, 2022. DOI: https://doi.org/10.1103/PhysRevB.106.014304, arXiv: https://arxiv.org/abs/2205.11067. Local PDF: `literature/gtu_tangent_2205.11067.pdf`.
- Alhambra, Cirac, "Locally Accurate Tensor Networks for Thermal States and Time Evolution", PRX Quantum 2, 040331, 2021. DOI: https://doi.org/10.1103/PRXQuantum.2.040331, arXiv: https://arxiv.org/abs/2106.00710. Local PDF: `literature/locally_accurate_tn_2106.00710.pdf`. This is mainly theoretical support for local time-evolution approximations, not a concrete backend.

## Triangular-Lattice Context

Xie et al., "Tensor renormalization of quantum many-body systems using projected entangled simplex states", PRX 4, 011025, 2014. DOI: https://doi.org/10.1103/PhysRevX.4.011025, arXiv: https://arxiv.org/abs/1307.5696. Local PDF: `literature/pess_simple_update_1307.5696.pdf`.

PESS is relevant because it explains why triangular and frustrated lattices can
benefit from simplex-aware entanglement structures. This repository should still
keep the current PEPS tensor convention with six virtual legs and one physical
leg unless a later project-local refactor explicitly chooses a different ansatz.

## Implementation Order

1. Keep dense local gates as the correctness oracle.
2. Maintain the current Simple Update path with explicit discarded-weight and
   norm diagnostics.
3. Strengthen the local star projection using SVD/HOSVD, ring Simple Update, or
   another PEPS-local factorization path.
4. Add NTU behind the same `apply_star_gate!`-style interface.
5. Add CTMRG only when an environment-aware observable or variational update
   needs it.
6. Treat fast Full Update, first-principles overlap optimization, and
   variational real-time PEPS methods as later accuracy upgrades, not as
   prerequisites for the first project-local PEPS dynamics workflow.
