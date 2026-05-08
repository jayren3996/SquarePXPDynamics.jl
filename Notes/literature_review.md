# Literature Review: Triangular PXP ScarFinder With PEPS

Prepared for the repository at `/Users/ren/Codex/PEPs`, based on the local design report and online primary literature reviewed on 2026-05-08.

## Scope

The target is not a general PEPS library. The algorithm should support ScarFinder for the 2D triangular-lattice PXP model:

- constrained 7-site star gates;
- real- and imaginary-time evolution;
- fixed-bond-dimension PEPS projection;
- blockade-violation diagnostics;
- candidate ranking for low-entanglement scar-like trajectories.

## 1. ScarFinder And Scar Diagnostics

### Core Reference

Ren, Hallam, Ying, Papic, "ScarFinder: a detector of optimal scar trajectories in quantum many-body dynamics", PRX Quantum 6, 040332, 2025. DOI: https://doi.org/10.1103/8g8w-nkwx, arXiv: https://arxiv.org/abs/2504.12383.

Local PDF: `literature/scarfinder_2504.12383.pdf`.

### Findings

ScarFinder does not first solve for scar eigenstates. It iterates real-time evolution and projection back to a low-entanglement variational manifold. The paper's operational loop is:

```text
initialize |psi> in manifold M
repeat:
  evolve |psi(dt)> = exp(-i H dt) |psi>
  project |psi'> = P_M |psi(dt)>
  correct energy / symmetries / constraints
  set |psi> = corrected |psi'>
return candidate |psi>
```

The key heuristic is that generic thermal components leave the low-entanglement manifold faster than scar components. Projection suppresses those thermal components, but it can also bias the state toward trivial low-entanglement states, so energy targeting and constraint enforcement are required.

ScarFinder distinguishes the projection interval `DeltaT` from the smaller microscopic time step used inside the time-evolution integrator. Larger `DeltaT` can improve thermal/scar separation but is limited by entanglement growth and truncation stability.

For PXP, the paper reports that truncation can generate forbidden adjacent excitations. Their 1D MPS implementation used an imaginary penalty term as a way to suppress constraint leakage. In this repository the more natural first-line strategy is already in the design report: local projected gates `P_blockade * U`, plus diagnostics that reject candidates with residual leakage.

### Related Scar References

- Ho, Choi, Pichler, Lukin, "Periodic Orbits, Entanglement, and Quantum Many-Body Scars in Constrained Models", PRL 122, 040603, 2019. DOI: https://doi.org/10.1103/PhysRevLett.122.040603, arXiv: https://arxiv.org/abs/1807.01815. Local PDF: `literature/periodic_orbits_scars_1807.01815.pdf`.
- Turner et al., "Quantum scarred eigenstates in a Rydberg atom chain", PRB 98, 155134, 2018. DOI: https://doi.org/10.1103/PhysRevB.98.155134, arXiv: https://arxiv.org/abs/1806.10933. Local PDF: `literature/weak_ergodicity_scars_1806.10933.pdf`.
- Choi et al., "Emergent SU(2) dynamics and perfect quantum many-body scars", PRL 122, 220603, 2019. DOI: https://doi.org/10.1103/PhysRevLett.122.220603, arXiv: https://arxiv.org/abs/1812.05561. Local PDF: `literature/emergent_su2_scars_1812.05561.pdf`.
- Khemani, Laumann, Chandran, "Signatures of integrability in the dynamics of Rydberg-blockaded chains", PRB 99, 161101(R), 2019. DOI: https://doi.org/10.1103/PhysRevB.99.161101, arXiv: https://arxiv.org/abs/1807.02108. Local PDF: `literature/rydberg_integrability_1807.02108.pdf`.

### Implementation Implications

ScarFinder-facing code should return a structured candidate record, not just a final PEPS:

- seed and unit cell;
- bond dimension `D`;
- projection interval `DeltaT`;
- microscopic step `dt`;
- energy density and target-energy error;
- blockade violation;
- discarded weight / truncation residual;
- revival score;
- local-observable oscillation score;
- entanglement proxy;
- reason for rejection or convergence.

The first ranking rule should be conservative:

```text
reject if blockade violation > tolerance
reject if energy correction diverges or drifts monotonically away from target
rank remaining candidates by low late-time entanglement proxy and revival/local-observable scores
```

## 2. PXP And Rydberg Literature

### Core 2D Reference

Lin, Calvera, Hsieh, "Quantum Many-Body Scar States in Two-Dimensional Rydberg Atom Arrays", PRB 101, 220304(R), 2020. DOI: https://doi.org/10.1103/PhysRevB.101.220304, arXiv: https://arxiv.org/abs/2003.04516.

Local PDF: `literature/2d_rydberg_scars_2003.04516.pdf`.

### Findings

The effective PXP model is obtained in the nearest-neighbor blockade regime by projecting the Rydberg Hamiltonian to configurations with no adjacent excitations and keeping locally constrained spin flips:

```text
H = sum_i X_i prod_{j in NN(i)} P_j
```

For the triangular lattice, a bulk site has six nearest neighbors, so the elementary local update for a central flip is naturally a 7-site star: center plus six triangular directions.

The 2D paper demonstrates that higher-dimensional PXP systems can host exact low-entanglement scar states and charge-density-wave quench oscillations. This justifies looking for low-entanglement trajectories in 2D, but it does not justify importing the 1D Z2/Z3 story as a fixed assumption. The triangular lattice is non-bipartite and frustrated, so the implementation should keep unit cells flexible and rank by measured dynamics.

### Triangular Star Blockade Detail

There are two related but distinct projectors:

1. The central PXP term projector, which checks whether all six neighbors of the center are unexcited before flipping the center.
2. The full 7-site local output blockade projector, which rejects any adjacent excitations inside the star.

For a triangular star, the full local output projector should check 12 local edges:

- 6 center-neighbor spokes;
- 6 neighbor-neighbor edges around the hexagonal ring.

This matters for `P_blockade * U`: after applying a dense 7-site gate, the local projection should remove forbidden output weight on all star edges, not only center-neighbor spokes. Tests should also compare `P_blockade * U` and `P_blockade * U * P_blockade` on constrained and unconstrained input vectors.

### Related Rydberg References

- Bernien et al., "Probing many-body dynamics on a 51-atom quantum simulator", Nature 551, 579, 2017. DOI: https://doi.org/10.1038/nature24622, arXiv: https://arxiv.org/abs/1707.04344. Local PDF: `literature/rydberg_51_atom_1707.04344.pdf`.
- Bluvstein et al., "Controlling quantum many-body dynamics in driven Rydberg atom arrays", Science 371, 1355, 2021. DOI: https://doi.org/10.1126/science.abg2530, arXiv: https://arxiv.org/abs/2012.12276. Local PDF: `literature/driven_rydberg_scars_2012.12276.pdf`.
- Ebadi et al., "Quantum phases of matter on a 256-atom programmable quantum simulator", Nature 595, 227, 2021. DOI: https://doi.org/10.1038/s41586-021-03582-4, arXiv: https://arxiv.org/abs/2012.12281. Local PDF: `literature/rydberg_256_atoms_2012.12281.pdf`.
- Scholl et al., "Quantum simulation of 2D antiferromagnets with hundreds of Rydberg atoms", Nature 595, 233, 2021. DOI: https://doi.org/10.1038/s41586-021-03585-1, arXiv: https://arxiv.org/abs/2012.12268. Local PDF: `literature/rydberg_2d_antiferromagnets_2012.12268.pdf`.
- Guo, Hu, Li, "Order-by-disorder and emergent Kosterlitz-Thouless phase in triangular Rydberg array", arXiv: https://arxiv.org/abs/2302.08963. Local PDF: `literature/triangular_rydberg_2302.08963.pdf`.
- Patil, "Quantum Monte Carlo simulations in the restricted Hilbert space of Rydberg atom arrays", SciPost Phys. 20, 022, 2026. DOI: https://doi.org/10.21468/SciPostPhys.20.1.022, arXiv: https://arxiv.org/abs/2309.00482. Local PDF: `literature/restricted_hilbert_rydberg_qmc_2309.00482.pdf`.

### Implementation Implications

Add diagnostics before broad model support:

- nearest-neighbor blockade violation density over all triangular edges;
- local 7-site invalid weight using the same edge list as the projector;
- center-flip eligibility probability;
- post-gate leakage before and after projection;
- 3-sublattice density imbalance for triangular density-wave candidates;
- density/revival contrast over a projection trajectory.

Keep one explicit convention for "excited/Rydberg" vs the repo basis. The project currently preserves `|up> = |0>`, `|down> = |1>`. The projector code should make the excitation convention visible so future agents cannot silently invert blockade logic.

## 3. PEPS/iPEPS Evolution Literature

### Core References

- Verstraete, Cirac, "Renormalization algorithms for Quantum-Many Body Systems in two and higher dimensions", arXiv: https://arxiv.org/abs/cond-mat/0407066. Local PDF: `literature/peps_origin_0407066.pdf`.
- Jordan, Orus, Vidal, Verstraete, Cirac, "Classical Simulation of Infinite-Size Quantum Lattice Systems in Two Spatial Dimensions", PRL 101, 250602, 2008. DOI: https://doi.org/10.1103/PhysRevLett.101.250602, arXiv: https://arxiv.org/abs/cond-mat/0703788. Local PDF: `literature/ipeps_jordan_orus_vidal_0703788.pdf`.
- Jiang, Weng, Xiang, "Accurate Determination of Tensor Network State of Quantum Lattice Models in Two Dimensions", PRL 101, 090603, 2008. DOI: https://doi.org/10.1103/PhysRevLett.101.090603, arXiv: https://arxiv.org/abs/0806.3719. Local PDF: `literature/simple_update_jiang_0806.3719.pdf`.
- Dziarmaga, "Time evolution of an infinite projected entangled pair state: a neighborhood tensor update", PRB 104, 094411, 2021. DOI: https://doi.org/10.1103/PhysRevB.104.094411, arXiv: https://arxiv.org/abs/2107.06635. Local PDF: `literature/ntu_ipeps_2107.06635.pdf`.

### Findings

All PEPS evolution algorithms face the same bottleneck: a local gate increases bond dimension, then the network must be projected back to fixed `D`.

Simple Update is the right first implementation:

- cheap and stable;
- stores diagonal bond spectra on virtual bonds;
- gives useful local truncation diagnostics;
- underestimates loop and long-range environment effects.

Full Update / fast Full Update improves accuracy by using an environment, but introduces CTMRG, gauge conditioning, regularization, and ALS stability issues. That is too much machinery for the first ScarFinder prototype.

Neighborhood Tensor Update is the best second implementation target. NTU optimizes with an exactly contractible local neighborhood metric rather than the full infinite environment. It is more accurate than Simple Update while remaining local enough to share infrastructure with star update diagnostics.

### Related PEPS Update References

- Wang, Verstraete, "Cluster update for tensor network states", arXiv: https://arxiv.org/abs/1110.4362. Local PDF: `literature/cluster_update_1110.4362.pdf`.
- Lubasch, Cirac, Banuls, "Algorithms for finite Projected Entangled Pair States", PRB 90, 064425, 2014. DOI: https://doi.org/10.1103/PhysRevB.90.064425, arXiv: https://arxiv.org/abs/1405.3259. Local PDF: `literature/finite_peps_algorithms_1405.3259.pdf`.
- Phien, Bengua, Tuan, Corboz, Orus, "The iPEPS algorithm, improved: fast full update and gauge fixing", PRB 92, 035142, 2015. DOI: https://doi.org/10.1103/PhysRevB.92.035142, arXiv: https://arxiv.org/abs/1503.05345. Local PDF: `literature/fast_full_update_1503.05345.pdf`.
- Czarnik, Dziarmaga, "Projected entangled pair states at finite temperature: Iterative self-consistent bond renormalization for exact imaginary time evolution", PRB 92, 035120, 2015. DOI: https://doi.org/10.1103/PhysRevB.92.035120, arXiv: https://arxiv.org/abs/1411.6778. Local PDF: `literature/self_consistent_bond_renorm_1411.6778.pdf`.
- Czarnik, Dziarmaga, "Time Evolution of an Infinite Projected Entangled Pair State: an Algorithm from First Principles", PRB 98, 045110, 2018. DOI: https://doi.org/10.1103/PhysRevB.98.045110, arXiv: https://arxiv.org/abs/1804.03872. Local PDF: `literature/first_principles_ipeps_1804.03872.pdf`.
- Czarnik, Dziarmaga, Corboz, "Time Evolution of an Infinite Projected Entangled Pair State: an Efficient Algorithm", PRB 99, 035115, 2019. DOI: https://doi.org/10.1103/PhysRevB.99.035115, arXiv: https://arxiv.org/abs/1811.05497. Local PDF: `literature/efficient_ipeps_evolution_1811.05497.pdf`.
- Hubig, Cirac, "Time-dependent study of disordered models with infinite projected entangled pair states", SciPost Phys. 6, 031, 2019. DOI: https://doi.org/10.21468/SciPostPhys.6.3.031, arXiv: https://arxiv.org/abs/1812.03801. Local PDF: `literature/realtime_disordered_ipeps_1812.03801.pdf`.
- Dziarmaga, "Time evolution of an infinite projected entangled pair state: a gradient tensor update in the tangent space", PRB 106, 014304, 2022. DOI: https://doi.org/10.1103/PhysRevB.106.014304, arXiv: https://arxiv.org/abs/2205.11067. Local PDF: `literature/gtu_tangent_2205.11067.pdf`.
- Alhambra, Cirac, "Locally Accurate Tensor Networks for Thermal States and Time Evolution", PRX Quantum 2, 040331, 2021. DOI: https://doi.org/10.1103/PRXQuantum.2.040331, arXiv: https://arxiv.org/abs/2106.00710. Local PDF: `literature/locally_accurate_tn_2106.00710.pdf`.

### Triangular-Lattice Tensor Network Note

Xie et al., "Tensor renormalization of quantum many-body systems using projected entangled simplex states", PRX 4, 011025, 2014. DOI: https://doi.org/10.1103/PhysRevX.4.011025, arXiv: https://arxiv.org/abs/1307.5696.

Local PDF: `literature/pess_simple_update_1307.5696.pdf`.

PESS is relevant because it explains why triangular/frustrated lattices often benefit from simplex-aware entanglement structures. This repo should still keep the current PEPS tensor convention with six virtual legs and one physical leg, because the project report explicitly targets triangular iPEPS tensors rather than introducing a new ansatz. PESS is context, not an immediate refactor target.

## Final Synthesis

The implementation should proceed in this order:

1. Correct dense 7-site triangular PXP gates and full-star blockade projectors.
2. Exact `D=1` star update oracle for product iPEPS.
3. General Simple Update star projection with per-spoke spectra and discarded-weight diagnostics.
4. PEPS-level local blockade and energy diagnostics.
5. ScarFinder driver using real-time projected evolution, fixed-D projection, optional imaginary-time target-energy correction, and candidate ranking.
6. NTU backend for the same `apply_gate!` interface.
7. Only later: full/fast-full update, CTMRG, and compressed 7-site gate/operator formats.
