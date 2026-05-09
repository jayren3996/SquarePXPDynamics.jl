# Project Goals

- Confirmed: Build internal Julia tooling for triangular-lattice PEPS/iPEPS dynamics that supports the 2D triangular PXP ScarFinder workflow.
  - Source: `AGENTS.md`
  - Source: `README.md`

- Confirmed: Near-term ScarFinder-facing priorities are constrained PXP gates, real- and imaginary-time evolution, fixed-bond-dimension evolve-project loops, blockade diagnostics, low-entanglement candidate ranking, and later ScarFinder orchestration.
  - Source: `AGENTS.md`
  - Source: `Notes/README.md`

- Confirmed: Avoid broad PEPS library work not needed for triangular PXP ScarFinder, including arbitrary graph PEPS, nested Julia packages, early CTMRG/full update, GPU/symmetry backends, or broad Rydberg Hamiltonian packaging.
  - Source: `AGENTS.md`
  - Source: `Notes/implementation_roadmap.md`

- Confirmed: The first accuracy ladder is dense-star correctness, Simple Update diagnostics, stronger local star projection such as SVD/HOSVD or ring Simple Update, then NTU, then CTMRG/full-environment methods if needed.
  - Source: `Notes/literature_review.md`
  - Source: `Notes/implementation_roadmap.md`
