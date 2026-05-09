# Physics Context

- Confirmed: The scientific target is ScarFinder for the 2D triangular-lattice PXP model. The PEPS/iPEPS code is internal tooling for that target, not a standalone tensor-network package.
  - Source: `AGENTS.md`
  - Source: `README.md`

- Confirmed: The triangular PXP local term is a 7-site star: a central spin flip dressed by projectors on the six nearest neighbors. A representative form is `h_c = P_1 P_2 P_3 P_4 P_5 P_6 X_c`.
  - Source: `docs/superpowers/specs/2026-05-08-triangular-peps-dynamics-design.md`
  - Source: `src/Models.jl`

- Confirmed: Constrained PXP evolution is enforced locally with projected gates `U_eff = P_blockade * U` for real time and `G_eff = P_blockade * G` for imaginary time. This removes forbidden local output support but does not eliminate possible constraint leakage introduced by approximate PEPS projection or truncation.
  - Source: `Notes/current_peps_evolution_solution.md`
  - Source: `src/Gates.jl`
  - Source: `src/Observables.jl`

- Confirmed: ScarFinder needs a repeated evolve-project map on a low-entanglement variational manifold, with ranking by truncation/discarded residuals, blockade violation, and entanglement proxies. Target-energy correction and full revival metrics are still future work in this repo.
  - Source: `Notes/implementation_roadmap.md`
  - Source: `src/ScarFinder.jl`

- Confirmed: The triangular-lattice cluster/stabilizer star model `K_i = X_i prod_{j in NN(i)} Z_j` is the preferred true-2D benchmark because the stabilizers commute and use the same 7-site star support as PXP.
  - Source: `docs/superpowers/specs/2026-05-08-triangular-peps-dynamics-design.md`
  - Source: `src/SolvableModels.jl`
