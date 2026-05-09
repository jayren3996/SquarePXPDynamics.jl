# Theoretical Constraints

- Confirmed: Same-color triangular star centers must have disjoint radius-1 stars. First-order evolution sweeps colors `1:7`; second-order evolution sweeps `1:7` then `7:1` with half-step layers.
  - Source: `Notes/implementation_roadmap.md`
  - Source: `src/Schedules.jl`
  - Source: `test/test_schedules.jl`

- Confirmed: Applying a local PEPS gate increases effective bond dimension, so fixed-`D` dynamics requires per-gate projection/truncation. The `dynamics_maxdim` / `scar_maxdim` split does not postpone all projection to the end of a projection interval.
  - Source: `Notes/current_peps_evolution_solution.md`
  - Source: `Notes/implemented_peps_algorithm_detail.md`

- Confirmed: Local one-site observables and nearest-neighbor blockade diagnostics are exact for `D=1` product states. For `D>1`, they are local screening diagnostics rather than contracted PEPS expectations.
  - Source: `Notes/current_peps_evolution_solution.md`
  - Source: `src/Observables.jl`

- Confirmed: Dense projected PXP gates satisfy the local relation that `P_blockade * U` removes forbidden output support, and on constrained inputs agrees with `P_blockade * U * P_blockade`.
  - Source: `Notes/implementation_roadmap.md`
  - Source: `test/test_gates.jl`

- Open question: Which diagnostics will be sufficient for production ScarFinder ranking before CTMRG or NTU is implemented remains unresolved.
  - Source: `Notes/implementation_roadmap.md`
  - Source: `Notes/current_peps_evolution_solution.md`
