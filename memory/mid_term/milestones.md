# Milestones

- Confirmed complete: Root Julia package scaffold with `ITensors`, `LinearAlgebra`, `Random`, and `Test`.
  - Source: `Project.toml`
  - Source: `git log`

- Confirmed complete: Triangular geometry, dense spin operators, PXP/blockade/star models, gates, schedules, and solvable benchmark helpers.
  - Source: `src/`
  - Source: `test/test_geometry.jl`
  - Source: `test/test_models.jl`
  - Source: `test/test_gates.jl`
  - Source: `test/test_schedules.jl`
  - Source: `test/test_solvable_models.jl`

- Confirmed complete: iPEPS state layer with one-site, three-site, and seven-site unit cells; product/random seeds; local observables; and hard state truncation.
  - Source: `src/States.jl`
  - Source: `src/Observables.jl`
  - Source: `test/test_states.jl`
  - Source: `test/test_observables.jl`

- Confirmed complete: Projected PXP step helpers, scheduled real/imaginary evolution, and per-step diagnostics.
  - Source: `src/Evolution.jl`
  - Source: `test/test_evolution.jl`

- Confirmed complete: Initial ScarFinder loop with deterministic seeds, repeated projected-PXP evolution, two-tier dimensions, blockade tolerance flagging, and deterministic ranking.
  - Source: `src/ScarFinder.jl`
  - Source: `test/test_scar_finder.jl`
  - Source: `git log`

- Inferred next milestone: Resolve current test/source drift around general `D>1` non-product star updates, then replace the local product-projection approximation with a true SVD/HOSVD, ring Simple Update, or NTU-compatible projection path.
  - Source: `Notes/current_peps_evolution_solution.md`
  - Source: `test/test_simple_update.jl`
