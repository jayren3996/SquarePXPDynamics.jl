# Active Tasks

- Resolve whether `test/test_simple_update.jl` should be updated to accept the current `D>1` local product-projection path or whether direct general dense-star updates should throw.
  - Source: `memory/mid_term/open_questions.md`

- Run the full root test suite after resolving test/source drift:
  - `julia --project=. -e 'using Pkg; Pkg.test()'`
  - Source: `AGENTS.md`

- Keep README/docs aligned with whether `D>1` behavior is a supported prototype path or an intentionally rejected direct API.
  - Source: `README.md`
  - Source: `Notes/current_peps_evolution_solution.md`
