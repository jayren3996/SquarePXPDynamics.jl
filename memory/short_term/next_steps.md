# Next Steps

1. Inspect and resolve the `D>1` Simple Update contradiction between `src/SimpleUpdate.jl`, `test/test_evolution.jl`, and `test/test_simple_update.jl`.
2. Run the full root Julia tests with `julia --project=. -e 'using Pkg; Pkg.test()'`.
3. If tests pass, update short-term handoff with the result.
4. Next implementation work should focus on a real local refactorization backend for dense non-product 7-site gates, unless ScarFinder diagnostics show a more urgent blocker.
