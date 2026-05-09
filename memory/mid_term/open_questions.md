# Open Questions

- Open question: `src/SimpleUpdate.jl`, `README.md`, and `Notes/current_peps_evolution_solution.md` describe a general `D>1` local product-projection path for dense non-product gates, and `test/test_evolution.jl` expects `D>1` projected PXP helpers to work. However, `test/test_simple_update.jl` still has a test named "general non-product star updates at D>1 fail explicitly" expecting an `ArgumentError`. Decide whether the test is stale or the implementation should reject that direct call.
  - Source: `src/SimpleUpdate.jl`
  - Source: `README.md`
  - Source: `Notes/current_peps_evolution_solution.md`
  - Source: `test/test_simple_update.jl`
  - Source: `test/test_evolution.jl`

- Open question: What exact local projection backend should replace the current `D>1` product-profile approximation: SVD/HOSVD, ring Simple Update, NTU, or another ScarFinder-specific projection?
  - Source: `Notes/implementation_roadmap.md`
  - Source: `Notes/literature_review.md`

- Open question: How should target-energy correction be implemented and validated with only local/star-energy estimates before environment contraction exists?
  - Source: `Notes/implementation_roadmap.md`
  - Source: `README.md`

- Open question: Which candidate-ranking metrics should become mandatory beyond discarded residual, blockade violation, and lambda entropy proxy: revival score, observable periodicity, energy drift, or three-sublattice density contrast?
  - Source: `Notes/implementation_roadmap.md`
  - Source: `src/ScarFinder.jl`

- Open question: When should CTMRG/full-environment observables be introduced, and which ScarFinder failure mode will justify that cost?
  - Source: `Notes/literature_review.md`
  - Source: `Notes/current_peps_evolution_solution.md`
