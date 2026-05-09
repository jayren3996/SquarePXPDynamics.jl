# Experiments And Results

- Confirmed: Current tests cover geometry/coloring, spin operators, model/projector construction, dense/projected gates, schedules, solvable models, state containers, observables, Simple Update paths, evolution, and ScarFinder.
  - Source: `test/runtests.jl`

- Confirmed: Dense projected PXP gate tests verify dense blockade projection behavior on allowed product inputs.
  - Source: `test/test_gates.jl`
  - Source: `test/test_simple_update.jl`

- Confirmed: Evolution tests check canonical color centers, identity preservation, Hamiltonian schedule weights, projected PXP smoke paths, imaginary projected PXP steps, `D>1` projected step helpers, and prebuilt-gate half-step behavior.
  - Source: `test/test_evolution.jl`

- Confirmed: ScarFinder tests check config validation, deterministic search output, `D=2` candidates over one-site/three-site/seven-site unit cells, hard truncation from dynamics to scar dimension, blockade tolerance flagging, and stable ranking.
  - Source: `test/test_scar_finder.jl`

- Not verified during memory creation: Full `Pkg.test()` result. Run it after resolving the open test/source drift noted in `open_questions.md`.
  - Source: current memory-curation session
