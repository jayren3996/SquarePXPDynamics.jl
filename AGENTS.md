# Agent Notes

## Project Direction

This repository is for 2D kagome-lattice PXP ScarFinder work. The PEPS/iPEPS code in this repo is internal tooling for that goal, not a standalone tensor-network package to be generalized independently.

Treat the PEPS layer as a subpackage-style module inside this project:

- Use the root `Project.toml` and root Julia environment.
- Do not create a nested Julia package, nested `Project.toml`, or separate PEPS environment.
- Keep PEPS code under the existing `src/` tree unless a later refactor introduces an internal namespace directory.
- Keep tests under the existing `test/` tree and run them with the root project:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

## Scientific Target

The eventual application is ScarFinder for the 2D kagome-lattice PXP model. PEPS features should be built only when they support that workflow:

- constrained PXP gates using `P_blockade * U`;
- real- and imaginary-time evolution needed by ScarFinder;
- fixed-bond-dimension evolve-project loops;
- blockade-violation diagnostics;
- low-entanglement candidate ranking and later ScarFinder orchestration.

Avoid broad PEPS library work that is not needed for kagome PXP ScarFinder.

## Current Implementation Boundary

The current PEPS implementation is an early internal tool. Prefer small, testable increments:

- preserve existing basis conventions: `|up> = |0>`, `|down> = |1>`;
- preserve dense 7-site star ordering: center first, then six neighbors in triangular direction order;
- keep projected PXP enforcement local and explicit;
- add correctness tests before performance work;
- use `ITensors.jl` through the root package dependency.

When adding new functionality, update README/docs so future agents understand whether it is ScarFinder-facing or only a temporary benchmark helper.

## Project Memory

Before substantial work, future agents should read `memory/README.md` and its recommended starter files. Run the project-memory-curator workflow only when explicitly requested. After substantial work, update short-term handoff memory only when asked; record important scientific, architectural, implementation, or workflow decisions in `memory/mid_term/decision_log.md`.
