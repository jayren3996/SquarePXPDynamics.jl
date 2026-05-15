# Physics Context

## Square-Lattice PXP Model

- Confirmed: The core model is square-lattice PXP dynamics, with one spin-1/2
  physical degree of freedom per square-lattice site.
- Confirmed: The local five-site square-star PXP term is
  `h_c = X_c P_down(right) P_down(up) P_down(left) P_down(down)`.
- Confirmed: The dense five-site Hamiltonian and blockade projector remain the
  source of truth for PXP local physics.
- Source: `README.md`
- Source: `notes/README.md`
- Source: `notes/2026-05-15-ipeps-literature-code-algorithm-notes.md`
- Source: `src/SquarePXP.jl`

## Blockade Constraint

- Confirmed: The local projected PXP gates remove forbidden local output support
  with an up center adjacent to an up neighbor.
- Confirmed: Approximate PEPS/iPEPS projection can still leak outside the
  constrained manifold, so blockade diagnostics are mandatory.
- Source: `README.md`
- Source: `notes/README.md`
- Source: `src/Observables.jl`

## iPEPS And Simple Update

- Confirmed: The implemented iPEPS state is a periodic Gamma-lambda
  simple-update style state, not a canonical iPEPS gauge.
- Confirmed: The five-site square-star update uses QR reduction of the four
  leaf tensors before applying the dense gate to the reduced core.
- Confirmed: Simple/local observables are development and regression
  diagnostics, not final CTMRG-quality physics measurements.
- Source: `README.md`
- Source: `notes/2026-05-15-chatgpt-pro-ipeps-review-plan.md`
- Source: `notes/2026-05-15-ipeps-literature-code-algorithm-notes.md`
- Source: `src/SquareIPEPS.jl`
- Source: `src/StarSimpleUpdate.jl`
- Source: `src/Observables.jl`

## ScarFinder Context

- Confirmed: ScarFinder is treated as an evolve-project workflow: evolve with
  local dynamics, project/truncate back into a low-entanglement variational
  manifold, measure diagnostics, then rank candidate trajectories.
- Confirmed: In this repository, ScarFinder should orchestrate evolution,
  measurement, ranking, and logging; it should not own low-level tensor index
  logic or CTMRG internals.
- Source: `notes/README.md`
- Source: `notes/2026-05-15-ipeps-literature-code-algorithm-notes.md`
- Source: `src/ScarFinder.jl`

## TFIM Benchmark Context

- Confirmed: A separate TFIM benchmark design uses the same five-site
  square-star machinery with local convention
  `h_c = -h X_c - (J/2) Z_c (Z_right + Z_up + Z_left + Z_down)`.
- Confirmed: The `J/2` factor is intended to avoid double-counting nearest
  neighbor bonds when summing star terms over every lattice site.
- Confirmed: The TFIM implementation was completed on branch
  `codex/infinite-tfim-benchmark` and then merged into `main` locally on
  2026-05-15.
- Source: `docs/superpowers/specs/2026-05-15-infinite-tfim-benchmark-design.md`
- Source: `docs/superpowers/notes/2026-05-15-current-work-infinite-tfim-benchmark.md`
