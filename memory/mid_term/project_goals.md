# Project Goals

## Primary Goal

- Confirmed: `SquarePXPDynamics` is Julia tooling for square-lattice PXP
  dynamics with PEPS/iPEPS methods.
- Confirmed: The long-term project goal is a ScarFinder workflow that generates
  low-entanglement candidate states, evolves and projects them under
  constrained PXP dynamics, diagnoses leakage/truncation/entropy, and ranks
  scar-like trajectories.
- Source: `README.md`
- Source: `notes/README.md`

## Active Scope

- Confirmed: Active scope is square-lattice PXP dynamics, PEPS/iPEPS tooling
  that directly supports ScarFinder, dense local gates, local diagnostics, and
  fixed-bond-dimension evolve-project loops.
- Confirmed: CTMRG/full-update machinery is later accuracy infrastructure, not
  a prerequisite for the current simple-update ScarFinder prototype.
- Source: `README.md`
- Source: `notes/README.md`

## Current Branch Goal

- Confirmed: A completed feature branch,
  `codex/infinite-tfim-benchmark`, adds a v1 infinite TFIM benchmark framework
  reusing the five-site square-star iPEPS update machinery while preserving PXP
  behavior.
- Confirmed: That branch is intended to produce simple-update benchmark
  records, not CTMRG-quality TFIM physics claims.
- Source: `docs/superpowers/specs/2026-05-15-infinite-tfim-benchmark-design.md`
- Source: `docs/superpowers/notes/2026-05-15-current-work-infinite-tfim-benchmark.md`

## Working Rule

- Confirmed: New features should answer a ScarFinder-facing or benchmark-facing
  question: evolve, project, diagnose, rank, serialize, or validate. Otherwise,
  leave them out until a concrete need appears.
- Source: `notes/README.md`
- Source: `docs/superpowers/specs/2026-05-15-infinite-tfim-benchmark-design.md`
