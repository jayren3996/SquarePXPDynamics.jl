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
  that directly supports ScarFinder, dense local gates, local and CTM-backed
  diagnostics, fixed-bond-dimension evolve-project loops, and benchmark/reference
  validation.
- Confirmed: CTMRG trust and S7b gauge-conditioning infrastructure are now
  available as accuracy/readiness layers. Production ScarFinder validation and
  physics-facing CTM workflows remain future work.
- Source: `README.md`
- Source: `notes/README.md`
- Source: `docs/superpowers/notes/2026-05-16-s0-s7-completion-audit.md`

## Current Integrated State

- Confirmed: Local `main` now includes the S0-S7 completion branch, including
  the v1 infinite TFIM benchmark framework, finite TFIM/MPS/PXP ED reference
  paths, S7a CTM trust, and S7b CTM gauge-readiness/conditioning APIs.
- Confirmed: Simple-update TFIM and ScarFinder outputs remain implementation
  diagnostics unless backed by appropriate CTM trust and finite-chi workflows.
- Source: `README.md`
- Source: `src/SquarePXPDynamics.jl`
- Source: `docs/superpowers/notes/2026-05-16-s0-s7-completion-audit.md`

## Working Rule

- Confirmed: New features should answer a ScarFinder-facing or benchmark-facing
  question: evolve, project, diagnose, rank, serialize, or validate. Otherwise,
  leave them out until a concrete need appears.
- Source: `notes/README.md`
- Source: `docs/superpowers/specs/2026-05-15-infinite-tfim-benchmark-design.md`
