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
- Confirmed: CTMRG trust, S7b gauge-conditioning infrastructure, trusted
  ScarFinder measurement backends, objective-based ranking, and PXP validation
  reports are now available as accuracy/readiness layers. Publication-grade
  ScarFinder campaigns and CTM-aware/full-update evolution remain future work.
- Source: `README.md`
- Source: `notes/README.md`
- Source: `docs/superpowers/notes/2026-05-16-s0-s7-completion-audit.md`
- Source: `docs/superpowers/notes/2026-05-17-gpt-roadmap-completion.md`

## Current Integrated State

- Confirmed: Local `main` now includes the S0-S7 completion work and the GPT
  PXP roadmap completion, including the v1 infinite TFIM benchmark framework,
  finite TFIM/MPS/PXP ED reference paths, S7a CTM trust, S7b CTM
  gauge-readiness/conditioning APIs, trusted ScarFinder backends, objective
  scoring, candidate metadata persistence, reproducible CTMRG initialization,
  PXP validation/convergence reports, and reverse-evolution diagnostics.
- Confirmed: Simple-update TFIM and ScarFinder outputs remain implementation
  diagnostics unless backed by appropriate CTM trust, finite-chi workflows, and
  convergence/error-budget reports.
- Source: `README.md`
- Source: `src/SquarePXPDynamics.jl`
- Source: `docs/superpowers/notes/2026-05-16-s0-s7-completion-audit.md`
- Source: `docs/superpowers/notes/2026-05-17-gpt-roadmap-completion.md`

## Working Rule

- Confirmed: New features should answer a ScarFinder-facing or benchmark-facing
  question: evolve, project, diagnose, rank, serialize, or validate. Otherwise,
  leave them out until a concrete need appears.
- Source: `notes/README.md`
- Source: `docs/superpowers/specs/2026-05-15-infinite-tfim-benchmark-design.md`
