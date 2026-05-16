# S0-S7 Slice 3 ScarFinder Energy Correction Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Complete the original S6 guarded energy-correction requirement while preserving S6-lite defaults.

**Architecture:** Add opt-in simple/local PXP energy correction to `ScarFinderParams` and `scarfinder!`. Correction uses existing imaginary-time `evolve!` attempts on copied states, accepts only strict diagnostic improvement toward `target_energy`, records accepted/rejected outcomes, and leaves CTM-trusted energy ranking deferred to S7b.

**Tech Stack:** Julia 1.12, existing `SquarePXPDynamics` iPEPS/evolution/observable modules, `Test`.

---

## Tasks

1. Add failing ScarFinder tests for correction parameters, default skip behavior,
   rejected non-improving correction, and non-worsening correction records.
2. Implement `target_energy`, `correction_time`, and `correction_attempts` in
   `ScarFinderParams` with validation and backward-compatible constructors.
3. Add correction fields to `ScarFinderIteration` and `ScarFinderCandidateScore`.
4. Implement correction attempts with `copy_state`, imaginary-time
   `TrotterParams`, and transactional state replacement only on strict
   improvement.
5. Extend CSV/JSON logs with correction fields appended after existing CTM
   fields so existing column order remains stable.
6. Update README and milestone memory.
7. Run focused ScarFinder/public-doc tests, then full package tests.

