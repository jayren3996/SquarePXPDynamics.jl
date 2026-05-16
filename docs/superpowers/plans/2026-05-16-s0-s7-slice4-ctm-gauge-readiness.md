# S0-S7 Slice 4 CTM Gauge Readiness Implementation Plan

> **Status:** Completed and merged locally. The unchecked boxes below are the
> original execution template, not current TODOs. See
> `docs/superpowers/notes/2026-05-16-s0-s7-completion-audit.md`.

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add the S7b readiness layer: CTM-backed bond norm diagnostics, a readiness predicate, and a transactional product/no-op `fix_bond_gauge!` path.

**Architecture:** Keep S7b gauge-readiness code in a new `CTMGaugeReadiness.jl` module so S7a simple-gauge diagnostics remain PEPSKit-free. Reuse the existing PEPSKit CTMRG context, PEPSKit `bondenv_fu` bond environment contraction, and CTM trust assessment. Do not mutate D>1 tensors in this slice; instead return a structured `:d_greater_than_one_not_implemented` result after all readiness checks pass.

**Tech Stack:** Julia 1.12, ITensors, TensorKit, PEPSKit, existing `SquareIPEPS`, `PEPSKitMeasurements`, `CTMTrust`, and `Test`.

---

## Tasks

1. Add failing `test/test_ctm_gauge_readiness.jl` coverage for:
   - exported public API names,
   - malformed norm-matrix diagnostics,
   - stale CTM context rejection,
   - D=1 product CTM bond-norm diagnostics for canonical `:right` and `:up`,
   - readiness rejection for untrusted CTM and malformed norm diagnostics,
   - product/no-op `fix_bond_gauge!` preserving simple observables and state version,
   - D>1 readiness returning an explicit not-implemented mutation result rather than silently mutating.
2. Add a public freshness helper in `PEPSKitMeasurements.jl`:
   `assert_fresh_pepskit_context(psi, ctx)`.
3. Implement `src/CTMGaugeReadiness.jl` with:
   - `CTMGaugePolicy`,
   - `CTMBondNormDiagnostic`,
   - `CTMGaugeReadiness`,
   - `BondGaugeFixInfo`,
   - `ctm_bond_norm_matrix`,
   - `ctm_bond_norm_diagnostic`,
   - `all_ctm_bond_norm_diagnostics`,
   - `ctm_ready_for_gauge_updates`,
   - `fix_bond_gauge!`.
4. Include/export the new module through `src/SquarePXPDynamics.jl`.
5. Update README, milestones, and decision log to record that Slice 4 adds
   readiness and a product/no-op gauge path while D>1 mutation remains Slice 5.
6. Run:
   `julia --project=. test/runtests.jl test_ctm_gauge_readiness.jl test_ctm_trust.jl test_pepskit_measurements.jl test_public_docs.jl`
7. Commit as `feat: add CTM gauge readiness diagnostics`.
