# S0-S7 Slice 5 D2 Gauge Conditioning Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Complete S7b by adding transactional D>1 `fix_bond_gauge!` mutation using validated CTM bond environments.

**Architecture:** Reuse the Slice 4 readiness gate. For D>1 bonds, use PEPSKit's `bondenv_fu`, `positive_approx`, `fixgauge_benv`, and `_fixgauge_benvXY` on the already converted PEPSKit tensors, then convert the two updated absorbed PEPSKit tensors back into the custom ITensors Gamma-lambda representation by deabsorbing the stored link weights. All factorization and conversion happens before mutating `psi`.

**Tech Stack:** Julia 1.12, ITensors, TensorKit, PEPSKit, existing `CTMGaugeReadiness`, `PEPSKitMeasurements`, and `Test`.

---

## Tasks

1. [x] Replace the Slice 4 D>1 explicit-nonmutation test with a RED test that expects `fix_bond_gauge!` to mutate a D=2 seeded positive-link fixture, increment `state_version`, keep simple observables finite, and stale the old CTM context.
2. [x] Add an extended CTM regression behind `SQUAREPXP_EXTENDED_TESTS=1` comparing `measure_ctm` before and after the D=2 gauge mutation with fresh contexts.
3. [x] Implement TensorKit-to-ITensors Gamma conversion for PEPSKit tensors with symmetric lambda deabsorption and zero-division protection.
4. [x] Implement horizontal and vertical D>1 gauge conditioning with PEPSKit bond-environment factorization.
5. [x] Update docs/memory to mark S7b/S7 complete once full verification passes.
6. [x] Run focused Slice 5 tests, extended CTM tests, full package tests, and then the final S0-S7 completion audit.
