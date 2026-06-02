# Aggressive-deletion execution plan (2026-06-02)

Grep-verified, ordered, one-slice-per-pass. Run the suite after each slice.
Net: ~2,100 src + ~1,400 test lines + 4 scripts + 1 dep. Two targets BLOCKED.

Execution order: 1 → 2 → 3 → 4 → 5. After Slice 5 specifically verify
`test_aqua`, `test_star_models`, `test_observables`, `test_star_simple_update`,
`test_public_docs`.

## Slice 1 — CTMGaugeReadiness (orphaned, −662 src/−290 test)
Only caller is `test/test_ctm_gauge_readiness.jl`; the `src/SquareIPEPS.jl:277`
hit is a comment. `git rm src/CTMGaugeReadiness.jl test/test_ctm_gauge_readiness.jl`.
In `src/SquarePXPDynamics.jl`: remove `include("CTMGaugeReadiness.jl")` (L26),
the `using .CTMGaugeReadinessModule:` block (L134-144), exports L309-312.
In `test/runtests.jl`: remove `"test_ctm_gauge_readiness.jl",` (L34).
PEPSKit/TensorKit stay (used by PEPSKitMeasurements).

## Slice 2 — GaugeDiagnostics (unused, −147 src/−150 test)
Only caller is `test/test_gauge_diagnostics.jl`. `git rm src/GaugeDiagnostics.jl
test/test_gauge_diagnostics.jl`. In main module: remove `include` L20, `using
.GaugeDiagnostics:` L98-100, exports L279-280. runtests: remove entry L24.

## Slice 3 — SquarePEPS (test-only, −127 src)
Repair FIRST: in `src/SquareIPEPS.jl` replace L7
`import ..SquarePEPS: physical_index, link_index` with
`function physical_index end` / `function link_index end`.
`git rm src/SquarePEPS.jl test/test_square_peps.jl`. In main module: remove
`include` L17, `using .SquarePEPS:` L79-80, export L267 — but ADD
`export physical_index, link_index` (L267 was their only export; they have
docstrings in SquareIPEPS so test_public_docs stays green). runtests: remove L6.

## Slice 4 — duplicate PXP audit campaign (−~210 src/−150 test/−40 script)
Keep `run_pxp_larger_d_benchmark`, `validate_pxp_ed_ipeps`,
`validate_pxp_reversibility`. KEEP shared helpers in PXPValidation.jl
(`_ctm_param_tuple`, `_maximum_or_zero`, `_minimum_or_zero`, `_finite_chi_max`,
`_audit_trust_status`) and PXPValidationSerialization.jl `_audit_csv_cell`.
`git rm scripts/pxp_audit_campaign.jl`. In `src/PXPValidation.jl` delete
high→low: `run_pxp_audit_campaign` 1311-1349, `_audit_summary` 1244-1309,
`_audit_validation_config` 1044-1063, `_audit_ctm_params` 1036-1042,
`PXPAuditReport` 500-512, `PXPAuditRun` 488-498, `PXPAuditSummary` 449-486,
`PXPAuditConfig` 350-447, exports 40-41. In `src/PXPValidationSerialization.jl`:
remove StructTypes loop L24-27, `write_pxp_audit_json` 143-151,
`PXP_AUDIT_CSV_HEADER` 163-192, `_audit_csv_row` 214-246, `write_pxp_audit_csv`
248-263 (KEEP `_audit_csv_cell` 194-212). Main module: remove
`using .PXPValidation:` L204-210, exports 337-338. In
`test/test_pxp_validation.jl`: delete audit testsets L463-EOF (last non-audit
`end` at 462; helpers at top stay). Keep the file in runtests.

## Slice 5 — TFIM axis (−948 src/−500 test/−381 script/−1 dep)
KEEP `test_tfim_schedule_reference.jl` (no TFIM symbols — Trotter schedule only)
and `StaticModel`/`AbstractModelProtocol`/`model_at` (PXP uses them).
`git rm` src/Benchmarks.jl, src/FiniteTFIMReference.jl,
src/FiniteMPSTFIMReference.jl, test/test_benchmarks.jl,
test/test_tfim_observables.jl, test/test_tfim_finite_reference.jl,
test/test_finite_mps_tfim_reference.jl, scripts/compare_tfim_tdvp_ipeps_t03.jl,
scripts/realistic_tfim_dynamics.jl, scripts/finite_mps_tfim_6x6.jl.
- `src/StarModels.jl`: delete `tfim_product_basis_energy` 216-230, `_tfim_z_value`
  205-214, `star_gate(::TFIMStarModel)` 134-139, `star_hamiltonian(::TFIMStarModel)`
  112-120, `_validate_finite_step`/`_validate_evolution` 93-102 (now dead),
  `tfim_pauli_convention` 86-91, `TFIMStarModel` 34-54, import L5; trim exports
  L12/14/15 to drop TFIMStarModel/tfim_pauli_convention/tfim_product_basis_energy.
- `src/Observables.jl`: delete `measure_tfim_simple` 555-606, `TFIMObservableSummary`
  + `_validate_finite_summary` 489-524, `tfim_energy_density_decomposed_simple`
  430-442, `tfim_energy_density_star_simple` 417-428, `_tfim_*` 403-415, import
  L9 `using ..StarModels: TFIMStarModel, star_hamiltonian`, exports 19-20.
- `src/SquarePXPDynamics.jl`: remove includes 29-31; trim `using .StarModels:`
  (drop TFIMStarModel L69, tfim_pauli_convention L73, tfim_product_basis_energy
  L78); remove `using .Observables:` tfim L109-110; remove `using .Benchmarks:`/
  FiniteTFIM/FiniteMPSTFIM L147-168; fix exports L263/265/266, remove 289-290,
  remove 315-322.
- `test/runtests.jl`: remove entries for test_tfim_observables,
  test_tfim_finite_reference, test_finite_mps_tfim_reference, test_benchmarks
  (L15,25,27,28,31). KEEP test_tfim_schedule_reference, test_star_models.
- Kept tests referencing TFIM: `test_star_models.jl` remove L35 + testsets
  51-67, 69-85, 87-106 and rewrite L109 `TFIMStarModel(1.0,2.0)`→`PXPStarModel(false)`;
  `test_observables.jl` delete 75-85; `test_star_simple_update.jl` remove TFIM tail
  573-577 in the 560-578 testset.
- `Project.toml`: remove `ITensorMPS` from [deps] L11 and [compat] L26 (Aqua
  test_stale_deps will FAIL otherwise — ITensorMPS only used by FiniteMPSTFIM).
- `test/test_public_docs.jl`: remove L25 (references deleted `run_benchmark`).

## BLOCKED (do NOT mechanically delete)
- **Dead-export sweep**: the remaining candidates (`PXPSquareSpaceGroupBasis`,
  `ScarFinderAuditRow/Stability`, `ScarFinderCompressionInfo`,
  `SQUARE_IPEPS_SNAPSHOT_FORMAT_VERSION`) are reachable public API (fields of
  kept exported types / accessor returns). Not dead.
- **ScarFinder placeholder objectives** (`TargetEnergyObjective`,
  `LowVarianceObjective`, `energy_variance_proxy`, `JSONCandidateStore`): fused
  into the KEPT `CompositeObjective` (fields/constructor/_score_value/serialization)
  and pinned by `test_scarfinder.jl` exact-string assertions; `JSONCandidateStore`
  is a working tested feature sharing `_write_candidate_metadata` with JLD2.
  This is a dedicated refactor (~80 src lines), not a delete slice — defer.
