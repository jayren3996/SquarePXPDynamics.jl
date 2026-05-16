# S0-S7 Slice 2 Diagnostics Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Improve S3/S4 auditability by recording all pre-update touched link minima in star-update diagnostics and adding model/protocol metadata to evolution logs.

**Architecture:** Preserve existing update math and logging behavior. Add one backward-compatible diagnostic field to `StarUpdateInfo` and one metadata field to `EvolutionLog`, with tests that prove existing PXP projected/unprojected choices are reconstructible from raw logs.

**Tech Stack:** Julia 1.12, existing `SquarePXPDynamics` modules, `Test`.

---

## Task 1: StarUpdateInfo Touched Lambda Diagnostics

**Files:**
- Modify: `test/test_star_simple_update.jl`
- Modify: `src/StarSimpleUpdate.jl`

- [ ] **Step 1: Add failing tests**

Add a test that sets distinct internal and external link weights before
`project_star!`, then asserts `info.touched_min_lambda` contains canonical
`BondKey` entries for both kinds of touched bonds.

- [ ] **Step 2: Verify RED**

Run:

```bash
julia --project=. test/runtests.jl test_star_simple_update.jl
```

Expected: fail because `StarUpdateInfo` has no `touched_min_lambda` field.

- [ ] **Step 3: Implement diagnostics**

Add `touched_min_lambda::Dict{BondKey,Float64}` to `StarUpdateInfo`. Compute it
transactionally before mutation by collecting the four center-leaf bonds and
the three external bonds around each leaf, canonicalizing with `bondkey`, and
recording `minimum(link_weight(psi, key.site, key.dir))`.

- [ ] **Step 4: Verify GREEN**

Run:

```bash
julia --project=. test/runtests.jl test_star_simple_update.jl
```

Expected: pass.

## Task 2: EvolutionLog Model Metadata

**Files:**
- Modify: `test/test_ipeps_evolution.jl`
- Modify: `src/IPEPSEvolution.jl`

- [ ] **Step 1: Add failing tests**

Add assertions that real logs expose model metadata:

```julia
@test explicit_log.model_metadata.model_type == "PXPStarModel"
@test explicit_log.model_metadata.pxp_projected === true
@test legacy_unprojected_log.model_metadata.pxp_projected === false
```

- [ ] **Step 2: Verify RED**

Run:

```bash
julia --project=. test/runtests.jl test_ipeps_evolution.jl
```

Expected: fail because `EvolutionLog` has no `model_metadata` field.

- [ ] **Step 3: Implement metadata**

Add `model_metadata::NamedTuple` to `EvolutionLog`. Build it from the
evolution protocol once per `evolve!` call:

```julia
model = model_at(protocol, 0.0, 1)
(
    protocol_type = string(nameof(typeof(protocol))),
    model_type = string(nameof(typeof(model))),
    pxp_projected = model isa PXPStarModel ? model.projected : nothing,
)
```

- [ ] **Step 4: Verify GREEN**

Run:

```bash
julia --project=. test/runtests.jl test_ipeps_evolution.jl test_public_docs.jl
```

Expected: pass.

## Task 3: Documentation And Verification

**Files:**
- Modify: `README.md`
- Modify: `memory/mid_term/milestones.md`

- [ ] **Step 1: Document Slice 2 diagnostics**

Mention `StarUpdateInfo.touched_min_lambda` and `EvolutionLog.model_metadata`
near the README status or shipped-feature list.

- [ ] **Step 2: Verify focused tests**

Run:

```bash
julia --project=. test/runtests.jl \
  test_star_simple_update.jl \
  test_ipeps_evolution.jl \
  test_public_docs.jl
```

Expected: pass.

- [ ] **Step 3: Commit Slice 2**

Run:

```bash
git add src/StarSimpleUpdate.jl src/IPEPSEvolution.jl test/test_star_simple_update.jl test/test_ipeps_evolution.jl README.md memory/mid_term/milestones.md
git commit -m "feat: strengthen S0-S7 update diagnostics"
```

