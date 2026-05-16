# S0-S7 Slice 1 Reconciliation Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Close the low-risk S0-S2/S5 reconciliation gaps by adding small iPEPS state helper APIs, link-weight normalization, and missing PXP `:x_plus` observable regression coverage.

**Architecture:** Preserve the current concrete `SquareIPEPSState` architecture. Add thin helper functions in `src/SquareIPEPS.jl`, export them through `src/SquarePXPDynamics.jl`, and strengthen existing test files without introducing speculative backend abstractions.

**Tech Stack:** Julia 1.12, ITensors, existing `SquarePXPDynamics` modules, `Test`.

---

## File Structure

- Modify `src/SquareIPEPS.jl`: add `unitcell_reps`, `physical_dim`, `simple_weight_dim`, `copy_state`, and `normalize_link_weights!`.
- Modify `src/SquarePXPDynamics.jl`: import and export the new public APIs.
- Modify `test/test_square_ipeps.jl`: cover state helper APIs and deep-copy behavior.
- Modify `test/test_square_ipeps_s2.jl`: cover `normalize_link_weights!`.
- Modify `test/test_observables_evolved.jl`: add PXP `:x_plus` dense five-site expectation coverage.
- Modify `README.md`: document the S0-S7 reconciliation status and helper APIs.
- Modify `memory/mid_term/milestones.md`: record Slice 1 completion evidence after implementation.

## Task 1: Add iPEPS Helper API Tests

**Files:**
- Modify: `test/test_square_ipeps.jl`

- [ ] **Step 1: Add failing tests for helper APIs**

Add this testset after `"product square iPEPS constructor"`:

```julia
@testset "square iPEPS public helper APIs" begin
    cell = PeriodicSquareUnitCell(4, 4)
    psi = product_square_ipeps(cell; state = :down, maxdim = 2)

    @test unitcell_reps(psi) == cell.reps
    @test unitcell_reps(psi) !== cell.reps
    @test physical_dim(psi, SquareCoord(1, 1)) == 2
    @test physical_dim(psi, SquareCoord(5, 5)) == 2
    @test simple_weight_dim(psi, SquareCoord(1, 1), :right) == 2
    @test simple_weight_dim(psi, SquareCoord(1, 1), :up) == 2
    @test simple_weight_dim(psi, SquareCoord(1, 1), :left) == 2
    @test simple_weight_dim(psi, SquareCoord(1, 1), :down) == 2
    @test_throws ArgumentError simple_weight_dim(psi, SquareCoord(1, 1), :diagonal)

    copied = copy_state(psi)
    @test copied !== psi
    @test copied.unitcell == psi.unitcell
    @test copied.maxdim == psi.maxdim
    @test copied.gauge == psi.gauge
    @test state_version(copied) == state_version(psi)
    @test log_norm(copied) == log_norm(psi)

    set_link_weight!(copied, SquareCoord(1, 1), :right, [0.6, 0.8])
    @test link_weight(copied, SquareCoord(1, 1), :right) == [0.6, 0.8]
    @test link_weight(psi, SquareCoord(1, 1), :right) == [1.0, 0.0]
end
```

- [ ] **Step 2: Run the focused test and verify RED**

Run:

```bash
julia --project=. test/runtests.jl test_square_ipeps.jl
```

Expected: fail with `UndefVarError` for `unitcell_reps`.

## Task 2: Implement iPEPS Helper APIs

**Files:**
- Modify: `src/SquareIPEPS.jl`
- Modify: `src/SquarePXPDynamics.jl`
- Test: `test/test_square_ipeps.jl`

- [ ] **Step 1: Export helper APIs from `SquareIPEPS`**

In `src/SquareIPEPS.jl`, update the export block:

```julia
export SquareIPEPSState
export product_square_ipeps, checkerboard_square_ipeps
export unitcell_reps, physical_dim, simple_weight_dim, copy_state
export physical_index, link_index
```

- [ ] **Step 2: Add helper implementations**

Insert after `log_norm(psi::SquareIPEPSState)::Float64 = psi.log_norm_value[]`:

```julia
"""
    unitcell_reps(psi)

Return a copy of the periodic unit-cell representatives used by `psi`.
"""
unitcell_reps(psi::SquareIPEPSState)::Vector{SquareCoord} = copy(psi.unitcell.reps)

"""
    physical_dim(psi, c)

Return the physical Hilbert-space dimension at coordinate `c`, after periodic
wrapping into the unit cell.
"""
physical_dim(psi::SquareIPEPSState, c::SquareCoord)::Int = dim(physical_index(psi, c))

"""
    simple_weight_dim(psi, c, dir)

Return the length of the simple-update link-weight vector on the nearest
neighbor bond from `c` in direction `dir`.
"""
function simple_weight_dim(psi::SquareIPEPSState, c::SquareCoord, dir::Symbol)::Int
    _validate_link_direction(dir)
    return length(_validated_link_weight(psi, c, dir))
end

"""
    copy_state(psi)

Return a deep mutable copy of `psi` with independent tensors, link-weight
vectors, mutation counter, and log-normalization ledger.
"""
function copy_state(psi::SquareIPEPSState)::SquareIPEPSState
    return SquareIPEPSState(
        psi.unitcell,
        Dict(c => copy(T) for (c, T) in psi.tensors),
        copy(psi.physical_indices),
        copy(psi.link_indices),
        Dict(key => copy(lambda) for (key, lambda) in psi.link_weights),
        psi.maxdim,
        psi.gauge,
        Ref(state_version(psi)),
        Ref(log_norm(psi)),
    )
end
```

- [ ] **Step 3: Import and export from top-level module**

In `src/SquarePXPDynamics.jl`, change:

```julia
using .SquareIPEPS: SquareIPEPSState, product_square_ipeps, checkerboard_square_ipeps
```

to:

```julia
using .SquareIPEPS: SquareIPEPSState, product_square_ipeps, checkerboard_square_ipeps
using .SquareIPEPS: unitcell_reps, physical_dim, simple_weight_dim, copy_state
```

Add to the export section:

```julia
export unitcell_reps, physical_dim, simple_weight_dim, copy_state
```

- [ ] **Step 4: Run focused helper tests and verify GREEN**

Run:

```bash
julia --project=. test/runtests.jl test_square_ipeps.jl test_public_docs.jl
```

Expected: pass.

- [ ] **Step 5: Commit helper APIs**

Run:

```bash
git add src/SquareIPEPS.jl src/SquarePXPDynamics.jl test/test_square_ipeps.jl
git commit -m "feat: add square iPEPS helper APIs"
```

## Task 3: Add Link-Weight Normalization

**Files:**
- Modify: `test/test_square_ipeps_s2.jl`
- Modify: `src/SquareIPEPS.jl`
- Modify: `src/SquarePXPDynamics.jl`

- [ ] **Step 1: Add failing tests**

Append to `test/test_square_ipeps_s2.jl` after `"link weight entropy diagnostics"`:

```julia
@testset "normalize all iPEPS link weights" begin
    cell = PeriodicSquareUnitCell(4, 4)
    psi = product_square_ipeps(cell; state = :down, maxdim = 2)
    c = SquareCoord(1, 1)

    set_link_weight!(psi, c, :right, [3.0, 4.0])
    set_link_weight!(psi, c, :up, [1.0, 1.0])
    version_before = state_version(psi)

    normalized = normalize_link_weights!(psi)

    @test normalized === psi
    @test state_version(psi) == version_before + 1
    @test link_weight(psi, c, :right) ≈ [0.6, 0.8]
    @test link_weight(psi, c, :up) ≈ [1 / sqrt(2), 1 / sqrt(2)]
    @test all(lambda -> isapprox(norm(lambda), 1; atol = 1e-12), values(psi.link_weights))

    stable_version = state_version(psi)
    normalize_link_weights!(psi)
    @test state_version(psi) == stable_version

    psi.link_weights[bondkey(cell, c, :right)] = [0.0, 0.0]
    @test_throws ArgumentError normalize_link_weights!(psi)
end
```

- [ ] **Step 2: Run focused S2 tests and verify RED**

Run:

```bash
julia --project=. test/runtests.jl test_square_ipeps_s2.jl
```

Expected: fail with `UndefVarError: normalize_link_weights! not defined`.

- [ ] **Step 3: Export `normalize_link_weights!`**

In `src/SquareIPEPS.jl`, change:

```julia
export weight_entropy, bond_entropy, all_bond_entropies
```

to:

```julia
export weight_entropy, bond_entropy, all_bond_entropies, normalize_link_weights!
```

- [ ] **Step 4: Implement normalization**

Insert after `all_bond_entropies`:

```julia
"""
    normalize_link_weights!(psi)

Normalize every stored simple-update link-weight vector in `psi` to unit
Euclidean norm. The state version is incremented only when at least one stored
vector changes. Invalid spectra, including all-zero vectors, throw
`ArgumentError`.
"""
function normalize_link_weights!(psi::SquareIPEPSState)::SquareIPEPSState
    changed = false
    normalized = Dict{BondKey,Vector{Float64}}()
    for (key, lambda) in psi.link_weights
        values = _validate_link_weight_values(link_index(psi, key.site, key.dir), lambda)
        scale = norm(values)
        isfinite(scale) && scale > 0 ||
            throw(ArgumentError("link weights must not all be zero"))
        next_values = values ./ scale
        normalized[key] = next_values
        changed |= !isapprox(values, next_values; atol = 1e-14, rtol = 1e-14)
    end
    if changed
        empty!(psi.link_weights)
        for (key, values) in normalized
            psi.link_weights[key] = values
        end
        _mark_mutated!(psi)
    end
    return psi
end
```

- [ ] **Step 5: Import and export top-level API**

In `src/SquarePXPDynamics.jl`, add `normalize_link_weights!` to the
`using .SquareIPEPS:` import line that currently imports entropy helpers, and
add it to the export section:

```julia
export weight_entropy, bond_entropy, all_bond_entropies, normalize_link_weights!
```

- [ ] **Step 6: Run focused tests and verify GREEN**

Run:

```bash
julia --project=. test/runtests.jl test_square_ipeps_s2.jl test_public_docs.jl
```

Expected: pass.

- [ ] **Step 7: Commit normalization**

Run:

```bash
git add src/SquareIPEPS.jl src/SquarePXPDynamics.jl test/test_square_ipeps_s2.jl
git commit -m "feat: normalize square iPEPS link weights"
```

## Task 4: Add PXP `:x_plus` Observable Regression

**Files:**
- Modify: `test/test_observables_evolved.jl`

- [ ] **Step 1: Add regression test**

Append after `"simple observables on product states"`:

```julia
@testset "PXP x-plus star expectation matches dense product reference" begin
    cell = PeriodicSquareUnitCell(10, 10)
    psi = product_square_ipeps(cell; state = :x_plus, maxdim = 1)
    c = SquareCoord(5, 5)
    Hstar = square_pxp_star_hamiltonian()

    plus = fill(inv(sqrt(2)), 2)
    dense = zeros(ComplexF64, 2^SQUARE_STAR_SITES)
    for values in Iterators.product((1:2 for _ = 1:SQUARE_STAR_SITES)...)
        dense[_dense_square_star_index_obs(values)] = prod(plus[value] for value in values)
    end
    expected = dot(dense, Hstar * dense) / dot(dense, dense)

    @test star_expectation_simple(psi, c, Hstar) ≈ expected atol = 1e-12
end
```

- [ ] **Step 2: Run focused observable tests**

Run:

```bash
julia --project=. test/runtests.jl test_observables_evolved.jl
```

Expected: pass. If it passes immediately, record it as regression coverage for
existing behavior rather than new behavior.

- [ ] **Step 3: Commit regression test**

Run:

```bash
git add test/test_observables_evolved.jl
git commit -m "test: cover PXP x-plus star observables"
```

## Task 5: Update User-Facing Milestone Docs

**Files:**
- Modify: `README.md`
- Modify: `memory/mid_term/milestones.md`

- [ ] **Step 1: Update README status**

In `README.md`, extend the status section with:

```markdown
The original S0-S7 implementation plan has been reconciled against the current
architecture in `docs/superpowers/specs/2026-05-16-s0-s7-completion-design.md`.
The current S0.5/S1 backend-facade items are superseded by the concrete custom
ITensors iPEPS stack unless a second update backend is introduced. Remaining
full S7 work is S7b: CTM local norm matrices, readiness checks, and
transactional gauge conditioning.
```

- [ ] **Step 2: Update milestone memory**

In `memory/mid_term/milestones.md`, add under completed/current work:

```markdown
- Confirmed: S0-S7 completion has been reconciled against the current
  architecture. Slice 1 adds public iPEPS helper APIs, link-weight
  normalization, and PXP `:x_plus` observable regression coverage while
  preserving the custom ITensors update plus PEPSKit measurement boundary.
- Source: `docs/superpowers/specs/2026-05-16-s0-s7-completion-design.md`
- Source: `src/SquareIPEPS.jl`
- Source: `test/test_square_ipeps.jl`
- Source: `test/test_square_ipeps_s2.jl`
- Source: `test/test_observables_evolved.jl`
```

- [ ] **Step 3: Run documentation checks**

Run:

```bash
git diff --check
julia --project=. test/runtests.jl test_public_docs.jl
```

Expected: both pass.

- [ ] **Step 4: Commit documentation**

Run:

```bash
git add README.md memory/mid_term/milestones.md
git commit -m "docs: record S0-S7 slice 1 reconciliation"
```

## Task 6: Final Slice 1 Verification

**Files:**
- Inspect: all changed files

- [ ] **Step 1: Run focused S0-S2/S5 checks**

Run:

```bash
julia --project=. test/runtests.jl \
  test_square_ipeps.jl \
  test_square_ipeps_s2.jl \
  test_observables_evolved.jl \
  test_public_docs.jl
```

Expected: pass.

- [ ] **Step 2: Run full package tests**

Run:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: pass.

- [ ] **Step 3: Inspect git status and recent commits**

Run:

```bash
git status --short --branch
git log --oneline -5
```

Expected: clean branch ahead of `main` by the Slice 1 commits.

