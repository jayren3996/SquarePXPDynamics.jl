# D2 Exact Finite Observables Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a tiny-cell exact finite iPEPS observable path and wire it into PXP validation/audit as an opt-in reference, without changing simple/local observable semantics.

**Architecture:** Keep `Observables.jl` as the cheap simple/local environment layer. Add `FiniteIPEPSObservables.jl` as a separate exact finite-contraction layer for small periodic cells, using the existing Gamma-lambda representation and absorbing each canonical periodic link weight exactly once. PXP validation gets nullable exact-finite density fields behind explicit `exact_finite_observables` and `exact_finite_max_sites` config fields; audit CSV/JSON expose these fields separately from simple and CTM fields.

**Tech Stack:** Julia 1.12, ITensors, existing `SquareIPEPSState`, existing PXP validation/audit structs, `Test`.

---

## Scope Boundaries

- Do not change `project_star!`.
- Do not change `local_density_simple`, `nearest_neighbor_density_simple`, or `star_expectation_simple`.
- Do not start CTM Stage 2.
- Do not design CTM-aware/full-update evolution.
- Do not add new CTM observables.
- Do not add tensor persistence.
- Do not change `log_norm` handling.
- Exact finite contractions are deliberately limited to tiny cells with a default `exact_finite_max_sites = 12`.
- “Exact finite” means exact contraction of the current `SquareIPEPSState`; it does not mean exact ED dynamics.

## File Structure

- Create: `src/FiniteIPEPSObservables.jl`
  - Owns dense finite contraction of `SquareIPEPSState`.
  - Owns exact one-site, nearest-neighbor, star, density, blockade, and PXP energy helpers.
- Modify: `src/SquarePXPDynamics.jl`
  - Includes the new module after `Observables.jl`.
  - Re-exports explicitly named `_finite` helper functions.
- Create: `test/test_finite_ipeps_observables.jl`
  - Verifies product limits and the known D=2 first-divergent state exact values.
- Modify: `test/runtests.jl`
  - Adds the new test file to `TEST_FILES`.
- Modify: `src/PXPValidation.jl`
  - Adds opt-in exact finite density measurement and nullable error fields.
  - Adds audit summary and CSV/JSON fields for exact finite density error.
- Modify: `test/test_pxp_validation.jl`
  - Verifies validation, JSON, audit summary, and CSV output include exact finite fields only when requested.
- Modify: `README.md`
  - Documents the opt-in exact finite tiny-cell validation path.
- Modify: `docs/superpowers/notes/2026-05-17-d2-measurement-localization.md`
  - Notes that the diagnostic helpers graduated into a scoped implementation.
- Modify: `memory/mid_term/decision_log.md`
  - Records the decision to keep simple observables unchanged and add opt-in exact finite references.

---

### Task 1: Red Test For Exact Finite Observable Module

**Files:**
- Create: `test/test_finite_ipeps_observables.jl`
- Modify: `test/runtests.jl`

- [ ] **Step 1: Add the new test file to the runner**

In `test/runtests.jl`, insert `"test_finite_ipeps_observables.jl"` immediately after `"test_observables_evolved.jl"`:

```julia
    "test_observables.jl",
    "test_observables_evolved.jl",
    "test_finite_ipeps_observables.jl",
    "test_gauge_diagnostics.jl",
```

- [ ] **Step 2: Write the failing exact finite tests**

Create `test/test_finite_ipeps_observables.jl`:

```julia
using LinearAlgebra

using SquarePXPDynamics.FiniteIPEPSObservables:
    dense_state_finite,
    exact_blockade_violation_finite,
    exact_density_finite,
    exact_nearest_neighbor_expectation_finite,
    exact_one_site_expectation_finite,
    exact_pxp_energy_density_finite,
    exact_star_expectation_finite

function _finite_obs_state_after_serial_stars(nstars; step = 0.02, maxdim = 2, cutoff = 1e-12)
    cell = PeriodicSquareUnitCell(3, 3)
    psi = product_square_ipeps(cell; state = :down, maxdim = 1)
    for center in cell.reps[1:nstars]
        project_star!(
            psi,
            center,
            step;
            evolution = :real,
            projected = true,
            maxdim,
            cutoff,
        )
    end
    return psi
end

@testset "exact finite iPEPS observables product limits" begin
    cell = PeriodicSquareUnitCell(3, 3)
    down = product_square_ipeps(cell; state = :down, maxdim = 1)
    up = product_square_ipeps(cell; state = :up, maxdim = 1)

    @test exact_density_finite(down) ≈ 0.0 atol = 1e-15
    @test exact_density_finite(up) ≈ 1.0 atol = 1e-15
    @test exact_blockade_violation_finite(down) ≈ 0.0 atol = 1e-15
    @test exact_blockade_violation_finite(up) ≈ 1.0 atol = 1e-15
    @test exact_pxp_energy_density_finite(down) ≈ 0.0 atol = 1e-15
end

@testset "exact finite dense state absorbs each canonical lambda once" begin
    cell = PeriodicSquareUnitCell(3, 3)
    psi = product_square_ipeps(cell; state = :down, maxdim = 2)
    for (i, key) in enumerate(sort(collect(keys(psi.link_weights)); by = k -> (k.site.y, k.site.x, String(k.dir))))
        set_link_weight!(psi, key.site, key.dir, [1.0 + i / 10, 0.25])
    end

    state = dense_state_finite(psi)
    expected = prod(link_weight(psi, key.site, key.dir)[1] for key in keys(psi.link_weights))

    @test count(x -> abs(x) > 1e-14, state) == 1
    @test state[end] ≈ expected atol = 1e-12 rtol = 1e-12
end

@testset "exact finite observables expose D2 simple measurement boundary" begin
    psi = _finite_obs_state_after_serial_stars(3)
    n = projector_up()
    nn = kron(n, n)
    zz = kron(pauli_z(), pauli_z())
    Hstar = square_pxp_star_hamiltonian()
    star_center_density = embed_one_site(n, 1, SQUARE_STAR_SITES)

    @test exact_density_finite(psi) ≈ 0.00013326224449912612 atol = 5e-15
    @test density_simple(psi) ≈ 0.000111054099003352 atol = 5e-15

    @test real(exact_one_site_expectation_finite(psi, SquareCoord(2, 1), n)) ≈
          0.0003997867121725804 atol = 5e-15
    @test local_density_simple(psi, SquareCoord(2, 1)) ≈
          0.0001999133387617944 atol = 5e-15

    @test real(exact_nearest_neighbor_expectation_finite(psi, SquareCoord(1, 1), :right, nn)) ≈
          nearest_neighbor_density_simple(psi, SquareCoord(1, 1), :right) atol = 5e-8
    @test real(exact_nearest_neighbor_expectation_finite(psi, SquareCoord(1, 1), :right, zz)) ≈
          0.9984005332366328 atol = 5e-15

    @test real(exact_star_expectation_finite(psi, SquareCoord(1, 1), star_center_density)) ≈
          0.00039994666951102603 atol = 5e-15
    @test star_expectation_simple(psi, SquareCoord(1, 1), Hstar) ≈
          exact_star_expectation_finite(psi, SquareCoord(1, 1), Hstar) atol = 5e-8
end

@testset "exact finite observables reject large cells by default" begin
    psi = product_square_ipeps(PeriodicSquareUnitCell(4, 4); state = :down, maxdim = 1)
    @test_throws ArgumentError exact_density_finite(psi)
end
```

- [ ] **Step 3: Run the test and verify RED**

Run:

```bash
julia --project=. test/runtests.jl test_finite_ipeps_observables.jl
```

Expected result:

```text
ERROR: LoadError: UndefVarError: `FiniteIPEPSObservables` not defined in `SquarePXPDynamics`
```

- [ ] **Step 4: Commit the red test**

```bash
git add test/runtests.jl test/test_finite_ipeps_observables.jl
git commit -m "test: specify exact finite iPEPS observables"
```

---

### Task 2: Implement Exact Finite Observable Module

**Files:**
- Create: `src/FiniteIPEPSObservables.jl`
- Modify: `src/SquarePXPDynamics.jl`
- Test: `test/test_finite_ipeps_observables.jl`

- [ ] **Step 1: Implement `FiniteIPEPSObservables.jl`**

Create `src/FiniteIPEPSObservables.jl`:

```julia
module FiniteIPEPSObservables

using ITensors
using LinearAlgebra

using ..SpinOps: projector_up
using ..SquareGeometry
using ..SquarePXP: SQUARE_STAR_SITES, square_pxp_star_hamiltonian
using ..SquareUnitCells
using ..SquareIPEPS

export dense_state_finite
export exact_one_site_expectation_finite, exact_nearest_neighbor_expectation_finite
export exact_star_expectation_finite, exact_density_finite
export exact_blockade_violation_finite, exact_pxp_energy_density_finite

const _DIRECTIONS = (:right, :up, :left, :down)

function _dense_index(values)
    idx = 1
    nsites = length(values)
    for (site, value) in enumerate(values)
        1 <= value <= 2 || throw(ArgumentError("basis values must be 1 or 2"))
        idx += (value - 1) * 2^(nsites - site)
    end
    return idx
end

function _check_tiny_finite_cell(psi::SquareIPEPSState, max_sites::Integer)
    nsites = length(psi.unitcell.reps)
    max_allowed = Int(max_sites)
    max_allowed >= 1 || throw(ArgumentError("max_sites must be positive"))
    nsites <= max_allowed || throw(
        ArgumentError(
            "exact finite contraction requested for $nsites sites; pass a larger max_sites explicitly",
        ),
    )
    all(c -> physical_dim(psi, c) == 2, psi.unitcell.reps) ||
        throw(ArgumentError("exact finite observables currently require physical dimension 2"))
    return nsites
end

function _absorb_all_weights_once(psi::SquareIPEPSState)
    tensors = Dict(c => copy(T) for (c, T) in psi.tensors)
    for key in keys(psi.link_weights)
        tensors[key.site] = absorb_link_weight(tensors[key.site], psi, key.site, key.dir)
    end
    return tensors
end

"""
    dense_state_finite(psi; max_sites = 12)

Return the dense `2^N` state vector obtained by exactly contracting the finite
periodic `SquareIPEPSState` in unit-cell representative order. Each canonical
periodic link weight is absorbed exactly once. This is an exact contraction of
the supplied iPEPS state, not an exact ED time-evolution reference.
"""
function dense_state_finite(psi::SquareIPEPSState; max_sites::Integer = 12)
    nsites = _check_tiny_finite_cell(psi, max_sites)
    tensors = _absorb_all_weights_once(psi)
    theta = ITensor()
    started = false
    for c in psi.unitcell.reps
        if started
            theta = @disable_warn_order theta * tensors[c]
        else
            theta = tensors[c]
            started = true
        end
    end

    phys = Tuple(physical_index(psi, c) for c in psi.unitcell.reps)
    for p in phys
        hasind(theta, p) || throw(ArgumentError("finite contraction lost physical index $p"))
    end

    state = zeros(ComplexF64, 2^nsites)
    for values in Iterators.product((1:2 for _ = 1:nsites)...)
        state[_dense_index(values)] = theta[(phys[i] => values[i] for i = 1:nsites)...]
    end
    return state
end

function _local_positions(cell::PeriodicSquareUnitCell, coords)
    sites = Tuple(wrap(cell, c) for c in coords)
    positions = Tuple(findfirst(==(site), cell.reps) for site in sites)
    all(!isnothing, positions) || throw(ArgumentError("observable sites must be in cell reps"))
    length(Set(positions)) == length(positions) ||
        throw(ArgumentError("observable sites must be distinct after wrapping"))
    return positions
end

function _local_expectation_from_state(state, nsites::Int, positions, O::AbstractMatrix)
    size(O) == (2^length(positions), 2^length(positions)) ||
        throw(ArgumentError("operator size does not match observable support"))
    normsq = sum(abs2, state)
    normsq > 0 || throw(ArgumentError("dense finite state has zero norm"))

    value = 0.0 + 0.0im
    for in_values in Iterators.product((1:2 for _ = 1:nsites)...)
        in_idx = _dense_index(in_values)
        amplitude = state[in_idx]
        iszero(amplitude) && continue
        local_in = ntuple(site -> in_values[positions[site]], length(positions))
        local_in_idx = _dense_index(local_in)
        for local_out in Iterators.product((1:2 for _ = 1:length(positions))...)
            local_out_idx = _dense_index(local_out)
            out_values = collect(in_values)
            for site = 1:length(positions)
                out_values[positions[site]] = local_out[site]
            end
            out_idx = _dense_index(Tuple(out_values))
            value += conj(state[out_idx]) * O[local_out_idx, local_in_idx] * amplitude
        end
    end
    return value / normsq
end

function _exact_local_expectation_finite(
    psi::SquareIPEPSState,
    coords,
    O::AbstractMatrix;
    max_sites::Integer = 12,
)
    nsites = _check_tiny_finite_cell(psi, max_sites)
    positions = _local_positions(psi.unitcell, coords)
    state = dense_state_finite(psi; max_sites)
    return _local_expectation_from_state(state, nsites, positions, O)
end

"""
    exact_one_site_expectation_finite(psi, c, O; max_sites = 12)

Return the exact finite contraction of one-site operator `O` at coordinate `c`
for the supplied `SquareIPEPSState`.
"""
function exact_one_site_expectation_finite(
    psi::SquareIPEPSState,
    c::SquareCoord,
    O::AbstractMatrix;
    max_sites::Integer = 12,
)
    size(O) == (2, 2) || throw(ArgumentError("one-site operator must be 2x2"))
    return _exact_local_expectation_finite(psi, (c,), O; max_sites)
end

"""
    exact_nearest_neighbor_expectation_finite(psi, c, dir, O; max_sites = 12)

Return the exact finite contraction of two-site nearest-neighbor operator `O`
on the bond from `c` in `dir`.
"""
function exact_nearest_neighbor_expectation_finite(
    psi::SquareIPEPSState,
    c::SquareCoord,
    dir::Symbol,
    O::AbstractMatrix;
    max_sites::Integer = 12,
)
    dir in _DIRECTIONS || throw(ArgumentError("direction must be :right, :up, :left, or :down"))
    size(O) == (4, 4) || throw(ArgumentError("two-site operator must be 4x4"))
    return _exact_local_expectation_finite(
        psi,
        (c, neighbor(psi.unitcell, c, dir)),
        O;
        max_sites,
    )
end

function _star_coords(cell::PeriodicSquareUnitCell, center::SquareCoord)
    c = wrap(cell, center)
    return (
        c,
        neighbor(cell, c, :right),
        neighbor(cell, c, :up),
        neighbor(cell, c, :left),
        neighbor(cell, c, :down),
    )
end

"""
    exact_star_expectation_finite(psi, center, O; max_sites = 12)

Return the exact finite contraction of a five-site square-star operator in
order `(center, right, up, left, down)`.
"""
function exact_star_expectation_finite(
    psi::SquareIPEPSState,
    center::SquareCoord,
    O::AbstractMatrix;
    max_sites::Integer = 12,
)
    size(O) == (2^SQUARE_STAR_SITES, 2^SQUARE_STAR_SITES) ||
        throw(ArgumentError("dense square-star operator must be 32x32"))
    return _exact_local_expectation_finite(psi, _star_coords(psi.unitcell, center), O; max_sites)
end

"""
    exact_density_finite(psi; max_sites = 12)

Return the average exact finite contraction of the Rydberg density over all
unit-cell representatives of the supplied iPEPS state.
"""
function exact_density_finite(psi::SquareIPEPSState; max_sites::Integer = 12)::Float64
    n = projector_up()
    values = [
        real(exact_one_site_expectation_finite(psi, c, n; max_sites)) for
        c in psi.unitcell.reps
    ]
    return sum(values) / length(values)
end

"""
    exact_blockade_violation_finite(psi; max_sites = 12)

Return the average exact finite contraction of nearest-neighbor `<n_i n_j>`
over canonical `:right` and `:up` bonds.
"""
function exact_blockade_violation_finite(
    psi::SquareIPEPSState;
    max_sites::Integer = 12,
)::Float64
    nn = kron(projector_up(), projector_up())
    total = 0.0
    count = 0
    for c in psi.unitcell.reps, dir in (:right, :up)
        total += real(exact_nearest_neighbor_expectation_finite(psi, c, dir, nn; max_sites))
        count += 1
    end
    return total / count
end

"""
    exact_pxp_energy_density_finite(psi; max_sites = 12)

Return the average exact finite contraction of the square-PXP star Hamiltonian
over all unit-cell representatives.
"""
function exact_pxp_energy_density_finite(
    psi::SquareIPEPSState;
    max_sites::Integer = 12,
)::Float64
    Hstar = square_pxp_star_hamiltonian()
    values = [
        real(exact_star_expectation_finite(psi, c, Hstar; max_sites)) for
        c in psi.unitcell.reps
    ]
    return sum(values) / length(values)
end

end
```

- [ ] **Step 2: Include and re-export the module**

In `src/SquarePXPDynamics.jl`, add after `include("Observables.jl")`:

```julia
include("FiniteIPEPSObservables.jl")
```

Add after the existing `.Observables` `using` block:

```julia
using .FiniteIPEPSObservables:
    dense_state_finite,
    exact_one_site_expectation_finite,
    exact_nearest_neighbor_expectation_finite,
    exact_star_expectation_finite,
    exact_density_finite,
    exact_blockade_violation_finite,
    exact_pxp_energy_density_finite
```

Add near the observable exports:

```julia
export dense_state_finite
export exact_one_site_expectation_finite, exact_nearest_neighbor_expectation_finite
export exact_star_expectation_finite, exact_density_finite
export exact_blockade_violation_finite, exact_pxp_energy_density_finite
```

- [ ] **Step 3: Run the exact finite tests and verify GREEN**

Run:

```bash
julia --project=. test/runtests.jl test_finite_ipeps_observables.jl
```

Expected result:

```text
Test Summary:              | Pass  Total
SquarePXPDynamics          |   16     16
```

- [ ] **Step 4: Run existing localization harness**

Run:

```bash
julia --project=. test/runtests.jl test_pxp_d2_localization.jl
```

Expected result:

```text
Test Summary:     | Pass  Broken  Total
SquarePXPDynamics |   50       5     55
```

- [ ] **Step 5: Commit exact finite module**

```bash
git add src/FiniteIPEPSObservables.jl src/SquarePXPDynamics.jl test/test_finite_ipeps_observables.jl test/runtests.jl
git commit -m "feat: add exact finite iPEPS observables"
```

---

### Task 3: Red Tests For Validation/Audit Exact-Finite Fields

**Files:**
- Modify: `test/test_pxp_validation.jl`

- [ ] **Step 1: Add validation report test**

Append this test near `"ED-vs-iPEPS validation report samples matched times"` in `test/test_pxp_validation.jl`:

```julia
@testset "PXP validation can attach exact finite tiny-cell density" begin
    config = PXPValidationConfig(
        3;
        total_time = 0.02,
        dt = 0.02,
        measure_every = 1,
        maxdim = 2,
        cutoff = 1e-12,
        exact_finite_observables = true,
        exact_finite_max_sites = 12,
    )
    report = validate_pxp_ed_ipeps(config; ctm_params = nothing)

    @test report.config.exact_finite_observables === true
    @test all(sample -> sample.exact_finite_density !== nothing, report.ipeps_samples)
    @test all(comparison -> comparison.ipeps_exact_finite_density !== nothing, report.comparisons)
    @test all(comparison -> comparison.density_error_exact_finite !== nothing, report.comparisons)
    @test report.comparisons[end].ipeps_exact_finite_density ≈
          0.0003996269892620211 atol = 5e-15
    @test abs(report.comparisons[end].density_error_exact_finite) < 1e-6
    @test abs(report.comparisons[end].density_error_simple) > 1e-4
end
```

Append this config-bound test near the same validation tests:

```julia
@testset "PXP validation rejects exact finite contraction above configured size" begin
    @test_throws ArgumentError PXPValidationConfig(
        4;
        total_time = 0.02,
        dt = 0.02,
        maxdim = 1,
        exact_finite_observables = true,
        exact_finite_max_sites = 12,
    )
end
```

- [ ] **Step 2: Add JSON artifact test assertions**

In `"PXP validation report writes JSON artifact"`, after `parsed.comparisons[1].ctm_reason === nothing`, add:

```julia
    @test haskey(parsed.config, :exact_finite_observables)
    @test parsed.config.exact_finite_observables === false
    @test parsed.config.exact_finite_max_sites == 12
    @test parsed.ipeps_samples[1].exact_finite_density === nothing
    @test parsed.comparisons[1].ipeps_exact_finite_density === nothing
    @test parsed.comparisons[1].density_error_exact_finite === nothing
```

Append this opt-in JSON round-trip test near `"PXP validation report writes JSON artifact"`:

```julia
@testset "PXP validation JSON preserves opt-in exact finite density" begin
    config = PXPValidationConfig(
        3;
        total_time = 0.02,
        dt = 0.02,
        measure_every = 1,
        maxdim = 2,
        cutoff = 1e-12,
        exact_finite_observables = true,
    )
    report = validate_pxp_ed_ipeps(config; ctm_params = nothing)
    path = tempname() * ".json"

    write_pxp_validation_json(report, path)
    parsed = JSON3.read(read(path, String))

    @test parsed.config.exact_finite_observables === true
    @test parsed.config.exact_finite_max_sites == 12
    @test parsed.ipeps_samples[end].exact_finite_density ≈ 0.0003996269892620211 atol =
        5e-15
    @test abs(parsed.comparisons[end].density_error_exact_finite) < 1e-6
end
```

- [ ] **Step 3: Add audit summary exact-field test**

In `"PXP audit campaign produces machine-readable summaries"`, after `@test summary.max_abs_density_error_ctm === nothing`, add:

```julia
    @test summary.max_abs_density_error_exact_finite === nothing
```

Then append a new test:

```julia
@testset "PXP audit campaign can opt into exact finite density summaries" begin
    config = PXPAuditConfig(;
        n_values = [3],
        total_time = 0.02,
        dt_values = [0.02],
        D_values = [2],
        cutoff_values = [1e-12],
        chi_values = Int[],
        exact_finite_observables = true,
        exact_finite_max_sites = 12,
    )
    report = run_pxp_audit_campaign(config; ctm_measure = _validation_fake_ctm_summary)
    summary = report.runs[1].summary

    @test report.config.exact_finite_observables === true
    @test summary.max_abs_density_error_exact_finite !== nothing
    @test summary.max_abs_density_error_exact_finite < 1e-6
    @test summary.max_abs_density_error_simple > 1e-4
end
```

- [ ] **Step 4: Add convergence and artifact serialization assertions**

Append this convergence test near `"PXP convergence report aggregates validation grid"`:

```julia
@testset "PXP convergence propagates exact finite density opt-in" begin
    base = PXPValidationConfig(
        3;
        total_time = 0.02,
        dt = 0.02,
        measure_every = 1,
        maxdim = 2,
        cutoff = 1e-12,
        exact_finite_observables = true,
    )
    sweep = PXPConvergenceConfig(
        base;
        dt_values = [0.02],
        D_values = [2],
        chi_values = Int[],
        cutoff_values = [1e-12],
    )

    report = validate_pxp_convergence(sweep)

    @test report.runs[1].config.exact_finite_observables === true
    @test report.runs[1].ipeps_samples[end].exact_finite_density !== nothing
    @test report.max_abs_density_error_exact_finite !== nothing
    @test report.max_abs_density_error_exact_finite < 1e-6
end
```

In `"PXP convergence report writes JSON artifact"`, after
`@test haskey(parsed.summary, :max_abs_density_error_simple)`, add:

```julia
    @test haskey(parsed.summary, :max_abs_density_error_exact_finite)
```

In `"PXP audit campaign writes JSON and CSV artifacts"`, after
`@test startswith(csv[1], "n,total_time,dt,D,cutoff")`, add:

```julia
    @test occursin("max_abs_density_error_exact_finite", csv[1])
```

Append this audit JSON/CSV opt-in round-trip test:

```julia
@testset "PXP audit JSON and CSV preserve opt-in exact finite summaries" begin
    config = PXPAuditConfig(;
        n_values = [3],
        total_time = 0.02,
        dt_values = [0.02],
        D_values = [2],
        cutoff_values = [1e-12],
        chi_values = Int[],
        exact_finite_observables = true,
    )
    report = run_pxp_audit_campaign(config)
    json_path = tempname() * ".json"
    csv_path = tempname() * ".csv"

    write_pxp_audit_json(report, json_path)
    write_pxp_audit_csv(report, csv_path)
    parsed = JSON3.read(read(json_path, String))
    csv = split(chomp(read(csv_path, String)), '\n')

    @test parsed.config.exact_finite_observables === true
    @test parsed.runs[1].summary.max_abs_density_error_exact_finite < 1e-6
    @test occursin("max_abs_density_error_exact_finite", csv[1])
    @test occursin("not_run,not_run", csv[2])
end
```

- [ ] **Step 5: Run validation tests and verify RED**

Run:

```bash
julia --project=. -e 'using Test, SquarePXPDynamics; include("test/test_pxp_validation.jl")'
```

Expected result:

```text
ERROR: MethodError: no method matching PXPValidationConfig(...; exact_finite_observables=true)
```

- [ ] **Step 6: Commit red validation tests**

```bash
git add test/test_pxp_validation.jl
git commit -m "test: require opt-in exact finite PXP validation fields"
```

---

### Task 4: Implement Validation/Audit Exact-Finite Integration

**Files:**
- Modify: `src/PXPValidation.jl`
- Test: `test/test_pxp_validation.jl`

- [ ] **Step 1: Import exact finite density**

At the top of `src/PXPValidation.jl`, add:

```julia
using ..FiniteIPEPSObservables: exact_density_finite
```

- [ ] **Step 2: Add config flags and early size validation**

Add these final fields to `PXPValidationConfig`:

```julia
    exact_finite_observables::Bool
    exact_finite_max_sites::Int
```

Add these keywords to its inner constructor:

```julia
exact_finite_observables::Bool = false,
exact_finite_max_sites::Integer = 12,
```

After schedule validation, add:

```julia
        exact_limit = _positive_int(exact_finite_max_sites, "exact_finite_max_sites")
        exact_finite_observables && n_int^2 > exact_limit &&
            throw(ArgumentError("exact finite observables require n^2 <= exact_finite_max_sites"))
```

Pass both values as final `new(...)` arguments:

```julia
            schedule,
            exact_finite_observables,
            exact_limit,
```

Update `_config_data(config::PXPValidationConfig)` to include:

```julia
        exact_finite_observables = config.exact_finite_observables,
        exact_finite_max_sites = config.exact_finite_max_sites,
```

Update `_copy_config` to preserve both fields:

```julia
        schedule = base.schedule,
        exact_finite_observables = base.exact_finite_observables,
        exact_finite_max_sites = base.exact_finite_max_sites,
```

- [ ] **Step 3: Add sample and comparison fields**

Add the final field to `PXPIPEPSSample`:

```julia
    exact_finite_density::Union{Nothing,Float64}
```

Add these fields to `PXPEDComparisonSample` after `ipeps_ctm_density`:

```julia
    ipeps_exact_finite_density::Union{Nothing,Float64}
```

and after `density_error_ctm`:

```julia
    density_error_exact_finite::Union{Nothing,Float64}
```

- [ ] **Step 4: Compute and compare exact finite density**

In `_comparison`, compute and pass the nullable exact value:

```julia
function _comparison(ed::PXPEEDSample, sample::PXPIPEPSSample)
    ctm_density = _ctm_density(sample.ctm)
    exact_density = sample.exact_finite_density
    return PXPEDComparisonSample(
        ed.step,
        ed.time,
        ed.return_probability,
        ed.excitation_density,
        sample.simple.density,
        ctm_density,
        exact_density,
        sample.simple.density - ed.excitation_density,
        ctm_density === nothing ? nothing : ctm_density - ed.excitation_density,
        exact_density === nothing ? nothing : exact_density - ed.excitation_density,
        sample.simple.blockade_violation,
        _ctm_blockade(sample.ctm),
        _ctm_trusted(sample.ctm),
        _ctm_reason(sample.ctm),
    )
end
```

In `validate_pxp_ed_ipeps`, after `simple = measure_simple(psi)`, add:

```julia
        exact_finite_density =
            config.exact_finite_observables ?
            exact_density_finite(psi; max_sites = config.exact_finite_max_sites) : nothing
```

Pass `exact_finite_density` as the final `PXPIPEPSSample` constructor argument:

```julia
                ctm,
                log_norm(psi),
                exact_finite_density,
```

- [ ] **Step 5: Serialize sample and comparison fields**

In `_ipeps_sample_data`, add:

```julia
        exact_finite_density = sample.exact_finite_density,
```

In `_comparison_data`, add:

```julia
        ipeps_exact_finite_density = sample.ipeps_exact_finite_density,
        density_error_exact_finite = sample.density_error_exact_finite,
```

- [ ] **Step 6: Add convergence summary field**

Add `max_abs_density_error_exact_finite::Union{Nothing,Float64}` to `PXPConvergenceReport` after `max_abs_density_error_ctm`.

In `validate_pxp_convergence`, add:

```julia
    exact_errors = [
        abs(c.density_error_exact_finite) for r in runs for c in r.comparisons
        if c.density_error_exact_finite !== nothing
    ]
```

and pass it in `PXPConvergenceReport(...)` after the CTM error field:

```julia
        _finite_max_or_nothing(exact_errors),
```

The resulting return block should be:

```julia
    return PXPConvergenceReport(
        config,
        runs,
        maximum(simple_errors),
        isempty(ctm_errors) ? nothing : maximum(ctm_errors),
        _finite_max_or_nothing(exact_errors),
        isempty(trust_flags) ? nothing : all(==(true), trust_flags),
    )
end
```

In `_convergence_report_data`, add:

```julia
            max_abs_density_error_exact_finite = report.max_abs_density_error_exact_finite,
```

- [ ] **Step 7: Add audit config and summary fields**

Add `exact_finite_observables::Bool` and `exact_finite_max_sites::Int` as the final fields of `PXPAuditConfig`, with constructor keywords:

```julia
exact_finite_observables::Bool = false,
exact_finite_max_sites::Integer = 12,
```

Add `exact_finite_max_sites::Int` as the final `PXPAuditConfig` field, validate it in the constructor, and pass it to validation configs in the constructor loop and `_audit_validation_config`:

```julia
        exact_limit = _positive_int(exact_finite_max_sites, "exact_finite_max_sites")
```

In the constructor validation loop:

```julia
                exact_finite_observables,
                exact_finite_max_sites = exact_limit,
```

Pass both values as final `PXPAuditConfig` constructor arguments:

```julia
            verbosity,
            seed,
            exact_finite_observables,
            exact_limit,
```

In `_audit_validation_config`:

```julia
        exact_finite_observables = config.exact_finite_observables,
        exact_finite_max_sites = config.exact_finite_max_sites,
```

Add both fields to `_audit_config_data`:

```julia
        exact_finite_observables = config.exact_finite_observables,
        exact_finite_max_sites = config.exact_finite_max_sites,
```

Add `max_abs_density_error_exact_finite::Union{Nothing,Float64}` to `PXPAuditSummary` after `max_abs_density_error_ctm`.

In `_audit_summary`, add:

```julia
    exact_density_errors = [
        abs(c.density_error_exact_finite) for c in validation.comparisons
        if c.density_error_exact_finite !== nothing
    ]
```

and pass:

```julia
        _finite_max_or_nothing(exact_density_errors),
```

after `_finite_max_or_nothing(ctm_density_errors)`.

- [ ] **Step 8: Serialize audit fields**

In `_audit_summary_data`, add:

```julia
        max_abs_density_error_exact_finite = summary.max_abs_density_error_exact_finite,
```

In `PXP_AUDIT_CSV_HEADER`, add:

```julia
    "max_abs_density_error_exact_finite",
```

immediately after `"max_abs_density_error_ctm"`.

In `_audit_csv_row`, add:

```julia
        summary.max_abs_density_error_exact_finite,
```

immediately after `summary.max_abs_density_error_ctm`.

- [ ] **Step 9: Run validation tests and verify GREEN**

Run:

```bash
julia --project=. -e 'using Test, SquarePXPDynamics; include("test/test_pxp_validation.jl")'
```

Expected result:

```text
Test Summary: ... | Pass  Total
```

- [ ] **Step 10: Commit validation integration**

```bash
git add src/PXPValidation.jl test/test_pxp_validation.jl
git commit -m "feat: add exact finite density to PXP validation"
```

---

### Task 5: Documentation And Decision Update

**Files:**
- Modify: `README.md`
- Modify: `docs/superpowers/notes/2026-05-17-d2-measurement-localization.md`
- Modify: `memory/mid_term/decision_log.md`

- [ ] **Step 1: Update README validation section**

In `README.md`, after the fast JSON artifact example, add:

```markdown
For tiny periodic cells, validation can attach an exact finite contraction
density alongside simple/local and CTM fields:

```julia
config = PXPValidationConfig(
    3;
    total_time = 0.02,
    dt = 0.02,
    maxdim = 2,
    exact_finite_observables = true,
    exact_finite_max_sites = 12,
)
report = validate_pxp_ed_ipeps(config; ctm_params = nothing)
```

The exact finite path is intentionally size-limited and uses dense `2^N`
contractions of the supplied `SquareIPEPSState`. It is a debugging and
tiny-cell validation reference, not exact ED dynamics and not a replacement for
CTM-backed thermodynamic measurements.
```

- [ ] **Step 2: Update D2 localization note**

In `docs/superpowers/notes/2026-05-17-d2-measurement-localization.md`, append:

```markdown
## Follow-Up Implementation Plan

Plan `docs/superpowers/plans/2026-05-17-d2-exact-finite-observables.md`
promotes the test-only exact finite helpers into a size-limited module and
wires exact finite density into PXP validation/audit as an opt-in field. The
plan keeps simple/local observable formulas unchanged.
```

- [ ] **Step 3: Update decision log**

In `memory/mid_term/decision_log.md`, add a new entry above older 2026-05-17 entries:

```markdown
## 2026-05-17 - Add Opt-In Exact Finite Observable References For Tiny Cells

Decision:

Add a size-limited exact finite iPEPS observable path for tiny periodic cells
and wire exact finite density into PXP validation/audit as an opt-in reference.
Keep simple/local observables unchanged.

Reason:

The D=2 PXP anomaly was localized to treating simple/local D>1 observables as
exact finite observables on a loopy periodic PEPS. Exact finite contraction is
useful for tiny debug cells, while simple/local and CTM measurements have
separate contracts.

Consequences:

No-CTM D>1 audit summaries can report both simple diagnostic error and exact
finite density error when requested. The exact path is dense and size-limited;
CTM Stage 2, CTM-aware/full-update design, new CTM observables, and tensor
persistence remain postponed.

Source:

`src/FiniteIPEPSObservables.jl`; `src/PXPValidation.jl`;
`test/test_finite_ipeps_observables.jl`; `test/test_pxp_validation.jl`

Status: active
```

- [ ] **Step 4: Run docs/source whitespace check**

Run:

```bash
git diff --check
```

Expected result: no output and exit code `0`.

- [ ] **Step 5: Commit docs**

```bash
git add README.md docs/superpowers/notes/2026-05-17-d2-measurement-localization.md memory/mid_term/decision_log.md
git commit -m "docs: document exact finite validation boundary"
```

---

### Task 6: Final Verification

**Files:**
- Read-only verification across changed modules and required regression paths.

- [ ] **Step 1: Run exact finite observable tests**

```bash
julia --project=. test/runtests.jl test_finite_ipeps_observables.jl
```

Expected result: all tests pass.

- [ ] **Step 2: Run D2 localization tests**

```bash
julia --project=. test/runtests.jl test_pxp_d2_localization.jl
```

Expected result:

```text
Test Summary:     | Pass  Broken  Total
SquarePXPDynamics |   50       5     55
```

- [ ] **Step 3: Run existing observable tests**

```bash
julia --project=. -e 'using Test, SquarePXPDynamics; include("test/test_observables.jl")'
```

Expected result: all tests pass.

- [ ] **Step 4: Run PXP validation tests**

```bash
julia --project=. -e 'using Test, SquarePXPDynamics; include("test/test_pxp_validation.jl")'
```

Expected result: all tests pass.

- [ ] **Step 5: Run full package tests**

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected result: full suite passes, with only the existing expected broken D2 exactness markers.

- [ ] **Step 6: Check final git state**

```bash
git status --short --branch
git log --oneline --decorate -5
```

Expected result: branch `codex/d2-anomaly-localization` contains the task commits, with no unrelated changes staged.
