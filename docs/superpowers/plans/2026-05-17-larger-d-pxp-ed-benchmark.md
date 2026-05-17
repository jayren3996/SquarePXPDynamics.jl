# M3 Larger-D PXP ED Benchmark Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build a documented larger-bond-dimension PXP dynamics benchmark campaign that compares iPEPS runs against exact finite PBC ED where that comparison is valid, with exact finite 3x3 iPEPS observables and honest global-observable labeling for 5x5 and 7x7 PBC ED.

**Architecture:** Keep the current EDKit-backed `FinitePXPEEDBenchmark.jl` path as the PBC symmetric ED reference and keep simple-update evolution unchanged. Extend the validation layer with a new M3 report that reuses one ED trajectory per `(n, dt, total_time, ED controls)` across a D/cutoff iPEPS sweep, records runtime and provenance fields, and labels observable modes explicitly. Exact finite iPEPS contraction is used only for tiny cells; larger odd PBC systems use symmetric PBC global density and ED return probability, not central-region observables.

**Tech Stack:** Julia 1.12, ITensors, EDKit, PEPSKit/TensorKit only through existing optional CTM hooks, JSON3, existing `test/runtests.jl` allowlist.

---

## Required Skill Flow For Implementation

- Before editing code, use `superpowers:test-driven-development`.
- Use `superpowers:systematic-debugging` for any test failure, benchmark mismatch, or unexpected D-sweep value.
- Use `superpowers:verification-before-completion` before claiming the implementation is complete.
- Do not run `project-memory-curator`.
- Do not update short-term handoff files unless the user explicitly asks.

## Multi-Agent Workstreams

- **Agent A: ED Capacity And Observable Semantics**
  - Owns `src/FinitePXPEEDBenchmark.jl`, ED observable-provenance helpers, and tests proving symmetric PBC ED cannot be labeled as central-region data.
  - Measures or scripts capacity for 5x5, 6x6, and 7x7 without putting huge jobs in default tests.
  - Recommends symmetric PBC global observables for M3; unreduced/open-boundary central-region ED is documented as a future path.

- **Agent B: iPEPS Larger-D Dynamics Harness**
  - Owns the reusable ED-plus-iPEPS validation loop in `src/PXPValidation.jl`.
  - Builds the D sweep for `D = 1, 2, 3, 4`, exact finite 3x3 density and return probability where possible, truncation/log-norm summaries, and reversibility diagnostics.

- **Agent C: Report Schema And CLI**
  - Owns M3 JSON/CSV writers and `scripts/pxp_larger_d_ed_benchmark.jl`.
  - Uses env-var parsing consistent with the existing scripts and writes nested JSON plus flat CSV.

- **Agent D: Tests And Docs**
  - Owns `test/test_pxp_larger_d_ed_benchmark.jl`, `test/runtests.jl`, README updates, the M3 note, and the decision-log entry.
  - Keeps 7x7 out of default tests and gates extended PXP ED smoke tests.

## Scope Boundaries

- Do not implement CTM-aware/full-update evolution.
- Do not add broad new CTM observables.
- Do not claim publication-grade physics.
- Do not treat `density_simple` as exact finite truth for D>1.
- Do not run huge 7x7 jobs in normal tests.
- Do not call a symmetric-sector global density a central 3x3 or local-region observable.
- Do not add an unreduced/open-boundary ED implementation in M3 unless the user explicitly expands the milestone after seeing this plan.

## File Map

- Modify `src/FiniteIPEPSObservables.jl`: add exact finite all-down return probability helper.
- Modify `src/FinitePXPEEDBenchmark.jl`: add observable-provenance helpers for the current symmetry-reduced PBC ED path and explicit rejecting methods for local/region observables on that basis.
- Modify `src/PXPValidation.jl`: add reusable ED-result validation helper plus M3 config, run, summary, report, JSON writer, and CSV writer.
- Modify `src/SquarePXPDynamics.jl`: include/re-export new public helpers and M3 benchmark API.
- Create `scripts/pxp_larger_d_ed_benchmark.jl`: env-var driven M3 benchmark runner.
- Create `test/test_pxp_larger_d_ed_benchmark.jl`: fast schema, semantics, and 3x3 exact-finite D-sweep tests.
- Modify `test/runtests.jl`: add the new test file to `TEST_FILES`.
- Modify `README.md`: document M3 benchmark modes, commands, and observable validity.
- Create `docs/superpowers/notes/2026-05-17-m3-larger-d-pxp-ed-benchmark.md`: record the benchmark contract and example campaign.
- Modify `memory/mid_term/decision_log.md`: record the comparison-boundary decision after implementation is real.

---

### Task 1: ED Observable Semantics And Exact Finite Return Probability

**Files:**
- Modify: `src/FiniteIPEPSObservables.jl`
- Modify: `src/FinitePXPEEDBenchmark.jl`
- Modify: `src/SquarePXPDynamics.jl`
- Test: `test/test_pxp_larger_d_ed_benchmark.jl`
- Test: `test/runtests.jl`

- [ ] **Step 1: Add the new test file to the runner**

In `test/runtests.jl`, insert the new test after `test_pxp_ed_benchmark.jl`:

```julia
    "test_pxp_ed_benchmark.jl",
    "test_pxp_larger_d_ed_benchmark.jl",
    "test_benchmarks.jl",
```

- [ ] **Step 2: Write failing ED provenance and exact return tests**

Create `test/test_pxp_larger_d_ed_benchmark.jl` with this first block:

```julia
using Test
using JSON3
using SquarePXPDynamics

const RUN_EXTENDED_PXP_ED_TESTS =
    get(ENV, "SQUAREPXP_EXTENDED_PXP_ED_TESTS", "") == "1" ||
    get(ENV, "SQUAREPXP_EXTENDED_TESTS", "") == "1"

function _csv_cell(lines, name)
    header = split(lines[1], ','; keepempty = true)
    row = split(lines[2], ','; keepempty = true)
    index = findfirst(==(name), header)
    index === nothing && error("missing CSV column $name")
    return row[index]
end

@testset "PXP ED observable provenance rejects central-region claims" begin
    basis = pxp_ed_space_group_basis(5)

    @test pxp_ed_boundary_condition(basis) === :periodic
    @test pxp_ed_symmetry_sector(basis) === :fully_symmetric_space_group
    @test pxp_ed_observable_scope(basis) === :pbc_global_site_average
    @test pxp_ed_reference_label(basis) == "finite_pbc_global_density"
    @test pxp_ed_group_order(basis) == 8 * 5^2

    @test_throws ArgumentError pxp_ed_site_density_operator(basis, 13)
    @test_throws ArgumentError pxp_ed_region_density_operator(basis, 1:9)
end

@testset "exact finite return probability is available only through tiny-cell contraction" begin
    cell = PeriodicSquareUnitCell(3, 3)
    down = product_square_ipeps(cell; state = :down, maxdim = 1)
    up = product_square_ipeps(cell; state = :up, maxdim = 1)

    @test exact_all_down_return_probability_finite(down; max_sites = 9) ≈ 1.0 atol = 1e-15
    @test exact_all_down_return_probability_finite(up; max_sites = 9) ≈ 0.0 atol = 1e-15

    large = product_square_ipeps(PeriodicSquareUnitCell(4, 4); state = :down, maxdim = 1)
    @test_throws ArgumentError exact_all_down_return_probability_finite(large; max_sites = 9)
end
```

- [ ] **Step 3: Run the test and verify RED**

Run:

```bash
julia --project=. test/runtests.jl test_pxp_larger_d_ed_benchmark.jl
```

Expected: FAIL with missing symbols such as:

```text
UndefVarError: `pxp_ed_boundary_condition` not defined
```

- [ ] **Step 4: Implement exact finite all-down return probability**

In `src/FiniteIPEPSObservables.jl`, add the export:

```julia
export exact_all_down_return_probability_finite
```

Then add this function after `exact_density_finite`:

```julia
"""
    exact_all_down_return_probability_finite(psi; max_sites = 12)

Return the normalized finite-contraction probability of the all-down product
state in the supplied periodic `SquareIPEPSState`. This is available only for
tiny cells accepted by [`dense_state_finite`](@ref).
"""
function exact_all_down_return_probability_finite(
    psi::SquareIPEPSState;
    max_sites::Integer = 12,
)::Float64
    nsites = _check_tiny_finite_cell(psi, max_sites)
    state = dense_state_finite(psi; max_sites)
    normsq = sum(abs2, state)
    normsq > 0 || throw(ArgumentError("dense finite state has zero norm"))
    all_down_index = 2^nsites
    return abs2(state[all_down_index]) / normsq
end
```

- [ ] **Step 5: Implement ED observable-provenance helpers**

In `src/FinitePXPEEDBenchmark.jl`, extend exports:

```julia
export pxp_ed_boundary_condition, pxp_ed_symmetry_sector, pxp_ed_observable_scope
export pxp_ed_reference_label, pxp_ed_site_density_operator, pxp_ed_region_density_operator
```

Add these methods after `pxp_ed_group_order`:

```julia
"""
    pxp_ed_boundary_condition(basis)

Return the boundary condition represented by a PXP ED basis.
"""
pxp_ed_boundary_condition(::PXPSquareSpaceGroupBasis) = :periodic

"""
    pxp_ed_symmetry_sector(basis)

Return the symmetry sector represented by a PXP ED basis.
"""
pxp_ed_symmetry_sector(basis::PXPSquareSpaceGroupBasis) =
    basis.point_group ? :fully_symmetric_space_group : :translation_symmetric

"""
    pxp_ed_observable_scope(basis)

Return the observable scope supported by the current PBC symmetry-reduced ED
basis. The value is global because local and central-region observables do not
preserve the selected symmetric sector.
"""
pxp_ed_observable_scope(::PXPSquareSpaceGroupBasis) = :pbc_global_site_average

"""
    pxp_ed_reference_label(basis)

Return a stable machine-readable label for the ED reference observable.
"""
pxp_ed_reference_label(::PXPSquareSpaceGroupBasis) = "finite_pbc_global_density"

"""
    pxp_ed_site_density_operator(basis, site)

Construct a site-density operator when the supplied basis supports local
observables. The current symmetry-reduced PBC basis rejects this request because
it would be a projected group average, not a site observable.
"""
function pxp_ed_site_density_operator(::PXPSquareSpaceGroupBasis, site::Integer)
    site >= 1 || throw(ArgumentError("site must be positive"))
    throw(ArgumentError("site density is not available in the symmetry-reduced PBC ED basis"))
end

"""
    pxp_ed_region_density_operator(basis, sites)

Construct a region-density operator when the supplied basis supports local
regions. The current symmetry-reduced PBC basis rejects this request because
there is no central region in a fully symmetric periodic basis.
"""
function pxp_ed_region_density_operator(::PXPSquareSpaceGroupBasis, sites)
    isempty(collect(sites)) && throw(ArgumentError("region sites must be nonempty"))
    throw(ArgumentError("region density is not available in the symmetry-reduced PBC ED basis"))
end
```

- [ ] **Step 6: Re-export the new helpers**

In `src/SquarePXPDynamics.jl`, add `exact_all_down_return_probability_finite` to the `using .FiniteIPEPSObservables` block and export list.

Add ED provenance helpers to the `using .FinitePXPEEDBenchmark` block and export list:

```julia
    pxp_ed_boundary_condition,
    pxp_ed_symmetry_sector,
    pxp_ed_observable_scope,
    pxp_ed_reference_label,
    pxp_ed_site_density_operator,
    pxp_ed_region_density_operator,
```

- [ ] **Step 7: Run focused tests**

Run:

```bash
julia --project=. test/runtests.jl test_finite_ipeps_observables.jl test_pxp_ed_benchmark.jl test_pxp_larger_d_ed_benchmark.jl test_public_docs.jl
```

Expected: PASS.

- [ ] **Step 8: Commit**

```bash
git add src/FiniteIPEPSObservables.jl src/FinitePXPEEDBenchmark.jl src/SquarePXPDynamics.jl test/runtests.jl test/test_pxp_larger_d_ed_benchmark.jl
git commit -m "feat: label PXP ED observable provenance"
```

---

### Task 2: Reusable Larger-D Benchmark Harness

**Files:**
- Modify: `src/PXPValidation.jl`
- Modify: `src/SquarePXPDynamics.jl`
- Test: `test/test_pxp_larger_d_ed_benchmark.jl`

- [ ] **Step 1: Write failing config and D-sweep tests**

Append to `test/test_pxp_larger_d_ed_benchmark.jl`:

```julia
@testset "larger-D PXP benchmark config validates controls" begin
    config = PXPLargerDBenchmarkConfig(;
        n_values = [3],
        total_time = 0.02,
        dt_values = [0.02],
        D_values = [1, 2, 3],
        cutoff_values = [1e-12],
        exact_finite_observables = true,
        exact_finite_max_sites = 9,
    )

    @test config.n_values == [3]
    @test config.D_values == [1, 2, 3]
    @test config.observable_mode === :auto
    @test config.ed_mode === :symmetric_pbc
    @test_throws ArgumentError PXPLargerDBenchmarkConfig(; n_values = Int[])
    @test_throws ArgumentError PXPLargerDBenchmarkConfig(; D_values = Int[])
    @test_throws ArgumentError PXPLargerDBenchmarkConfig(; observable_mode = :central_region)
    @test_throws ArgumentError PXPLargerDBenchmarkConfig(; ed_mode = :open_boundary)
end

@testset "larger-D PXP benchmark separates exact finite and simple diagnostics" begin
    config = PXPLargerDBenchmarkConfig(;
        n_values = [3],
        total_time = 0.02,
        dt_values = [0.02],
        D_values = [1, 2, 3],
        cutoff_values = [1e-12],
        exact_finite_observables = true,
        exact_finite_max_sites = 9,
    )

    report = run_pxp_larger_d_benchmark(config)

    @test length(report.runs) == 3
    @test all(run -> run.summary.observable_mode === :exact_finite, report.runs)
    @test all(run -> run.summary.ed_observable_scope === :pbc_global_site_average, report.runs)
    @test all(run -> run.summary.density_error_exact_finite !== nothing, report.runs)
    @test all(run -> run.summary.return_probability_error !== nothing, report.runs)
    @test all(run -> run.summary.density_error_simple !== nothing, report.runs)
    @test all(run -> run.summary.ed_runtime_seconds >= 0, report.runs)
    @test all(run -> run.summary.ipeps_runtime_seconds >= 0, report.runs)
    @test all(run -> run.summary.max_truncerr >= 0, report.runs)
    @test all(run -> run.summary.log_norm_delta_abs >= 0, report.runs)
    @test all(run -> run.summary.reversibility_density_drift >= 0, report.runs)

    d2 = only(run for run in report.runs if run.summary.D == 2)
    @test abs(d2.summary.density_error_exact_finite) < 1e-6
    @test abs(d2.summary.density_error_simple) > 1e-4
end
```

- [ ] **Step 2: Run the test and verify RED**

Run:

```bash
julia --project=. test/runtests.jl test_pxp_larger_d_ed_benchmark.jl
```

Expected: FAIL with:

```text
UndefVarError: `PXPLargerDBenchmarkConfig` not defined
```

- [ ] **Step 3: Add public M3 types**

In `src/PXPValidation.jl`, add these exports near the audit exports:

```julia
export PXPLargerDBenchmarkConfig, PXPLargerDBenchmarkSummary
export PXPLargerDBenchmarkRun, PXPLargerDBenchmarkReport
export run_pxp_larger_d_benchmark
```

Add these type definitions after `PXPAuditReport`:

```julia
"""
    PXPLargerDBenchmarkConfig(; kwargs...)

Controls the M3 larger-D PXP dynamics benchmark. `ed_mode = :symmetric_pbc`
uses the current finite PBC ED path, and `observable_mode = :auto` selects
`:exact_finite` for tiny cells with exact finite observables enabled and
`:symmetric_pbc_ed_global` otherwise.
"""
struct PXPLargerDBenchmarkConfig
    n_values::Vector{Int}
    total_time::Float64
    dt_values::Vector{Float64}
    D_values::Vector{Int}
    cutoff_values::Vector{Float64}
    measure_every::Int
    order::Int
    schedule::Symbol
    initial_state::Symbol
    point_group::Bool
    use_sparse::Bool
    ed_tol::Float64
    ed_m_init::Int
    ed_m_max::Int
    ed_extend_step::Int
    ed_mode::Symbol
    observable_mode::Symbol
    chi_values::Vector{Int}
    ctm_tol::Float64
    ctm_maxiter::Int
    ctm_verbosity::Int
    ctm_seed::Union{Nothing,Int}
    exact_finite_observables::Bool
    exact_finite_max_sites::Int
end
```

Add the inner constructor with the same validators already used by `PXPAuditConfig`:

```julia
function PXPLargerDBenchmarkConfig(;
    n_values = [3],
    total_time::Real = 0.02,
    dt_values = [0.02],
    D_values = [1, 2, 3, 4],
    cutoff_values = [1e-12],
    measure_every::Integer = 1,
    order::Integer = 1,
    schedule::Symbol = :serial,
    initial_state::Symbol = :down,
    point_group::Bool = true,
    use_sparse::Bool = true,
    ed_tol::Real = 1e-10,
    ed_m_init::Integer = 30,
    ed_m_max::Integer = 60,
    ed_extend_step::Integer = 10,
    ed_mode::Symbol = :symmetric_pbc,
    observable_mode::Symbol = :auto,
    chi_values = Int[],
    ctm_tol::Real = 1e-8,
    ctm_maxiter::Integer = 100,
    ctm_verbosity::Integer = 0,
    ctm_seed = 0,
    exact_finite_observables::Bool = false,
    exact_finite_max_sites::Integer = 12,
)
    n_grid = [_positive_int(n, "n_values entry") for n in n_values]
    dt_grid = [_finite_positive(dt, "dt_values entry") for dt in dt_values]
    D_grid = [_positive_int(D, "D_values entry") for D in D_values]
    cutoff_grid = [_finite_positive(cutoff, "cutoff_values entry") for cutoff in cutoff_values]
    chi_grid = [_positive_int(chi, "chi_values entry") for chi in chi_values]
    isempty(n_grid) && throw(ArgumentError("n_values must be nonempty"))
    isempty(dt_grid) && throw(ArgumentError("dt_values must be nonempty"))
    isempty(D_grid) && throw(ArgumentError("D_values must be nonempty"))
    isempty(cutoff_grid) && throw(ArgumentError("cutoff_values must be nonempty"))
    ed_mode === :symmetric_pbc ||
        throw(ArgumentError("ed_mode must be :symmetric_pbc for M3"))
    observable_mode in (:auto, :exact_finite, :symmetric_pbc_ed_global, :ctm_trusted) ||
        throw(ArgumentError("observable_mode must be :auto, :exact_finite, :symmetric_pbc_ed_global, or :ctm_trusted"))
    initial_state in (:down, :all_down) ||
        throw(ArgumentError("initial_state must be :down or :all_down"))

    cadence = _positive_int(measure_every, "measure_every")
    ord = _positive_int(order, "order")
    ord in (1, 2) || throw(ArgumentError("order must be 1 or 2"))
    schedule in (:serial, :five_color) ||
        throw(ArgumentError("schedule must be :serial or :five_color"))
    total = _finite_nonnegative(total_time, "total_time")
    tol = _finite_positive(ed_tol, "ed_tol")
    m_init = _positive_int(ed_m_init, "ed_m_init")
    m_max = _positive_int(ed_m_max, "ed_m_max")
    m_max >= m_init || throw(ArgumentError("ed_m_max must be at least ed_m_init"))
    extend = _positive_int(ed_extend_step, "ed_extend_step")
    ctm_tol_f = _finite_positive(ctm_tol, "ctm_tol")
    ctm_maxiter_i = _positive_int(ctm_maxiter, "ctm_maxiter")
    ctm_verbosity_i = _nonnegative_int(ctm_verbosity, "ctm_verbosity")
    seed = ctm_seed === nothing ? nothing : _nonnegative_int(ctm_seed, "ctm_seed")
    exact_limit = _positive_int(exact_finite_max_sites, "exact_finite_max_sites")

    return PXPLargerDBenchmarkConfig(
        n_grid,
        total,
        dt_grid,
        D_grid,
        cutoff_grid,
        cadence,
        ord,
        schedule,
        initial_state,
        point_group,
        use_sparse,
        tol,
        m_init,
        m_max,
        extend,
        ed_mode,
        observable_mode,
        chi_grid,
        ctm_tol_f,
        ctm_maxiter_i,
        ctm_verbosity_i,
        seed,
        exact_finite_observables,
        exact_limit,
    )
end
```

Add summary/run/report structs:

```julia
"""
    PXPLargerDBenchmarkSummary

Flat per-run M3 benchmark row. ED fields describe the finite PBC symmetric ED
reference; exact finite fields are nullable and are populated only for tiny
iPEPS cells where exact finite contraction was enabled.
"""
struct PXPLargerDBenchmarkSummary
    n::Int
    D::Int
    dt::Float64
    cutoff::Float64
    total_time::Float64
    ed_mode::Symbol
    observable_mode::Symbol
    ed_boundary_condition::Symbol
    ed_symmetry_sector::Symbol
    ed_observable_scope::Symbol
    ed_reference_label::String
    ed_basis_dimension::Int
    ed_constrained_dimension::Int
    ed_group_order::Int
    ed_hamiltonian_nnz::Union{Nothing,Int}
    ed_runtime_seconds::Float64
    ipeps_runtime_seconds::Float64
    reversibility_runtime_seconds::Float64
    density_error_simple::Float64
    density_error_exact_finite::Union{Nothing,Float64}
    density_error_ctm::Union{Nothing,Float64}
    return_probability_error::Union{Nothing,Float64}
    ed_return_probability::Float64
    ed_excitation_density::Float64
    ipeps_simple_density::Float64
    ipeps_exact_finite_density::Union{Nothing,Float64}
    ipeps_ctm_density::Union{Nothing,Float64}
    max_truncerr::Float64
    log_norm_initial::Float64
    log_norm_final::Float64
    log_norm_delta_abs::Float64
    reversibility_density_drift::Float64
    ctm_trust_status::Symbol
    ctm_trust_reason::Symbol
    notes::Vector{String}
    warnings::Vector{String}
end

"""One M3 benchmark run containing the full validation report and flat summary."""
struct PXPLargerDBenchmarkRun
    validation::PXPValidationReport
    reversibility::PXPReversibilityReport
    summary::PXPLargerDBenchmarkSummary
end

"""Complete M3 benchmark campaign report."""
struct PXPLargerDBenchmarkReport
    config::PXPLargerDBenchmarkConfig
    runs::Vector{PXPLargerDBenchmarkRun}
    metadata::PXPValidationMetadata
end
```

- [ ] **Step 4: Factor validation so ED can be reused across D**

In `src/PXPValidation.jl`, extract the iPEPS half of `validate_pxp_ed_ipeps` into a private helper:

```julia
function _validate_pxp_ipeps_against_ed(
    config::PXPValidationConfig,
    ed_result::PXPEEDBenchmarkResult;
    ctm_params = nothing,
    trust_policy::CTMTrustPolicy = CTMTrustPolicy(),
    ctm_measure = measure_ctm,
)::PXPValidationReport
    psi = _validation_initial_state(config)
    trotter = _validation_trotter(config)
    samples = PXPIPEPSSample[]
    last_time = 0.0

    for ed_sample in ed_result.samples
        interval = ed_sample.time - last_time
        evolution = if iszero(interval)
            nothing
        else
            evolve!(psi, interval; params = trotter)
        end
        last_time = ed_sample.time

        simple = measure_simple(psi)
        exact_finite_density = config.exact_finite_observables ?
            exact_density_finite(psi; max_sites = config.exact_finite_max_sites) : nothing
        ctm = _maybe_trusted_ctm(psi, ctm_params, trust_policy, ctm_measure)
        push!(
            samples,
            PXPIPEPSSample(
                ed_sample.step,
                ed_sample.time,
                simple,
                evolution,
                ctm,
                log_norm(psi),
                exact_finite_density,
            ),
        )
    end

    comparisons = [
        _comparison(ed_sample, ipeps_sample) for
        (ed_sample, ipeps_sample) in zip(ed_result.samples, samples)
    ]
    return PXPValidationReport(config, ed_result, samples, comparisons, _validation_metadata())
end
```

Then replace the body of `validate_pxp_ed_ipeps` with:

```julia
ed_result = run_pxp_ed_benchmark(_validation_ed_config(config))
return _validate_pxp_ipeps_against_ed(
    config,
    ed_result;
    ctm_params,
    trust_policy,
    ctm_measure,
)
```

- [ ] **Step 5: Implement M3 run helpers**

Still in `src/PXPValidation.jl`, add helper functions:

First extend the existing `FiniteIPEPSObservables` import:

```julia
using ..FiniteIPEPSObservables: exact_density_finite, exact_all_down_return_probability_finite
```

Then add:

```julia
function _larger_d_ctm_params(config::PXPLargerDBenchmarkConfig)
    isempty(config.chi_values) && return nothing
    if config.ctm_seed === nothing
        return Tuple(
            PEPSKitCTMRGParams(chi, config.ctm_tol, config.ctm_maxiter, config.ctm_verbosity)
            for chi in config.chi_values
        )
    else
        return Tuple(
            PEPSKitCTMRGParams(
                chi,
                config.ctm_tol,
                config.ctm_maxiter,
                config.ctm_verbosity;
                seed = config.ctm_seed,
            ) for chi in config.chi_values
        )
    end
end

function _larger_d_ed_config(config::PXPLargerDBenchmarkConfig, n::Int, dt::Float64)
    return PXPEEDBenchmarkConfig(
        n;
        total_time = config.total_time,
        dt,
        measure_every = config.measure_every,
        initial_state = config.initial_state,
        point_group = config.point_group,
        use_sparse = config.use_sparse,
        tol = config.ed_tol,
        m_init = config.ed_m_init,
        m_max = config.ed_m_max,
        extend_step = config.ed_extend_step,
    )
end

function _larger_d_validation_config(
    config::PXPLargerDBenchmarkConfig,
    n::Int,
    dt::Float64,
    D::Int,
    cutoff::Float64,
)
    use_exact = config.exact_finite_observables && n^2 <= config.exact_finite_max_sites
    return PXPValidationConfig(
        n;
        total_time = config.total_time,
        dt,
        measure_every = config.measure_every,
        initial_state = config.initial_state,
        point_group = config.point_group,
        use_sparse = config.use_sparse,
        ed_tol = config.ed_tol,
        ed_m_init = config.ed_m_init,
        ed_m_max = config.ed_m_max,
        ed_extend_step = config.ed_extend_step,
        order = config.order,
        maxdim = D,
        cutoff,
        schedule = config.schedule,
        exact_finite_observables = use_exact,
        exact_finite_max_sites = config.exact_finite_max_sites,
    )
end
```

Add mode and summary helpers:

```julia
function _larger_d_observable_mode(config::PXPLargerDBenchmarkConfig, run_config::PXPValidationConfig)
    config.observable_mode !== :auto && return config.observable_mode
    run_config.exact_finite_observables && return :exact_finite
    !isempty(config.chi_values) && return :ctm_trusted
    return :symmetric_pbc_ed_global
end

function _last_evolution_max_truncerr(samples)
    evolutions = [sample.evolution for sample in samples if sample.evolution !== nothing]
    return _maximum_or_zero([e.max_truncerr for e in evolutions])
end

function _exact_return_probability_or_nothing(sample::PXPIPEPSSample, config::PXPValidationConfig)
    config.exact_finite_observables || return nothing
    psi = _validation_initial_state(config)
    evolve!(psi, sample.time; params = _validation_trotter(config))
    return exact_all_down_return_probability_finite(psi; max_sites = config.exact_finite_max_sites)
end
```

If recomputing the final exact return probability from scratch is too slow during implementation, use it only for `n^2 <= exact_finite_max_sites`; for 3x3 this is acceptable. Keep it out of 5x5 and 7x7 paths.

- [ ] **Step 6: Implement `run_pxp_larger_d_benchmark`**

Add:

```julia
"""
    run_pxp_larger_d_benchmark(config = PXPLargerDBenchmarkConfig(); kwargs...)

Run the M3 larger-D PXP ED benchmark campaign. ED is run once for each
`(n, dt)` pair and reused for every requested `D` and cutoff at that pair.
"""
function run_pxp_larger_d_benchmark(
    config::PXPLargerDBenchmarkConfig = PXPLargerDBenchmarkConfig();
    trust_policy::CTMTrustPolicy = CTMTrustPolicy(),
    ctm_measure = measure_ctm,
)::PXPLargerDBenchmarkReport
    ctm_params = _larger_d_ctm_params(config)
    runs = PXPLargerDBenchmarkRun[]

    for n in config.n_values, dt in config.dt_values
        ed_result = nothing
        ed_seconds = @elapsed begin
            ed_result = run_pxp_ed_benchmark(_larger_d_ed_config(config, n, dt))
        end
        for D in config.D_values, cutoff in config.cutoff_values
            run_config = _larger_d_validation_config(config, n, dt, D, cutoff)
            validation = nothing
            ipeps_seconds = @elapsed begin
                validation = _validate_pxp_ipeps_against_ed(
                    run_config,
                    ed_result;
                    ctm_params,
                    trust_policy,
                    ctm_measure,
                )
            end
            psi = _validation_initial_state(run_config)
            reversibility = nothing
            reversibility_seconds = @elapsed begin
                reversibility = validate_pxp_reversibility(
                    psi,
                    run_config.total_time;
                    params = _validation_trotter(run_config),
                )
            end
            push!(
                runs,
                PXPLargerDBenchmarkRun(
                    validation,
                    reversibility,
                    _larger_d_summary(
                        config,
                        run_config,
                        validation,
                        reversibility,
                        ed_seconds,
                        ipeps_seconds,
                        reversibility_seconds,
                    ),
                ),
            )
        end
    end

    return PXPLargerDBenchmarkReport(config, runs, _validation_metadata())
end
```

Implement `_larger_d_summary` by taking the final comparison and final ED sample:

```julia
function _larger_d_summary(
    config::PXPLargerDBenchmarkConfig,
    run_config::PXPValidationConfig,
    validation::PXPValidationReport,
    reversibility::PXPReversibilityReport,
    ed_seconds::Float64,
    ipeps_seconds::Float64,
    reversibility_seconds::Float64,
)::PXPLargerDBenchmarkSummary
    final_comparison = validation.comparisons[end]
    final_sample = validation.ipeps_samples[end]
    final_ed = validation.ed_result.samples[end]
    log_norms = [sample.log_norm for sample in validation.ipeps_samples]
    trust_status, trust_reason = _audit_trust_status(validation.ipeps_samples)
    mode = _larger_d_observable_mode(config, run_config)
    exact_return = final_sample.exact_finite_density === nothing ?
        nothing : _exact_return_probability_or_nothing(final_sample, run_config)

    warnings = String[]
    run_config.maxdim > 1 && push!(
        warnings,
        "density_simple is a diagnostic for D>1 loopy PEPS, not exact finite truth",
    )
    run_config.exact_finite_observables || push!(
        warnings,
        "exact finite iPEPS observables were not available for this run",
    )
    push!(warnings, "symmetric PBC ED observables are global site averages")

    return PXPLargerDBenchmarkSummary(
        run_config.n,
        run_config.maxdim,
        run_config.dt,
        run_config.cutoff,
        run_config.total_time,
        config.ed_mode,
        mode,
        :periodic,
        validation.ed_result.point_group ? :fully_symmetric_space_group : :translation_symmetric,
        :pbc_global_site_average,
        "finite_pbc_global_density",
        validation.ed_result.basis_dimension,
        validation.ed_result.constrained_dimension,
        validation.ed_result.group_order,
        validation.ed_result.hamiltonian_nnz,
        ed_seconds,
        ipeps_seconds,
        reversibility_seconds,
        final_comparison.density_error_simple,
        final_comparison.density_error_exact_finite,
        final_comparison.density_error_ctm,
        exact_return === nothing ? nothing : exact_return - final_ed.return_probability,
        final_ed.return_probability,
        final_ed.excitation_density,
        final_comparison.ipeps_simple_density,
        final_comparison.ipeps_exact_finite_density,
        final_comparison.ipeps_ctm_density,
        _last_evolution_max_truncerr(validation.ipeps_samples),
        isempty(log_norms) ? 0.0 : first(log_norms),
        isempty(log_norms) ? 0.0 : last(log_norms),
        isempty(log_norms) ? 0.0 : abs(last(log_norms) - first(log_norms)),
        reversibility.density_drift,
        trust_status,
        trust_reason,
        ["M3 larger-D PXP ED benchmark"],
        warnings,
    )
end
```

- [ ] **Step 7: Re-export M3 types/functions**

In `src/SquarePXPDynamics.jl`, add the new `PXPValidation` imports and exports:

```julia
    PXPLargerDBenchmarkConfig,
    PXPLargerDBenchmarkSummary,
    PXPLargerDBenchmarkRun,
    PXPLargerDBenchmarkReport,
    run_pxp_larger_d_benchmark,
```

- [ ] **Step 8: Run focused tests**

Run:

```bash
julia --project=. test/runtests.jl test_pxp_larger_d_ed_benchmark.jl test_pxp_validation.jl test_public_docs.jl
```

Expected: PASS.

- [ ] **Step 9: Commit**

```bash
git add src/PXPValidation.jl src/SquarePXPDynamics.jl test/test_pxp_larger_d_ed_benchmark.jl
git commit -m "feat: add larger-D PXP ED benchmark harness"
```

---

### Task 3: JSON And CSV Report Writers

**Files:**
- Modify: `src/PXPValidation.jl`
- Modify: `src/SquarePXPDynamics.jl`
- Test: `test/test_pxp_larger_d_ed_benchmark.jl`

- [ ] **Step 1: Write failing schema tests**

Append:

```julia
@testset "larger-D PXP benchmark JSON and CSV preserve required schema" begin
    config = PXPLargerDBenchmarkConfig(;
        n_values = [3],
        total_time = 0.02,
        dt_values = [0.02],
        D_values = [1, 2],
        cutoff_values = [1e-12],
        exact_finite_observables = true,
        exact_finite_max_sites = 9,
    )
    report = run_pxp_larger_d_benchmark(config)
    json_path = tempname() * ".json"
    csv_path = tempname() * ".csv"

    @test write_pxp_larger_d_benchmark_json(report, json_path) == json_path
    @test write_pxp_larger_d_benchmark_csv(report, csv_path) == csv_path

    parsed = JSON3.read(read(json_path, String))
    @test parsed.schema_version == 1
    @test parsed.config.ed_mode == "symmetric_pbc"
    @test parsed.config.observable_mode == "auto"
    @test length(parsed.runs) == 2
    @test parsed.runs[1].summary.ed_observable_scope == "pbc_global_site_average"
    @test parsed.runs[1].summary.observable_mode == "exact_finite"
    @test parsed.runs[2].summary.D == 2
    @test !any(k -> occursin("central", lowercase(String(k))), keys(parsed.runs[1].summary))

    csv = split(chomp(read(csv_path, String)), '\n')
    header = split(csv[1], ','; keepempty = true)
    required = [
        "n",
        "D",
        "dt",
        "cutoff",
        "total_time",
        "ed_basis_dimension",
        "ed_constrained_dimension",
        "ed_group_order",
        "ed_runtime_seconds",
        "ipeps_runtime_seconds",
        "observable_mode",
        "density_error_simple",
        "density_error_exact_finite",
        "density_error_ctm",
        "return_probability_error",
        "max_truncerr",
        "log_norm_initial",
        "log_norm_final",
        "log_norm_delta_abs",
        "reversibility_density_drift",
        "ctm_trust_status",
        "ctm_trust_reason",
        "notes",
        "warnings",
    ]
    @test all(name -> name in header, required)
    @test !any(name -> occursin("central", lowercase(name)), header)
    @test length(csv) == 3
    @test parse(Float64, _csv_cell(csv, "density_error_exact_finite")) < 1e-5
end
```

- [ ] **Step 2: Run and verify RED**

Run:

```bash
julia --project=. test/runtests.jl test_pxp_larger_d_ed_benchmark.jl
```

Expected: FAIL with:

```text
UndefVarError: `write_pxp_larger_d_benchmark_json` not defined
```

- [ ] **Step 3: Add serializers and writers**

In `src/PXPValidation.jl`, export:

```julia
export write_pxp_larger_d_benchmark_json, write_pxp_larger_d_benchmark_csv
```

Add serializers:

```julia
function _larger_d_config_data(config::PXPLargerDBenchmarkConfig)
    return (;
        n_values = config.n_values,
        total_time = config.total_time,
        dt_values = config.dt_values,
        D_values = config.D_values,
        cutoff_values = config.cutoff_values,
        measure_every = config.measure_every,
        order = config.order,
        schedule = String(config.schedule),
        initial_state = String(config.initial_state),
        point_group = config.point_group,
        use_sparse = config.use_sparse,
        ed_tol = config.ed_tol,
        ed_m_init = config.ed_m_init,
        ed_m_max = config.ed_m_max,
        ed_extend_step = config.ed_extend_step,
        ed_mode = String(config.ed_mode),
        observable_mode = String(config.observable_mode),
        chi_values = config.chi_values,
        ctm_tol = config.ctm_tol,
        ctm_maxiter = config.ctm_maxiter,
        ctm_verbosity = config.ctm_verbosity,
        ctm_seed = config.ctm_seed,
        exact_finite_observables = config.exact_finite_observables,
        exact_finite_max_sites = config.exact_finite_max_sites,
    )
end

function _larger_d_summary_data(summary::PXPLargerDBenchmarkSummary)
    return (;
        n = summary.n,
        D = summary.D,
        dt = summary.dt,
        cutoff = summary.cutoff,
        total_time = summary.total_time,
        ed_mode = String(summary.ed_mode),
        observable_mode = String(summary.observable_mode),
        ed_boundary_condition = String(summary.ed_boundary_condition),
        ed_symmetry_sector = String(summary.ed_symmetry_sector),
        ed_observable_scope = String(summary.ed_observable_scope),
        ed_reference_label = summary.ed_reference_label,
        ed_basis_dimension = summary.ed_basis_dimension,
        ed_constrained_dimension = summary.ed_constrained_dimension,
        ed_group_order = summary.ed_group_order,
        ed_hamiltonian_nnz = summary.ed_hamiltonian_nnz,
        ed_runtime_seconds = summary.ed_runtime_seconds,
        ipeps_runtime_seconds = summary.ipeps_runtime_seconds,
        reversibility_runtime_seconds = summary.reversibility_runtime_seconds,
        density_error_simple = summary.density_error_simple,
        density_error_exact_finite = summary.density_error_exact_finite,
        density_error_ctm = summary.density_error_ctm,
        return_probability_error = summary.return_probability_error,
        ed_return_probability = summary.ed_return_probability,
        ed_excitation_density = summary.ed_excitation_density,
        ipeps_simple_density = summary.ipeps_simple_density,
        ipeps_exact_finite_density = summary.ipeps_exact_finite_density,
        ipeps_ctm_density = summary.ipeps_ctm_density,
        max_truncerr = summary.max_truncerr,
        log_norm_initial = summary.log_norm_initial,
        log_norm_final = summary.log_norm_final,
        log_norm_delta_abs = summary.log_norm_delta_abs,
        reversibility_density_drift = summary.reversibility_density_drift,
        ctm_trust_status = String(summary.ctm_trust_status),
        ctm_trust_reason = String(summary.ctm_trust_reason),
        notes = summary.notes,
        warnings = summary.warnings,
    )
end

function _larger_d_run_data(run::PXPLargerDBenchmarkRun)
    return (;
        summary = _larger_d_summary_data(run.summary),
        validation = _report_data(run.validation),
        reversibility = _reversibility_data(run.reversibility),
    )
end

function _larger_d_report_data(report::PXPLargerDBenchmarkReport)
    return (;
        schema_version = 1,
        config = _larger_d_config_data(report.config),
        metadata = _metadata_data(report.metadata),
        runs = _larger_d_run_data.(report.runs),
    )
end

"""
    write_pxp_larger_d_benchmark_json(report, path)

Write a nested M3 larger-D PXP ED benchmark report as JSON.
"""
function write_pxp_larger_d_benchmark_json(
    report::PXPLargerDBenchmarkReport,
    path::AbstractString,
)
    open(path, "w") do io
        JSON3.write(io, _larger_d_report_data(report))
        write(io, '\n')
    end
    return path
end
```

Add CSV header and writer:

```julia
const PXP_LARGER_D_CSV_HEADER = [
    "n",
    "D",
    "dt",
    "cutoff",
    "total_time",
    "ed_mode",
    "observable_mode",
    "ed_boundary_condition",
    "ed_symmetry_sector",
    "ed_observable_scope",
    "ed_reference_label",
    "ed_basis_dimension",
    "ed_constrained_dimension",
    "ed_group_order",
    "ed_hamiltonian_nnz",
    "ed_runtime_seconds",
    "ipeps_runtime_seconds",
    "reversibility_runtime_seconds",
    "density_error_simple",
    "density_error_exact_finite",
    "density_error_ctm",
    "return_probability_error",
    "ed_return_probability",
    "ed_excitation_density",
    "ipeps_simple_density",
    "ipeps_exact_finite_density",
    "ipeps_ctm_density",
    "max_truncerr",
    "log_norm_initial",
    "log_norm_final",
    "log_norm_delta_abs",
    "reversibility_density_drift",
    "ctm_trust_status",
    "ctm_trust_reason",
    "notes",
    "warnings",
]

function _larger_d_csv_cell(xs::Vector{String})
    return _audit_csv_cell(join(xs, ";"))
end

_larger_d_csv_cell(x) = _audit_csv_cell(x)

function _larger_d_csv_row(summary::PXPLargerDBenchmarkSummary)
    values = (
        summary.n,
        summary.D,
        summary.dt,
        summary.cutoff,
        summary.total_time,
        summary.ed_mode,
        summary.observable_mode,
        summary.ed_boundary_condition,
        summary.ed_symmetry_sector,
        summary.ed_observable_scope,
        summary.ed_reference_label,
        summary.ed_basis_dimension,
        summary.ed_constrained_dimension,
        summary.ed_group_order,
        summary.ed_hamiltonian_nnz,
        summary.ed_runtime_seconds,
        summary.ipeps_runtime_seconds,
        summary.reversibility_runtime_seconds,
        summary.density_error_simple,
        summary.density_error_exact_finite,
        summary.density_error_ctm,
        summary.return_probability_error,
        summary.ed_return_probability,
        summary.ed_excitation_density,
        summary.ipeps_simple_density,
        summary.ipeps_exact_finite_density,
        summary.ipeps_ctm_density,
        summary.max_truncerr,
        summary.log_norm_initial,
        summary.log_norm_final,
        summary.log_norm_delta_abs,
        summary.reversibility_density_drift,
        summary.ctm_trust_status,
        summary.ctm_trust_reason,
        summary.notes,
        summary.warnings,
    )
    return join(_larger_d_csv_cell.(values), ",")
end

"""
    write_pxp_larger_d_benchmark_csv(report, path)

Write the flat M3 benchmark summary rows as CSV.
"""
function write_pxp_larger_d_benchmark_csv(
    report::PXPLargerDBenchmarkReport,
    path::AbstractString,
)
    open(path, "w") do io
        println(io, join(PXP_LARGER_D_CSV_HEADER, ","))
        for run in report.runs
            println(io, _larger_d_csv_row(run.summary))
        end
    end
    return path
end
```

- [ ] **Step 4: Re-export writers**

In `src/SquarePXPDynamics.jl`, add both writer names to the `using .PXPValidation` block and export list.

- [ ] **Step 5: Run focused tests**

Run:

```bash
julia --project=. test/runtests.jl test_pxp_larger_d_ed_benchmark.jl test_pxp_validation.jl test_public_docs.jl
```

Expected: PASS.

- [ ] **Step 6: Commit**

```bash
git add src/PXPValidation.jl src/SquarePXPDynamics.jl test/test_pxp_larger_d_ed_benchmark.jl
git commit -m "feat: write larger-D PXP ED benchmark reports"
```

---

### Task 4: Env-Driven Benchmark Script

**Files:**
- Create: `scripts/pxp_larger_d_ed_benchmark.jl`
- Test: `test/test_pxp_larger_d_ed_benchmark.jl`

- [ ] **Step 1: Write failing script smoke test**

Append this test:

```julia
@testset "larger-D PXP benchmark script exists" begin
    script = joinpath(dirname(@__DIR__), "scripts", "pxp_larger_d_ed_benchmark.jl")
    @test isfile(script)
    text = read(script, String)
    @test occursin("SQUAREPXP_LARGERD_N", text)
    @test occursin("SQUAREPXP_LARGERD_EXACT_FINITE", text)
    @test occursin("write_pxp_larger_d_benchmark_json", text)
    @test occursin("write_pxp_larger_d_benchmark_csv", text)
end
```

- [ ] **Step 2: Run and verify RED**

Run:

```bash
julia --project=. test/runtests.jl test_pxp_larger_d_ed_benchmark.jl
```

Expected: FAIL because `scripts/pxp_larger_d_ed_benchmark.jl` does not exist.

- [ ] **Step 3: Create the script**

Create `scripts/pxp_larger_d_ed_benchmark.jl`:

```julia
#!/usr/bin/env julia

using Pkg

project_root = dirname(@__DIR__)
Pkg.activate(project_root; io = devnull)

using SquarePXPDynamics

function _env_value(name::String, default::AbstractString)
    value = get(ENV, name, "")
    return isempty(value) ? String(default) : value
end

function _env_bool(name::String, default::Bool)
    value = lowercase(strip(_env_value(name, string(default))))
    value in ("1", "true", "yes", "on") && return true
    value in ("0", "false", "no", "off") && return false
    throw(ArgumentError("$name must be one of 1,true,yes,on,0,false,no,off"))
end

function _env_int(name::String, default::Int)
    return parse(Int, _env_value(name, string(default)))
end

function _env_float(name::String, default::Float64)
    return parse(Float64, _env_value(name, string(default)))
end

function _env_symbol(name::String, default::Symbol)
    return Symbol(_env_value(name, String(default)))
end

function _env_list(::Type{T}, name::String, default::String) where {T}
    raw = strip(_env_value(name, default))
    isempty(raw) && return T[]
    return parse.(T, split(raw, ","))
end

function _env_optional_int(name::String, default::String)
    raw = strip(_env_value(name, default))
    return isempty(raw) || lowercase(raw) == "nothing" ? nothing : parse(Int, raw)
end

config = PXPLargerDBenchmarkConfig(;
    n_values = _env_list(Int, "SQUAREPXP_LARGERD_N", "3"),
    total_time = _env_float("SQUAREPXP_LARGERD_TOTAL_TIME", 0.02),
    dt_values = _env_list(Float64, "SQUAREPXP_LARGERD_DT", "0.02"),
    D_values = _env_list(Int, "SQUAREPXP_LARGERD_D", "1,2,3,4"),
    cutoff_values = _env_list(Float64, "SQUAREPXP_LARGERD_CUTOFF", "1e-12"),
    measure_every = _env_int("SQUAREPXP_LARGERD_MEASURE_EVERY", 1),
    order = _env_int("SQUAREPXP_LARGERD_ORDER", 1),
    schedule = _env_symbol("SQUAREPXP_LARGERD_SCHEDULE", :serial),
    initial_state = _env_symbol("SQUAREPXP_LARGERD_INITIAL_STATE", :down),
    point_group = _env_bool("SQUAREPXP_LARGERD_POINT_GROUP", true),
    use_sparse = _env_bool("SQUAREPXP_LARGERD_USE_SPARSE", true),
    ed_tol = _env_float("SQUAREPXP_LARGERD_ED_TOL", 1e-10),
    ed_m_init = _env_int("SQUAREPXP_LARGERD_ED_M_INIT", 30),
    ed_m_max = _env_int("SQUAREPXP_LARGERD_ED_M_MAX", 60),
    ed_extend_step = _env_int("SQUAREPXP_LARGERD_ED_EXTEND_STEP", 10),
    ed_mode = _env_symbol("SQUAREPXP_LARGERD_ED_MODE", :symmetric_pbc),
    observable_mode = _env_symbol("SQUAREPXP_LARGERD_OBSERVABLE_MODE", :auto),
    chi_values = _env_list(Int, "SQUAREPXP_LARGERD_CHI", ""),
    ctm_tol = _env_float("SQUAREPXP_LARGERD_CTM_TOL", 1e-8),
    ctm_maxiter = _env_int("SQUAREPXP_LARGERD_CTM_MAXITER", 100),
    ctm_verbosity = _env_int("SQUAREPXP_LARGERD_CTM_VERBOSITY", 0),
    ctm_seed = _env_optional_int("SQUAREPXP_LARGERD_CTM_SEED", "0"),
    exact_finite_observables = _env_bool("SQUAREPXP_LARGERD_EXACT_FINITE", false),
    exact_finite_max_sites = _env_int("SQUAREPXP_LARGERD_EXACT_FINITE_MAX_SITES", 12),
)

json_out = _env_value(
    "SQUAREPXP_LARGERD_JSON",
    joinpath(project_root, "artifacts", "pxp_larger_d_ed_benchmark.json"),
)
csv_out = _env_value(
    "SQUAREPXP_LARGERD_CSV",
    joinpath(project_root, "artifacts", "pxp_larger_d_ed_benchmark.csv"),
)

mkpath(dirname(json_out))
mkpath(dirname(csv_out))

report = run_pxp_larger_d_benchmark(config)
write_pxp_larger_d_benchmark_json(report, json_out)
write_pxp_larger_d_benchmark_csv(report, csv_out)

println(json_out)
println(csv_out)
```

- [ ] **Step 4: Run focused script test**

Run:

```bash
julia --project=. test/runtests.jl test_pxp_larger_d_ed_benchmark.jl
```

Expected: PASS.

- [ ] **Step 5: Manual smoke command**

Run:

```bash
SQUAREPXP_LARGERD_N=3 \
SQUAREPXP_LARGERD_DT=0.02 \
SQUAREPXP_LARGERD_D=1,2 \
SQUAREPXP_LARGERD_CUTOFF=1e-12 \
SQUAREPXP_LARGERD_TOTAL_TIME=0.02 \
SQUAREPXP_LARGERD_EXACT_FINITE=true \
SQUAREPXP_LARGERD_EXACT_FINITE_MAX_SITES=9 \
SQUAREPXP_LARGERD_JSON=artifacts/m3-smoke.json \
SQUAREPXP_LARGERD_CSV=artifacts/m3-smoke.csv \
julia --project=. scripts/pxp_larger_d_ed_benchmark.jl
```

Expected output:

```text
artifacts/m3-smoke.json
artifacts/m3-smoke.csv
```

Inspect:

```bash
head -n 2 artifacts/m3-smoke.csv
```

Expected: header includes `observable_mode`, `density_error_exact_finite`, `return_probability_error`, and `ed_observable_scope`.

- [ ] **Step 6: Commit**

```bash
git add scripts/pxp_larger_d_ed_benchmark.jl test/test_pxp_larger_d_ed_benchmark.jl
git commit -m "feat: add larger-D PXP ED benchmark script"
```

---

### Task 5: ED Capacity Smoke And Extended Gates

**Files:**
- Modify: `test/test_pxp_larger_d_ed_benchmark.jl`
- Modify: `README.md`

- [ ] **Step 1: Add fast 5x5 capacity semantics test**

Append:

```julia
@testset "PXP symmetric PBC ED capacity includes 5x5 smoke metadata" begin
    config = PXPLargerDBenchmarkConfig(;
        n_values = [5],
        total_time = 0.0,
        dt_values = [0.01],
        D_values = [1],
        cutoff_values = [1e-12],
        exact_finite_observables = true,
        exact_finite_max_sites = 9,
    )

    report = run_pxp_larger_d_benchmark(config)
    summary = only(report.runs).summary

    @test summary.n == 5
    @test summary.ed_basis_dimension == 188
    @test summary.ed_constrained_dimension == 25_531
    @test summary.ed_group_order == 8 * 5^2
    @test summary.observable_mode === :symmetric_pbc_ed_global
    @test summary.density_error_exact_finite === nothing
    @test summary.return_probability_error === nothing
    @test any(w -> occursin("exact finite", w), summary.warnings)
end
```

This run has `total_time = 0.0`, so it builds the 5x5 ED reference and one D=1 initial iPEPS sample without long evolution.

- [ ] **Step 2: Add extended non-default smoke**

Append:

```julia
if RUN_EXTENDED_PXP_ED_TESTS
    @testset "extended PXP larger-D benchmark smoke" begin
        config = PXPLargerDBenchmarkConfig(;
            n_values = [3],
            total_time = 0.02,
            dt_values = [0.02],
            D_values = [1, 2, 3, 4],
            cutoff_values = [1e-12],
            exact_finite_observables = true,
            exact_finite_max_sites = 9,
        )
        report = run_pxp_larger_d_benchmark(config)
        @test length(report.runs) == 4
        @test all(run -> run.summary.density_error_exact_finite !== nothing, report.runs)
        @test all(run -> isfinite(run.summary.log_norm_delta_abs), report.runs)
    end
end
```

- [ ] **Step 3: Run default focused tests**

Run:

```bash
julia --project=. test/runtests.jl test_pxp_larger_d_ed_benchmark.jl
```

Expected: PASS without extended test output.

- [ ] **Step 4: Run extended focused tests**

Run:

```bash
SQUAREPXP_EXTENDED_PXP_ED_TESTS=1 julia --project=. test/runtests.jl test_pxp_larger_d_ed_benchmark.jl
```

Expected: PASS with the extended D=1,2,3,4 testset. No 7x7 job runs.

- [ ] **Step 5: Document manual 7x7 boundary command**

In `README.md`, add this command to the M3 section:

```bash
SQUAREPXP_LARGERD_N=7 \
SQUAREPXP_LARGERD_DT=0.01 \
SQUAREPXP_LARGERD_D=1 \
SQUAREPXP_LARGERD_TOTAL_TIME=0.0 \
SQUAREPXP_LARGERD_USE_SPARSE=false \
SQUAREPXP_LARGERD_JSON=artifacts/m3-7x7-capacity.json \
SQUAREPXP_LARGERD_CSV=artifacts/m3-7x7-capacity.csv \
julia --project=. scripts/pxp_larger_d_ed_benchmark.jl
```

Document expected interpretation:

- `7x7` is the largest square PBC size supported by the current UInt64 basis.
- This command is a capacity boundary probe, not a default test.
- If it is too slow or memory-heavy, use the existing `scripts/pxp_ed_7x7_benchmark.jl` with `PXP_ED_USE_SPARSE=false` for ED-only diagnosis and record the runtime boundary in the M3 note.

- [ ] **Step 6: Commit**

```bash
git add test/test_pxp_larger_d_ed_benchmark.jl README.md
git commit -m "test: add M3 PXP ED capacity smoke"
```

---

### Task 6: README, M3 Note, And Decision Log

**Files:**
- Modify: `README.md`
- Create: `docs/superpowers/notes/2026-05-17-m3-larger-d-pxp-ed-benchmark.md`
- Modify: `memory/mid_term/decision_log.md`

- [ ] **Step 1: Update README feature and validation sections**

Add a bullet near the existing PXP validation bullets:

```markdown
- M3 larger-D PXP ED benchmark reports via `run_pxp_larger_d_benchmark`, with
  exact finite 3x3 iPEPS observables when enabled and symmetric PBC ED global
  density/return-probability metadata for larger odd cells (`src/PXPValidation.jl`).
```

Add a new subsection under PXP validation reports:

```markdown
### M3 larger-D PXP ED benchmark

Use `run_pxp_larger_d_benchmark` for larger-D sweeps against finite PBC ED:

```julia
config = PXPLargerDBenchmarkConfig(;
    n_values = [3],
    total_time = 0.02,
    dt_values = [0.02, 0.01],
    D_values = [1, 2, 3, 4],
    cutoff_values = [1e-12],
    exact_finite_observables = true,
    exact_finite_max_sites = 9,
)
report = run_pxp_larger_d_benchmark(config)
write_pxp_larger_d_benchmark_json(report, "artifacts/m3-larger-d.json")
write_pxp_larger_d_benchmark_csv(report, "artifacts/m3-larger-d.csv")
```

or from the shell:

```bash
SQUAREPXP_LARGERD_N=3 \
SQUAREPXP_LARGERD_DT=0.02,0.01 \
SQUAREPXP_LARGERD_D=1,2,3,4 \
SQUAREPXP_LARGERD_CUTOFF=1e-12 \
SQUAREPXP_LARGERD_TOTAL_TIME=0.02 \
SQUAREPXP_LARGERD_EXACT_FINITE=true \
SQUAREPXP_LARGERD_EXACT_FINITE_MAX_SITES=9 \
julia --project=. scripts/pxp_larger_d_ed_benchmark.jl
```

The ED reference is finite periodic and symmetry-reduced. Its density is a
global site average in the selected symmetric sector. It is not a central 3x3
or local-window observable. For D>1, `density_error_simple` remains a simple
environment diagnostic; use `density_error_exact_finite` for 3x3 finite
validation and CTM-trusted fields only when finite-chi trust sweeps were run.
```

- [ ] **Step 2: Document extended test gate**

In the development section, add:

```markdown
Run the optional PXP ED benchmark smoke:

```bash
SQUAREPXP_EXTENDED_PXP_ED_TESTS=1 julia --project=. test/runtests.jl test_pxp_larger_d_ed_benchmark.jl
```

This does not run 7x7. Large 7x7 probes are manual benchmark commands.
```

- [ ] **Step 3: Create the M3 note**

Create `docs/superpowers/notes/2026-05-17-m3-larger-d-pxp-ed-benchmark.md`:

```markdown
# M3 Larger-D PXP ED Benchmark

Date: 2026-05-17

## Contract

The M3 benchmark compares iPEPS dynamics against finite PBC ED through
observable modes that are explicitly labeled in JSON and CSV artifacts.

- `exact_finite`: tiny-cell exact finite contraction of the current iPEPS state.
- `simple_diagnostic`: simple/local environment diagnostic, not exact for D>1.
- `ctm_trusted`: CTM-backed value only when finite-chi trust is attached.
- `symmetric_pbc_ed_global`: finite PBC ED global site-averaged density and ED
  return probability in the selected symmetric sector.

The symmetric PBC ED path does not provide central 3x3 observables. PBC has no
physical center, and the current ED basis is fully symmetry reduced by default.

## Default Fast Campaign

```bash
SQUAREPXP_LARGERD_N=3 \
SQUAREPXP_LARGERD_DT=0.02 \
SQUAREPXP_LARGERD_D=1,2,3,4 \
SQUAREPXP_LARGERD_CUTOFF=1e-12 \
SQUAREPXP_LARGERD_TOTAL_TIME=0.02 \
SQUAREPXP_LARGERD_EXACT_FINITE=true \
SQUAREPXP_LARGERD_EXACT_FINITE_MAX_SITES=9 \
julia --project=. scripts/pxp_larger_d_ed_benchmark.jl
```

## Larger Odd PBC Campaign

Use `n = 5` for the first larger odd benchmark. Use `n = 7` only as a manual
capacity/runtime boundary probe.

## Postponed

Unreduced PBC 5x5 ED can support local operators, but a central region is still
not physically privileged under PBC. Open-boundary ED is the cleaner route for
literal central 3x3 observables and should be planned as a separate milestone.
```

- [ ] **Step 4: Add decision-log entry**

Append to `memory/mid_term/decision_log.md`:

```markdown
## 2026-05-17 - Use Symmetric PBC ED Only For Global M3 Observables

Decision:

M3 larger-D PXP dynamics benchmarks compare against the current symmetry-reduced
finite PBC ED path only through global sector-preserving observables: ED return
probability and global site-averaged excitation density. Exact finite iPEPS
contraction is used for tiny 3x3 validation when enabled. Central-region
observables are not claimed for symmetric PBC ED.

Reason:

The current ED basis is reduced by translations and, by default, the square
point group. A local or central-region operator does not preserve that basis as
a literal local observable; after projection it becomes a group average. PBC
also has no physical center.

Consequences:

5x5 and 7x7 PBC ED benchmarks are scientifically honest global comparisons.
Literal central 3x3 comparisons require a future unreduced/open-boundary ED
path and are outside M3.

Source:

`src/FinitePXPEEDBenchmark.jl`; `src/PXPValidation.jl`;
`docs/superpowers/notes/2026-05-17-m3-larger-d-pxp-ed-benchmark.md`

Status: active
```

- [ ] **Step 5: Run docs-oriented checks**

Run:

```bash
rg -n "central|symmetric_pbc_ed_global|run_pxp_larger_d_benchmark|SQUAREPXP_LARGERD" README.md docs/superpowers/notes/2026-05-17-m3-larger-d-pxp-ed-benchmark.md memory/mid_term/decision_log.md
julia --project=. test/runtests.jl test_public_docs.jl
```

Expected: `rg` shows the new M3 docs and the public doc test passes.

- [ ] **Step 6: Commit**

```bash
git add README.md docs/superpowers/notes/2026-05-17-m3-larger-d-pxp-ed-benchmark.md memory/mid_term/decision_log.md
git commit -m "docs: document M3 PXP ED benchmark semantics"
```

---

### Task 7: Verification And Benchmark Commands

**Files:**
- No new files unless failures require fixes.

- [ ] **Step 1: Run focused exact finite and validation tests**

Run:

```bash
julia --project=. test/runtests.jl test_finite_ipeps_observables.jl test_pxp_validation.jl
```

Expected: PASS.

- [ ] **Step 2: Run M3 tests**

Run:

```bash
julia --project=. test/runtests.jl test_pxp_ed_benchmark.jl test_pxp_larger_d_ed_benchmark.jl
```

Expected: PASS. The 5x5 `total_time = 0.0` smoke should complete without running long 7x7 jobs.

- [ ] **Step 3: Run docs/public API tests**

Run:

```bash
julia --project=. test/runtests.jl test_public_docs.jl test_aqua.jl
```

Expected: PASS.

- [ ] **Step 4: Run extended PXP ED smoke only if runtime is acceptable**

Run:

```bash
SQUAREPXP_EXTENDED_PXP_ED_TESTS=1 julia --project=. test/runtests.jl test_pxp_larger_d_ed_benchmark.jl
```

Expected: PASS. No 7x7 job runs.

- [ ] **Step 5: Run manual 3x3 D=1,2,3,4 campaign**

Run:

```bash
SQUAREPXP_LARGERD_N=3 \
SQUAREPXP_LARGERD_DT=0.02,0.01 \
SQUAREPXP_LARGERD_D=1,2,3,4 \
SQUAREPXP_LARGERD_CUTOFF=1e-12 \
SQUAREPXP_LARGERD_TOTAL_TIME=0.02 \
SQUAREPXP_LARGERD_EXACT_FINITE=true \
SQUAREPXP_LARGERD_EXACT_FINITE_MAX_SITES=9 \
SQUAREPXP_LARGERD_JSON=artifacts/m3-larger-d-3x3.json \
SQUAREPXP_LARGERD_CSV=artifacts/m3-larger-d-3x3.csv \
julia --project=. scripts/pxp_larger_d_ed_benchmark.jl
```

Expected:

```text
artifacts/m3-larger-d-3x3.json
artifacts/m3-larger-d-3x3.csv
```

Then inspect:

```bash
head -n 5 artifacts/m3-larger-d-3x3.csv
```

Expected: rows for D=1,2,3,4 at each dt, exact finite fields populated, simple diagnostic fields present separately.

- [ ] **Step 6: Run manual 5x5 PBC global smoke**

Run:

```bash
SQUAREPXP_LARGERD_N=5 \
SQUAREPXP_LARGERD_DT=0.01 \
SQUAREPXP_LARGERD_D=1,2 \
SQUAREPXP_LARGERD_CUTOFF=1e-12 \
SQUAREPXP_LARGERD_TOTAL_TIME=0.0 \
SQUAREPXP_LARGERD_EXACT_FINITE=true \
SQUAREPXP_LARGERD_EXACT_FINITE_MAX_SITES=9 \
SQUAREPXP_LARGERD_JSON=artifacts/m3-larger-d-5x5-capacity.json \
SQUAREPXP_LARGERD_CSV=artifacts/m3-larger-d-5x5-capacity.csv \
julia --project=. scripts/pxp_larger_d_ed_benchmark.jl
```

Expected: CSV has `ed_basis_dimension = 188`, `ed_constrained_dimension = 25531`, `observable_mode = symmetric_pbc_ed_global`, and blank exact finite fields.

- [ ] **Step 7: Decide whether to run full suite**

If focused tests complete quickly, run:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: full default suite passes. If runtime is too high, state that focused verification passed and full suite was deferred because of runtime.

- [ ] **Step 8: Final commit or PR prep**

After verification:

```bash
git status --short
git log --oneline -5
```

Expected: clean worktree after commits, with M3 commits visible at the top.

---

## Acceptance Mapping

- Documented larger-D benchmark script/API: Tasks 2, 3, 4, and 6.
- 3x3 exact-finite benchmark over D=1,2,3 and preferably D=4: Tasks 2, 5, and 7.
- Feasible odd PBC ED path at minimum 5x5: Task 5 uses 5x5 symmetric PBC global mode.
- 7x7 ED-only smoke or documented runtime boundary: Tasks 5 and 7 document manual 7x7 boundary commands.
- JSON and CSV reports: Task 3.
- Clear observable modes: Tasks 1, 2, 3, and 6.
- Tests for schema and comparison semantics: Tasks 1, 2, 3, and 5.
- README instructions with example commands: Task 6.
- Decision note stating valid and invalid conclusions: Task 6.

## Postponed By Design

- Unreduced constrained 5x5 PBC ED for local observables: feasible, but not needed for M3 global PBC comparison and still not a physical central-region observable under PBC.
- Open-boundary ED central 3x3: scientifically cleaner for literal central-region comparisons, but not PBC and should be a separate milestone.
- CTM-aware/full-update evolution: outside this milestone.
- Publication-grade finite-chi, D, dt, cutoff, and unit-cell sweeps: future audit campaign after M3 tooling exists.

## Self-Review

- Spec coverage: The plan covers ED capacity/semantics, larger-D iPEPS harness, JSON/CSV schema, CLI/env controls, tests/docs, exact finite 3x3, 5x5 PBC global mode, and 7x7 boundary documentation.
- Placeholder scan: The plan intentionally avoids placeholder tokens and unscoped "add tests" steps. Every task has files, concrete tests, commands, expected outcomes, and commit points.
- Type consistency: Public names are consistently `PXPLargerDBenchmarkConfig`, `PXPLargerDBenchmarkSummary`, `PXPLargerDBenchmarkRun`, `PXPLargerDBenchmarkReport`, `run_pxp_larger_d_benchmark`, `write_pxp_larger_d_benchmark_json`, and `write_pxp_larger_d_benchmark_csv`.
