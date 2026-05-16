# PXP Validation Report And ED Harness Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build the first production-facing PXP validation layer: CTM-trusted measurements, an ED-vs-iPEPS short-time comparison harness, and machine-readable validation artifacts.

**Architecture:** Add a focused `PXPValidation` module that composes existing `measure_simple`, `validate_ctm_sweep`, `assess_ctm_trust`, `evolve!`, and `run_pxp_ed_benchmark` APIs. Keep ScarFinder unchanged in this slice; the output report becomes the trusted measurement object ScarFinder can consume in the next slice.

**Tech Stack:** Julia 1.12, ITensors-backed `SquareIPEPSState`, PEPSKit CTMRG measurement adapter, EDKit finite PXP reference, JSON3 artifact writing, existing package test runner.

---

## File Structure

- Create: `src/PXPValidation.jl`
  - Owns `TrustedCTMMeasurement`, `PXPValidationConfig`, iPEPS trajectory samples, ED comparison samples, validation report metadata, `measure_ctm_trusted`, `validate_pxp_ed_ipeps`, and `write_pxp_validation_json`.
- Modify: `src/SquarePXPDynamics.jl`
  - Include `PXPValidation.jl` after `FinitePXPEEDBenchmark.jl`.
  - Re-export the new public validation types and functions.
- Create: `test/test_pxp_validation.jl`
  - Unit tests for trusted CTM composition, validation config checks, ED-vs-iPEPS trajectory reporting, JSON writing, and fake-CTM trusted comparisons.
- Modify: `test/runtests.jl`
  - Add `test_pxp_validation.jl` before `test_scarfinder.jl`.
- Create: `scripts/validate_pxp_ed_ipeps.jl`
  - Fast CLI-style validation script that writes a JSON artifact without running CTMRG by default.
- Modify: `README.md`
  - Document the new validation report path and make clear it is the first trusted validation harness, not a ScarFinder production ranking change.
- Modify: `memory/mid_term/decision_log.md`
  - Record the architectural decision to make CTM-trusted validation reports the bridge between prototype evolution and future ScarFinder ranking.

---

### Task 1: CTM-Trusted Measurement API

**Files:**
- Create: `src/PXPValidation.jl`
- Create: `test/test_pxp_validation.jl`
- Modify: `src/SquarePXPDynamics.jl`
- Modify: `test/runtests.jl`
- Test: `test/test_pxp_validation.jl`
- Test: `test/test_public_docs.jl`

- [ ] **Step 1: Write the failing CTM-trust tests**

Add this file:

```julia
# test/test_pxp_validation.jl
using Test
using SquarePXPDynamics

function _validation_fake_ctm_summary(params; density, blockade, energy, accepted = true)
    return CTMObservableSummary(
        density,
        density,
        density,
        blockade,
        energy,
        CTMRGDiagnostics(
            params.chi,
            params.tol,
            params.maxiter,
            params.maxiter,
            params.tol / 10,
            true,
            accepted,
        ),
    )
end

@testset "trusted CTM measurement composes finite chi sweep and trust policy" begin
    cell = PeriodicSquareUnitCell(3, 3)
    psi = product_square_ipeps(cell; state = :down, maxdim = 1)
    params = (
        PEPSKitCTMRGParams(2, 1e-5, 4, 0),
        PEPSKitCTMRGParams(4, 1e-6, 4, 0),
    )

    trusted = measure_ctm_trusted(
        psi;
        params,
        measure = (state; params) -> _validation_fake_ctm_summary(
            params;
            density = 0.2 + params.chi * 1e-5,
            blockade = params.chi * 1e-6,
            energy = -0.1 - params.chi * 1e-5,
        ),
    )

    @test trusted isa TrustedCTMMeasurement
    @test length(trusted.points) == 2
    @test trusted.measurement === trusted.points[end].measurement
    @test trusted.trust.trusted === true
    @test trusted.trust.reason === :trusted
    @test trusted.measurement.diagnostics.chi == 4
end

@testset "trusted CTM measurement records rejected finite chi drift" begin
    cell = PeriodicSquareUnitCell(3, 3)
    psi = product_square_ipeps(cell; state = :down, maxdim = 1)
    params = (
        PEPSKitCTMRGParams(2, 1e-5, 4, 0),
        PEPSKitCTMRGParams(4, 1e-6, 4, 0),
    )

    rejected = measure_ctm_trusted(
        psi;
        params,
        measure = (state; params) -> _validation_fake_ctm_summary(
            params;
            density = 0.2,
            blockade = 0.0,
            energy = params.chi == 2 ? -0.1 : -0.2,
        ),
    )

    @test rejected.trust.trusted === false
    @test rejected.trust.reason === :energy_delta_too_large
    @test rejected.trust.finite_chi_energy_delta > CTMTrustPolicy().max_energy_delta
end
```

- [ ] **Step 2: Run the focused tests and verify they fail for missing symbols**

Run:

```bash
julia --project=. test/runtests.jl test_pxp_validation.jl
```

Expected: FAIL with `UndefVarError` for `measure_ctm_trusted` or `TrustedCTMMeasurement`.

- [ ] **Step 3: Add the initial `PXPValidation` module**

Create `src/PXPValidation.jl`:

```julia
module PXPValidation

using ..SquareIPEPS: SquareIPEPSState
using ..Observables: measure_simple
using ..PEPSKitMeasurements:
    CTMObservableSummary,
    CTMValidationPoint,
    measure_ctm,
    validate_ctm_sweep
using ..CTMTrust: CTMTrustAssessment, CTMTrustPolicy, assess_ctm_trust

export TrustedCTMMeasurement, measure_ctm_trusted

"""
    TrustedCTMMeasurement(measurement, points, trust)

Finite-`chi` CTMRG measurement bundle for one iPEPS state. `points` stores the
full validation sweep, `measurement` is the last sweep point's CTM observable
summary, and `trust` is the finite-`chi` assessment returned by
[`assess_ctm_trust`](@ref).
"""
struct TrustedCTMMeasurement
    measurement::CTMObservableSummary
    points::Vector{CTMValidationPoint}
    trust::CTMTrustAssessment

    function TrustedCTMMeasurement(
        measurement::CTMObservableSummary,
        points::Vector{CTMValidationPoint},
        trust::CTMTrustAssessment,
    )
        isempty(points) && throw(ArgumentError("trusted CTM measurement requires at least one sweep point"))
        points[end].measurement == measurement ||
            throw(ArgumentError("measurement must match the final CTM validation point"))
        return new(measurement, points, trust)
    end
end

"""
    measure_ctm_trusted(psi; params, policy = CTMTrustPolicy(),
                        reference = measure_simple(psi), measure = measure_ctm)

Run a CTMRG validation sweep for `psi`, assess finite-`chi` trust, and return a
[`TrustedCTMMeasurement`](@ref). The `measure` keyword exists so tests and
benchmark scripts can supply deterministic synthetic CTM summaries without
running PEPSKit CTMRG.
"""
function measure_ctm_trusted(
    psi::SquareIPEPSState;
    params,
    policy::CTMTrustPolicy = CTMTrustPolicy(),
    reference = measure_simple(psi),
    measure = measure_ctm,
)::TrustedCTMMeasurement
    points = validate_ctm_sweep(psi; params, reference, measure)
    assessment = assess_ctm_trust(points; policy)
    return TrustedCTMMeasurement(points[end].measurement, points, assessment)
end

end
```

- [ ] **Step 4: Wire the new module into the top-level package**

Modify `src/SquarePXPDynamics.jl`:

```julia
include("FinitePXPEEDBenchmark.jl")
include("PXPValidation.jl")
include("ScarFinder.jl")
```

Add this `using` block after the existing `FinitePXPEEDBenchmark` imports:

```julia
using .PXPValidation:
    TrustedCTMMeasurement,
    measure_ctm_trusted
```

Add these exports near the finite PXP ED exports:

```julia
export TrustedCTMMeasurement, measure_ctm_trusted
```

- [ ] **Step 5: Register the new test file**

Modify `test/runtests.jl` and insert the new test before `test_scarfinder.jl`:

```julia
    "test_ctm_gauge_readiness.jl",
    "test_pxp_validation.jl",
    "test_scarfinder.jl",
```

- [ ] **Step 6: Run the focused tests and public-doc check**

Run:

```bash
julia --project=. test/runtests.jl test_pxp_validation.jl test_public_docs.jl
```

Expected: PASS.

- [ ] **Step 7: Commit Task 1**

Run:

```bash
git add src/PXPValidation.jl src/SquarePXPDynamics.jl test/runtests.jl test/test_pxp_validation.jl
git commit -m "feat: add trusted CTM validation measurement"
```

---

### Task 2: ED-vs-iPEPS Validation Report

**Files:**
- Modify: `src/PXPValidation.jl`
- Modify: `src/SquarePXPDynamics.jl`
- Modify: `test/test_pxp_validation.jl`
- Test: `test/test_pxp_validation.jl`
- Test: `test/test_public_docs.jl`

- [ ] **Step 1: Add failing report and trajectory tests**

Append to `test/test_pxp_validation.jl`:

```julia
@testset "PXP validation config validates ED and iPEPS controls" begin
    config = PXPValidationConfig(3; total_time = 0.02, dt = 0.01, measure_every = 1)

    @test config.n == 3
    @test config.total_time ≈ 0.02
    @test config.dt ≈ 0.01
    @test config.measure_every == 1
    @test config.initial_state === :down
    @test config.order == 1
    @test config.maxdim == 1
    @test config.schedule === :serial

    @test_throws ArgumentError PXPValidationConfig(2; total_time = 0.02, dt = 0.01)
    @test_throws ArgumentError PXPValidationConfig(3; total_time = 0.025, dt = 0.01)
    @test_throws ArgumentError PXPValidationConfig(3; total_time = 0.02, dt = 0.0)
    @test_throws ArgumentError PXPValidationConfig(3; total_time = 0.02, dt = 0.01, order = 3)
    @test_throws ArgumentError PXPValidationConfig(3; total_time = 0.02, dt = 0.01, schedule = :five_color)
end

@testset "ED-vs-iPEPS validation report samples matched times" begin
    config = PXPValidationConfig(3; total_time = 0.02, dt = 0.01, measure_every = 1)
    report = validate_pxp_ed_ipeps(config; ctm_params = nothing)

    @test report isa PXPValidationReport
    @test report.config === config
    @test report.ed_result.lattice_size == (3, 3)
    @test length(report.ed_result.samples) == 3
    @test length(report.ipeps_samples) == 3
    @test length(report.comparisons) == 3
    @test [s.step for s in report.ipeps_samples] == [0, 1, 2]
    @test [s.time for s in report.ipeps_samples] ≈ [0.0, 0.01, 0.02] atol = 1e-12
    @test report.ipeps_samples[1].evolution === nothing
    @test report.ipeps_samples[2].evolution isa EvolutionLog
    @test report.ipeps_samples[1].ctm === nothing
    @test report.comparisons[1].density_error_simple ≈ 0.0 atol = 1e-12
    @test all(c -> isfinite(c.density_error_simple), report.comparisons)
    @test all(c -> c.ipeps_ctm_density === nothing, report.comparisons)
    @test report.metadata.julia_version == string(VERSION)
end

@testset "ED-vs-iPEPS validation can attach trusted fake CTM" begin
    config = PXPValidationConfig(3; total_time = 0.01, dt = 0.01, measure_every = 1)
    params = (
        PEPSKitCTMRGParams(2, 1e-5, 4, 0),
        PEPSKitCTMRGParams(4, 1e-6, 4, 0),
    )

    report = validate_pxp_ed_ipeps(
        config;
        ctm_params = params,
        ctm_measure = (state; params) -> begin
            simple = measure_simple(state)
            return _validation_fake_ctm_summary(
                params;
                density = simple.density + params.chi * 1e-5,
                blockade = simple.blockade_violation + params.chi * 1e-6,
                energy = simple.pxp_energy_density - params.chi * 1e-5,
            )
        end,
    )

    @test all(sample -> sample.ctm isa TrustedCTMMeasurement, report.ipeps_samples)
    @test all(comparison -> comparison.ctm_trusted === true, report.comparisons)
    @test all(comparison -> comparison.ctm_reason === :trusted, report.comparisons)
    @test all(comparison -> comparison.ipeps_ctm_density !== nothing, report.comparisons)
end
```

- [ ] **Step 2: Run the focused test and verify the new report symbols fail**

Run:

```bash
julia --project=. test/runtests.jl test_pxp_validation.jl
```

Expected: FAIL with missing `PXPValidationConfig`, `PXPValidationReport`, or `validate_pxp_ed_ipeps`.

- [ ] **Step 3: Extend imports and exports in `src/PXPValidation.jl`**

At the top of `src/PXPValidation.jl`, replace the imports and export list with:

```julia
using ..SquareIPEPS: SquareIPEPSState, log_norm, product_square_ipeps
using ..SquareUnitCells: PeriodicSquareUnitCell
using ..Observables: SimpleObservableSummary, measure_simple
using ..PEPSKitMeasurements:
    CTMObservableSummary,
    CTMValidationPoint,
    measure_ctm,
    validate_ctm_sweep
using ..CTMTrust: CTMTrustAssessment, CTMTrustPolicy, assess_ctm_trust
using ..IPEPSEvolution: EvolutionLog, TrotterParams, evolve!
using ..FinitePXPEEDBenchmark:
    PXPEEDBenchmarkConfig,
    PXPEEDBenchmarkResult,
    PXPEEDSample,
    run_pxp_ed_benchmark

export TrustedCTMMeasurement, measure_ctm_trusted
export PXPValidationConfig, PXPValidationMetadata, PXPIPEPSSample
export PXPEDComparisonSample, PXPValidationReport, validate_pxp_ed_ipeps
```

- [ ] **Step 4: Add validation report types**

Append these definitions after `measure_ctm_trusted` in `src/PXPValidation.jl`:

```julia
function _finite_nonnegative(value::Real, label::String)
    converted = Float64(value)
    isfinite(converted) && converted >= 0 ||
        throw(ArgumentError("$label must be finite and nonnegative"))
    return converted
end

function _finite_positive(value::Real, label::String)
    converted = Float64(value)
    isfinite(converted) && converted > 0 ||
        throw(ArgumentError("$label must be finite and positive"))
    return converted
end

function _nonnegative_int(value::Integer, label::String)
    converted = Int(value)
    converted >= 0 || throw(ArgumentError("$label must be nonnegative"))
    return converted
end

function _positive_int(value::Integer, label::String)
    converted = Int(value)
    converted >= 1 || throw(ArgumentError("$label must be at least 1"))
    return converted
end

"""
    PXPValidationConfig(n; kwargs...)

Controls for a short-time finite-ED versus iPEPS PXP validation run. The ED
side uses [`PXPEEDBenchmarkConfig`](@ref); the iPEPS side uses a periodic
`n x n` unit cell, all-down product initialization, real-time PXP evolution,
and a serial star schedule by default so `3 x 3` smoke validation is supported.
"""
struct PXPValidationConfig
    n::Int
    total_time::Float64
    dt::Float64
    measure_every::Int
    initial_state::Symbol
    point_group::Bool
    use_sparse::Bool
    ed_tol::Float64
    ed_m_init::Int
    ed_m_max::Int
    ed_extend_step::Int
    order::Int
    maxdim::Int
    cutoff::Float64
    schedule::Symbol

    function PXPValidationConfig(
        n::Integer;
        total_time::Real = 0.02,
        dt::Real = 0.01,
        measure_every::Integer = 1,
        initial_state::Symbol = :down,
        point_group::Bool = true,
        use_sparse::Bool = true,
        ed_tol::Real = 1e-10,
        ed_m_init::Integer = 30,
        ed_m_max::Integer = 60,
        ed_extend_step::Integer = 10,
        order::Integer = 1,
        maxdim::Integer = 1,
        cutoff::Real = 1e-12,
        schedule::Symbol = :serial,
    )
        n_int = _positive_int(n, "n")
        total = _finite_nonnegative(total_time, "total_time")
        step = _finite_positive(dt, "dt")
        cadence = _positive_int(measure_every, "measure_every")
        initial_state in (:down, :all_down) ||
            throw(ArgumentError("initial_state must be :down or :all_down"))
        tol = _finite_positive(ed_tol, "ed_tol")
        m_init = _positive_int(ed_m_init, "ed_m_init")
        m_max = _positive_int(ed_m_max, "ed_m_max")
        m_max >= m_init || throw(ArgumentError("ed_m_max must be at least ed_m_init"))
        extend = _positive_int(ed_extend_step, "ed_extend_step")
        ord = _positive_int(order, "order")
        ord in (1, 2) || throw(ArgumentError("order must be 1 or 2"))
        dim = _positive_int(maxdim, "maxdim")
        trunc_cutoff = _finite_nonnegative(cutoff, "cutoff")
        schedule in (:serial, :five_color) ||
            throw(ArgumentError("schedule must be :serial or :five_color"))
        schedule === :five_color && n_int % 5 != 0 &&
            throw(ArgumentError("five_color validation requires n divisible by 5"))

        PXPEEDBenchmarkConfig(
            n_int;
            total_time = total,
            dt = step,
            measure_every = cadence,
            initial_state,
            point_group,
            use_sparse,
            tol,
            m_init,
            m_max,
            extend_step = extend,
        )
        TrotterParams(step, ord, :real, dim, trunc_cutoff; schedule)

        return new(
            n_int,
            total,
            step,
            cadence,
            initial_state,
            point_group,
            use_sparse,
            tol,
            m_init,
            m_max,
            extend,
            ord,
            dim,
            trunc_cutoff,
            schedule,
        )
    end
end

"""
    PXPValidationMetadata

Reproducibility metadata attached to a PXP validation report. `git_commit` is
`nothing` when the repository commit cannot be read from the current process.
"""
struct PXPValidationMetadata
    git_commit::Union{Nothing,String}
    julia_version::String
    project_path::Union{Nothing,String}
end

"""
    PXPIPEPSSample

One iPEPS validation sample at an ED measurement time. The first sample has
`evolution === nothing`; later samples store the incremental [`EvolutionLog`](@ref)
for the time interval since the previous sample. `ctm` is present only when a
finite-`chi` CTM sweep was requested.
"""
struct PXPIPEPSSample
    step::Int
    time::Float64
    simple::SimpleObservableSummary
    evolution::Union{Nothing,EvolutionLog}
    ctm::Union{Nothing,TrustedCTMMeasurement}
    log_norm::Float64
end

"""
    PXPEDComparisonSample

Observable comparison at one matched ED/iPEPS sample time. Density errors are
reported as `iPEPS - ED`; CTM fields are `nothing` when no CTM sweep was run.
"""
struct PXPEDComparisonSample
    step::Int
    time::Float64
    ed_return_probability::Float64
    ed_excitation_density::Float64
    ipeps_simple_density::Float64
    ipeps_ctm_density::Union{Nothing,Float64}
    density_error_simple::Float64
    density_error_ctm::Union{Nothing,Float64}
    simple_blockade_violation::Float64
    ctm_blockade_violation::Union{Nothing,Float64}
    ctm_trusted::Union{Nothing,Bool}
    ctm_reason::Union{Nothing,Symbol}
end

"""
    PXPValidationReport

Complete short-time finite PXP validation artifact: input config, finite ED
result, matched iPEPS trajectory samples, observable comparisons, and
reproducibility metadata.
"""
struct PXPValidationReport
    config::PXPValidationConfig
    ed_result::PXPEEDBenchmarkResult
    ipeps_samples::Vector{PXPIPEPSSample}
    comparisons::Vector{PXPEDComparisonSample}
    metadata::PXPValidationMetadata
end
```

- [ ] **Step 5: Add trajectory and comparison functions**

Append these functions after the report types:

```julia
function _validation_ed_config(config::PXPValidationConfig)
    return PXPEEDBenchmarkConfig(
        config.n;
        total_time = config.total_time,
        dt = config.dt,
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

function _validation_trotter(config::PXPValidationConfig)
    return TrotterParams(
        config.dt,
        config.order,
        :real,
        config.maxdim,
        config.cutoff;
        schedule = config.schedule,
    )
end

function _validation_initial_state(config::PXPValidationConfig)
    cell = PeriodicSquareUnitCell(config.n, config.n)
    state = config.initial_state === :all_down ? :down : config.initial_state
    return product_square_ipeps(cell; state, maxdim = config.maxdim)
end

function _git_commit()
    try
        commit = chomp(read(`git rev-parse HEAD`, String))
        isempty(commit) && return nothing
        return commit
    catch
        return nothing
    end
end

function _validation_metadata()
    project = Base.active_project()
    return PXPValidationMetadata(_git_commit(), string(VERSION), project)
end

function _maybe_trusted_ctm(
    psi::SquareIPEPSState,
    ctm_params,
    trust_policy::CTMTrustPolicy,
    ctm_measure,
)
    ctm_params === nothing && return nothing
    return measure_ctm_trusted(
        psi;
        params = ctm_params,
        policy = trust_policy,
        measure = ctm_measure,
    )
end

function _ctm_density(ctm::Nothing)
    return nothing
end

function _ctm_density(ctm::TrustedCTMMeasurement)
    return ctm.measurement.density
end

function _ctm_blockade(ctm::Nothing)
    return nothing
end

function _ctm_blockade(ctm::TrustedCTMMeasurement)
    return ctm.measurement.blockade_violation
end

function _ctm_trusted(ctm::Nothing)
    return nothing
end

function _ctm_trusted(ctm::TrustedCTMMeasurement)
    return ctm.trust.trusted
end

function _ctm_reason(ctm::Nothing)
    return nothing
end

function _ctm_reason(ctm::TrustedCTMMeasurement)
    return ctm.trust.reason
end

function _comparison(ed::PXPEEDSample, sample::PXPIPEPSSample)
    ctm_density = _ctm_density(sample.ctm)
    return PXPEDComparisonSample(
        ed.step,
        ed.time,
        ed.return_probability,
        ed.excitation_density,
        sample.simple.density,
        ctm_density,
        sample.simple.density - ed.excitation_density,
        ctm_density === nothing ? nothing : ctm_density - ed.excitation_density,
        sample.simple.blockade_violation,
        _ctm_blockade(sample.ctm),
        _ctm_trusted(sample.ctm),
        _ctm_reason(sample.ctm),
    )
end

"""
    validate_pxp_ed_ipeps(config; ctm_params = nothing,
                          trust_policy = CTMTrustPolicy(),
                          ctm_measure = measure_ctm)

Run a finite periodic PXP ED trajectory and a matched all-down iPEPS trajectory
on an `n x n` unit cell. iPEPS samples are measured at the same times as ED.
When `ctm_params` is supplied, every iPEPS sample also receives a trusted CTM
measurement bundle from [`measure_ctm_trusted`](@ref).
"""
function validate_pxp_ed_ipeps(
    config::PXPValidationConfig;
    ctm_params = nothing,
    trust_policy::CTMTrustPolicy = CTMTrustPolicy(),
    ctm_measure = measure_ctm,
)::PXPValidationReport
    ed_result = run_pxp_ed_benchmark(_validation_ed_config(config))
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
            ),
        )
    end

    comparisons = [
        _comparison(ed_sample, ipeps_sample) for
        (ed_sample, ipeps_sample) in zip(ed_result.samples, samples)
    ]
    return PXPValidationReport(
        config,
        ed_result,
        samples,
        comparisons,
        _validation_metadata(),
    )
end
```

- [ ] **Step 6: Re-export the report API**

Modify `src/SquarePXPDynamics.jl` and replace the `PXPValidation` import block with:

```julia
using .PXPValidation:
    TrustedCTMMeasurement,
    measure_ctm_trusted,
    PXPValidationConfig,
    PXPValidationMetadata,
    PXPIPEPSSample,
    PXPEDComparisonSample,
    PXPValidationReport,
    validate_pxp_ed_ipeps
```

Replace the validation export line with:

```julia
export TrustedCTMMeasurement, measure_ctm_trusted
export PXPValidationConfig, PXPValidationMetadata, PXPIPEPSSample
export PXPEDComparisonSample, PXPValidationReport, validate_pxp_ed_ipeps
```

- [ ] **Step 7: Run focused tests and docs check**

Run:

```bash
julia --project=. test/runtests.jl test_pxp_validation.jl test_public_docs.jl
```

Expected: PASS.

- [ ] **Step 8: Commit Task 2**

Run:

```bash
git add src/PXPValidation.jl src/SquarePXPDynamics.jl test/test_pxp_validation.jl
git commit -m "feat: add PXP ED iPEPS validation report"
```

---

### Task 3: Machine-Readable JSON Artifacts And Script

**Files:**
- Modify: `src/PXPValidation.jl`
- Modify: `src/SquarePXPDynamics.jl`
- Modify: `test/test_pxp_validation.jl`
- Create: `scripts/validate_pxp_ed_ipeps.jl`
- Test: `test/test_pxp_validation.jl`

- [ ] **Step 1: Add failing JSON writer test**

Append to `test/test_pxp_validation.jl`:

```julia
@testset "PXP validation report writes JSON artifact" begin
    config = PXPValidationConfig(3; total_time = 0.01, dt = 0.01, measure_every = 1)
    report = validate_pxp_ed_ipeps(config; ctm_params = nothing)
    path = tempname() * ".json"

    written = write_pxp_validation_json(report, path)
    data = read(path, String)

    @test written == path
    @test occursin("\"config\"", data)
    @test occursin("\"comparisons\"", data)
    @test occursin("\"density_error_simple\"", data)
    @test occursin("\"metadata\"", data)
end
```

- [ ] **Step 2: Run the focused test and verify the writer symbol fails**

Run:

```bash
julia --project=. test/runtests.jl test_pxp_validation.jl
```

Expected: FAIL with missing `write_pxp_validation_json`.

- [ ] **Step 3: Add JSON3 import and writer export**

At the top of `src/PXPValidation.jl`, add:

```julia
using JSON3
```

Update the export list:

```julia
export PXPEDComparisonSample, PXPValidationReport, validate_pxp_ed_ipeps
export write_pxp_validation_json
```

- [ ] **Step 4: Add JSON conversion helpers and writer**

Append to `src/PXPValidation.jl`:

```julia
function _json_value(value::Nothing)
    return nothing
end

function _json_value(value::Symbol)
    return String(value)
end

function _config_data(config::PXPValidationConfig)
    return (;
        n = config.n,
        total_time = config.total_time,
        dt = config.dt,
        measure_every = config.measure_every,
        initial_state = String(config.initial_state),
        point_group = config.point_group,
        use_sparse = config.use_sparse,
        ed_tol = config.ed_tol,
        ed_m_init = config.ed_m_init,
        ed_m_max = config.ed_m_max,
        ed_extend_step = config.ed_extend_step,
        order = config.order,
        maxdim = config.maxdim,
        cutoff = config.cutoff,
        schedule = String(config.schedule),
    )
end

function _metadata_data(metadata::PXPValidationMetadata)
    return (;
        git_commit = metadata.git_commit,
        julia_version = metadata.julia_version,
        project_path = metadata.project_path,
    )
end

function _simple_data(summary::SimpleObservableSummary)
    return (;
        density = summary.density,
        density_even = summary.density_even,
        density_odd = summary.density_odd,
        blockade_violation = summary.blockade_violation,
        pxp_energy_density = summary.pxp_energy_density,
        mean_bond_entropy = summary.mean_bond_entropy,
        max_bond_entropy = summary.max_bond_entropy,
    )
end

function _ctm_diagnostics_data(diagnostics::Nothing)
    return nothing
end

function _ctm_diagnostics_data(diagnostics)
    return (;
        chi = diagnostics.chi,
        tol = diagnostics.tol,
        maxiter = diagnostics.maxiter,
        iterations = diagnostics.iterations,
        residual = diagnostics.residual,
        converged = diagnostics.converged,
        accepted = diagnostics.accepted,
    )
end

function _ctm_summary_data(summary::CTMObservableSummary)
    return (;
        density = summary.density,
        density_even = summary.density_even,
        density_odd = summary.density_odd,
        blockade_violation = summary.blockade_violation,
        pxp_energy_density = summary.pxp_energy_density,
        diagnostics = _ctm_diagnostics_data(summary.diagnostics),
    )
end

function _ctm_point_data(point::CTMValidationPoint)
    return (;
        chi = point.params.chi,
        tol = point.params.tol,
        maxiter = point.params.maxiter,
        verbosity = point.params.verbosity,
        measurement = _ctm_summary_data(point.measurement),
        delta_density = point.delta_density,
        delta_density_even = point.delta_density_even,
        delta_density_odd = point.delta_density_odd,
        delta_blockade_violation = point.delta_blockade_violation,
        delta_pxp_energy_density = point.delta_pxp_energy_density,
    )
end

function _trust_data(trust::CTMTrustAssessment)
    return (;
        trusted = trust.trusted,
        reason = String(trust.reason),
        message = trust.message,
        compared_points = trust.compared_points,
        finite_chi_density_delta = trust.finite_chi_density_delta,
        finite_chi_blockade_delta = trust.finite_chi_blockade_delta,
        finite_chi_energy_delta = trust.finite_chi_energy_delta,
        observed_max_residual = trust.observed_max_residual,
    )
end

function _trusted_ctm_data(ctm::Nothing)
    return nothing
end

function _trusted_ctm_data(ctm::TrustedCTMMeasurement)
    return (;
        measurement = _ctm_summary_data(ctm.measurement),
        points = _ctm_point_data.(ctm.points),
        trust = _trust_data(ctm.trust),
    )
end

function _evolution_data(evolution::Nothing)
    return nothing
end

function _evolution_data(evolution::EvolutionLog)
    return (;
        total_time = evolution.total_time,
        nsteps = evolution.nsteps,
        dt = evolution.params.dt,
        order = evolution.params.order,
        evolution = String(evolution.params.evolution),
        maxdim = evolution.params.maxdim,
        cutoff = evolution.params.cutoff,
        schedule = String(evolution.params.schedule),
        max_truncerr = evolution.max_truncerr,
        max_bond_entropy = evolution.max_bond_entropy,
        mean_bond_entropy = evolution.mean_bond_entropy,
        log_norm_before = evolution.log_norm_before,
        log_norm_after = evolution.log_norm_after,
        log_norm_delta = evolution.log_norm_delta,
        model_metadata = evolution.model_metadata,
    )
end

function _ed_sample_data(sample::PXPEEDSample)
    return (;
        step = sample.step,
        time = sample.time,
        norm = sample.norm,
        return_probability = sample.return_probability,
        excitation_density = sample.excitation_density,
    )
end

function _ed_result_data(result::PXPEEDBenchmarkResult)
    return (;
        lattice_size = result.lattice_size,
        basis_dimension = result.basis_dimension,
        constrained_dimension = result.constrained_dimension,
        group_order = result.group_order,
        point_group = result.point_group,
        hamiltonian_nnz = result.hamiltonian_nnz,
        samples = _ed_sample_data.(result.samples),
    )
end

function _ipeps_sample_data(sample::PXPIPEPSSample)
    return (;
        step = sample.step,
        time = sample.time,
        simple = _simple_data(sample.simple),
        evolution = _evolution_data(sample.evolution),
        ctm = _trusted_ctm_data(sample.ctm),
        log_norm = sample.log_norm,
    )
end

function _comparison_data(sample::PXPEDComparisonSample)
    return (;
        step = sample.step,
        time = sample.time,
        ed_return_probability = sample.ed_return_probability,
        ed_excitation_density = sample.ed_excitation_density,
        ipeps_simple_density = sample.ipeps_simple_density,
        ipeps_ctm_density = sample.ipeps_ctm_density,
        density_error_simple = sample.density_error_simple,
        density_error_ctm = sample.density_error_ctm,
        simple_blockade_violation = sample.simple_blockade_violation,
        ctm_blockade_violation = sample.ctm_blockade_violation,
        ctm_trusted = sample.ctm_trusted,
        ctm_reason = _json_value(sample.ctm_reason),
    )
end

function _report_data(report::PXPValidationReport)
    return (;
        config = _config_data(report.config),
        metadata = _metadata_data(report.metadata),
        ed_result = _ed_result_data(report.ed_result),
        ipeps_samples = _ipeps_sample_data.(report.ipeps_samples),
        comparisons = _comparison_data.(report.comparisons),
    )
end

"""
    write_pxp_validation_json(report, path)

Write a [`PXPValidationReport`](@ref) to `path` as a JSON artifact containing
configuration, metadata, ED samples, iPEPS samples, CTM trust data when present,
and matched observable comparisons.
"""
function write_pxp_validation_json(report::PXPValidationReport, path::AbstractString)
    open(path, "w") do io
        JSON3.write(io, _report_data(report))
        write(io, '\n')
    end
    return path
end
```

- [ ] **Step 5: Re-export the JSON writer**

Modify `src/SquarePXPDynamics.jl` and add `write_pxp_validation_json` to the `using .PXPValidation` block and export list:

```julia
using .PXPValidation:
    TrustedCTMMeasurement,
    measure_ctm_trusted,
    PXPValidationConfig,
    PXPValidationMetadata,
    PXPIPEPSSample,
    PXPEDComparisonSample,
    PXPValidationReport,
    validate_pxp_ed_ipeps,
    write_pxp_validation_json

export TrustedCTMMeasurement, measure_ctm_trusted
export PXPValidationConfig, PXPValidationMetadata, PXPIPEPSSample
export PXPEDComparisonSample, PXPValidationReport, validate_pxp_ed_ipeps
export write_pxp_validation_json
```

- [ ] **Step 6: Add the fast validation script**

Create `scripts/validate_pxp_ed_ipeps.jl`:

```julia
using Pkg

project_root = dirname(@__DIR__)
Pkg.activate(project_root; io = devnull)

using SquarePXPDynamics

function _env_int(name::String, default::Int)
    return parse(Int, get(ENV, name, string(default)))
end

function _env_float(name::String, default::Float64)
    return parse(Float64, get(ENV, name, string(default)))
end

function _env_symbol(name::String, default::Symbol)
    return Symbol(get(ENV, name, String(default)))
end

out = get(
    ENV,
    "SQUAREPXP_PXP_VALIDATION_OUT",
    joinpath(project_root, "artifacts", "pxp_validation_report.json"),
)

config = PXPValidationConfig(
    _env_int("SQUAREPXP_PXP_VALIDATION_N", 3);
    total_time = _env_float("SQUAREPXP_PXP_VALIDATION_TOTAL_TIME", 0.02),
    dt = _env_float("SQUAREPXP_PXP_VALIDATION_DT", 0.01),
    measure_every = _env_int("SQUAREPXP_PXP_VALIDATION_MEASURE_EVERY", 1),
    order = _env_int("SQUAREPXP_PXP_VALIDATION_ORDER", 1),
    maxdim = _env_int("SQUAREPXP_PXP_VALIDATION_MAXDIM", 1),
    cutoff = _env_float("SQUAREPXP_PXP_VALIDATION_CUTOFF", 1e-12),
    schedule = _env_symbol("SQUAREPXP_PXP_VALIDATION_SCHEDULE", :serial),
)

report = validate_pxp_ed_ipeps(config; ctm_params = nothing)
mkpath(dirname(out))
write_pxp_validation_json(report, out)
println(out)
```

- [ ] **Step 7: Run focused tests and the script**

Run:

```bash
julia --project=. test/runtests.jl test_pxp_validation.jl test_public_docs.jl
julia --project=. scripts/validate_pxp_ed_ipeps.jl
```

Expected:

- First command: PASS.
- Second command: prints `.../artifacts/pxp_validation_report.json`.

- [ ] **Step 8: Commit Task 3**

Run:

```bash
git add src/PXPValidation.jl src/SquarePXPDynamics.jl test/test_pxp_validation.jl scripts/validate_pxp_ed_ipeps.jl
git commit -m "feat: write PXP validation report artifacts"
```

---

### Task 4: Documentation And Decision Log

**Files:**
- Modify: `README.md`
- Modify: `memory/mid_term/decision_log.md`
- Test: documentation smoke through `rg`

- [ ] **Step 1: Document the validation report path in README**

Add a bullet near the existing CTM/ED/ScarFinder feature list:

```markdown
- PXP validation reports that compare finite ED all-down trajectories against
  matched iPEPS trajectories and optionally attach trusted finite-`chi` CTM
  measurement sweeps via `validate_pxp_ed_ipeps` and
  `write_pxp_validation_json` (`src/PXPValidation.jl`).
```

Add a short section after the existing CTM guidance:

```markdown
### PXP validation reports

The first production-facing validation path is
`validate_pxp_ed_ipeps(PXPValidationConfig(...))`. It runs a finite periodic
PXP ED trajectory, evolves a matched all-down iPEPS trajectory on the same
unit cell, and reports density differences at shared sample times. Passing
`ctm_params = (...)` attaches `measure_ctm_trusted` output at every sample:
the final CTM measurement, the finite-`chi` sweep points, and the
`assess_ctm_trust` result.

For a fast JSON artifact without CTMRG:

```julia
config = PXPValidationConfig(3; total_time = 0.02, dt = 0.01)
report = validate_pxp_ed_ipeps(config; ctm_params = nothing)
write_pxp_validation_json(report, "artifacts/pxp_validation_report.json")
```

or from the shell:

```bash
julia --project=. scripts/validate_pxp_ed_ipeps.jl
```

This is a validation harness, not a ScarFinder ranking change. ScarFinder still
uses its existing simple/local default until a later slice wires validation
reports into candidate objectives.
```

- [ ] **Step 2: Record the architectural decision**

Append this entry to `memory/mid_term/decision_log.md`:

```markdown
## 2026-05-16 - Promote CTM-Trusted PXP Validation Reports

Decision:

Add a focused `PXPValidation` layer that composes existing CTM sweep, CTM
trust, iPEPS evolution, and finite PXP ED APIs into a machine-readable
validation report.

Reason:

The next reliability step is to make trusted measurement and ED comparison a
normal workflow before changing ScarFinder ranking. A narrow report layer
preserves the existing explicit `measure_simple` / `measure_ctm` boundary and
avoids a broad measurement facade before a second production backend exists.

Consequences:

Short-time PXP runs can now produce reproducible JSON artifacts with ED
density references, iPEPS diagnostics, optional finite-`chi` CTM trust, and
run metadata. ScarFinder can consume this report shape in a follow-on slice.

Source:

`src/PXPValidation.jl`; `test/test_pxp_validation.jl`;
`scripts/validate_pxp_ed_ipeps.jl`

Status: active
```

- [ ] **Step 3: Run documentation smoke checks**

Run:

```bash
rg -n "PXP validation reports|validate_pxp_ed_ipeps|measure_ctm_trusted|PXPValidation" README.md memory/mid_term/decision_log.md src test scripts
git diff --check
```

Expected:

- `rg` finds the new README, source, test, script, and decision-log references.
- `git diff --check` exits successfully.

- [ ] **Step 4: Commit Task 4**

Run:

```bash
git add README.md memory/mid_term/decision_log.md
git commit -m "docs: document PXP validation reports"
```

---

### Task 5: Final Verification

**Files:**
- No file changes required unless verification finds a defect.
- Test: full package suite.

- [ ] **Step 1: Run focused validation and public API checks**

Run:

```bash
julia --project=. test/runtests.jl test_pxp_validation.jl test_public_docs.jl
```

Expected: PASS.

- [ ] **Step 2: Run the fast artifact script**

Run:

```bash
rm -f artifacts/pxp_validation_report.json
julia --project=. scripts/validate_pxp_ed_ipeps.jl
test -s artifacts/pxp_validation_report.json
```

Expected: script prints `artifacts/pxp_validation_report.json`, and the file exists with nonzero size.

- [ ] **Step 3: Run the full suite**

Run:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: PASS.

- [ ] **Step 4: Inspect final status**

Run:

```bash
git status --short --branch
```

Expected: clean working tree on the implementation branch, or only expected generated artifacts under `artifacts/`.

---

## Self-Review Checklist

- Spec coverage: The plan implements CTM-trusted measurement, ED-vs-iPEPS validation, JSON artifacts, docs, and the decision-log entry needed for the P0 reliability slice.
- Placeholder scan: The plan contains concrete file paths, commands, type names, function signatures, and code blocks for each code change.
- Type consistency: Public symbols are consistently named `TrustedCTMMeasurement`, `measure_ctm_trusted`, `PXPValidationConfig`, `PXPIPEPSSample`, `PXPEDComparisonSample`, `PXPValidationReport`, `validate_pxp_ed_ipeps`, and `write_pxp_validation_json`.
- Scope check: ScarFinder objectives, candidate persistence, deterministic PEPSKit CTMRG initialization, new observables, and full-update algorithms are intentionally outside this implementation slice.

Plan complete and saved to `docs/superpowers/plans/2026-05-16-pxp-validation-report-ed-harness.md`. Two execution options:

1. **Subagent-Driven (recommended)** - Dispatch a fresh subagent per task, review between tasks, fast iteration.
2. **Inline Execution** - Execute tasks in this session using executing-plans, batch execution with checkpoints.
