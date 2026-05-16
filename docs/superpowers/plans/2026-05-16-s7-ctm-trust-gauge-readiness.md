# S7 CTM Trust And Gauge Readiness Implementation Plan

> **Status:** Completed and merged locally as part of S0-S7. The unchecked boxes
> below are the original execution template, not current TODOs. See
> `docs/superpowers/notes/2026-05-16-s0-s7-completion-audit.md`.

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add S7a CTMRG trust assessment, finite-chi audit output, and read-only simple-gauge diagnostics without implementing mutating gauge fixing.

**Architecture:** Add `CTMTrust.jl` as a leaf module over existing `PEPSKitMeasurements` records and `GaugeDiagnostics.jl` as a leaf module over `SquareIPEPS`/`SquareUnitCells`. Export only stable public APIs through `SquarePXPDynamics.jl`, keep PEPSKit-specific measurement code in `PEPSKitMeasurements.jl`, and keep S7b gauge mutation as a documented handoff.

**Tech Stack:** Julia 1.12, ITensors.jl, LinearAlgebra, existing `SquarePXPDynamics` modules, `Test`.

---

## File Structure

- Create `src/CTMTrust.jl`: trust policy, trust assessment, finite-chi drift calculation, and trust CSV writer.
- Create `src/GaugeDiagnostics.jl`: read-only simple-gauge local bond norm diagnostics.
- Modify `src/SquarePXPDynamics.jl`: include the two new modules, import and export stable public symbols.
- Create `test/test_ctm_trust.jl`: synthetic trust-policy, finite-chi, and CSV tests.
- Create `test/test_gauge_diagnostics.jl`: D=1, D=2 diagonal, D=2 off-diagonal, alias-direction, and validation tests.
- Modify `test/test_pepskit_measurements.jl`: add stale aggregate CTM energy check.
- Modify `test/test_public_docs.jl`: add explicit docstring checks for new exported symbols.
- Modify `test/runtests.jl`: include the two new focused test files.
- Modify `README.md`: document the CTM trust workflow and clarify S7a vs S7b.
- Create `docs/superpowers/notes/2026-05-16-s7b-gauge-fixing-handoff.md`: follow-up constraints for mutating gauge fixing.

## Task 1: Add CTM Trust Tests

**Files:**
- Create: `test/test_ctm_trust.jl`
- Modify: `test/runtests.jl`

- [ ] **Step 1: Add the test file to the runner**

Patch `test/runtests.jl` so the new CTM trust tests run before PEPSKit measurement integration tests:

```julia
const TEST_FILES = [
    "test_spinops.jl",
    "test_square_geometry.jl",
    "test_square_pxp.jl",
    "test_star_models.jl",
    "test_square_peps.jl",
    "test_square_unitcells.jl",
    "test_square_ipeps.jl",
    "test_square_ipeps_s2.jl",
    "test_star_simple_update.jl",
    "test_ipeps_evolution.jl",
    "test_observables_evolved.jl",
    "test_tfim_observables.jl",
    "test_tfim_schedule_reference.jl",
    "test_benchmarks.jl",
    "test_ctm_trust.jl",
    "test_pepskit_measurements.jl",
    "test_scarfinder.jl",
    "test_public_docs.jl",
    "test_aqua.jl",
]
```

- [ ] **Step 2: Write the failing CTM trust tests**

Create `test/test_ctm_trust.jl`:

```julia
function _trust_diag_test(chi; tol = 1e-8, residual = tol / 10, accepted = true)
    return CTMRGDiagnostics(chi, tol, 100, 12, residual, true, accepted)
end

function _trust_summary_test(;
    density,
    blockade,
    energy,
    diag = _trust_diag_test(4),
)
    return CTMObservableSummary(density, density, density, blockade, energy, diag)
end

function _trust_point_test(;
    chi,
    tol,
    density,
    blockade,
    energy,
    reference_density = 0.0,
    reference_blockade = 0.0,
    reference_energy = 0.0,
    diag = _trust_diag_test(chi; tol),
)
    reference = CTMObservableSummary(
        reference_density,
        reference_density,
        reference_density,
        reference_blockade,
        reference_energy,
        CTMRGDiagnostics(chi, tol, 100, 12, tol / 10, true, true),
    )
    measurement = _trust_summary_test(
        density = density,
        blockade = blockade,
        energy = energy,
        diag = diag,
    )
    return CTMValidationPoint(PEPSKitCTMRGParams(chi, tol, 100, 0), reference, measurement)
end

@testset "CTM trust policy validation" begin
    policy = CTMTrustPolicy()

    @test policy.min_points == 2
    @test policy.require_accepted_diagnostics === true
    @test policy.max_density_delta ≈ 1e-3
    @test policy.max_blockade_delta ≈ 1e-4
    @test policy.max_energy_delta ≈ 1e-3
    @test policy.max_residual === nothing

    @test_throws ArgumentError CTMTrustPolicy(1, true, 1e-3, 1e-4, 1e-3, nothing)
    @test_throws ArgumentError CTMTrustPolicy(2, true, -1.0, 1e-4, 1e-3, nothing)
    @test_throws ArgumentError CTMTrustPolicy(2, true, 1e-3, NaN, 1e-3, nothing)
    @test_throws ArgumentError CTMTrustPolicy(2, true, 1e-3, 1e-4, Inf, nothing)
    @test_throws ArgumentError CTMTrustPolicy(2, true, 1e-3, 1e-4, 1e-3, -1e-8)
end

@testset "CTM trust assessment validation" begin
    assessment = CTMTrustAssessment(true, :trusted, "trusted", 2, 1e-5, 1e-6, 1e-5, 1e-9)

    @test assessment.trusted === true
    @test assessment.reason === :trusted
    @test assessment.compared_points == 2
    @test assessment.finite_chi_density_delta ≈ 1e-5
    @test assessment.observed_max_residual ≈ 1e-9
    @test_throws ArgumentError CTMTrustAssessment(true, :density_delta_too_large, "bad", 2, 1e-5, 1e-6, 1e-5, 1e-9)
    @test_throws ArgumentError CTMTrustAssessment(false, :bad_reason, "bad", 2, nothing, nothing, nothing, nothing)
    @test_throws ArgumentError CTMTrustAssessment(false, :too_few_points, "bad", -1, nothing, nothing, nothing, nothing)
    @test_throws ArgumentError CTMTrustAssessment(false, :too_few_points, "bad", 0, Inf, nothing, nothing, nothing)
end

@testset "CTM trust accepts stable final chi window" begin
    points = [
        _trust_point_test(chi = 2, tol = 1e-5, density = 9.0, blockade = 3.0, energy = -4.0),
        _trust_point_test(
            chi = 4,
            tol = 1e-6,
            density = 10.0,
            blockade = 5.0,
            energy = -7.0,
            reference_density = -100.0,
            reference_blockade = -100.0,
            reference_energy = 100.0,
        ),
        _trust_point_test(
            chi = 8,
            tol = 1e-8,
            density = 10.0002,
            blockade = 5.00002,
            energy = -7.0002,
            reference_density = -100.0,
            reference_blockade = -100.0,
            reference_energy = 100.0,
        ),
    ]

    assessment = assess_ctm_trust(points)

    @test assessment.trusted === true
    @test assessment.reason === :trusted
    @test assessment.compared_points == 2
    @test assessment.finite_chi_density_delta ≈ 2e-4 atol = 1e-12
    @test assessment.finite_chi_blockade_delta ≈ 2e-5 atol = 1e-12
    @test assessment.finite_chi_energy_delta ≈ 2e-4 atol = 1e-12
    @test assessment.observed_max_residual !== nothing
    @test assessment.observed_max_residual <= 1e-6
end

@testset "CTM trust rejects expected trust failures" begin
    stable4 = _trust_point_test(chi = 4, tol = 1e-6, density = 0.1, blockade = 0.01, energy = -0.2)
    stable8 = _trust_point_test(chi = 8, tol = 1e-8, density = 0.1001, blockade = 0.01001, energy = -0.2001)

    too_few = assess_ctm_trust([stable4])
    @test too_few.trusted === false
    @test too_few.reason === :too_few_points
    @test too_few.compared_points == 1

    nonmonotonic = assess_ctm_trust([stable8, stable4])
    @test nonmonotonic.trusted === false
    @test nonmonotonic.reason === :nonmonotonic_sweep

    loose_tol = assess_ctm_trust([
        _trust_point_test(chi = 4, tol = 1e-8, density = 0.1, blockade = 0.01, energy = -0.2),
        _trust_point_test(chi = 8, tol = 1e-6, density = 0.1001, blockade = 0.01001, energy = -0.2001),
    ])
    @test loose_tol.trusted === false
    @test loose_tol.reason === :nonmonotonic_sweep

    missing_diag = assess_ctm_trust([
        _trust_point_test(chi = 4, tol = 1e-6, density = 0.1, blockade = 0.01, energy = -0.2),
        _trust_point_test(
            chi = 8,
            tol = 1e-8,
            density = 0.1001,
            blockade = 0.01001,
            energy = -0.2001,
            diag = nothing,
        ),
    ])
    @test missing_diag.trusted === false
    @test missing_diag.reason === :missing_diagnostics

    unaccepted = assess_ctm_trust([
        stable4,
        _trust_point_test(
            chi = 8,
            tol = 1e-8,
            density = 0.1001,
            blockade = 0.01001,
            energy = -0.2001,
            diag = _trust_diag_test(8; tol = 1e-8, accepted = false),
        ),
    ])
    @test unaccepted.trusted === false
    @test unaccepted.reason === :unaccepted_diagnostics

    missing_residual = assess_ctm_trust(
        [
            stable4,
            _trust_point_test(
                chi = 8,
                tol = 1e-8,
                density = 0.1001,
                blockade = 0.01001,
                energy = -0.2001,
                diag = CTMRGDiagnostics(8, 1e-8, 100, 12, nothing, true, true),
            ),
        ];
        policy = CTMTrustPolicy(2, true, 1e-3, 1e-4, 1e-3, 1e-6),
    )
    @test missing_residual.trusted === false
    @test missing_residual.reason === :missing_residual

    residual_large = assess_ctm_trust(
        [
            stable4,
            _trust_point_test(
                chi = 8,
                tol = 1e-8,
                density = 0.1001,
                blockade = 0.01001,
                energy = -0.2001,
                diag = _trust_diag_test(8; tol = 1e-8, residual = 1e-4),
            ),
        ];
        policy = CTMTrustPolicy(2, true, 1e-3, 1e-4, 1e-3, 1e-6),
    )
    @test residual_large.trusted === false
    @test residual_large.reason === :residual_too_large
end

@testset "CTM trust rejects finite-chi observable drift" begin
    base = _trust_point_test(chi = 4, tol = 1e-6, density = 0.2, blockade = 0.01, energy = -0.3)

    density_drift = assess_ctm_trust([
        base,
        _trust_point_test(chi = 8, tol = 1e-8, density = 0.202, blockade = 0.01001, energy = -0.3001),
    ])
    @test density_drift.trusted === false
    @test density_drift.reason === :density_delta_too_large

    blockade_drift = assess_ctm_trust([
        base,
        _trust_point_test(chi = 8, tol = 1e-8, density = 0.2001, blockade = 0.0102, energy = -0.3001),
    ])
    @test blockade_drift.trusted === false
    @test blockade_drift.reason === :blockade_delta_too_large

    energy_drift = assess_ctm_trust([
        base,
        _trust_point_test(chi = 8, tol = 1e-8, density = 0.2001, blockade = 0.01001, energy = -0.302),
    ])
    @test energy_drift.trusted === false
    @test energy_drift.reason === :energy_delta_too_large
end

@testset "CTM trust CSV audit output" begin
    points = [
        _trust_point_test(chi = 4, tol = 1e-6, density = 0.2, blockade = 0.01, energy = -0.3),
        _trust_point_test(chi = 8, tol = 1e-8, density = 0.2001, blockade = 0.01001, energy = -0.3001),
    ]
    path = tempname() * ".csv"

    write_ctm_trust_csv(points, path)
    csv = read(path, String)

    @test occursin("trust_policy_min_points", csv)
    @test occursin("trust_finite_chi_density_delta", csv)
    @test occursin("trust_observed_max_residual", csv)
    @test occursin("trusted", lowercase(csv))
    @test occursin("8,1.0e-8,100", csv)
    @test !occursin("message", lowercase(csv))
end

@testset "CTM trust malformed inputs throw" begin
    @test_throws ArgumentError assess_ctm_trust(Any["not a point"])
end
```

- [ ] **Step 3: Run the new CTM trust tests and verify they fail**

Run:

```bash
julia --project=. test/runtests.jl test_ctm_trust.jl
```

Expected: FAIL during loading or execution because `CTMTrustPolicy`, `CTMTrustAssessment`, `assess_ctm_trust`, and `write_ctm_trust_csv` are not defined.

- [ ] **Step 4: Commit the failing tests**

```bash
git add test/runtests.jl test/test_ctm_trust.jl
git commit -m "test: add CTM trust policy coverage"
```

## Task 2: Implement CTMTrust Module

**Files:**
- Create: `src/CTMTrust.jl`
- Modify: `src/SquarePXPDynamics.jl`
- Test: `test/test_ctm_trust.jl`
- Test: `test/test_public_docs.jl`

- [ ] **Step 1: Create `src/CTMTrust.jl`**

Create `src/CTMTrust.jl`:

```julia
module CTMTrust

using ..PEPSKitMeasurements:
    PEPSKitCTMRGParams, CTMRGDiagnostics, CTMObservableSummary, CTMValidationPoint

export CTMTrustPolicy, CTMTrustAssessment, assess_ctm_trust, write_ctm_trust_csv

const CTM_TRUST_REASONS = (
    :trusted,
    :too_few_points,
    :nonmonotonic_sweep,
    :missing_diagnostics,
    :unaccepted_diagnostics,
    :missing_residual,
    :residual_too_large,
    :density_delta_too_large,
    :blockade_delta_too_large,
    :energy_delta_too_large,
)

"""
    CTMTrustPolicy([min_points, require_accepted_diagnostics, max_density_delta, max_blockade_delta, max_energy_delta, max_residual])

Software-level finite-chi acceptance thresholds for CTMRG validation sweeps.
The default policy requires the final two sweep points to have accepted CTMRG
diagnostics and small CTM-to-CTM observable drift. These thresholds are
regression trust defaults, not universal physics-quality criteria.
"""
struct CTMTrustPolicy
    min_points::Int
    require_accepted_diagnostics::Bool
    max_density_delta::Float64
    max_blockade_delta::Float64
    max_energy_delta::Float64
    max_residual::Union{Float64,Nothing}

    function CTMTrustPolicy(
        min_points::Integer,
        require_accepted_diagnostics::Bool,
        max_density_delta::Real,
        max_blockade_delta::Real,
        max_energy_delta::Real,
        max_residual::Union{Real,Nothing},
    )
        min_points >= 2 || throw(ArgumentError("min_points must be at least 2"))
        density = _finite_nonnegative(max_density_delta, "max_density_delta")
        blockade = _finite_nonnegative(max_blockade_delta, "max_blockade_delta")
        energy = _finite_nonnegative(max_energy_delta, "max_energy_delta")
        residual = if max_residual === nothing
            nothing
        else
            _finite_nonnegative(max_residual, "max_residual")
        end
        return new(
            Int(min_points),
            require_accepted_diagnostics,
            density,
            blockade,
            energy,
            residual,
        )
    end
end

CTMTrustPolicy() = CTMTrustPolicy(2, true, 1e-3, 1e-4, 1e-3, nothing)

"""
    CTMTrustAssessment

Structured result from [`assess_ctm_trust`](@ref). `finite_chi_*_delta`
fields are maximum adjacent CTM-to-CTM drifts over the assessed final sweep
window. `observed_max_residual` is the largest residual reported in that
window when residual metadata is available.
"""
struct CTMTrustAssessment
    trusted::Bool
    reason::Symbol
    message::String
    compared_points::Int
    finite_chi_density_delta::Union{Float64,Nothing}
    finite_chi_blockade_delta::Union{Float64,Nothing}
    finite_chi_energy_delta::Union{Float64,Nothing}
    observed_max_residual::Union{Float64,Nothing}

    function CTMTrustAssessment(
        trusted::Bool,
        reason::Symbol,
        message::AbstractString,
        compared_points::Integer,
        finite_chi_density_delta::Union{Real,Nothing},
        finite_chi_blockade_delta::Union{Real,Nothing},
        finite_chi_energy_delta::Union{Real,Nothing},
        observed_max_residual::Union{Real,Nothing},
    )
        reason in CTM_TRUST_REASONS ||
            throw(ArgumentError("unknown CTM trust reason: $reason"))
        compared_points >= 0 ||
            throw(ArgumentError("compared_points must be nonnegative"))
        (!trusted || reason === :trusted) ||
            throw(ArgumentError("trusted assessments must use reason :trusted"))
        density = _maybe_finite_nonnegative(finite_chi_density_delta, "finite_chi_density_delta")
        blockade =
            _maybe_finite_nonnegative(finite_chi_blockade_delta, "finite_chi_blockade_delta")
        energy = _maybe_finite_nonnegative(finite_chi_energy_delta, "finite_chi_energy_delta")
        residual = _maybe_finite_nonnegative(observed_max_residual, "observed_max_residual")
        return new(
            trusted,
            reason,
            String(message),
            Int(compared_points),
            density,
            blockade,
            energy,
            residual,
        )
    end
end

function _finite_nonnegative(value::Real, label::AbstractString)::Float64
    converted = Float64(value)
    isfinite(converted) && converted >= 0 ||
        throw(ArgumentError("$label must be finite and nonnegative"))
    return converted
end

function _maybe_finite_nonnegative(value::Union{Real,Nothing}, label::AbstractString)
    value === nothing && return nothing
    return _finite_nonnegative(value, label)
end

function _trust_assessment(
    trusted::Bool,
    reason::Symbol,
    message::AbstractString,
    compared_points::Integer,
    density_delta,
    blockade_delta,
    energy_delta,
    residual,
)
    return CTMTrustAssessment(
        trusted,
        reason,
        message,
        compared_points,
        density_delta,
        blockade_delta,
        energy_delta,
        residual,
    )
end

function _as_points(points)
    collected = collect(points)
    all(point -> point isa CTMValidationPoint, collected) ||
        throw(ArgumentError("points must contain only CTMValidationPoint values"))
    return CTMValidationPoint[point for point in collected]
end

function _assert_finite_summary(summary::CTMObservableSummary, label::AbstractString)
    all(
        isfinite,
        (
            summary.density,
            summary.density_even,
            summary.density_odd,
            summary.blockade_violation,
            summary.pxp_energy_density,
        ),
    ) || throw(ArgumentError("$label CTM summary must be finite"))
    return summary
end

function _assert_finite_points(points::Vector{CTMValidationPoint})
    for (i, point) in pairs(points)
        _assert_finite_summary(point.reference, "reference point $i")
        _assert_finite_summary(point.measurement, "measurement point $i")
    end
    return points
end

function _observed_max_residual(window)
    residuals = Float64[]
    for point in window
        diagnostics = point.diagnostics
        diagnostics === nothing && continue
        diagnostics.residual === nothing && continue
        push!(residuals, diagnostics.residual)
    end
    return isempty(residuals) ? nothing : maximum(residuals)
end

function _check_monotonic_window(window)
    for i in 2:length(window)
        previous = window[i - 1].params
        current = window[i].params
        current.chi > previous.chi || return false
        current.tol <= previous.tol || return false
    end
    return true
end

function _max_adjacent_delta(window, field::Symbol)
    values = Float64[]
    for i in 2:length(window)
        current = getfield(window[i].measurement, field)
        previous = getfield(window[i - 1].measurement, field)
        push!(values, abs(current - previous))
    end
    return maximum(values)
end

function _diagnostics_missing(window)
    return any(point -> point.diagnostics === nothing, window)
end

function _diagnostics_unaccepted(window)
    return any(point -> point.diagnostics !== nothing && point.diagnostics.accepted === false, window)
end

function _residual_missing(window)
    return any(
        point -> point.diagnostics === nothing || point.diagnostics.residual === nothing,
        window,
    )
end

"""
    assess_ctm_trust(points; policy = CTMTrustPolicy())

Assess whether the final `policy.min_points` entries in a CTMRG validation
sweep are stable enough for downstream measurement trust. Trust deltas compare
adjacent CTM measurements in the final sweep window and never use
`CTMValidationPoint.delta_*`, which are CTM-minus-simple-reference deltas.
"""
function assess_ctm_trust(
    points;
    policy::CTMTrustPolicy = CTMTrustPolicy(),
)::CTMTrustAssessment
    collected = _assert_finite_points(_as_points(points))
    npoints = length(collected)
    if npoints < policy.min_points
        return _trust_assessment(
            false,
            :too_few_points,
            "need at least $(policy.min_points) CTM sweep points, got $npoints",
            npoints,
            nothing,
            nothing,
            nothing,
            _observed_max_residual(collected),
        )
    end

    window = collected[(end - policy.min_points + 1):end]
    compared = length(window)
    observed_residual = _observed_max_residual(window)

    if !_check_monotonic_window(window)
        return _trust_assessment(
            false,
            :nonmonotonic_sweep,
            "final CTM sweep window must have strictly increasing chi and nonincreasing tol",
            compared,
            nothing,
            nothing,
            nothing,
            observed_residual,
        )
    end

    if policy.require_accepted_diagnostics && _diagnostics_missing(window)
        return _trust_assessment(false, :missing_diagnostics, "CTM diagnostics missing", compared, nothing, nothing, nothing, observed_residual)
    end
    if policy.require_accepted_diagnostics && _diagnostics_unaccepted(window)
        return _trust_assessment(false, :unaccepted_diagnostics, "CTM diagnostics were not accepted", compared, nothing, nothing, nothing, observed_residual)
    end
    if policy.max_residual !== nothing && _residual_missing(window)
        return _trust_assessment(false, :missing_residual, "CTM residual required but missing", compared, nothing, nothing, nothing, observed_residual)
    end
    if policy.max_residual !== nothing && observed_residual !== nothing &&
       observed_residual > policy.max_residual
        return _trust_assessment(false, :residual_too_large, "CTM residual exceeds policy threshold", compared, nothing, nothing, nothing, observed_residual)
    end

    density_delta = _max_adjacent_delta(window, :density)
    blockade_delta = _max_adjacent_delta(window, :blockade_violation)
    energy_delta = _max_adjacent_delta(window, :pxp_energy_density)

    if density_delta > policy.max_density_delta
        return _trust_assessment(false, :density_delta_too_large, "finite-chi density drift exceeds policy threshold", compared, density_delta, blockade_delta, energy_delta, observed_residual)
    elseif blockade_delta > policy.max_blockade_delta
        return _trust_assessment(false, :blockade_delta_too_large, "finite-chi blockade drift exceeds policy threshold", compared, density_delta, blockade_delta, energy_delta, observed_residual)
    elseif energy_delta > policy.max_energy_delta
        return _trust_assessment(false, :energy_delta_too_large, "finite-chi energy drift exceeds policy threshold", compared, density_delta, blockade_delta, energy_delta, observed_residual)
    end

    return _trust_assessment(
        true,
        :trusted,
        "final CTM sweep window satisfies trust policy",
        compared,
        density_delta,
        blockade_delta,
        energy_delta,
        observed_residual,
    )
end

function _csv_value(value::Nothing)
    return ""
end

function _csv_value(value::Bool)
    return string(value)
end

function _csv_value(value::Real)
    return string(value)
end

function _csv_value(value::Symbol)
    return String(value)
end

function _diagnostic_field(::Nothing, ::Symbol)
    return nothing
end

function _diagnostic_field(diagnostics::CTMRGDiagnostics, field::Symbol)
    return getfield(diagnostics, field)
end

"""
    write_ctm_trust_csv(points, path; policy = CTMTrustPolicy())

Write CTM sweep points and the internally computed trust assessment to `path`
as a flat CSV. The trust assessment fields are repeated on each row so the file
is self-contained for audit and spreadsheet use.
"""
function write_ctm_trust_csv(
    points,
    path::AbstractString;
    policy::CTMTrustPolicy = CTMTrustPolicy(),
)
    collected = _as_points(points)
    assessment = assess_ctm_trust(collected; policy)
    header = (
        "chi",
        "tol",
        "maxiter",
        "verbosity",
        "ctm_density",
        "ctm_blockade_violation",
        "ctm_pxp_energy_density",
        "ctm_iterations",
        "ctm_residual",
        "ctm_converged",
        "ctm_accepted",
        "trust_policy_min_points",
        "trust_policy_require_accepted_diagnostics",
        "trust_policy_max_density_delta",
        "trust_policy_max_blockade_delta",
        "trust_policy_max_energy_delta",
        "trust_policy_max_residual",
        "trust_trusted",
        "trust_reason",
        "trust_compared_points",
        "trust_finite_chi_density_delta",
        "trust_finite_chi_blockade_delta",
        "trust_finite_chi_energy_delta",
        "trust_observed_max_residual",
    )
    open(path, "w") do io
        println(io, join(header, ","))
        for point in collected
            row = (
                point.params.chi,
                point.params.tol,
                point.params.maxiter,
                point.params.verbosity,
                point.measurement.density,
                point.measurement.blockade_violation,
                point.measurement.pxp_energy_density,
                _diagnostic_field(point.diagnostics, :iterations),
                _diagnostic_field(point.diagnostics, :residual),
                _diagnostic_field(point.diagnostics, :converged),
                _diagnostic_field(point.diagnostics, :accepted),
                policy.min_points,
                policy.require_accepted_diagnostics,
                policy.max_density_delta,
                policy.max_blockade_delta,
                policy.max_energy_delta,
                policy.max_residual,
                assessment.trusted,
                assessment.reason,
                assessment.compared_points,
                assessment.finite_chi_density_delta,
                assessment.finite_chi_blockade_delta,
                assessment.finite_chi_energy_delta,
                assessment.observed_max_residual,
            )
            println(io, join(map(_csv_value, row), ","))
        end
    end
    return path
end

end
```

- [ ] **Step 2: Wire CTMTrust into the top-level module**

Patch `src/SquarePXPDynamics.jl`:

```julia
include("PEPSKitMeasurements.jl")
include("CTMTrust.jl")
include("StarSimpleUpdate.jl")
```

Add imports after the `PEPSKitMeasurements` imports:

```julia
using .CTMTrust:
    CTMTrustPolicy, CTMTrustAssessment, assess_ctm_trust, write_ctm_trust_csv
```

Add exports next to the CTM exports:

```julia
export CTMTrustPolicy, CTMTrustAssessment, assess_ctm_trust, write_ctm_trust_csv
```

- [ ] **Step 3: Run the CTM trust tests**

Run:

```bash
julia --project=. test/runtests.jl test_ctm_trust.jl
```

Expected: PASS.

- [ ] **Step 4: Run public doc tests**

Run:

```bash
julia --project=. test/runtests.jl test_public_docs.jl
```

Expected: PASS. If it fails, add or fix docstrings in `src/CTMTrust.jl`; do not weaken `test/test_public_docs.jl`.

- [ ] **Step 5: Commit CTMTrust implementation**

```bash
git add src/CTMTrust.jl src/SquarePXPDynamics.jl test/test_ctm_trust.jl test/runtests.jl
git commit -m "feat: add CTM trust assessment"
```

## Task 3: Add Gauge Diagnostics Tests

**Files:**
- Create: `test/test_gauge_diagnostics.jl`
- Modify: `test/runtests.jl`

- [ ] **Step 1: Add the gauge diagnostics test file to the runner**

Patch `test/runtests.jl` so the new test runs near the other iPEPS/simple-observable tests:

```julia
const TEST_FILES = [
    "test_spinops.jl",
    "test_square_geometry.jl",
    "test_square_pxp.jl",
    "test_star_models.jl",
    "test_square_peps.jl",
    "test_square_unitcells.jl",
    "test_square_ipeps.jl",
    "test_square_ipeps_s2.jl",
    "test_star_simple_update.jl",
    "test_ipeps_evolution.jl",
    "test_observables_evolved.jl",
    "test_gauge_diagnostics.jl",
    "test_tfim_observables.jl",
    "test_tfim_schedule_reference.jl",
    "test_benchmarks.jl",
    "test_ctm_trust.jl",
    "test_pepskit_measurements.jl",
    "test_scarfinder.jl",
    "test_public_docs.jl",
    "test_aqua.jl",
]
```

- [ ] **Step 2: Write the failing gauge diagnostics tests**

Create `test/test_gauge_diagnostics.jl`:

```julia
using ITensors
using LinearAlgebra

function _seeded_offdiagonal_d2_gauge_state_test(cell)
    psi = product_square_ipeps(cell; state = :down, maxdim = 2)
    for c in cell.reps
        p = physical_index(psi, c)
        left = link_index(psi, c, :left)
        right = link_index(psi, c, :right)
        up = link_index(psi, c, :up)
        down = link_index(psi, c, :down)
        T = ITensor(ComplexF64, p, left, right, up, down)
        for pv = 1:dim(p), lv = 1:dim(left), rv = 1:dim(right), uv = 1:dim(up), dv = 1:dim(down)
            re = 0.01 * (11c.x + 7c.y + 5pv + 3lv + 2rv + uv + dv)
            im = 0.005 * (13c.x - 3c.y + 2pv + lv - rv + uv - dv)
            T[p=>pv, left=>lv, right=>rv, up=>uv, down=>dv] = complex(re, im)
        end
        psi.tensors[c] = T
        set_link_weight!(psi, c, :right, [0.8, 0.6])
        set_link_weight!(psi, c, :up, [0.6, 0.8])
    end
    return psi
end

@testset "simple gauge diagnostics vanish on product states" begin
    cell = PeriodicSquareUnitCell(4, 4)
    psi = product_square_ipeps(cell; state = :down, maxdim = 1)
    c = SquareCoord(2, 2)

    diagnostic = gauge_diagnostic_simple(psi, c, :right)
    alias = gauge_diagnostic_simple(psi, neighbor(cell, c, :right), :left)
    deviations = all_gauge_deviations_simple(psi)

    @test diagnostic isa SimpleGaugeDiagnostic
    @test diagnostic.bond == bondkey(cell, c, :right)
    @test alias.bond == diagnostic.bond
    @test alias.deviation ≈ diagnostic.deviation atol = 1e-14
    @test diagnostic.deviation ≈ 0 atol = 1e-14
    @test diagnostic.frobenius_norm > 0
    @test diagnostic.diagonal_min >= 0
    @test diagnostic.diagonal_max >= diagnostic.diagonal_min
    @test diagnostic.diagonal_condition >= 1 || isinf(diagnostic.diagonal_condition)
    @test length(deviations) == length(psi.link_weights)
    @test all(value -> isapprox(value, 0; atol = 1e-14), values(deviations))
end

@testset "simple gauge diagnostics handle diagonal D=2 states" begin
    cell = PeriodicSquareUnitCell(4, 4)
    psi = product_square_ipeps(cell; state = :down, maxdim = 2)

    for c in cell.reps, dir in (:right, :up)
        diagnostic = gauge_diagnostic_simple(psi, c, dir)
        @test diagnostic.deviation ≈ 0 atol = 1e-14
        @test isfinite(diagnostic.frobenius_norm)
        @test diagnostic.diagonal_min >= 0
        @test diagnostic.diagonal_max >= diagnostic.diagonal_min
    end
end

@testset "simple gauge diagnostics detect off-diagonal D=2 fixture" begin
    cell = PeriodicSquareUnitCell(4, 4)
    psi = _seeded_offdiagonal_d2_gauge_state_test(cell)
    c = SquareCoord(2, 2)

    diagnostic = gauge_diagnostic_simple(psi, c, :right)
    deviations = all_gauge_deviations_simple(psi)

    @test diagnostic.deviation > 1e-10
    @test isfinite(diagnostic.deviation)
    @test isfinite(diagnostic.frobenius_norm)
    @test diagnostic.frobenius_norm > 0
    @test diagnostic.diagonal_min >= 0
    @test diagnostic.diagonal_max >= diagnostic.diagonal_min
    @test all(isfinite, values(deviations))
    @test maximum(values(deviations)) > 1e-10
end

@testset "simple gauge diagnostics validate directions and gauge" begin
    psi = product_square_ipeps(PeriodicSquareUnitCell(4, 4); state = :down, maxdim = 1)
    @test_throws ArgumentError gauge_diagnostic_simple(psi, SquareCoord(2, 2), :bad)

    bad = SquareIPEPSState(
        psi.unitcell,
        psi.tensors,
        psi.physical_indices,
        psi.link_indices,
        psi.link_weights,
        psi.maxdim,
        :not_simple,
        Ref(state_version(psi)),
        Ref(log_norm(psi)),
    )
    @test_throws ArgumentError gauge_diagnostic_simple(bad, SquareCoord(2, 2), :right)
end
```

- [ ] **Step 3: Run the gauge diagnostics tests and verify they fail**

Run:

```bash
julia --project=. test/runtests.jl test_gauge_diagnostics.jl
```

Expected: FAIL during loading or execution because `SimpleGaugeDiagnostic`, `gauge_diagnostic_simple`, `gauge_deviation_simple`, and `all_gauge_deviations_simple` are not defined.

- [ ] **Step 4: Commit the failing gauge diagnostics tests**

```bash
git add test/runtests.jl test/test_gauge_diagnostics.jl
git commit -m "test: add simple gauge diagnostics coverage"
```

## Task 4: Implement GaugeDiagnostics Module

**Files:**
- Create: `src/GaugeDiagnostics.jl`
- Modify: `src/SquarePXPDynamics.jl`
- Test: `test/test_gauge_diagnostics.jl`
- Test: `test/test_public_docs.jl`

- [ ] **Step 1: Create `src/GaugeDiagnostics.jl`**

Create `src/GaugeDiagnostics.jl`:

```julia
module GaugeDiagnostics

using ITensors
using LinearAlgebra
using ..SquareGeometry: SquareCoord
using ..SquareUnitCells: BondKey, bondkey, neighbor
using ..SquareIPEPS:
    SquareIPEPSState,
    physical_index,
    link_index,
    link_weight,
    absorb_link_weight

export SimpleGaugeDiagnostic
export gauge_diagnostic_simple, gauge_deviation_simple, all_gauge_deviations_simple

const _DIRECTIONS = (:right, :up, :left, :down)

"""
    SimpleGaugeDiagnostic

Read-only diagnostic for one canonical simple-update bond. `deviation` is the
relative Frobenius norm of the off-diagonal part of the local two-site bond
norm matrix in the current lambda basis. This is not a canonical-gauge
certificate; it is a local warning signal for poor simple-update gauges.
"""
struct SimpleGaugeDiagnostic
    bond::BondKey
    deviation::Float64
    frobenius_norm::Float64
    diagonal_min::Float64
    diagonal_max::Float64
    diagonal_condition::Float64

    function SimpleGaugeDiagnostic(
        bond::BondKey,
        deviation::Real,
        frobenius_norm::Real,
        diagonal_min::Real,
        diagonal_max::Real,
        diagonal_condition::Real,
    )
        dev = Float64(deviation)
        frob = Float64(frobenius_norm)
        dmin = Float64(diagonal_min)
        dmax = Float64(diagonal_max)
        cond = Float64(diagonal_condition)
        isfinite(dev) && dev >= 0 ||
            throw(ArgumentError("gauge deviation must be finite and nonnegative"))
        isfinite(frob) && frob > 0 ||
            throw(ArgumentError("gauge norm must be finite and positive"))
        isfinite(dmin) || throw(ArgumentError("diagonal_min must be finite"))
        isfinite(dmax) || throw(ArgumentError("diagonal_max must be finite"))
        dmax >= dmin || throw(ArgumentError("diagonal_max must be at least diagonal_min"))
        (isfinite(cond) && cond >= 0) || isinf(cond) ||
            throw(ArgumentError("diagonal_condition must be nonnegative or Inf"))
        return new(bond, dev, frob, dmin, dmax, cond)
    end
end

function _validate_direction(dir::Symbol)
    dir in _DIRECTIONS ||
        throw(ArgumentError("direction must be :right, :up, :left, or :down"))
    return dir
end

function _opposite_direction(dir::Symbol)
    dir === :right && return :left
    dir === :left && return :right
    dir === :up && return :down
    dir === :down && return :up
    throw(ArgumentError("direction must be :right, :up, :left, or :down"))
end

function _external_dirs(target_dir::Symbol)
    return Tuple(dir for dir in _DIRECTIONS if dir !== target_dir)
end

function _assert_simple_gauge(psi::SquareIPEPSState)
    psi.gauge === :simple || throw(ArgumentError("gauge diagnostics require psi.gauge === :simple"))
    return nothing
end

function _canonical_bond_and_endpoint(psi::SquareIPEPSState, c::SquareCoord, dir::Symbol)
    _validate_direction(dir)
    key = bondkey(psi.unitcell, c, dir)
    return key, key.site, key.dir, neighbor(psi.unitcell, key.site, key.dir)
end

function _external_absorbed_endpoint(
    psi::SquareIPEPSState,
    c::SquareCoord,
    target_dir::Symbol,
)
    T = copy(psi.tensors[c])
    for dir in _external_dirs(target_dir)
        T = absorb_link_weight(T, psi, c, dir)
    end
    return T
end

function _endpoint_gram_matrix(
    psi::SquareIPEPSState,
    c::SquareCoord,
    target_dir::Symbol,
)
    T = _external_absorbed_endpoint(psi, c, target_dir)
    p = physical_index(psi, c)
    target = link_index(psi, c, target_dir)
    external = map(dir -> link_index(psi, c, dir), _external_dirs(target_dir))
    data = Array(T, p, target, external...)
    dim_target = dim(target)
    gram = zeros(ComplexF64, dim_target, dim_target)
    for i = 1:dim_target, j = 1:dim_target
        total = 0.0 + 0.0im
        for s in axes(data, 1), e1 in axes(data, 3), e2 in axes(data, 4), e3 in axes(data, 5)
            total += conj(data[s, i, e1, e2, e3]) * data[s, j, e1, e2, e3]
        end
        gram[i, j] = total
    end
    return gram
end

function _bond_norm_matrix(
    psi::SquareIPEPSState,
    a::SquareCoord,
    dir::Symbol,
    b::SquareCoord,
)
    left_gram = _endpoint_gram_matrix(psi, a, dir)
    right_gram = _endpoint_gram_matrix(psi, b, _opposite_direction(dir))
    size(left_gram) == size(right_gram) ||
        throw(ArgumentError("endpoint Gram matrices have incompatible dimensions"))
    norm_matrix = left_gram .* right_gram
    return (norm_matrix + norm_matrix') / 2
end

"""
    gauge_diagnostic_simple(psi, c, dir)::SimpleGaugeDiagnostic

Return read-only local simple-gauge diagnostics for the periodic nearest
neighbor bond from `c` in `dir`. All four directions are accepted and mapped to
the canonical stored `BondKey`.
"""
function gauge_diagnostic_simple(
    psi::SquareIPEPSState,
    c::SquareCoord,
    dir::Symbol,
)::SimpleGaugeDiagnostic
    _assert_simple_gauge(psi)
    key, a, canonical_dir, b = _canonical_bond_and_endpoint(psi, c, dir)
    N = _bond_norm_matrix(psi, a, canonical_dir, b)
    normN = norm(N)
    isfinite(normN) && normN > 0 ||
        throw(ArgumentError("local bond norm matrix must have finite nonzero norm"))
    diagonal = real.(diag(N))
    offdiag = N - Diagonal(diag(N))
    deviation = norm(offdiag) / normN
    diagonal_min = minimum(diagonal)
    diagonal_max = maximum(diagonal)
    diagonal_condition = diagonal_min > 0 ? diagonal_max / diagonal_min : Inf
    return SimpleGaugeDiagnostic(
        key,
        deviation,
        normN,
        diagonal_min,
        diagonal_max,
        diagonal_condition,
    )
end

"""
    gauge_deviation_simple(psi, c, dir)::Float64

Return only the off-diagonal local simple-gauge deviation for the periodic
nearest-neighbor bond from `c` in `dir`.
"""
function gauge_deviation_simple(
    psi::SquareIPEPSState,
    c::SquareCoord,
    dir::Symbol,
)::Float64
    return gauge_diagnostic_simple(psi, c, dir).deviation
end

"""
    all_gauge_deviations_simple(psi)::Dict{BondKey,Float64}

Return read-only simple-gauge deviations for every canonical stored link
weight in `psi`.
"""
function all_gauge_deviations_simple(psi::SquareIPEPSState)
    _assert_simple_gauge(psi)
    return Dict(
        key => gauge_deviation_simple(psi, key.site, key.dir) for
        key in keys(psi.link_weights)
    )
end

end
```

- [ ] **Step 2: Wire GaugeDiagnostics into the top-level module**

Patch `src/SquarePXPDynamics.jl`:

```julia
include("SquareIPEPS.jl")
include("GaugeDiagnostics.jl")
include("StarModels.jl")
```

Add imports after the `SquareIPEPS` imports:

```julia
using .GaugeDiagnostics:
    SimpleGaugeDiagnostic,
    gauge_diagnostic_simple,
    gauge_deviation_simple,
    all_gauge_deviations_simple
```

Add exports near the other iPEPS diagnostics:

```julia
export SimpleGaugeDiagnostic
export gauge_diagnostic_simple, gauge_deviation_simple, all_gauge_deviations_simple
```

- [ ] **Step 3: Run the gauge diagnostics tests**

Run:

```bash
julia --project=. test/runtests.jl test_gauge_diagnostics.jl
```

Expected: PASS.

- [ ] **Step 4: Run public doc tests**

Run:

```bash
julia --project=. test/runtests.jl test_public_docs.jl
```

Expected: PASS. If it fails, add docstrings to every exported symbol in `src/GaugeDiagnostics.jl`; do not weaken doc coverage.

- [ ] **Step 5: Commit GaugeDiagnostics implementation**

```bash
git add src/GaugeDiagnostics.jl src/SquarePXPDynamics.jl test/test_gauge_diagnostics.jl test/runtests.jl
git commit -m "feat: add simple gauge diagnostics"
```

## Task 5: Add CTM Stale-Context Regression For Aggregate Energy

**Files:**
- Modify: `test/test_pepskit_measurements.jl`

- [ ] **Step 1: Add the stale aggregate energy assertion**

In `test/test_pepskit_measurements.jl`, inside the existing `"five-site CTM star expectation smoke test"` testset after the current stale `local_density_ctm` assertion, add:

```julia
        @test_throws ArgumentError pxp_energy_density_ctm(psi_stale, ctx_stale)
```

The surrounding block should end like this:

```julia
        psi_stale = product_square_ipeps(cell; state = :down, maxdim = 1)
        ctx_stale = pepskit_ctmrg_context(psi_stale; params)
        project_star!(psi_stale, c, dt; evolution = :real, projected = true, maxdim = 1)
        @test_throws ArgumentError local_density_ctm(psi_stale, c, ctx_stale)
        @test_throws ArgumentError pxp_energy_density_ctm(psi_stale, ctx_stale)
```

- [ ] **Step 2: Run the PEPSKit measurement tests**

Run:

```bash
julia --project=. test/runtests.jl test_pepskit_measurements.jl
```

Expected: PASS.

- [ ] **Step 3: Commit the stale-context regression**

```bash
git add test/test_pepskit_measurements.jl
git commit -m "test: cover stale CTM energy contexts"
```

## Task 6: Update README CTM Trust Workflow

**Files:**
- Modify: `README.md`

- [ ] **Step 1: Replace the finite-chi validation example**

In the CTM guidance section that currently shows `validate_ctm_sweep`, replace the example with:

```julia
points = validate_ctm_sweep(
    psi;
    params = [
        PEPSKitCTMRGParams(4, 1e-6, 50, 0),
        PEPSKitCTMRGParams(8, 1e-8, 100, 0),
    ],
)
assessment = assess_ctm_trust(points)
write_ctm_validation_csv(points, "ctm-validation.csv")
write_ctm_trust_csv(points, "ctm-trust.csv")
```

Immediately after the example, add:

```markdown
`assess_ctm_trust` compares the final finite-`chi` CTM measurements against
each other; it does not use the simple/local reference deltas stored in
`CTMValidationPoint`. A trusted assessment is a measurement-validation signal,
not permission to run gauge-changing updates. Mutating full-update-style gauge
conditioning is deferred to S7b.
```

- [ ] **Step 2: Run a README smoke check**

Run:

```bash
rg -n "assess_ctm_trust|write_ctm_trust_csv|S7b|CTMValidationPoint" README.md
```

Expected: output includes all four terms.

- [ ] **Step 3: Commit README update**

```bash
git add README.md
git commit -m "docs: document CTM trust workflow"
```

## Task 7: Add S7b Gauge-Fixing Handoff Note

**Files:**
- Create: `docs/superpowers/notes/2026-05-16-s7b-gauge-fixing-handoff.md`

- [ ] **Step 1: Create the S7b handoff note**

Create `docs/superpowers/notes/2026-05-16-s7b-gauge-fixing-handoff.md`:

```markdown
# S7b Gauge-Fixing Handoff

Date: 2026-05-16

## Prerequisites From S7a

- `assess_ctm_trust(points)` must be available for finite-chi measurement
  validation.
- `gauge_diagnostic_simple(psi, c, dir)` must be available for read-only local
  simple-gauge diagnostics.
- CTM contexts must reject stale state usage through state id/version checks.
- S7a trust remains a measurement-validation signal. It is not sufficient by
  itself to authorize mutating gauge updates.

## Local CTM Norm-Matrix Requirements

Future `fix_bond_gauge!` work must define and test environment-backed local
bond norm matrices before mutating tensors. Required diagnostics:

- finite entries,
- Hermiticity residual,
- positive-semidefinite eigenvalue floor,
- condition number or reciprocal condition number,
- bond-direction coverage over canonical `:right` and `:up` links,
- clear failure behavior when the local environment is singular or indefinite.

## Mutation Contract

Future gauge-fixing code must:

- accept a fresh `PEPSKitMeasurementContext` or build one explicitly,
- verify measurement trust and local norm-matrix quality before mutation,
- perform all local factorization steps before writing into `psi`,
- increment `state_version(psi)` after mutation,
- invalidate old CTM contexts by relying on the existing state-version guard,
- preserve Gamma-lambda invariants or document any representation change.

## Required Gauge-Invariant Tests

S7b tests must compare observables and CTM summaries, not raw tensor entries.
Required checks:

- D=1 product states remain unchanged under any no-op gauge path.
- D=2 seeded states keep finite simple observables after gauge conditioning.
- Fresh CTM context measurements before/after a pure gauge change agree within
  documented tolerance.
- Old CTM contexts throw after gauge mutation.
- Singular or ill-conditioned norm matrices fail without partially mutating
  the state.

## Known Blockers For Mutating `fix_bond_gauge!`

- S7a does not compute CTM local bond norm matrices.
- S7a does not define whitening or ALS/full-update truncation factors.
- PEPSKit environment internals needed for local norm matrices still need a
  focused feasibility check.
- The project has not chosen whether `fix_bond_gauge!` should update only the
  current Gamma tensors or introduce a richer gauge state.
```

- [ ] **Step 2: Commit the handoff note**

```bash
git add docs/superpowers/notes/2026-05-16-s7b-gauge-fixing-handoff.md
git commit -m "docs: add S7b gauge-fixing handoff"
```

## Task 8: Run Focused Verification

**Files:**
- No source edits expected.

- [ ] **Step 1: Run focused new tests**

Run:

```bash
julia --project=. test/runtests.jl test_ctm_trust.jl test_gauge_diagnostics.jl test_public_docs.jl
```

Expected: PASS.

- [ ] **Step 2: Run CTM adapter regression tests**

Run:

```bash
julia --project=. test/runtests.jl test_pepskit_measurements.jl
```

Expected: PASS.

- [ ] **Step 3: Run full package tests**

Run:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: PASS.

- [ ] **Step 4: Run extended CTM tests manually**

Run:

```bash
SQUAREPXP_EXTENDED_TESTS=1 julia --project=. test/runtests.jl test_pepskit_measurements.jl
```

Expected: PASS. If this is too slow for local iteration, capture the runtime and the last successful focused test output in the final handoff instead of claiming extended verification passed.

- [ ] **Step 5: Run whitespace verification**

Run:

```bash
git diff --check
```

Expected: no output.

## Task 9: Final Integration Commit And Handoff

**Files:**
- Modify if needed: `memory/mid_term/decision_log.md`

- [ ] **Step 1: Inspect the final diff**

Run:

```bash
git status --short
git diff --stat
```

Expected: only S7a implementation, tests, README, and S7b handoff files are changed. Do not include unrelated local files such as experimental TFIM scripts unless the user explicitly asks.

- [ ] **Step 2: Record any new S7a implementation decision**

If implementation discovers a meaningful decision beyond the approved spec, append it to `memory/mid_term/decision_log.md` with source references. Use this shape:

```markdown
## 2026-05-16 - S7a Uses Separate CTM Trust And Gauge Diagnostics Modules

Decision:

Keep CTM finite-chi trust policy in `src/CTMTrust.jl` and read-only simple-gauge
diagnostics in `src/GaugeDiagnostics.jl`, with no mutating gauge-fixing API in
S7a.

Reason:

CTM measurement trust and gauge-update readiness have different correctness
requirements. S7a validates measurement stability and local simple-gauge
diagnostics, while S7b must separately validate environment norm matrices
before mutating tensors.

Source:

`docs/superpowers/specs/2026-05-16-s7-ctm-trust-gauge-readiness-design.md`;
`src/CTMTrust.jl`; `src/GaugeDiagnostics.jl`

Status: active
```

- [ ] **Step 3: Commit final integration**

If there are remaining S7a changes after the earlier task commits, run:

```bash
git add src/CTMTrust.jl src/GaugeDiagnostics.jl src/SquarePXPDynamics.jl \
    test/test_ctm_trust.jl test/test_gauge_diagnostics.jl \
    test/test_pepskit_measurements.jl test/test_public_docs.jl test/runtests.jl \
    README.md docs/superpowers/notes/2026-05-16-s7b-gauge-fixing-handoff.md \
    memory/mid_term/decision_log.md
git commit -m "feat: complete S7 CTM trust readiness"
```

Expected: commit succeeds. If every task already committed its changes and the tree has no remaining S7a edits, skip this commit and record that no final integration commit was needed.

- [ ] **Step 4: Final handoff summary**

The final response should report:

```text
Implemented S7a CTM trust and gauge readiness.

Verification:
- test_ctm_trust.jl: PASS
- test_gauge_diagnostics.jl: PASS
- test_public_docs.jl: PASS
- test_pepskit_measurements.jl: PASS
- full Pkg.test(): PASS
- extended CTM tests: PASS or not run, with reason

Notes:
- S7a does not implement mutating fix_bond_gauge!.
- S7b handoff is in docs/superpowers/notes/2026-05-16-s7b-gauge-fixing-handoff.md.
```
