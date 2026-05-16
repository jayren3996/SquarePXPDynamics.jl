module CTMTrust

using ..PEPSKitMeasurements: CTMRGDiagnostics, CTMObservableSummary, CTMValidationPoint

export CTMTrustPolicy, CTMTrustAssessment, assess_ctm_trust, write_ctm_trust_csv

const _TRUST_REASONS = (
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

const _TRUST_CSV_HEADER = (
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

function _finite_nonnegative(value::Real, label::String)
    isfinite(value) && value >= 0 || throw(ArgumentError("$label must be finite and nonnegative"))
    return Float64(value)
end

function _optional_finite_nonnegative(value::Nothing, label::String)
    return nothing
end

function _optional_finite_nonnegative(value::Real, label::String)
    return _finite_nonnegative(value, label)
end

"""
    CTMTrustPolicy(min_points=2, require_accepted_diagnostics=true,
                   max_density_delta=1e-3, max_blockade_delta=1e-4,
                   max_energy_delta=1e-3, max_residual=nothing)

Thresholds used by [`assess_ctm_trust`](@ref) to decide whether the final
finite-`chi` CTMRG validation window is stable enough to trust as a measurement
signal. The policy validates that at least two points are compared and that all
numeric thresholds are finite and nonnegative.
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
        return new(
            Int(min_points),
            require_accepted_diagnostics,
            _finite_nonnegative(max_density_delta, "max_density_delta"),
            _finite_nonnegative(max_blockade_delta, "max_blockade_delta"),
            _finite_nonnegative(max_energy_delta, "max_energy_delta"),
            _optional_finite_nonnegative(max_residual, "max_residual"),
        )
    end
end

CTMTrustPolicy() = CTMTrustPolicy(2, true, 1e-3, 1e-4, 1e-3, nothing)

"""
    CTMTrustAssessment

Result returned by [`assess_ctm_trust`](@ref). `trusted` is true only when
`reason == :trusted`; otherwise `reason` records the first failed policy check,
the finite-`chi` drift fields report maximum adjacent CTM-to-CTM observable
changes when available, and `observed_max_residual` reports the largest
available CTMRG residual in the assessed records.
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
        reason in _TRUST_REASONS || throw(ArgumentError("invalid CTM trust reason $reason"))
        trusted == (reason === :trusted) ||
            throw(ArgumentError("trusted assessments must use reason :trusted"))
        compared_points >= 0 || throw(ArgumentError("compared_points must be nonnegative"))
        return new(
            trusted,
            reason,
            String(message),
            Int(compared_points),
            _optional_finite_nonnegative(finite_chi_density_delta, "finite_chi_density_delta"),
            _optional_finite_nonnegative(
                finite_chi_blockade_delta,
                "finite_chi_blockade_delta",
            ),
            _optional_finite_nonnegative(finite_chi_energy_delta, "finite_chi_energy_delta"),
            _optional_finite_nonnegative(observed_max_residual, "observed_max_residual"),
        )
    end
end

function _diagnostic_field(::Nothing, ::Symbol)
    return nothing
end

function _diagnostic_field(diagnostics::CTMRGDiagnostics, field::Symbol)
    return getfield(diagnostics, field)
end

function _assert_finite_summary(summary::CTMObservableSummary, label::String)
    all(
        isfinite,
        (
            summary.density,
            summary.density_even,
            summary.density_odd,
            summary.blockade_violation,
            summary.pxp_energy_density,
            summary.sublattice_imbalance,
            summary.checkerboard_structure_factor,
        ),
    ) || throw(ArgumentError("$label CTM summary fields must be finite"))
    return summary
end

function _collect_points(points)
    collected = collect(points)
    for (idx, point) in pairs(collected)
        point isa CTMValidationPoint ||
            throw(ArgumentError("CTM trust point $idx is not a CTMValidationPoint"))
        _assert_finite_summary(point.reference, "reference")
        _assert_finite_summary(point.measurement, "measurement")
    end
    return collected
end

function _observed_max_residual(points)
    residuals = Float64[]
    for point in points
        residual = _diagnostic_field(point.diagnostics, :residual)
        residual === nothing && continue
        push!(residuals, residual)
    end
    return isempty(residuals) ? nothing : maximum(residuals)
end

function _finite_chi_deltas(points)
    length(points) >= 2 || return (nothing, nothing, nothing)
    density_delta = 0.0
    blockade_delta = 0.0
    energy_delta = 0.0
    for idx in 2:length(points)
        previous = points[idx - 1].measurement
        current = points[idx].measurement
        density_delta = max(density_delta, abs(current.density - previous.density))
        blockade_delta = max(
            blockade_delta,
            abs(current.blockade_violation - previous.blockade_violation),
        )
        energy_delta = max(
            energy_delta,
            abs(current.pxp_energy_density - previous.pxp_energy_density),
        )
    end
    return (density_delta, blockade_delta, energy_delta)
end

function _assessment(
    reason::Symbol,
    message::String,
    compared_points::Int,
    density_delta,
    blockade_delta,
    energy_delta,
    observed_max_residual,
)
    return CTMTrustAssessment(
        reason === :trusted,
        reason,
        message,
        compared_points,
        density_delta,
        blockade_delta,
        energy_delta,
        observed_max_residual,
    )
end

"""
    assess_ctm_trust(points; policy=CTMTrustPolicy())

Assess a CTMRG validation sweep as a measurement-validation signal. Only the
final `policy.min_points` records are compared after checking that their
environment bond dimensions strictly increase and tolerances do not increase.
Finite-`chi` drifts are maximum adjacent CTM-to-CTM changes in density,
blockade violation, and PXP energy density; they do not use reference deltas.
"""
function assess_ctm_trust(points; policy::CTMTrustPolicy = CTMTrustPolicy())
    collected = _collect_points(points)
    compared_points = min(length(collected), policy.min_points)

    if length(collected) < policy.min_points
        return _assessment(
            :too_few_points,
            "fewer CTM validation points than required by policy",
            length(collected),
            nothing,
            nothing,
            nothing,
            _observed_max_residual(collected),
        )
    end

    window = collected[(end - policy.min_points + 1):end]
    observed_max_residual = _observed_max_residual(window)
    density_delta, blockade_delta, energy_delta = _finite_chi_deltas(window)

    for idx in 2:length(window)
        previous = window[idx - 1]
        current = window[idx]
        if current.params.chi <= previous.params.chi || current.params.tol > previous.params.tol
            return _assessment(
                :nonmonotonic_sweep,
                "final CTM validation window must have strictly increasing chi and nonincreasing tol",
                compared_points,
                density_delta,
                blockade_delta,
                energy_delta,
                observed_max_residual,
            )
        end
    end

    if policy.require_accepted_diagnostics
        any(point -> point.diagnostics === nothing, window) && return _assessment(
            :missing_diagnostics,
            "final CTM validation window is missing diagnostics",
            compared_points,
            density_delta,
            blockade_delta,
            energy_delta,
            observed_max_residual,
        )
        any(point -> !point.diagnostics.accepted, window) && return _assessment(
            :unaccepted_diagnostics,
            "final CTM validation window contains unaccepted diagnostics",
            compared_points,
            density_delta,
            blockade_delta,
            energy_delta,
            observed_max_residual,
        )
    end

    if policy.max_residual !== nothing
        any(point -> _diagnostic_field(point.diagnostics, :residual) === nothing, window) &&
            return _assessment(
                :missing_residual,
                "final CTM validation window is missing residual diagnostics",
                compared_points,
                density_delta,
                blockade_delta,
                energy_delta,
                observed_max_residual,
            )
        if observed_max_residual > policy.max_residual
            return _assessment(
                :residual_too_large,
                "final CTM validation window residual exceeds policy threshold",
                compared_points,
                density_delta,
                blockade_delta,
                energy_delta,
                observed_max_residual,
            )
        end
    end

    density_delta > policy.max_density_delta && return _assessment(
        :density_delta_too_large,
        "final CTM validation window density drift exceeds policy threshold",
        compared_points,
        density_delta,
        blockade_delta,
        energy_delta,
        observed_max_residual,
    )
    blockade_delta > policy.max_blockade_delta && return _assessment(
        :blockade_delta_too_large,
        "final CTM validation window blockade drift exceeds policy threshold",
        compared_points,
        density_delta,
        blockade_delta,
        energy_delta,
        observed_max_residual,
    )
    energy_delta > policy.max_energy_delta && return _assessment(
        :energy_delta_too_large,
        "final CTM validation window energy drift exceeds policy threshold",
        compared_points,
        density_delta,
        blockade_delta,
        energy_delta,
        observed_max_residual,
    )

    return _assessment(
        :trusted,
        "final CTM validation window satisfies trust policy",
        compared_points,
        density_delta,
        blockade_delta,
        energy_delta,
        observed_max_residual,
    )
end

function _csv_value(::Nothing)
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

"""
    write_ctm_trust_csv(points, path; policy=CTMTrustPolicy())

Write CTM validation points plus repeated trust policy and assessment fields to
`path` as CSV. The trust assessment is computed internally from `points`; the
free-form assessment message is intentionally omitted so the CSV remains stable
for audit and regression use.
"""
function write_ctm_trust_csv(
    points,
    path::AbstractString;
    policy::CTMTrustPolicy = CTMTrustPolicy(),
)
    collected = _collect_points(points)
    assessment = assess_ctm_trust(collected; policy)
    open(path, "w") do io
        println(io, join(_TRUST_CSV_HEADER, ","))
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
            println(io, join(_csv_value.(row), ","))
        end
    end
    return path
end

end
