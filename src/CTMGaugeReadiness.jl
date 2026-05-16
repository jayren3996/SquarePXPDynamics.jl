module CTMGaugeReadinessModule

using LinearAlgebra
using PEPSKit
using TensorKit

import ..PEPSKitMeasurements
using ..SquareGeometry: SquareCoord
using ..SquareUnitCells: BondKey, bondkey
using ..SquareIPEPS: SquareIPEPSState, link_weight
using ..PEPSKitMeasurements: PEPSKitMeasurementContext, assert_fresh_pepskit_context
using ..CTMTrust: CTMTrustAssessment

export CTMGaugePolicy, CTMBondNormDiagnostic, CTMGaugeReadiness, BondGaugeFixInfo
export ctm_bond_norm_matrix, ctm_bond_norm_diagnostic
export all_ctm_bond_norm_diagnostics, ctm_ready_for_gauge_updates
export fix_bond_gauge!

const _DIAG_REASONS = (
    :accepted,
    :nonfinite_entries,
    :zero_norm,
    :nonhermitian,
    :indefinite,
    :ill_conditioned,
)

const _READINESS_REASONS = (
    :ready,
    :stale_context,
    :untrusted_ctm,
    :missing_bond_norm,
    :bad_bond_norm,
)

"""
    CTMGaugePolicy(; max_hermiticity_residual=1e-8,
                   min_psd_eigenvalue=-1e-10,
                   min_rcond=1e-12,
                   require_all_bonds=true)

Thresholds for deciding whether CTM local bond norm matrices are suitable for
gauge-changing updates. `require_all_bonds` asks readiness checks to see every
canonical bond in the iPEPS unit cell before approving mutation.
"""
struct CTMGaugePolicy
    max_hermiticity_residual::Float64
    min_psd_eigenvalue::Float64
    min_rcond::Float64
    require_all_bonds::Bool

    function CTMGaugePolicy(
        max_hermiticity_residual::Real,
        min_psd_eigenvalue::Real,
        min_rcond::Real,
        require_all_bonds::Bool,
    )
        herm = Float64(max_hermiticity_residual)
        floor = Float64(min_psd_eigenvalue)
        rcond = Float64(min_rcond)
        isfinite(herm) && herm >= 0 ||
            throw(ArgumentError("max_hermiticity_residual must be finite and nonnegative"))
        isfinite(floor) ||
            throw(ArgumentError("min_psd_eigenvalue must be finite"))
        isfinite(rcond) && 0 <= rcond <= 1 ||
            throw(ArgumentError("min_rcond must be finite and in [0, 1]"))
        return new(herm, floor, rcond, require_all_bonds)
    end
end

function CTMGaugePolicy(;
    max_hermiticity_residual::Real = 1e-8,
    min_psd_eigenvalue::Real = -1e-10,
    min_rcond::Real = 1e-12,
    require_all_bonds::Bool = true,
)
    return CTMGaugePolicy(
        max_hermiticity_residual,
        min_psd_eigenvalue,
        min_rcond,
        require_all_bonds,
    )
end

"""
    CTMBondNormDiagnostic(bond, matrix; policy=CTMGaugePolicy())

Diagnostics for one CTM-backed local bond norm matrix. `matrix` is stored after
global phase alignment and Hermitian symmetrization; `hermiticity_residual`
records the relative residual before symmetrization. `accepted` is true only
when the matrix is finite, approximately Hermitian, positive semidefinite
within policy tolerance, and well-conditioned enough for gauge-changing work.
"""
struct CTMBondNormDiagnostic
    bond::BondKey
    matrix::Matrix{ComplexF64}
    hermiticity_residual::Float64
    frobenius_norm::Float64
    eigen_min::Float64
    eigen_max::Float64
    rcond::Float64
    accepted::Bool
    reject_reason::Union{Nothing,Symbol}

    function CTMBondNormDiagnostic(
        bond::BondKey,
        matrix::AbstractMatrix;
        policy::CTMGaugePolicy = CTMGaugePolicy(),
    )
        size(matrix, 1) == size(matrix, 2) ||
            throw(ArgumentError("CTM bond norm matrix must be square"))
        size(matrix, 1) >= 1 ||
            throw(ArgumentError("CTM bond norm matrix must be nonempty"))
        raw = Matrix{ComplexF64}(matrix)
        if !all(_finite_complex, raw)
            return new(
                bond,
                raw,
                Inf,
                Inf,
                NaN,
                NaN,
                0.0,
                false,
                :nonfinite_entries,
            )
        end

        phased = _phase_align(raw)
        frob = Float64(norm(phased))
        if !(isfinite(frob) && frob > 0)
            return new(bond, phased, Inf, frob, NaN, NaN, 0.0, false, :zero_norm)
        end

        hermiticity_residual = Float64(norm(phased - adjoint(phased)) / frob)
        sym = Matrix{ComplexF64}((phased + adjoint(phased)) / 2)
        eigs = eigvals(Hermitian(sym))
        eigen_min = Float64(minimum(eigs))
        eigen_max = Float64(maximum(eigs))
        reciprocal_condition = eigen_max > 0 ? max(eigen_min, 0.0) / eigen_max : 0.0
        reason = _reject_reason(
            hermiticity_residual,
            eigen_min,
            reciprocal_condition,
            policy,
        )
        return new(
            bond,
            sym,
            hermiticity_residual,
            frob,
            eigen_min,
            eigen_max,
            reciprocal_condition,
            reason === nothing,
            reason,
        )
    end
end

"""
    CTMGaugeReadiness

Structured readiness result combining CTM context freshness, finite-chi CTM
trust, and local CTM bond norm diagnostics.
"""
struct CTMGaugeReadiness
    ready::Bool
    reason::Symbol
    message::String
    trust::CTMTrustAssessment
    diagnostics::Dict{BondKey,CTMBondNormDiagnostic}

    function CTMGaugeReadiness(
        ready::Bool,
        reason::Symbol,
        message::AbstractString,
        trust::CTMTrustAssessment,
        diagnostics::Dict{BondKey,CTMBondNormDiagnostic},
    )
        reason in _READINESS_REASONS ||
            throw(ArgumentError("invalid CTM gauge readiness reason $reason"))
        ready == (reason === :ready) ||
            throw(ArgumentError("ready readiness results must use reason :ready"))
        return new(ready, reason, String(message), trust, diagnostics)
    end
end

"""
    BondGaugeFixInfo

Result returned by [`fix_bond_gauge!`](@ref). Slice 4 supports a transactional
D=1 product/no-op path and reports D>1 as an explicit nonmutation.
"""
struct BondGaugeFixInfo
    bond::BondKey
    mutated::Bool
    reason::Symbol
    readiness::CTMGaugeReadiness
end

_finite_complex(z) = isfinite(real(z)) && isfinite(imag(z))

function _phase_align(matrix::Matrix{ComplexF64})
    trm = tr(matrix)
    if _finite_complex(trm) && abs(trm) > 0
        return matrix .* cis(-angle(trm))
    else
        return matrix
    end
end

function _reject_reason(
    hermiticity_residual::Float64,
    eigen_min::Float64,
    rcond::Float64,
    policy::CTMGaugePolicy,
)
    hermiticity_residual <= policy.max_hermiticity_residual || return :nonhermitian
    eigen_min >= policy.min_psd_eigenvalue || return :indefinite
    rcond >= policy.min_rcond || return :ill_conditioned
    return nothing
end

function _matrix_from_bondenv(benv)
    rows = TensorKit.dim(TensorKit.codomain(benv))
    cols = TensorKit.dim(TensorKit.domain(benv))
    return reshape(ComplexF64.(collect(benv.data)), rows, cols)
end

function _horizontal_bondenv_matrix(peps, env, row::Int, col::Int)
    nc = size(peps, 2)
    next_col = mod1(col + 1, nc)
    X, _, _, Y = PEPSKit._qr_bond(peps[row, col], peps[row, next_col])
    return _matrix_from_bondenv(PEPSKit.bondenv_fu(row, col, X, Y, env))
end

function _ctm_bond_norm_raw_matrix(
    psi::SquareIPEPSState,
    bond::BondKey,
    ctx::PEPSKitMeasurementContext,
)
    site = PEPSKitMeasurements._squarecoord_to_cartesianindex(psi.unitcell, bond.site)
    row, col = Tuple(site)
    if bond.dir === :right
        return _horizontal_bondenv_matrix(ctx.peps, ctx.env, row, col)
    elseif bond.dir === :up
        rotated_peps = rotr90(ctx.peps)
        rotated_env = rotr90(ctx.env)
        rotated_site = PEPSKit.siterotr90(site, size(ctx.peps))
        rotated_row, rotated_col = Tuple(rotated_site)
        return _horizontal_bondenv_matrix(
            rotated_peps,
            rotated_env,
            rotated_row,
            rotated_col,
        )
    else
        throw(ArgumentError("BondKey direction must be :right or :up"))
    end
end

"""
    ctm_bond_norm_diagnostic(psi, c, dir, ctx; policy=CTMGaugePolicy())

Compute CTM-backed bond norm diagnostics for the canonical nearest-neighbor
bond selected by `(c, dir)`. The CTM context must be fresh for `psi`.
"""
function ctm_bond_norm_diagnostic(
    psi::SquareIPEPSState,
    c::SquareCoord,
    dir::Symbol,
    ctx::PEPSKitMeasurementContext;
    policy::CTMGaugePolicy = CTMGaugePolicy(),
)::CTMBondNormDiagnostic
    assert_fresh_pepskit_context(psi, ctx)
    bond = bondkey(psi.unitcell, c, dir)
    raw = _ctm_bond_norm_raw_matrix(psi, bond, ctx)
    return CTMBondNormDiagnostic(bond, raw; policy)
end

"""
    ctm_bond_norm_matrix(psi, c, dir, ctx; policy=CTMGaugePolicy())

Return the phase-aligned, Hermitian CTM local bond norm matrix used in
[`ctm_bond_norm_diagnostic`](@ref).
"""
function ctm_bond_norm_matrix(
    psi::SquareIPEPSState,
    c::SquareCoord,
    dir::Symbol,
    ctx::PEPSKitMeasurementContext;
    policy::CTMGaugePolicy = CTMGaugePolicy(),
)::Matrix{ComplexF64}
    return ctm_bond_norm_diagnostic(psi, c, dir, ctx; policy).matrix
end

"""
    all_ctm_bond_norm_diagnostics(psi, ctx; policy=CTMGaugePolicy())

Compute CTM bond norm diagnostics for every canonical link-weight bond in the
iPEPS unit cell.
"""
function all_ctm_bond_norm_diagnostics(
    psi::SquareIPEPSState,
    ctx::PEPSKitMeasurementContext;
    policy::CTMGaugePolicy = CTMGaugePolicy(),
)::Dict{BondKey,CTMBondNormDiagnostic}
    assert_fresh_pepskit_context(psi, ctx)
    return Dict{BondKey,CTMBondNormDiagnostic}(
        bond => ctm_bond_norm_diagnostic(psi, bond.site, bond.dir, ctx; policy) for
        bond in keys(psi.link_weights)
    )
end

function _collect_diagnostics(
    diagnostics::AbstractDict,
)::Dict{BondKey,CTMBondNormDiagnostic}
    collected = Dict{BondKey,CTMBondNormDiagnostic}()
    for (bond, diagnostic) in diagnostics
        bond isa BondKey ||
            throw(ArgumentError("CTM gauge diagnostic keys must be BondKey values"))
        diagnostic isa CTMBondNormDiagnostic ||
            throw(ArgumentError("CTM gauge diagnostic values must be CTMBondNormDiagnostic"))
        collected[bond] = diagnostic
    end
    return collected
end

function _readiness(
    reason::Symbol,
    message::String,
    trust::CTMTrustAssessment,
    diagnostics::Dict{BondKey,CTMBondNormDiagnostic},
)
    return CTMGaugeReadiness(reason === :ready, reason, message, trust, diagnostics)
end

"""
    ctm_ready_for_gauge_updates(psi, ctx, trust; diagnostics=nothing,
                                policy=CTMGaugePolicy())

Return a structured readiness decision for gauge-changing updates. Readiness
requires a fresh CTM context, trusted finite-`chi` assessment, required bond
coverage, and accepted local CTM bond norm diagnostics.
"""
function ctm_ready_for_gauge_updates(
    psi::SquareIPEPSState,
    ctx::PEPSKitMeasurementContext,
    trust::CTMTrustAssessment;
    diagnostics = nothing,
    policy::CTMGaugePolicy = CTMGaugePolicy(),
)::CTMGaugeReadiness
    try
        assert_fresh_pepskit_context(psi, ctx)
    catch err
        err isa ArgumentError || rethrow()
        return _readiness(
            :stale_context,
            sprint(showerror, err),
            trust,
            Dict{BondKey,CTMBondNormDiagnostic}(),
        )
    end

    if !trust.trusted
        return _readiness(
            :untrusted_ctm,
            trust.message,
            trust,
            Dict{BondKey,CTMBondNormDiagnostic}(),
        )
    end

    collected = diagnostics === nothing ?
        all_ctm_bond_norm_diagnostics(psi, ctx; policy) :
        _collect_diagnostics(diagnostics)

    if policy.require_all_bonds
        required = Set(keys(psi.link_weights))
        observed = Set(keys(collected))
        if !issubset(required, observed)
            return _readiness(
                :missing_bond_norm,
                "CTM bond norm diagnostics do not cover every canonical unit-cell bond",
                trust,
                collected,
            )
        end
    end

    for diagnostic in values(collected)
        if !diagnostic.accepted
            return _readiness(
                :bad_bond_norm,
                "CTM bond norm diagnostic $(diagnostic.bond) rejected: $(diagnostic.reject_reason)",
                trust,
                collected,
            )
        end
    end

    return _readiness(
        :ready,
        "CTM context, trust assessment, and local bond norm diagnostics are ready",
        trust,
        collected,
    )
end

"""
    fix_bond_gauge!(psi, c, dir, ctx, trust; diagnostics=nothing,
                    policy=CTMGaugePolicy())

Transactional S7b gauge-fix entry point. Slice 4 only performs readiness
checks and the D=1 product/no-op path; D>1 mutation is reported explicitly for
the next gauge-conditioning slice.
"""
function fix_bond_gauge!(
    psi::SquareIPEPSState,
    c::SquareCoord,
    dir::Symbol,
    ctx::PEPSKitMeasurementContext,
    trust::CTMTrustAssessment;
    diagnostics = nothing,
    policy::CTMGaugePolicy = CTMGaugePolicy(),
)::BondGaugeFixInfo
    bond = bondkey(psi.unitcell, c, dir)
    readiness = ctm_ready_for_gauge_updates(
        psi,
        ctx,
        trust;
        diagnostics,
        policy,
    )
    readiness.ready || return BondGaugeFixInfo(bond, false, readiness.reason, readiness)
    if length(link_weight(psi, bond.site, bond.dir)) == 1
        return BondGaugeFixInfo(bond, false, :product_noop, readiness)
    end
    return BondGaugeFixInfo(bond, false, :d_greater_than_one_not_implemented, readiness)
end

end
