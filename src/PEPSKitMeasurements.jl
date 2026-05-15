module PEPSKitMeasurements

using ITensors
using PEPSKit
using TensorKit
using ..SquareGeometry: SquareCoord
using ..SquarePXP: SQUARE_STAR_SITES, square_pxp_star_hamiltonian
using ..SquareUnitCells: PeriodicSquareUnitCell, wrap, neighbor
using ..SquareIPEPS:
    SquareIPEPSState, physical_index, link_index, link_weight, state_version

export PEPSKitCTMRGParams, PEPSKitMeasurementContext, CTMRGDiagnostics
export CTMObservableSummary
export to_pepskit_infinitepeps, pepskit_ctmrg_context
export local_density_ctm, nearest_neighbor_density_ctm, blockade_violation_ctm
export star_expectation_ctm, pxp_energy_density_ctm, measure_ctm, ctm_diagnostics

const _DIRECTIONS = (:right, :up, :left, :down)

"""
    PEPSKitCTMRGParams(chi, tol, maxiter, verbosity)

Validated CTMRG controls for the PEPSKit measurement adapter. `chi` is the
environment bond dimension, `tol` and `maxiter` are forwarded to
`PEPSKit.leading_boundary`, and `verbosity` follows PEPSKit's CTMRG logging
levels. This struct configures measurements only; it is not used by the custom
ITensors simple-update dynamics.
"""
struct PEPSKitCTMRGParams
    chi::Int
    tol::Float64
    maxiter::Int
    verbosity::Int

    function PEPSKitCTMRGParams(
        chi::Integer,
        tol::Real,
        maxiter::Integer,
        verbosity::Integer,
    )
        chi >= 1 || throw(ArgumentError("chi must be at least 1"))
        isfinite(tol) && tol > 0 || throw(ArgumentError("tol must be finite and positive"))
        maxiter >= 1 || throw(ArgumentError("maxiter must be at least 1"))
        verbosity >= 0 || throw(ArgumentError("verbosity must be nonnegative"))
        return new(Int(chi), Float64(tol), Int(maxiter), Int(verbosity))
    end
end

"""
    CTMRGDiagnostics

Structured summary of CTMRG convergence metadata exposed by the measurement
adapter. `iterations`, `residual`, and `converged` are `nothing` when PEPSKit
does not expose a recognized field. `accepted` is the adapter's downstream
trust flag for ranking/logging.
"""
struct CTMRGDiagnostics
    chi::Int
    tol::Float64
    maxiter::Int
    iterations::Union{Int,Nothing}
    residual::Union{Float64,Nothing}
    converged::Union{Bool,Nothing}
    accepted::Bool

    function CTMRGDiagnostics(
        chi::Integer,
        tol::Real,
        maxiter::Integer,
        iterations::Union{Integer,Nothing},
        residual::Union{Real,Nothing},
        converged::Union{Bool,Nothing},
        accepted::Bool,
    )
        chi >= 1 || throw(ArgumentError("chi must be at least 1"))
        isfinite(tol) && tol > 0 || throw(ArgumentError("tol must be finite and positive"))
        maxiter >= 1 || throw(ArgumentError("maxiter must be at least 1"))
        if iterations !== nothing
            iterations >= 0 || throw(ArgumentError("iterations must be nonnegative"))
        end
        if residual !== nothing
            isfinite(residual) && residual >= 0 ||
                throw(ArgumentError("residual must be finite and nonnegative"))
        end
        return new(
            Int(chi),
            Float64(tol),
            Int(maxiter),
            iterations === nothing ? nothing : Int(iterations),
            residual === nothing ? nothing : Float64(residual),
            converged,
            accepted,
        )
    end
end

"""
    PEPSKitMeasurementContext

Reusable PEPSKit measurement cache containing the converted `InfinitePEPS`, the
CTMRG environment returned by PEPSKit, the raw CTMRG info object, and the
parameters used to construct it. Local operators are cached by unit-cell
coordinate and operator identity so repeated measurements reuse the same
PEPSKit `LocalOperator` without rerunning CTMRG.
"""
struct PEPSKitMeasurementContext
    peps::Any
    env::Any
    info::Any
    params::PEPSKitCTMRGParams
    diagnostics::CTMRGDiagnostics
    source_state_id::UInt
    source_state_version::Int
    operator_cache::Dict{Any,Any}
end

PEPSKitMeasurementContext(
    peps,
    env,
    info,
    params::PEPSKitCTMRGParams,
    diagnostics::CTMRGDiagnostics,
    source_state_id::UInt,
    source_state_version::Integer,
) = PEPSKitMeasurementContext(
    peps,
    env,
    info,
    params,
    diagnostics,
    source_state_id,
    Int(source_state_version),
    Dict{Any,Any}(),
)

"""
    CTMObservableSummary

Summary of PEPSKit CTMRG observables for a custom ITensors square iPEPS state:
total density, even/odd sublattice densities, nearest-neighbor blockade
violation, and the five-site PXP energy density.
"""
struct CTMObservableSummary
    density::Float64
    density_even::Float64
    density_odd::Float64
    blockade_violation::Float64
    pxp_energy_density::Float64
    diagnostics::Union{CTMRGDiagnostics,Nothing}
end

CTMObservableSummary(
    density::Float64,
    density_even::Float64,
    density_odd::Float64,
    blockade_violation::Float64,
    pxp_energy_density::Float64,
    diagnostics::CTMRGDiagnostics,
) = invoke(
    CTMObservableSummary,
    Tuple{
        Float64,
        Float64,
        Float64,
        Float64,
        Float64,
        Union{CTMRGDiagnostics,Nothing},
    },
    density,
    density_even,
    density_odd,
    blockade_violation,
    pxp_energy_density,
    diagnostics,
)

CTMObservableSummary(
    density::Real,
    density_even::Real,
    density_odd::Real,
    blockade_violation::Real,
    pxp_energy_density::Real,
) = CTMObservableSummary(
    Float64(density),
    Float64(density_even),
    Float64(density_odd),
    Float64(blockade_violation),
    Float64(pxp_energy_density),
    nothing,
)

CTMObservableSummary(
    density::Real,
    density_even::Real,
    density_odd::Real,
    blockade_violation::Real,
    pxp_energy_density::Real,
    diagnostics::CTMRGDiagnostics,
) = CTMObservableSummary(
    Float64(density),
    Float64(density_even),
    Float64(density_odd),
    Float64(blockade_violation),
    Float64(pxp_energy_density),
    diagnostics,
)

"""
    ctm_diagnostics(ctx)

Return the structured CTMRG diagnostics attached to a PEPSKit measurement
context.
"""
ctm_diagnostics(ctx::PEPSKitMeasurementContext)::CTMRGDiagnostics = ctx.diagnostics

function _maybe_info_value(info, names)
    for name in names
        if info isa NamedTuple && haskey(info, name)
            return getfield(info, name)
        elseif hasproperty(info, name)
            return getproperty(info, name)
        elseif info isa AbstractDict && haskey(info, name)
            return info[name]
        elseif info isa AbstractDict && haskey(info, String(name))
            return info[String(name)]
        end
    end
    return nothing
end

function _maybe_int(value)
    value === nothing && return nothing
    value isa Integer || return nothing
    return Int(value)
end

function _maybe_float(value)
    value === nothing && return nothing
    value isa Real || return nothing
    converted = Float64(value)
    isfinite(converted) && converted >= 0 || return nothing
    return converted
end

function _maybe_bool(value)
    value === nothing && return nothing
    value isa Bool || return nothing
    return value
end

function _ctmrg_diagnostics(params::PEPSKitCTMRGParams, info)::CTMRGDiagnostics
    iterations = _maybe_int(
        _maybe_info_value(info, (:iterations, :iteration, :iter, :niter, :numiter)),
    )
    residual = _maybe_float(
        _maybe_info_value(
            info,
            (:residual, :err, :error, :normres, :conv_error, :truncation_error),
        ),
    )
    converged = _maybe_bool(_maybe_info_value(info, (:converged, :conv, :isconverged)))
    accepted =
        converged === true ||
        (converged === nothing && residual !== nothing && residual <= params.tol)
    return CTMRGDiagnostics(
        params.chi,
        params.tol,
        params.maxiter,
        iterations,
        residual,
        converged,
        accepted,
    )
end

function _assert_fresh_context(psi::SquareIPEPSState, ctx::PEPSKitMeasurementContext)
    objectid(psi) == ctx.source_state_id ||
        throw(ArgumentError("PEPSKit measurement context belongs to a different iPEPS state"))
    state_version(psi) == ctx.source_state_version || throw(
        ArgumentError(
            "PEPSKit measurement context is stale because the iPEPS state was mutated",
        ),
    )
    return nothing
end

function _validate_direction(dir::Symbol)
    dir in _DIRECTIONS ||
        throw(ArgumentError("direction must be :right, :up, :left, or :down"))
    return dir
end

function _dense_index(values)
    idx = 1
    nsites = length(values)
    for (site, value) in enumerate(values)
        1 <= value <= 2 || throw(ArgumentError("basis values must be 1 or 2"))
        idx += (value - 1) * 2^(nsites - site)
    end
    return idx
end

function _real_expectation(value; atol = 1e-8)
    z = ComplexF64(value)
    isfinite(real(z)) && isfinite(imag(z)) ||
        throw(ArgumentError("expectation value must be finite"))
    abs(imag(z)) <= atol || throw(
        ArgumentError(
            "Hermitian CTMRG observable produced non-negligible imaginary part $(imag(z))",
        ),
    )
    return Float64(real(z))
end

function _physical_lattice(cell::PeriodicSquareUnitCell)
    return fill(TensorKit.ComplexSpace(2), cell.Ly, cell.Lx)
end

"""
    _squarecoord_to_cartesianindex(cell, c)

Map `SquareCoord(x, y)` to PEPSKit's matrix coordinate `CartesianIndex(row,
col)`. Columns follow `x`; rows are flipped as `row = Ly - y + 1`, so the
custom `:up` neighbor maps to PEPSKit north (`row - 1`) and `:right` maps to
PEPSKit east (`col + 1`).
"""
function _squarecoord_to_cartesianindex(cell::PeriodicSquareUnitCell, c::SquareCoord)
    site = wrap(cell, c)
    return CartesianIndex(cell.Ly - site.y + 1, site.x)
end

function _local_neighbor_cartesianindex(
    cell::PeriodicSquareUnitCell,
    c::SquareCoord,
    dir::Symbol,
)
    _validate_direction(dir)
    center = _squarecoord_to_cartesianindex(cell, c)
    dir === :right && return center + CartesianIndex(0, 1)
    dir === :up && return center + CartesianIndex(-1, 0)
    dir === :left && return center + CartesianIndex(0, -1)
    dir === :down && return center + CartesianIndex(1, 0)
    throw(ArgumentError("direction must be :right, :up, :left, or :down"))
end

function _coord_from_row_col(cell::PeriodicSquareUnitCell, row::Int, col::Int)
    return SquareCoord(col, cell.Ly - row + 1)
end

function _star_sites_cartesian(cell::PeriodicSquareUnitCell, center::SquareCoord)
    c = wrap(cell, center)
    wrapped_sites = (
        c,
        neighbor(cell, c, :right),
        neighbor(cell, c, :up),
        neighbor(cell, c, :left),
        neighbor(cell, c, :down),
    )
    length(Set(wrapped_sites)) == 5 || throw(
        ArgumentError(
            "wrapped square star must contain five distinct unit-cell representatives",
        ),
    )
    pepskit_center = _squarecoord_to_cartesianindex(cell, c)
    return (
        pepskit_center,
        pepskit_center + CartesianIndex(0, 1),
        pepskit_center + CartesianIndex(-1, 0),
        pepskit_center + CartesianIndex(0, -1),
        pepskit_center + CartesianIndex(1, 0),
    )
end

_pepskit_star_sites(cell::PeriodicSquareUnitCell, center::SquareCoord) =
    _star_sites_cartesian(cell, center)

function _dense_operator_tensor_map(O::AbstractMatrix, nsites::Int)
    size(O) == (2^nsites, 2^nsites) ||
        throw(ArgumentError("dense operator must be $(2^nsites)x$(2^nsites)"))
    p = TensorKit.ComplexSpace(2)
    data = zeros(ComplexF64, ntuple(Returns(2), 2nsites))
    for out_values in Iterators.product((1:2 for _ = 1:nsites)...)
        out_idx = _dense_index(out_values)
        for in_values in Iterators.product((1:2 for _ = 1:nsites)...)
            in_idx = _dense_index(in_values)
            data[out_values..., in_values...] = O[out_idx, in_idx]
        end
    end
    spaces = ntuple(_ -> p, nsites)
    product_space = reduce(⊗, spaces)
    return TensorKit.TensorMap(data, product_space ← product_space)
end

function _pepskit_pxp_star_tensormap()
    return _dense_operator_tensor_map(square_pxp_star_hamiltonian(), SQUARE_STAR_SITES)
end

function _pepskit_density_operator(cell::PeriodicSquareUnitCell, c::SquareCoord)
    n = ComplexF64[1 0; 0 0]
    op = _dense_operator_tensor_map(n, 1)
    return PEPSKit.LocalOperator(
        _physical_lattice(cell),
        (_squarecoord_to_cartesianindex(cell, c),) => op,
    )
end

function _pepskit_twosite_nn_operator(
    cell::PeriodicSquareUnitCell,
    c::SquareCoord,
    dir::Symbol,
)
    _validate_direction(dir)
    n = ComplexF64[1 0; 0 0]
    op = _dense_operator_tensor_map(kron(n, n), 2)
    sites = (
        _squarecoord_to_cartesianindex(cell, c),
        _local_neighbor_cartesianindex(cell, c, dir),
    )
    return PEPSKit.LocalOperator(_physical_lattice(cell), sites => op)
end

function _pepskit_pxp_star_operator(cell::PeriodicSquareUnitCell, center::SquareCoord)
    op = _pepskit_pxp_star_tensormap()
    return PEPSKit.LocalOperator(
        _physical_lattice(cell),
        _pepskit_star_sites(cell, center) => op,
    )
end

function _pepskit_pxp_energy_operator(cell::PeriodicSquareUnitCell)
    op = _pepskit_pxp_star_tensormap()
    terms = [_pepskit_star_sites(cell, c) => op for c in cell.reps]
    return PEPSKit.LocalOperator(_physical_lattice(cell), terms...)
end

function _require_dense_index(index::Index, label)
    ITensors.hasqns(index) && throw(
        ArgumentError(
            "PEPSKit measurement adapter supports only dense, QN-free ITensors indices; $label has QNs",
        ),
    )
    return index
end

function _site_spaces(psi::SquareIPEPSState, c::SquareCoord)
    p = TensorKit.ComplexSpace(2)
    north = TensorKit.ComplexSpace(ITensors.dim(link_index(psi, c, :up)))
    east = TensorKit.ComplexSpace(ITensors.dim(link_index(psi, c, :right)))
    south = TensorKit.ComplexSpace(ITensors.dim(link_index(psi, c, :down)))'
    west = TensorKit.ComplexSpace(ITensors.dim(link_index(psi, c, :left)))'
    return p, north, east, south, west
end

function _absorbed_site_array(psi::SquareIPEPSState, c::SquareCoord)
    T = psi.tensors[c]
    pidx = _require_dense_index(physical_index(psi, c), "physical index at $c")
    ITensors.dim(pidx) == 2 ||
        throw(ArgumentError("physical index at $c must have dimension 2"))

    nidx = _require_dense_index(link_index(psi, c, :up), "up link at $c")
    eidx = _require_dense_index(link_index(psi, c, :right), "right link at $c")
    sidx = _require_dense_index(link_index(psi, c, :down), "down link at $c")
    widx = _require_dense_index(link_index(psi, c, :left), "left link at $c")
    for idx in (pidx, nidx, eidx, sidx, widx)
        ITensors.hasind(T, idx) ||
            throw(ArgumentError("site tensor at $c is missing expected index $idx"))
    end

    data = Array(T, pidx, nidx, eidx, sidx, widx)
    lambdas = (
        sqrt.(link_weight(psi, c, :up)),
        sqrt.(link_weight(psi, c, :right)),
        sqrt.(link_weight(psi, c, :down)),
        sqrt.(link_weight(psi, c, :left)),
    )
    size(data) ==
    (2, length(lambdas[1]), length(lambdas[2]), length(lambdas[3]), length(lambdas[4])) ||
        throw(ArgumentError("site tensor at $c has incompatible PEPS leg dimensions"))

    for n in axes(data, 2), e in axes(data, 3), s in axes(data, 4), w in axes(data, 5)
        scale = lambdas[1][n] * lambdas[2][e] * lambdas[3][s] * lambdas[4][w]
        @views data[:, n, e, s, w] .*= scale
    end
    return ComplexF64.(data)
end

function _pepskit_site_tensor(psi::SquareIPEPSState, c::SquareCoord)
    p, north, east, south, west = _site_spaces(psi, c)
    data = _absorbed_site_array(psi, c)
    return TensorKit.TensorMap(data, p ← north ⊗ east ⊗ south ⊗ west)
end

"""
    to_pepskit_infinitepeps(psi::SquareIPEPSState)

Convert the custom ITensors Γ-λ `SquareIPEPSState` into a PEPSKit
`InfinitePEPS` for CTMRG measurement. This adapter preserves the physical basis
`1 = :up`/Rydberg and `2 = :down`/vacancy. `SquareCoord(x, y)` maps to PEPSKit
matrix coordinates `(row = Ly - y + 1, col = x)`, so `:up` is PEPSKit north and
`:right` is PEPSKit east. PEPSKit virtual legs are ordered `N,E,S,W`, mapped
from custom `up,right,down,left`. Because the source state is in Γ-λ
simple-update form, each λ spectrum is absorbed symmetrically as `sqrt(λ)` into
both endpoint tensors before constructing the ordinary PEPS tensors. The input
state is not mutated. Only dense, QN-free ITensors indices are supported.
"""
function to_pepskit_infinitepeps(psi::SquareIPEPSState)
    cell = psi.unitcell
    tensors = [
        _pepskit_site_tensor(psi, _coord_from_row_col(cell, row, col)) for
        row = 1:cell.Ly, col = 1:cell.Lx
    ]
    return PEPSKit.InfinitePEPS(tensors)
end

"""
    pepskit_ctmrg_context(psi; params = PEPSKitCTMRGParams(8, 1e-8, 100, 0))

Convert `psi` with [`to_pepskit_infinitepeps`](@ref), construct a PEPSKit
`CTMRGEnv`, run `PEPSKit.leading_boundary`, and return a
`PEPSKitMeasurementContext`. The raw PEPSKit convergence `info` object is
preserved in the context; this function does not require or hide perfect CTMRG
convergence.
"""
function pepskit_ctmrg_context(
    psi::SquareIPEPSState;
    params::PEPSKitCTMRGParams = PEPSKitCTMRGParams(8, 1e-8, 100, 0),
)::PEPSKitMeasurementContext
    peps = to_pepskit_infinitepeps(psi)
    chi_space = TensorKit.ComplexSpace(params.chi)
    env0 = PEPSKit.CTMRGEnv(randn, ComplexF64, peps, chi_space)
    env, info = PEPSKit.leading_boundary(
        env0,
        peps;
        alg = :simultaneous,
        tol = params.tol,
        miniter = 1,
        maxiter = params.maxiter,
        verbosity = params.verbosity,
        trunc = PEPSKit.truncrank(params.chi),
    )
    diagnostics = _ctmrg_diagnostics(params, info)
    return PEPSKitMeasurementContext(
        peps,
        env,
        info,
        params,
        diagnostics,
        objectid(psi),
        state_version(psi),
    )
end

function _expectation(ctx::PEPSKitMeasurementContext, op)
    return _real_expectation(PEPSKit.expectation_value(ctx.peps, op, ctx.env))
end

function _operator_cache_key(
    tag::Symbol,
    cell::PeriodicSquareUnitCell,
    center::SquareCoord,
    O::AbstractMatrix,
)
    c = wrap(cell, center)
    return (tag, cell.Lx, cell.Ly, c.x, c.y, size(O), objectid(O))
end

function _pepskit_star_localoperator(
    cell::PeriodicSquareUnitCell,
    center::SquareCoord,
    O::AbstractMatrix,
)
    size(O) == (2^SQUARE_STAR_SITES, 2^SQUARE_STAR_SITES) ||
        throw(ArgumentError("dense square-star operator must be 32x32"))
    op = _dense_operator_tensor_map(O, SQUARE_STAR_SITES)
    return PEPSKit.LocalOperator(
        _physical_lattice(cell),
        _pepskit_star_sites(cell, center) => op,
    )
end

function _cached_star_localoperator(
    psi::SquareIPEPSState,
    center::SquareCoord,
    O::AbstractMatrix,
    ctx::PEPSKitMeasurementContext,
)
    key = _operator_cache_key(:star, psi.unitcell, center, O)
    return get!(ctx.operator_cache, key) do
        _pepskit_star_localoperator(psi.unitcell, center, O)
    end
end

function _cached_pxp_energy_operator(
    psi::SquareIPEPSState,
    ctx::PEPSKitMeasurementContext,
)
    cell = psi.unitcell
    key = (:pxp_energy_density, cell.Lx, cell.Ly)
    return get!(ctx.operator_cache, key) do
        _pepskit_pxp_energy_operator(cell)
    end
end

"""
    local_density_ctm(psi, c, ctx)::Float64

Measure the one-site Rydberg density `<n_c>` at coordinate `c` using PEPSKit
`expectation_value` and a precomputed CTMRG context. Basis index `1` is the
Rydberg/up state.
"""
function local_density_ctm(
    psi::SquareIPEPSState,
    c::SquareCoord,
    ctx::PEPSKitMeasurementContext,
)::Float64
    _assert_fresh_context(psi, ctx)
    return _expectation(ctx, _pepskit_density_operator(psi.unitcell, c))
end

"""
    nearest_neighbor_density_ctm(psi, c, dir, ctx)::Float64

Measure the nearest-neighbor pair density `<n_c n_neighbor>` for `dir` in
`:right`, `:up`, `:left`, or `:down` using PEPSKit CTMRG.
"""
function nearest_neighbor_density_ctm(
    psi::SquareIPEPSState,
    c::SquareCoord,
    dir::Symbol,
    ctx::PEPSKitMeasurementContext,
)::Float64
    _assert_fresh_context(psi, ctx)
    _validate_direction(dir)
    return _expectation(ctx, _pepskit_twosite_nn_operator(psi.unitcell, c, dir))
end

"""
    star_expectation_ctm(psi, center, O, ctx)::ComplexF64

Measure the normalized five-site square-star expectation of dense `32x32`
operator `O` using PEPSKit `expectation_value` and the precomputed CTMRG
context `ctx`. The dense site order is `(center, right, up, left, down)`, with
physical basis `1 = :up`/Rydberg and `2 = :down`/vacancy. PEPSKit site keys are
literal tuples in the adapter coordinate convention
`SquareCoord(x,y) -> CartesianIndex(row = Ly - y + 1, col = x)`.

This function reuses `ctx`; it does not run CTMRG. The corresponding
`LocalOperator` is cached in `ctx.operator_cache`.
"""
function star_expectation_ctm(
    psi::SquareIPEPSState,
    center::SquareCoord,
    O::AbstractMatrix,
    ctx::PEPSKitMeasurementContext,
)::ComplexF64
    _assert_fresh_context(psi, ctx)
    ctx.peps !== nothing || throw(ArgumentError("PEPSKit measurement context has no PEPS"))
    op = _cached_star_localoperator(psi, center, O, ctx)
    value = PEPSKit.expectation_value(ctx.peps, op, ctx.env)
    _real_expectation(value)
    return ComplexF64(value)
end

"""
    blockade_violation_ctm(psi, ctx)::Float64

Return the average nearest-neighbor blockade violation over canonical
periodic `:right` and `:up` bonds using PEPSKit CTMRG.
"""
function blockade_violation_ctm(
    psi::SquareIPEPSState,
    ctx::PEPSKitMeasurementContext,
)::Float64
    _assert_fresh_context(psi, ctx)
    total = 0.0
    count = 0
    for c in psi.unitcell.reps, dir in (:right, :up)
        total += nearest_neighbor_density_ctm(psi, c, dir, ctx)
        count += 1
    end
    return total / count
end

"""
    pxp_energy_density_ctm(psi, ctx)::Float64

Return the unit-cell average five-site square-star PXP energy density using
PEPSKit CTMRG. The dense `square_pxp_star_hamiltonian()` from `SquarePXP.jl` is
the source of truth, in star order `(center, right, up, left, down)` and
physical basis `1 = :up`/Rydberg, `2 = :down`/vacancy. The supplied context is
reused for every center; this function does not rerun CTMRG.
"""
function pxp_energy_density_ctm(
    psi::SquareIPEPSState,
    ctx::PEPSKitMeasurementContext,
)::Float64
    _assert_fresh_context(psi, ctx)
    ctx.peps !== nothing || throw(ArgumentError("PEPSKit measurement context has no PEPS"))
    total = _expectation(ctx, _cached_pxp_energy_operator(psi, ctx))
    return total / length(psi.unitcell.reps)
end

function _density_ctm(
    psi::SquareIPEPSState,
    ctx::PEPSKitMeasurementContext;
    sublattice = nothing,
)
    reps = if sublattice === nothing
        psi.unitcell.reps
    elseif sublattice === :even
        [c for c in psi.unitcell.reps if iseven(c.x + c.y)]
    elseif sublattice === :odd
        [c for c in psi.unitcell.reps if isodd(c.x + c.y)]
    else
        throw(ArgumentError("sublattice must be nothing, :even, or :odd"))
    end
    isempty(reps) && throw(ArgumentError("selected sublattice is empty"))
    return sum(local_density_ctm(psi, c, ctx) for c in reps) / length(reps)
end

"""
    measure_ctm(psi; params = PEPSKitCTMRGParams(8, 1e-8, 100, 0))::CTMObservableSummary

Build one PEPSKit CTMRG measurement context for `psi` and reuse it to compute
density, even/odd densities, nearest-neighbor blockade violation, and
five-site PXP energy density. Density, blockade, and PXP energy are all
CTMRG-backed PEPSKit `expectation_value` measurements. This is a measurement
backend for states produced by the custom ITensors simple-update engine; it
does not update or evolve the state.
"""
function measure_ctm(
    psi::SquareIPEPSState;
    params::PEPSKitCTMRGParams = PEPSKitCTMRGParams(8, 1e-8, 100, 0),
)::CTMObservableSummary
    ctx = pepskit_ctmrg_context(psi; params)
    return CTMObservableSummary(
        _density_ctm(psi, ctx),
        _density_ctm(psi, ctx; sublattice = :even),
        _density_ctm(psi, ctx; sublattice = :odd),
        blockade_violation_ctm(psi, ctx),
        pxp_energy_density_ctm(psi, ctx),
        ctx.diagnostics,
    )
end

end
