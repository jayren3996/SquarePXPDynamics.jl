module PEPSKitObservables

# PEPSKit-native PXP observables on the native `PXPIPEPSState`. Two layers live
# here:
#   - CTMRG observables: PEPSKit `LocalOperator`s (density, blockade, PXP star,
#     energy) measured via `CTMRGEnv` / `leading_boundary` / `expectation_value`
#     (no ITensors conversion).
#   - Simple/local observables: computed directly on the convention-"B"
#     `InfinitePEPS` (each stored Γ carries `sqrt(λ)`) by contracting small open
#     patches with their conjugates, no CTMRG. For a product (D = 1) state these
#     are exact; for D > 1 they agree with the legacy simple observables in the
#     same Vidal gauge. This layer also provides the `measure_simple` /
#     `evolve!` / `reverse_evolve!` adapters that let the ScarFinder
#     orchestration (M6) run on either tensor engine through ordinary multiple
#     dispatch, plus the `pxpipeps_from_square_ipeps` converter.

using TensorKit
using PEPSKit
using PEPSKit: SUWeight
using Random: MersenneTwister

using ..Internals: _validate_direction
using ..SquareGeometry: SquareCoord
using ..SquareUnitCells: PeriodicSquareUnitCell, wrap, neighbor
using ..SquarePXP: SQUARE_STAR_SITES, square_pxp_star_hamiltonian
using ..SquareIPEPS: SquareIPEPSState, weight_entropy, link_weight
using ..PEPSKitMeasurements:
    PEPSKitCTMRGParams, default_ctmrg_params, to_pepskit_infinitepeps
using ..PEPSKitBackend:
    PXPIPEPSState,
    squarecoord_to_cartesianindex,
    cartesianindex_to_squarecoord,
    measurement_peps,
    state_version,
    log_norm,
    _dense_operator_tensor_map
using ..PEPSKitEvolution: evolve_pepskit!
using ..StarSimpleUpdate: StarUpdateInfo

import ..Observables: measure_simple
import ..IPEPSEvolution: evolve!, reverse_evolve!
using ..Observables: SimpleObservableSummary
using ..IPEPSEvolution: TrotterParams, EvolutionLog

export density_operator, blockade_operator, pxp_star_operator, pxp_energy_operator
export PEPSKitNativeContext, pepskit_native_context
export measure_ctm_pepskit, density_ctm_pepskit
export blockade_violation_ctm_pepskit, pxp_energy_density_ctm_pepskit
export PEPSKitNativeSummary
export pxpipeps_from_square_ipeps

# === CTMRG observable layer (was PEPSKitPXPObservables) ===

# dense -> TensorMap (`_dense_operator_tensor_map`) is shared from PEPSKitBackend.

_physical_lattice(cell::PeriodicSquareUnitCell) = fill(ComplexSpace(2), cell.Ly, cell.Lx)

# Rydberg/excited number operator: |excited><excited| = |1><1|.
_density_matrix() = ComplexF64[1 0; 0 0]

# `_validate_direction` is shared from Internals.

function _neighbor_cartesianindex(cell::PeriodicSquareUnitCell, c::SquareCoord, dir::Symbol)
    _validate_direction(dir)
    center = squarecoord_to_cartesianindex(cell, c)
    dir === :right && return center + CartesianIndex(0, 1)
    dir === :up && return center + CartesianIndex(-1, 0)
    dir === :left && return center + CartesianIndex(0, -1)
    return center + CartesianIndex(1, 0)
end

# Star sites in (center, right, up, left, down) PEPSKit CartesianIndex order.
function _star_sites_cartesian(cell::PeriodicSquareUnitCell, center::SquareCoord)
    c = wrap(cell, center)
    wrapped = (
        c,
        neighbor(cell, c, :right),
        neighbor(cell, c, :up),
        neighbor(cell, c, :left),
        neighbor(cell, c, :down),
    )
    length(Set(wrapped)) == SQUARE_STAR_SITES || throw(
        ArgumentError(
            "wrapped square star must contain five distinct unit-cell representatives",
        ),
    )
    pc = squarecoord_to_cartesianindex(cell, c)
    return (
        pc,
        pc + CartesianIndex(0, 1),
        pc + CartesianIndex(-1, 0),
        pc + CartesianIndex(0, -1),
        pc + CartesianIndex(1, 0),
    )
end

# --- LocalOperators ---

"""
    density_operator(cell, site)

PEPSKit `LocalOperator` for the Rydberg density `n = |excited><excited|` at
`site`.
"""
function density_operator(cell::PeriodicSquareUnitCell, site::SquareCoord)
    op = _dense_operator_tensor_map(_density_matrix(), 1)
    return PEPSKit.LocalOperator(
        _physical_lattice(cell),
        (squarecoord_to_cartesianindex(cell, site),) => op,
    )
end

"""
    blockade_operator(cell)

PEPSKit `LocalOperator` summing `n_c n_neighbor` over every `(rep, :right)` and
`(rep, :up)` nearest-neighbour bond: the total blockade-violation operator.
"""
function blockade_operator(cell::PeriodicSquareUnitCell)
    n = _density_matrix()
    op = _dense_operator_tensor_map(kron(n, n), 2)
    terms = Pair{NTuple{2, CartesianIndex{2}}, typeof(op)}[]
    for c in cell.reps, dir in (:right, :up)
        sites = (
            squarecoord_to_cartesianindex(cell, c),
            _neighbor_cartesianindex(cell, c, dir),
        )
        push!(terms, sites => op)
    end
    return PEPSKit.LocalOperator(_physical_lattice(cell), terms...)
end

"""
    pxp_star_operator(cell, center)

PEPSKit `LocalOperator` for the single five-site PXP star term at `center`.
"""
function pxp_star_operator(cell::PeriodicSquareUnitCell, center::SquareCoord)
    op = _dense_operator_tensor_map(square_pxp_star_hamiltonian(), SQUARE_STAR_SITES)
    return PEPSKit.LocalOperator(
        _physical_lattice(cell),
        _star_sites_cartesian(cell, center) => op,
    )
end

"""
    pxp_energy_operator(cell)

PEPSKit `LocalOperator` summing the PXP star term over every unit-cell center.
"""
function pxp_energy_operator(cell::PeriodicSquareUnitCell)
    op = _dense_operator_tensor_map(square_pxp_star_hamiltonian(), SQUARE_STAR_SITES)
    terms = [_star_sites_cartesian(cell, c) => op for c in cell.reps]
    return PEPSKit.LocalOperator(_physical_lattice(cell), terms...)
end

# --- CTMRG context + measurements ---

_ctm_initializer(::Nothing) = randn
function _ctm_initializer(seed::Int)
    rng = MersenneTwister(seed)
    return (T, dims...) -> randn(rng, T, dims...)
end

"""
    PEPSKitNativeContext

Cached CTMRG context for a PEPSKit-native state: the absorbed measurement PEPS,
the converged `CTMRGEnv`, the raw convergence `info`, the params, and the source
state version for staleness checks.
"""
struct PEPSKitNativeContext
    cell::PeriodicSquareUnitCell
    peps::PEPSKit.InfinitePEPS
    env::Any
    info::Any
    params::PEPSKitCTMRGParams
    source_version::Int
end

"""
    pepskit_native_context(state; params = default_ctmrg_params())

Absorb the SUWeight into a measurement PEPS, build a `CTMRGEnv`, run
`leading_boundary`, and return a reusable `PEPSKitNativeContext`.
"""
function pepskit_native_context(
    state::PXPIPEPSState;
    params::PEPSKitCTMRGParams = default_ctmrg_params(),
)::PEPSKitNativeContext
    peps = measurement_peps(state)
    chi_space = ComplexSpace(params.chi)
    env0 = PEPSKit.CTMRGEnv(_ctm_initializer(params.seed), ComplexF64, peps, chi_space)
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
    return PEPSKitNativeContext(
        state.unitcell, peps, env, info, params, state_version(state),
    )
end

function _real_expectation(value; atol = 1e-8)
    z = ComplexF64(value)
    abs(imag(z)) <= atol * max(1, abs(real(z))) || throw(
        ArgumentError(
            "Hermitian CTMRG observable produced non-negligible imaginary part $(imag(z))",
        ),
    )
    return Float64(real(z))
end

_expectation(ctx::PEPSKitNativeContext, op) =
    _real_expectation(PEPSKit.expectation_value(ctx.peps, op, ctx.env))

"""
    density_ctm_pepskit(ctx, site)

Rydberg density at `site` from the CTMRG context `ctx`.
"""
function density_ctm_pepskit(ctx::PEPSKitNativeContext, site::SquareCoord)
    return _expectation(ctx, density_operator(ctx.cell, site))
end

"""
    blockade_violation_ctm_pepskit(ctx)

Total nearest-neighbour blockade violation `sum <n_c n_neighbor>` from `ctx`.
"""
function blockade_violation_ctm_pepskit(ctx::PEPSKitNativeContext)
    return _expectation(ctx, blockade_operator(ctx.cell))
end

"""
    pxp_energy_density_ctm_pepskit(ctx)

PXP energy per unit-cell site: `<H_star summed over centers> / nsites`.
"""
function pxp_energy_density_ctm_pepskit(ctx::PEPSKitNativeContext)
    nsites = length(ctx.cell.reps)
    return _expectation(ctx, pxp_energy_operator(ctx.cell)) / nsites
end

"""
    PEPSKitNativeSummary

Bundle of CTMRG-measured PXP observables on a PEPSKit-native state.
"""
struct PEPSKitNativeSummary
    chi::Int
    mean_density::Float64
    blockade_violation::Float64
    pxp_energy_density::Float64
end

"""
    measure_ctm_pepskit(state; params = default_ctmrg_params())

Run CTMRG once on the PEPSKit-native `state` and return a
`PEPSKitNativeSummary` of the mean Rydberg density, total blockade violation,
and PXP energy density.
"""
function measure_ctm_pepskit(
    state::PXPIPEPSState;
    params::PEPSKitCTMRGParams = default_ctmrg_params(),
)::PEPSKitNativeSummary
    ctx = pepskit_native_context(state; params)
    reps = ctx.cell.reps
    mean_density = sum(c -> density_ctm_pepskit(ctx, c), reps) / length(reps)
    return PEPSKitNativeSummary(
        params.chi,
        mean_density,
        blockade_violation_ctm_pepskit(ctx),
        pxp_energy_density_ctm_pepskit(ctx),
    )
end

# === Simple/local observable layer + adapters (was PEPSKitSimpleObservables) ===

# --- dense convention-B patches -------------------------------------------
#
# Local "simple" expectations on the Vidal/convention-"B" state. The canonical
# simple environment places the squared Schmidt weight λ² on every *cut*
# (external) bond of the patch and the bare bond weight λ once on every *internal*
# (shared) bond. Convention-B tensors already carry sqrt(λ) on each leg, so:
#   - an internal bond contracts sqrt(λ)·sqrt(λ) = λ once — left as-is;
#   - an external bond would give only sqrt(λ)² = λ in the bra–ket, so each
#     external leg is multiplied by an extra sqrt(λ) (from the SUWeight ledger)
#     to reach the full λ², matching the legacy simple-update convention.

# Array leg index (in [physical, N, E, S, W]) carrying the bond in `dir`.
_leg_dim(dir::Symbol) = dir === :up ? 2 : dir === :right ? 3 : dir === :down ? 4 :
                        dir === :left ? 5 : throw(ArgumentError("bad direction $dir"))

# Schmidt weights λ (Vector) on the bond in `dir` from site `c`. East/north
# bonds are stored at the site; west/south are the east/north bonds of the
# left/down neighbour.
function _leg_weight(state::PXPIPEPSState, c::SquareCoord, dir::Symbol)
    cell = state.unitcell
    if dir === :right
        idx = squarecoord_to_cartesianindex(cell, c)
        return state.weights.data[1, idx[1], idx[2]].data
    elseif dir === :up
        idx = squarecoord_to_cartesianindex(cell, c)
        return state.weights.data[2, idx[1], idx[2]].data
    elseif dir === :left
        return _leg_weight(state, neighbor(cell, c, :left), :right)
    elseif dir === :down
        return _leg_weight(state, neighbor(cell, c, :down), :up)
    else
        throw(ArgumentError("bad direction $dir"))
    end
end

# Dense site array [physical, N, E, S, W] with each external direction's leg
# multiplied by sqrt(λ) so cut bonds carry the full λ² in the bra–ket.
function _site_array(state::PXPIPEPSState, c::SquareCoord, externals)
    idx = squarecoord_to_cartesianindex(state.unitcell, c)
    A = convert(Array, state.peps[idx[1], idx[2]])
    for dir in externals
        d = _leg_dim(dir)
        v = sqrt.(_leg_weight(state, c, dir))
        shape = ntuple(i -> i == d ? length(v) : 1, ndims(A))
        A = A .* reshape(v, shape)
    end
    return A
end

const _ALL_DIRS = (:up, :right, :down, :left)

# Single-site reduced density matrix ρ[p, p'] = Σ_virtual A[p, …] conj(A[p', …]).
function _site_rdm(A::AbstractArray)
    M = reshape(A, size(A, 1), :)
    return M * M'
end

function _local_density(state::PXPIPEPSState, c::SquareCoord)
    ρ = _site_rdm(_site_array(state, c, _ALL_DIRS))
    num = real(ρ[1, 1])
    den = real(ρ[1, 1] + ρ[2, 2])
    den > 0 || throw(ArgumentError("non-positive single-site norm at $c"))
    return num / den
end

# Two-site reduced density matrix on the bond from `c` in `dir` (`:right`/`:up`),
# legs (p1, p2, p1', p2'). The shared bond stays sqrt(λ) on each tensor (→ λ once);
# the six external legs are weighted to the full λ².
function _bond_rdm(state::PXPIPEPSState, c::SquareCoord, dir::Symbol)
    other = neighbor(state.unitcell, c, dir)
    if dir === :right            # center E (leg 3) ↔ neighbour W (leg 5)
        A0 = _site_array(state, c, (:up, :down, :left))
        A1 = _site_array(state, other, (:up, :right, :down))
        θ = ncon((A0, A1), ([-1, -3, 1, -4, -5], [-2, -6, -7, -8, 1]))
    elseif dir === :up           # center N (leg 2) ↔ neighbour S (leg 4)
        A0 = _site_array(state, c, (:right, :down, :left))
        A1 = _site_array(state, other, (:up, :right, :left))
        θ = ncon((A0, A1), ([-1, 1, -3, -4, -5], [-2, -6, -7, 1, -8]))
    else
        throw(ArgumentError("bond direction must be :right or :up"))
    end
    # trace the six external legs against the conjugate patch
    return ncon((θ, conj(θ)), ([-1, -2, 1, 2, 3, 4, 5, 6], [-3, -4, 1, 2, 3, 4, 5, 6]))
end

function _nn_excited(state::PXPIPEPSState, c::SquareCoord, dir::Symbol)
    ρ = _bond_rdm(state, c, dir)
    num = real(ρ[1, 1, 1, 1])
    den = 0.0
    for a in 1:size(ρ, 1), b in 1:size(ρ, 2)
        den += real(ρ[a, b, a, b])
    end
    den > 0 || throw(ArgumentError("non-positive two-site norm at $c $dir"))
    return num / den
end

# Big-endian 5-site linear index in (center, right, up, left, down) order,
# matching `square_pxp_star_hamiltonian()`.
_didx5(a, b, c, d, e) = 1 + (a - 1) * 16 + (b - 1) * 8 + (c - 1) * 4 + (d - 1) * 2 + (e - 1)

function _star_hamiltonian_tensor()
    H = square_pxp_star_hamiltonian()
    T = zeros(ComplexF64, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2)
    for cp in 1:2, rp in 1:2, up in 1:2, lp in 1:2, dp in 1:2
        R = _didx5(cp, rp, up, lp, dp)
        for c in 1:2, r in 1:2, u in 1:2, l in 1:2, d in 1:2
            T[cp, rp, up, lp, dp, c, r, u, l, d] = H[R, _didx5(c, r, u, l, d)]
        end
    end
    return T
end

# Five-site star patch θ5 with legs (center,right,up,left,down) + 12 externals.
# The center's four legs are all internal (shared with the leaves); each leaf
# shares one leg with the center and exposes three external legs.
function _star_patch(state::PXPIPEPSState, center::SquareCoord)
    cell = state.unitcell
    A0 = _site_array(state, center, ())                                  # all internal
    Ar = _site_array(state, neighbor(cell, center, :right), (:up, :right, :down))
    Au = _site_array(state, neighbor(cell, center, :up), (:up, :right, :left))
    Al = _site_array(state, neighbor(cell, center, :left), (:up, :down, :left))
    Ad = _site_array(state, neighbor(cell, center, :down), (:right, :down, :left))
    # bonds: 101 = N0↔up.S, 102 = E0↔right.W, 103 = S0↔down.N, 104 = W0↔left.E
    return ncon(
        (A0, Ar, Au, Al, Ad),
        (
            [-1, 101, 102, 103, 104],   # center [p, N, E, S, W]
            [-2, -6, -7, -8, 102],      # right  [p, N, E, S, W=bond]
            [-3, -9, -10, 101, -11],    # up     [p, N, E, S=bond, W]
            [-4, -12, 104, -13, -14],   # left   [p, N, E=bond, S, W]
            [-5, 103, -15, -16, -17],   # down   [p, N=bond, E, S, W]
        ),
    )
end

function _star_energy(state::PXPIPEPSState, center::SquareCoord, Hten)
    θ = _star_patch(state, center)
    ext = collect(6:17)
    num = ncon(
        (conj(θ), Hten, θ),
        (
            vcat([1, 2, 3, 4, 5], ext),
            [1, 2, 3, 4, 5, 18, 19, 20, 21, 22],
            vcat([18, 19, 20, 21, 22], ext),
        ),
    )
    den = ncon((conj(θ), θ), (vcat([1, 2, 3, 4, 5], ext), vcat([1, 2, 3, 4, 5], ext)))
    return real(num) / real(den)
end

# --- bond entropies from the SUWeight ledger ------------------------------

function _bond_entropies(state::PXPIPEPSState)
    entropies = Float64[]
    for w in state.weights.data
        push!(entropies, weight_entropy(w.data))
    end
    isempty(entropies) && throw(ArgumentError("state has no bond entropies"))
    return entropies
end

# --- public adapters -------------------------------------------------------

"""
    measure_simple(state::PXPIPEPSState)::SimpleObservableSummary

Compute cheap deterministic simple-update diagnostics for a PEPSKit-native
square iPEPS state directly from its convention-"B" `InfinitePEPS` and `SUWeight`
ledger. This does not run CTMRG. The conventions match the legacy
[`measure_simple`](@ref) so the two tensor engines are comparable: `density` is
the unit-cell mean Rydberg density, `blockade_violation` is the mean
nearest-neighbour excited-pair density over canonical `:right`/`:up` bonds, and
`pxp_energy_density` is the mean five-site star Hamiltonian expectation.
"""
function measure_simple(state::PXPIPEPSState)::SimpleObservableSummary
    cell = state.unitcell
    reps = cell.reps
    isempty(reps) && throw(ArgumentError("unit cell is empty"))
    even = [c for c in reps if iseven(c.x + c.y)]
    odd = [c for c in reps if isodd(c.x + c.y)]

    dens = Dict(c => _local_density(state, c) for c in reps)
    density = sum(values(dens)) / length(reps)
    density_even = isempty(even) ? density : sum(dens[c] for c in even) / length(even)
    density_odd = isempty(odd) ? density : sum(dens[c] for c in odd) / length(odd)

    blockade = 0.0
    for c in reps, dir in (:right, :up)
        blockade += _nn_excited(state, c, dir)
    end
    blockade /= (2 * length(reps))

    Hten = _star_hamiltonian_tensor()
    energy = sum(_star_energy(state, c, Hten) for c in reps) / length(reps)

    ent = _bond_entropies(state)
    return SimpleObservableSummary(
        density,
        density_even,
        density_odd,
        blockade,
        energy,
        sum(ent) / length(ent),
        maximum(ent),
    )
end

"""
    evolve!(state::PXPIPEPSState, total_time; params::TrotterParams,
            protocol = nothing)::EvolutionLog

PEPSKit-native evolution adapter exposing the legacy [`evolve!`](@ref)
interface. Maps the [`TrotterParams`](@ref) onto `evolve_pepskit!` and returns a
legacy [`EvolutionLog`](@ref) so the ScarFinder orchestration can consume either
tensor engine. Bond-entropy diagnostics come from the `SUWeight` ledger;
`layer_infos` is empty (the per-star records live in the native
`StarUpdateInfoPK`). Custom evolution protocols are not yet supported on this
engine.
"""
function evolve!(
    state::PXPIPEPSState,
    total_time::Real;
    params::TrotterParams,
    protocol = nothing,
)::EvolutionLog
    return _evolve_native!(state, total_time, params, protocol; reverse = false)
end

"""
    reverse_evolve!(state::PXPIPEPSState, total_time; params, protocol = nothing)

Evolve the PEPSKit-native `state` backward in real time, applying the inverse of
each forward Trotter step (substeps in reverse order with negated layer times).
Mirrors the legacy [`reverse_evolve!`](@ref); `params.evolution` must be `:real`
and custom protocols are not supported on `:pepskit_simple` yet. Used by the
native multi-step reversibility audit.
"""
function reverse_evolve!(
    state::PXPIPEPSState,
    total_time::Real;
    params::TrotterParams,
    protocol = nothing,
)::EvolutionLog
    params.evolution === :real ||
        throw(ArgumentError("reverse_evolve! is defined only for real-time evolution"))
    return _evolve_native!(state, total_time, params, protocol; reverse = true)
end

function _evolve_native!(
    state::PXPIPEPSState,
    total_time::Real,
    params::TrotterParams,
    protocol;
    reverse::Bool,
)::EvolutionLog
    protocol === nothing || throw(ArgumentError(
        "custom evolution protocols are not supported on the PEPSKit-native " *
        "backend (:pepskit_simple) yet",
    ))
    before = log_norm(state)
    logpk = evolve_pepskit!(
        state,
        total_time;
        dt = params.dt,
        order = params.order,
        schedule = params.schedule,
        maxdim = params.maxdim,
        cutoff = params.cutoff,
        rel_floor = params.rel_floor,
        evolution = params.evolution,
        projected = true,
        reverse = reverse,
    )
    ent = _bond_entropies(state)
    after = logpk.log_norm_after
    return EvolutionLog(
        Float64(total_time),
        params,
        (; tensor_engine = :pepskit_simple),
        logpk.nsteps,
        Vector{Vector{StarUpdateInfo}}(),
        logpk.max_truncerr,
        maximum(ent),
        sum(ent) / length(ent),
        before,
        after,
        after - before,
    )
end

"""
    pxpipeps_from_square_ipeps(psi::SquareIPEPSState)::PXPIPEPSState

Convert a legacy ITensors Γ-λ `SquareIPEPSState` into a PEPSKit-native
`PXPIPEPSState` in the shared convention-"B" gauge: the `InfinitePEPS` comes
from [`to_pepskit_infinitepeps`](@ref) (each λ absorbed symmetrically as
`sqrt(λ)`) and the `SUWeight` ledger carries the canonical `:right` (east) and
`:up` (north) link spectra. Used to validate that the native and legacy
simple-observable conventions agree in the same gauge for `D > 1`.
"""
function pxpipeps_from_square_ipeps(psi::SquareIPEPSState)::PXPIPEPSState
    cell = psi.unitcell
    peps = to_pepskit_infinitepeps(psi)
    Nr, Nc = cell.Ly, cell.Lx
    W = DiagonalTensorMap{Float64,ComplexSpace,Vector{Float64}}
    data = Array{W,3}(undef, 2, Nr, Nc)
    for row in 1:Nr, col in 1:Nc
        c = cartesianindex_to_squarecoord(cell, CartesianIndex(row, col))
        λe = Float64.(link_weight(psi, c, :right))
        λn = Float64.(link_weight(psi, c, :up))
        data[1, row, col] = DiagonalTensorMap(λe, ComplexSpace(length(λe)))
        data[2, row, col] = DiagonalTensorMap(λn, ComplexSpace(length(λn)))
    end
    return PXPIPEPSState(peps, SUWeight(data), cell, Ref(state_version(psi)), Ref(log_norm(psi)))
end

end # module
