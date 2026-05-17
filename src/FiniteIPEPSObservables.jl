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
