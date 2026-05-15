module Observables

using ITensors
using ..SpinOps: projector_up
using ..SquareGeometry
using ..SquarePXP: SQUARE_STAR_SITES, square_pxp_star_hamiltonian
using ..SquareUnitCells
using ..SquareIPEPS

export local_density_simple, density_simple, sublattice_densities
export nearest_neighbor_density_simple, blockade_violation_simple
export star_expectation_simple, pxp_energy_density_simple
export mean_bond_entropy, max_bond_entropy
export SimpleObservableSummary, measure_simple

const _DIRECTIONS = (:right, :up, :left, :down)

function _validate_direction(dir::Symbol)
    dir in _DIRECTIONS ||
        throw(ArgumentError("direction must be :right, :up, :left, or :down"))
    return dir
end

function _opposite_dir(dir::Symbol)
    dir === :right && return :left
    dir === :up && return :down
    dir === :left && return :right
    dir === :down && return :up
    throw(ArgumentError("direction must be :right, :up, :left, or :down"))
end

function _external_dirs_for_leaf(dir::Symbol)
    dir === :right && return (:right, :up, :down)
    dir === :up && return (:left, :right, :up)
    dir === :left && return (:left, :up, :down)
    dir === :down && return (:left, :right, :down)
    throw(ArgumentError("direction must be :right, :up, :left, or :down"))
end

function _real_expectation(value; atol = 1e-10)
    z = ComplexF64(value)
    isfinite(real(z)) && isfinite(imag(z)) ||
        throw(ArgumentError("expectation value must be finite"))
    abs(imag(z)) <= atol || throw(
        ArgumentError(
            "Hermitian observable produced non-negligible imaginary part $(imag(z))",
        ),
    )
    return Float64(real(z))
end

function _positive_norm(value; atol = 1e-12)
    norm_value = _real_expectation(value; atol)
    norm_value > atol || throw(ArgumentError("local simple observable patch has zero norm"))
    return norm_value
end

function _site_tensor_with_weights(
    psi::SquareIPEPSState,
    c::SquareCoord;
    dirs = _DIRECTIONS,
)
    site = wrap(psi.unitcell, c)
    T = copy(psi.tensors[site])
    for dir in dirs
        _validate_direction(dir)
        T = absorb_link_weight(T, psi, site, dir)
    end
    return T
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

function _dense_operator_itensor(O::AbstractMatrix, phys::NTuple{N,Index}) where {N}
    size(O) == (2^N, 2^N) || throw(ArgumentError("dense operator must be $(2^N)x$(2^N)"))
    all(p -> dim(p) == 2, phys) ||
        throw(ArgumentError("physical indices must have dimension 2"))

    out = prime.(phys)
    data = zeros(ComplexF64, ntuple(Returns(2), 2N))
    for out_values in Iterators.product((1:2 for _ = 1:N)...)
        out_idx = _dense_index(out_values)
        for in_values in Iterators.product((1:2 for _ = 1:N)...)
            in_idx = _dense_index(in_values)
            data[out_values..., in_values...] = O[out_idx, in_idx]
        end
    end
    return ITensor(data, out..., phys...)
end

function _prime_physical(T::ITensor, phys)
    out = T
    for p in phys
        out = prime(out, p)
    end
    return out
end

function _expectation_from_patch(theta::ITensor, O::AbstractMatrix, phys)
    phys_tuple = Tuple(phys)
    op = _dense_operator_itensor(O, phys_tuple)
    numerator =
        @disable_warn_order scalar(dag(_prime_physical(theta, phys_tuple)) * (op * theta))
    denominator = @disable_warn_order scalar(dag(theta) * theta)
    return ComplexF64(numerator / _positive_norm(denominator))
end

function _selected_reps(cell::PeriodicSquareUnitCell, sublattice)
    if sublattice === nothing
        return cell.reps
    elseif sublattice === :even
        return [c for c in cell.reps if iseven(c.x + c.y)]
    elseif sublattice === :odd
        return [c for c in cell.reps if isodd(c.x + c.y)]
    else
        throw(ArgumentError("sublattice must be nothing, :even, or :odd"))
    end
end

function _star_coords(psi::SquareIPEPSState, center::SquareCoord)
    cell = psi.unitcell
    c = wrap(cell, center)
    return (
        center = c,
        right = neighbor(cell, c, :right),
        up = neighbor(cell, c, :up),
        left = neighbor(cell, c, :left),
        down = neighbor(cell, c, :down),
    )
end

function _validate_distinct_star(psi::SquareIPEPSState, center::SquareCoord)
    coords = _star_coords(psi, center)
    sites = (coords.center, coords.right, coords.up, coords.left, coords.down)
    length(Set(sites)) == SQUARE_STAR_SITES || throw(
        ArgumentError(
            "wrapped square star must contain five distinct unit-cell representatives",
        ),
    )
    return coords
end

function _star_patch_tensor(psi::SquareIPEPSState, center::SquareCoord)
    coords = _validate_distinct_star(psi, center)
    # Internal center-to-leaf lambdas are absorbed into the center tensor only.
    # Each leaf carries just its three external lambdas, so every star-patch
    # bond weight is counted once in this simple-update local environment.
    theta = _site_tensor_with_weights(psi, coords.center)
    for dir in _DIRECTIONS
        leaf = getproperty(coords, dir)
        theta = @disable_warn_order theta * _site_tensor_with_weights(
            psi,
            leaf;
            dirs = _external_dirs_for_leaf(dir),
        )
    end
    phys = (
        physical_index(psi, coords.center),
        physical_index(psi, coords.right),
        physical_index(psi, coords.up),
        physical_index(psi, coords.left),
        physical_index(psi, coords.down),
    )
    return theta, phys
end

"""
    local_density_simple(psi, c)::Float64

Return the simple-update local-environment Rydberg density `<P_up>` at site
`c`. The local Γ tensor is copied, all four incident λ weights are absorbed
into its virtual legs, and the one-site bra-ket expectation is normalized by
the same local patch norm. Basis index `1` is `:up`/Rydberg.
"""
function local_density_simple(psi::SquareIPEPSState, c::SquareCoord)::Float64
    site = wrap(psi.unitcell, c)
    A = _site_tensor_with_weights(psi, site)
    p = physical_index(psi, site)
    value = _expectation_from_patch(A, projector_up(), (p,))
    return _real_expectation(value)
end

"""
    density_simple(psi; sublattice = nothing)::Float64

Return the average simple-update local Rydberg density. Optional `sublattice`
may be `nothing`, `:even`, or `:odd`. These are cheap local-environment
diagnostics for Γ-λ iPEPS states, not CTMRG-quality measurements.
"""
function density_simple(psi::SquareIPEPSState; sublattice = nothing)::Float64
    reps = _selected_reps(psi.unitcell, sublattice)
    isempty(reps) && throw(ArgumentError("selected sublattice is empty"))
    return sum(local_density_simple(psi, c) for c in reps) / length(reps)
end

"""
    sublattice_densities(psi)

Return `(even = ..., odd = ...)` simple-update Rydberg densities for the even
and odd parity sublattices of a square iPEPS state.
"""
function sublattice_densities(psi::SquareIPEPSState)
    return (
        even = density_simple(psi; sublattice = :even),
        odd = density_simple(psi; sublattice = :odd),
    )
end

"""
    nearest_neighbor_density_simple(psi, c, dir)::Float64

Return the simple-update two-site expectation `<n_c n_neighbor>` for the
nearest-neighbor bond from `c` in `dir`. External λ weights are absorbed into
their carrying tensors, while the internal bond λ is absorbed into the first
site only before the shared bond is contracted, so the internal λ is counted
exactly once.
"""
function nearest_neighbor_density_simple(
    psi::SquareIPEPSState,
    c::SquareCoord,
    dir::Symbol,
)::Float64
    _validate_direction(dir)
    site = wrap(psi.unitcell, c)
    other = neighbor(psi.unitcell, site, dir)
    opposite = _opposite_dir(dir)

    left_dirs = Tuple(d for d in _DIRECTIONS if d !== dir)
    right_dirs = Tuple(d for d in _DIRECTIONS if d !== opposite)
    A = _site_tensor_with_weights(psi, site; dirs = left_dirs)
    A = absorb_link_weight(A, psi, site, dir)
    B = _site_tensor_with_weights(psi, other; dirs = right_dirs)
    theta = A * B

    p1 = physical_index(psi, site)
    p2 = physical_index(psi, other)
    value = _expectation_from_patch(theta, kron(projector_up(), projector_up()), (p1, p2))
    return _real_expectation(value)
end

"""
    blockade_violation_simple(psi)::Float64

Return the average nearest-neighbor excited-pair density over canonical
periodic `:right` and `:up` bonds using the simple-update two-site local
environment.
"""
function blockade_violation_simple(psi::SquareIPEPSState)::Float64
    total = 0.0
    count = 0
    for c in psi.unitcell.reps
        for dir in (:right, :up)
            total += nearest_neighbor_density_simple(psi, c, dir)
            count += 1
        end
    end
    return total / count
end

"""
    star_expectation_simple(psi, center, O)::ComplexF64

Return the normalized five-site simple-update star expectation for dense
`32x32` operator `O` in square-star order `(center, right, up, left, down)`.
The four internal center-leaf λ weights and all external star-patch λ weights
are each included exactly once. This is a local diagnostic, not a CTMRG
environment measurement.
"""
function star_expectation_simple(
    psi::SquareIPEPSState,
    center::SquareCoord,
    O::AbstractMatrix,
)::ComplexF64
    size(O) == (2^SQUARE_STAR_SITES, 2^SQUARE_STAR_SITES) ||
        throw(ArgumentError("dense square-star operator must be 32x32"))
    theta, phys = _star_patch_tensor(psi, center)
    return _expectation_from_patch(theta, O, phys)
end

"""
    pxp_energy_density_simple(psi)::Float64

Return the unit-cell average of the local square-star PXP Hamiltonian
expectation using [`star_expectation_simple`](@ref). The dense
`square_pxp_star_hamiltonian()` convention is the source of truth.
"""
function pxp_energy_density_simple(psi::SquareIPEPSState)::Float64
    Hstar = square_pxp_star_hamiltonian()
    value =
        sum(star_expectation_simple(psi, c, Hstar) for c in psi.unitcell.reps) /
        length(psi.unitcell.reps)
    return _real_expectation(value)
end

"""
    mean_bond_entropy(psi)::Float64

Return the mean entropy of all canonical simple-update link-weight spectra.
Product states have zero mean bond entropy.
"""
function mean_bond_entropy(psi::SquareIPEPSState)::Float64
    entropies = collect(values(all_bond_entropies(psi)))
    isempty(entropies) && throw(ArgumentError("state has no bond entropies"))
    result = sum(entropies) / length(entropies)
    isfinite(result) || throw(ArgumentError("mean bond entropy must be finite"))
    return Float64(result)
end

"""
    max_bond_entropy(psi)::Float64

Return the maximum entropy of all canonical simple-update link-weight spectra.
Product states have zero maximum bond entropy.
"""
function max_bond_entropy(psi::SquareIPEPSState)::Float64
    entropies = collect(values(all_bond_entropies(psi)))
    isempty(entropies) && throw(ArgumentError("state has no bond entropies"))
    result = maximum(entropies)
    isfinite(result) || throw(ArgumentError("maximum bond entropy must be finite"))
    return Float64(result)
end

"""
    SimpleObservableSummary

Compact deterministic diagnostics from simple-update/local-environment
observables: total density, even/odd sublattice densities, nearest-neighbor
blockade violation, five-site PXP energy density, and link-entropy summaries.
"""
struct SimpleObservableSummary
    density::Float64
    density_even::Float64
    density_odd::Float64
    blockade_violation::Float64
    pxp_energy_density::Float64
    mean_bond_entropy::Float64
    max_bond_entropy::Float64
end

"""
    measure_simple(psi)::SimpleObservableSummary

Compute cheap deterministic simple-update diagnostics for a custom ITensors
square iPEPS state. This does not run CTMRG or refresh any environment.
"""
function measure_simple(psi::SquareIPEPSState)::SimpleObservableSummary
    densities = sublattice_densities(psi)
    return SimpleObservableSummary(
        density_simple(psi),
        densities.even,
        densities.odd,
        blockade_violation_simple(psi),
        pxp_energy_density_simple(psi),
        mean_bond_entropy(psi),
        max_bond_entropy(psi),
    )
end

end
