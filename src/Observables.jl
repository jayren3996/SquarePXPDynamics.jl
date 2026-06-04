module Observables

using ITensors
using ..SpinOps: pauli_x, pauli_y, pauli_z, projector_up
using ..SquareGeometry
using ..SquarePXP: SQUARE_STAR_SITES, square_pxp_star_hamiltonian
using ..SquareUnitCells
using ..SquareIPEPS
using ..Internals: _DIRECTIONS, _validate_direction, _opposite_direction

export local_density_simple, density_simple, sublattice_densities
export sublattice_imbalance_simple, checkerboard_structure_factor_simple
export nearest_neighbor_density_simple, blockade_violation_simple
export star_expectation_simple, pxp_energy_density_simple
export mean_bond_entropy, max_bond_entropy
export SimpleObservableSummary, measure_simple
export local_x_simple, local_y_simple, local_z_simple, nearest_neighbor_zz_simple

const _opposite_dir = _opposite_direction

function _external_dirs_for_leaf(dir::Symbol)
    dir === :right && return (:right, :up, :down)
    dir === :up && return (:left, :right, :up)
    dir === :left && return (:left, :up, :down)
    dir === :down && return (:left, :right, :down)
    throw(ArgumentError("direction must be :right, :up, :left, or :down"))
end

function _validate_canonical_bond_dir(dir::Symbol)
    (dir === :right || dir === :up) ||
        throw(ArgumentError("TFIM canonical bond direction must be :right or :up"))
    return dir
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

function _one_site_expectation_simple(psi::SquareIPEPSState, c::SquareCoord, O::AbstractMatrix)
    size(O) == (2, 2) || throw(ArgumentError("one-site operator must be 2x2"))
    site = wrap(psi.unitcell, c)
    A = _site_tensor_with_weights(psi, site)
    p = physical_index(psi, site)
    return _expectation_from_patch(A, O, (p,))
end

function _nearest_neighbor_expectation_simple(
    psi::SquareIPEPSState,
    c::SquareCoord,
    dir::Symbol,
    O::AbstractMatrix,
)
    _validate_direction(dir)
    size(O) == (4, 4) || throw(ArgumentError("two-site operator must be 4x4"))
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
    return _expectation_from_patch(theta, O, (p1, p2))
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
    return _real_expectation(_one_site_expectation_simple(psi, c, projector_up()))
end

"""
    local_x_simple(psi, c)::Float64

Return the simple-update local-environment Pauli-X expectation at site `c`.
All four incident link weights are absorbed into the one-site local patch.
"""
function local_x_simple(psi::SquareIPEPSState, c::SquareCoord)::Float64
    return _real_expectation(_one_site_expectation_simple(psi, c, pauli_x()))
end

"""
    local_y_simple(psi, c)::Float64

Return the simple-update local-environment Pauli-Y expectation at site `c`.
All four incident link weights are absorbed into the one-site local patch.
"""
function local_y_simple(psi::SquareIPEPSState, c::SquareCoord)::Float64
    return _real_expectation(_one_site_expectation_simple(psi, c, pauli_y()))
end

"""
    local_z_simple(psi, c)::Float64

Return the simple-update local-environment Pauli-Z expectation at site `c`.
The TFIM convention is `Z |up> = +|up>` and `Z |down> = -|down>`.
"""
function local_z_simple(psi::SquareIPEPSState, c::SquareCoord)::Float64
    return _real_expectation(_one_site_expectation_simple(psi, c, pauli_z()))
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
    sublattice_imbalance_simple(psi)::Float64

Return the even-minus-odd Rydberg density contrast from simple/local
measurements.
"""
function sublattice_imbalance_simple(psi::SquareIPEPSState)::Float64
    even, odd = sublattice_densities(psi)
    return even - odd
end

"""
    checkerboard_structure_factor_simple(psi)::Float64

Return the squared checkerboard density contrast from simple/local
measurements.
"""
function checkerboard_structure_factor_simple(psi::SquareIPEPSState)::Float64
    return sublattice_imbalance_simple(psi)^2
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
    value =
        _nearest_neighbor_expectation_simple(psi, c, dir, kron(projector_up(), projector_up()))
    return _real_expectation(value)
end

"""
    nearest_neighbor_zz_simple(psi, c, dir)::Float64

Return the simple-update two-site Pauli-ZZ expectation for the canonical TFIM
bond from `c` in `dir`. Only `:right` and `:up` are accepted, so each periodic
nearest-neighbor bond appears once in unit-cell averages.
"""
function nearest_neighbor_zz_simple(
    psi::SquareIPEPSState,
    c::SquareCoord,
    dir::Symbol,
)::Float64
    _validate_canonical_bond_dir(dir)
    value = _nearest_neighbor_expectation_simple(psi, c, dir, kron(pauli_z(), pauli_z()))
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

# Example

```julia
using SquarePXPDynamics
cell = PeriodicSquareUnitCell(10, 10)
psi = product_square_ipeps(cell; state = :down, maxdim = 1)
summary = measure_simple(psi)
# summary.density == 0.0 for the all-down product state
```
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



module FiniteIPEPSObservables

using ITensors
using LinearAlgebra

using ..SpinOps: projector_up
using ..SquareGeometry
using ..SquarePXP: SQUARE_STAR_SITES, square_pxp_star_hamiltonian
using ..SquareUnitCells
using ..SquareIPEPS
using ..Internals: _DIRECTIONS

export dense_state_finite
export exact_one_site_expectation_finite, exact_nearest_neighbor_expectation_finite
export exact_star_expectation_finite, exact_density_finite
export exact_all_down_return_probability_finite
export exact_blockade_violation_finite, exact_pxp_energy_density_finite

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

# ---------------------------------------------------------------------------
# Memory-efficient EXACT oracle: double-layer (boundary) contraction
# ---------------------------------------------------------------------------
# `dense_state_finite` builds the single-layer 2^N state vector; its intermediate
# carries the accumulating physical legs (2^k) AND the full torus vertical
# boundary (~2*Lx dangling D-bonds) at once, peaking at ~2^(N/2)*D^(2Lx) -> OOM
# for D>=5 on 4x4. Here we contract the DOUBLE layer (ket * conj(bra) with the
# physical legs closed locally) so there is NO 2^k factor; densities are scalars
# <n_c> = Z_c / Z. The construction below (E-tensors + per-bond combiners) is the
# correctness foundation, validated by a naive full contraction against the dense
# path; the memory-efficient boundary sweep builds on these same E-tensors.
# See notes/stage2-truncation/improvement-roadmap.md (priority #1).

_identity_2x2() = ComplexF64[1.0 0.0; 0.0 1.0]

# One combiner per undirected bond, shared by BOTH incident sites so the combined
# D^2 index is the same Index object across neighbors (and across the periodic
# wrap), which lets neighbor E-tensors auto-contract on a single doubled leg.
function _bond_combiners(psi::SquareIPEPSState)
    combs = Dict{Index,ITensor}()
    for c in psi.unitcell.reps
        for dir in (:right, :up)            # canonical owners cover every bond once
            li = link_index(psi, c, dir)
            haskey(combs, li) && continue
            combs[li] = combiner(li, prime(li))
        end
    end
    return combs
end

# Double-layer site tensor for rep `c`: ket * op * conj(bra), physical legs closed,
# four doubled virtual legs combined to dim-D^2 indices shared with neighbors.
# `op` is a 2x2 single-site operator inserted between ket and bra (identity -> norm).
function _double_layer_site(
    psi::SquareIPEPSState,
    folded::Dict{SquareCoord,ITensor},
    c::SquareCoord,
    combs::Dict{Index,ITensor};
    op::AbstractMatrix = _identity_2x2(),
)
    K = folded[c]
    p = physical_index(psi, c)
    Kp = prime(dag(K))                       # conj(bra), all indices at plev 1
    opT = ITensor(ComplexF64, prime(p), p)   # op[out,in]: p(plev0,ket) -> p'(plev1,bra)
    for i = 1:2, j = 1:2
        opT[prime(p)=>i, p=>j] = op[i, j]
    end
    raw = @disable_warn_order (K * opT) * Kp
    E = raw
    for dir in (:left, :right, :up, :down)
        E = @disable_warn_order E * combs[link_index(psi, c, dir)]
    end
    return E
end

# Naive full contraction of the double-layer network to a scalar (Z or Z_c).
# Correctness reference only: it forms dense intermediates whose boundary scales
# as D^(2*Lx) per layer, so it is tractable on tiny cells (3x3) but not the memory
# win. The boundary sweep replaces this for large cells while matching it exactly.
# Combined (D^2) index for the bond owning `li` (the shared combiner output).
_cidx(combs::Dict{Index,ITensor}, li::Index) = combinedind(combs[li])

# Exact two-site SVD recompression sweep around the column ring. With a tiny
# cutoff this is lossless (it only refactors), but it caps each horizontal bond
# at the true Schmidt rank so the boundary stays factorized instead of growing
# D^2 per absorbed row.
function _recompress_ring!(B::Vector{ITensor}; cutoff::Float64 = 1e-14)
    Lx = length(B)
    Lx >= 2 || return B
    sweep(order) = for x in order
        xn = x == Lx ? 1 : x + 1
        bond = commoninds(B[x], B[xn])
        isempty(bond) && continue
        T = @disable_warn_order B[x] * B[xn]
        linds = uniqueinds(B[x], B[xn])
        U, S, V, _, _, _ = svd(T, linds...; cutoff = cutoff)
        B[x] = U
        B[xn] = @disable_warn_order S * V
    end
    sweep(1:Lx)              # forward (recompresses the seam bond at x=Lx)
    sweep(reverse(1:Lx))     # backward, for a tighter gauge
    return B
end

# Build the column-ring boundary for the double layer and contract it to a scalar
# (Z if `op_site===nothing`, else Z_c with `projector_up()` inserted at `op_site`).
# Vertical sweep over rows: the boundary is a ring of Lx tensors; each carries its
# frontier vertical leg and a RENAMED copy of the row-1 bottom leg (the vertical
# wrap), closed by a delta at the end. The horizontal wrap closes via the ring.
function _boundary_scalar(
    psi::SquareIPEPSState,
    folded::Dict{SquareCoord,ITensor},
    combs::Dict{Index,ITensor};
    op_site::Union{Nothing,SquareCoord} = nothing,
    cutoff::Float64 = 1e-14,
)
    cell = psi.unitcell
    Lx, Ly = cell.Lx, cell.Ly
    siteop(c) = (op_site !== nothing && c == op_site) ? projector_up() : _identity_2x2()

    # Row 1: rename each bottom (wrap) leg to a fresh index so it does not clash
    # with the row-Ly frontier (which is the same torus Index) at closing time.
    B = Vector{ITensor}(undef, Lx)
    widx = Vector{Index}(undef, Lx)
    for x = 1:Lx
        c = SquareCoord(x, 1)
        E = _double_layer_site(psi, folded, c, combs; op = siteop(c))
        cd = _cidx(combs, link_index(psi, c, :down))
        w = Index(dim(cd), "wrap,$x")
        B[x] = replaceind(E, cd, w)
        widx[x] = w
    end

    # Absorb rows 2..Ly with a left-to-right ZIP so the running horizontal bond is
    # SVD-recompressed after every column. Forming the whole row at once would hold
    # Lx rank-6 tensors of size ~(D^2)^6 = D^12 simultaneously (OOM at D>=5); the
    # zip keeps only one such tensor live at a time (the seam column), then peels
    # an orthogonalized site off and carries a small bond onward.
    for y = 2:Ly
        Erow = [_double_layer_site(psi, folded, SquareCoord(x, y), combs; op = siteop(SquareCoord(x, y)))
                for x = 1:Lx]
        newB = Vector{ITensor}(undef, Lx)
        carry = nothing
        for x = 1:Lx
            # (carry * B[x]) * Erow[x]: carry already absorbed the left bonds into a
            # small index, so only the seam column (x=1, no carry) forms the big tensor.
            T = carry === nothing ? (@disable_warn_order B[x] * Erow[x]) :
                (@disable_warn_order (carry * B[x]) * Erow[x])
            if x < Lx
                xn = x + 1
                rights = unique(vcat(commoninds(T, B[xn]), commoninds(T, Erow[xn])))
                lefts = setdiff(collect(inds(T)), rights)
                U, S, V, _, _, _ = svd(T, lefts...; cutoff = cutoff)
                newB[x] = U
                carry = @disable_warn_order S * V
            else
                newB[Lx] = T            # keeps the seam bonds (shared with newB[1]) open
            end
        end
        B = newB
        _recompress_ring!(B; cutoff = cutoff)   # merges/recompresses all bonds incl. the seam
    end

    # Close the vertical wrap (row-Ly frontier == row-1 bottom == renamed widx).
    for x = 1:Lx
        fcu = _cidx(combs, link_index(psi, SquareCoord(x, Ly), :up))
        B[x] = @disable_warn_order B[x] * delta(fcu, widx[x])
    end
    # Contract the ring (the horizontal wrap closes on the final product).
    acc = B[1]
    for x = 2:Lx
        acc = @disable_warn_order acc * B[x]
    end
    return scalar(acc)
end

function _boundary_density_finite(psi::SquareIPEPSState; max_sites::Integer = 16)::Float64
    _check_tiny_finite_cell(psi, max_sites)
    folded = _absorb_all_weights_once(psi)
    combs = _bond_combiners(psi)
    reps = psi.unitcell.reps
    Z = _boundary_scalar(psi, folded, combs; op_site = nothing)
    abs(Z) > 1e-300 || throw(ArgumentError("double-layer norm is zero"))
    # Each site's Z_c is an independent boundary sweep (it builds its own tensors
    # and shares only the read-only `folded`/`combs`), so the N sweeps run in
    # parallel across Julia threads when available — the dominant cost is N+1
    # sweeps, and this turns it into ~wall-time of one. (Falls back to serial with
    # a single thread.) Caching per-row environments would cut it to ~2 sweeps; see
    # notes/stage2-truncation/improvement-roadmap.md.
    contribs = Vector{Float64}(undef, length(reps))
    Threads.@threads for i in eachindex(reps)
        contribs[i] = real(_boundary_scalar(psi, folded, combs; op_site = reps[i]) / Z)
    end
    return sum(contribs) / length(reps)
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
    exact_density_finite(psi; max_sites = 12, method = :auto)

Return the average exact finite contraction of the Rydberg density over all
unit-cell representatives of the supplied iPEPS state.

`method` selects the contraction backend; both are exact and agree to machine
precision:

- `:dense` — single-layer `dense_state_finite` (the original path). Builds the
  `2^N` state, so it is limited to small cells; its intermediate peaks at
  ~`2^(N/2)·D^(2Lx)`. On a large-memory host this comfortably reaches `D = 6` on a
  16-site (4×4) torus (a few seconds/call); on memory-constrained hosts it can OOM
  at `D ≥ 5` — pass `:boundary` there.
- `:boundary` — double-layer column-ring boundary contraction. Closes physical
  legs locally (no `2^N` factor) and keeps the boundary factorized, so its memory
  is bounded regardless of `D`. This is the path for cells too large for `2^N`
  (`> 16` sites) and the memory-safe fallback when dense OOMs. It is slower than
  dense on full-rank `D ≥ 5` bonds (per-site sweeps; environment caching is the
  documented next optimization).
- `:auto` (default) — `:dense` for `nsites ≤ 16` (fast where the `2^N` state fits),
  `:boundary` for larger cells. Never changes the result.
"""
function exact_density_finite(
    psi::SquareIPEPSState;
    max_sites::Integer = 12,
    method::Symbol = :auto,
)::Float64
    nsites = _check_tiny_finite_cell(psi, max_sites)
    chosen = method === :auto ? (nsites > 16 ? :boundary : :dense) : method
    chosen === :boundary && return _boundary_density_finite(psi; max_sites)
    chosen === :dense || throw(ArgumentError("method must be :auto, :dense, or :boundary"))
    # Build the exact dense state ONCE and read every site's density from it.
    # (The per-site exact_one_site_expectation_finite rebuilds the dense
    # contraction each call — 16x redundant on a 4x4 cell; this makes the
    # exact, CTM-free oracle tractable at 16 sites.)
    state = dense_state_finite(psi; max_sites)
    n = projector_up()
    reps = psi.unitcell.reps
    total = 0.0
    for c in reps
        positions = _local_positions(psi.unitcell, (c,))
        total += real(_local_expectation_from_state(state, nsites, positions, n))
    end
    return total / length(reps)
end

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
