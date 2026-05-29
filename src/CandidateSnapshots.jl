module CandidateSnapshots

using ITensors: ITensor, Index, dim
using JLD2
using ..SquareGeometry: SquareCoord
using ..SquareUnitCells: PeriodicSquareUnitCell, neighbor
using ..SquareIPEPS: SquareIPEPSState, BondKey, log_norm, state_version

export write_square_ipeps_snapshot, load_square_ipeps_snapshot
export SQUARE_IPEPS_SNAPSHOT_FORMAT_VERSION

"""
    SQUARE_IPEPS_SNAPSHOT_FORMAT_VERSION

Integer constant identifying the on-disk format used by
[`write_square_ipeps_snapshot`](@ref) and required by
[`load_square_ipeps_snapshot`](@ref). Loaders reject files whose stored
`format_version` does not match this constant — bump the value when the
serialized layout changes.
"""
const SQUARE_IPEPS_SNAPSHOT_FORMAT_VERSION = 1
const _PRIMARY_LINK_DIRECTIONS = (:right, :up)

"""
    write_square_ipeps_snapshot(path, psi::SquareIPEPSState)

Write a SquareIPEPSState to `path` as a JLD2 snapshot. Site tensors are saved
in canonical index order `(physical, left, right, up, down)`; link weights and
per-link dimensions are saved for the `:right` and `:up` representatives, with
`:left`/`:down` reconstructed as wrapped mirrors at load time. ITensors `Index`
objects are not stored — fresh indices are created on load using the same tag
scheme as `product_square_ipeps`. Format version
[`SQUARE_IPEPS_SNAPSHOT_FORMAT_VERSION`](@ref) is written into the file.
"""
function write_square_ipeps_snapshot(path::AbstractString, psi::SquareIPEPSState)
    cell = psi.unitcell
    site_arrays = Dict{Tuple{Int,Int},Array{ComplexF64,5}}()
    for c in cell.reps
        T = psi.tensors[c]
        p = psi.physical_indices[c]
        left = psi.link_indices[(c, :left)]
        right = psi.link_indices[(c, :right)]
        up = psi.link_indices[(c, :up)]
        down = psi.link_indices[(c, :down)]
        site_arrays[(c.x, c.y)] = ComplexF64.(Array(T, p, left, right, up, down))
    end

    link_weights = Dict{Tuple{Int,Int,Symbol},Vector{Float64}}()
    link_dims = Dict{Tuple{Int,Int,Symbol},Int}()
    for c in cell.reps
        for direction in _PRIMARY_LINK_DIRECTIONS
            key = BondKey(c, direction)
            link_weights[(c.x, c.y, direction)] = copy(psi.link_weights[key])
            link_dims[(c.x, c.y, direction)] = dim(psi.link_indices[(c, direction)])
        end
    end

    JLD2.jldsave(
        String(path);
        format_version = SQUARE_IPEPS_SNAPSHOT_FORMAT_VERSION,
        Lx = cell.Lx,
        Ly = cell.Ly,
        maxdim = psi.maxdim,
        gauge = psi.gauge,
        log_norm_value = log_norm(psi),
        mutation_version = state_version(psi),
        site_arrays = site_arrays,
        link_weights = link_weights,
        link_dims = link_dims,
    )
    return String(path)
end

"""
    load_square_ipeps_snapshot(path)::SquareIPEPSState

Reconstruct a SquareIPEPSState from a JLD2 snapshot written by
[`write_square_ipeps_snapshot`](@ref). Fresh ITensors `Index` objects are
created with the canonical tag scheme; indices are shared between neighboring
representatives so left/down link indices reuse the right/up indices of their
neighbors.
"""
function load_square_ipeps_snapshot(path::AbstractString)::SquareIPEPSState
    data = JLD2.load(String(path))

    format_version = data["format_version"]
    format_version == SQUARE_IPEPS_SNAPSHOT_FORMAT_VERSION || throw(
        ArgumentError(
            "unsupported SquareIPEPSState snapshot format_version $format_version " *
            "(expected $SQUARE_IPEPS_SNAPSHOT_FORMAT_VERSION)",
        ),
    )

    cell = PeriodicSquareUnitCell(Int(data["Lx"]), Int(data["Ly"]))

    saved_dims = data["link_dims"]
    saved_weights = data["link_weights"]
    link_indices = Dict{Tuple{SquareCoord,Symbol},Index}()
    link_weights = Dict{BondKey,Vector{Float64}}()
    for c in cell.reps
        right_dim = Int(saved_dims[(c.x, c.y, :right)])
        right_index = Index(right_dim, "link,$(c.x),$(c.y),right")
        link_indices[(c, :right)] = right_index
        link_indices[(neighbor(cell, c, :right), :left)] = right_index
        link_weights[BondKey(c, :right)] =
            Float64.(saved_weights[(c.x, c.y, :right)])

        up_dim = Int(saved_dims[(c.x, c.y, :up)])
        up_index = Index(up_dim, "link,$(c.x),$(c.y),up")
        link_indices[(c, :up)] = up_index
        link_indices[(neighbor(cell, c, :up), :down)] = up_index
        link_weights[BondKey(c, :up)] = Float64.(saved_weights[(c.x, c.y, :up)])
    end

    physical_indices = Dict{SquareCoord,Index}()
    tensors = Dict{SquareCoord,ITensor}()
    saved_site = data["site_arrays"]
    for c in cell.reps
        p = Index(2, "phys,$(c.x),$(c.y)")
        physical_indices[c] = p
        left = link_indices[(c, :left)]
        right = link_indices[(c, :right)]
        up = link_indices[(c, :up)]
        down = link_indices[(c, :down)]
        arr = saved_site[(c.x, c.y)]
        tensors[c] = ITensor(arr, p, left, right, up, down)
    end

    return SquareIPEPSState(
        cell,
        tensors,
        physical_indices,
        link_indices,
        link_weights,
        Int(data["maxdim"]),
        Symbol(data["gauge"]),
        Ref(Int(data["mutation_version"])),
        Ref(Float64(data["log_norm_value"])),
    )
end

end # module
