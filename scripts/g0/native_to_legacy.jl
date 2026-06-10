# Native (PXPIPEPSState / PEPSKit) -> legacy (SquareIPEPSState / ITensors) state
# adapter, for running the legacy EXACT finite-torus observables — in particular
# the memory-bounded double-layer `:boundary` oracle (Observables.jl, 2026-06-03)
# — on natively-evolved states with cells > 16 sites (e.g. the 6x6 torus).
#
# Construction: the contractable state is `measurement_peps(state)` (convention
# B, each Γ carries sqrt(λ) per leg), so the legacy state gets those tensors
# verbatim with ALL-ONES link weights — `_absorb_all_weights_once` is then a
# no-op and the legacy folded network is bit-identical to the native one.
#
# Leg/coordinate conventions (PEPSKitBackend.squarecoord_to_cartesianindex):
#   native peps[row, col], row = Ly - y + 1, col = x; legs [P, N, E, S, W];
#   legacy SquareCoord(x, y): N=up, E=right, S=down, W=left.
# Link sharing: (c,:right) is the same Index as (right-neighbor,:left); (c,:up)
# same as (up-neighbor,:down) — matching native E==east-neighbor W, N==north S.
#
# Gate (legacy_boundary_gate.jl): on 4x4 the adapter must reproduce the native
# exact_density_finite through BOTH legacy paths (:dense and :boundary) to ~1e-9.

using SquarePXPDynamics
using ITensors
const _idim = ITensors.dim   # `dim` clashes with TensorKit's when both are loaded

function native_to_legacy(state::PXPIPEPSState)::SquareIPEPSState
    cell = state.unitcell
    peps = measurement_peps(state)

    arrs = Dict{SquareCoord,Array{ComplexF64,5}}()
    for c in cell.reps
        idx = squarecoord_to_cartesianindex(cell, c)
        arrs[c] = convert(Array{ComplexF64,5}, convert(Array, peps[idx[1], idx[2]]))
    end

    links = Dict{Tuple{SquareCoord,Symbol},Index}()
    weights = Dict{BondKey,Vector{Float64}}()
    for c in cell.reps
        A = arrs[c]
        dN, dE = size(A, 2), size(A, 3)
        ri = Index(dE, "link,$(c.x),$(c.y),right")
        links[(c, :right)] = ri
        links[(neighbor(cell, c, :right), :left)] = ri
        weights[BondKey(c, :right)] = ones(Float64, dE)
        ui = Index(dN, "link,$(c.x),$(c.y),up")
        links[(c, :up)] = ui
        links[(neighbor(cell, c, :up), :down)] = ui
        weights[BondKey(c, :up)] = ones(Float64, dN)
    end

    physical = Dict(c => Index(2, "phys,$(c.x),$(c.y)") for c in cell.reps)
    tensors = Dict{SquareCoord,ITensor}()
    maxd = 1
    for c in cell.reps
        A = arrs[c]
        p = physical[c]
        l, r = links[(c, :left)], links[(c, :right)]
        u, d = links[(c, :up)], links[(c, :down)]
        size(A) == (2, _idim(u), _idim(r), _idim(d), _idim(l)) || throw(ArgumentError(
            "bond-dimension mismatch at $c: array $(size(A)) vs links " *
            "(2, $(_idim(u)), $(_idim(r)), $(_idim(d)), $(_idim(l))) — native bond dims " *
            "must agree across shared bonds"))
        tensors[c] = ITensor(A, p, u, r, d, l)   # A axes [P,N,E,S,W] = (p,u,r,d,l)
        maxd = max(maxd, _idim(u), _idim(r))
    end

    return SquareIPEPSState(
        cell, tensors, physical, links, weights, maxd, :simple, Ref(0), Ref(0.0),
    )
end
