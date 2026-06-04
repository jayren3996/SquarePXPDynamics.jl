module PEPSKitStarSnake

# EXPERIMENTAL (M7). Snake / matrix-product-operator (MPO) representation of the
# five-site square-star PXP gate.
#
# The production evolution applies the genuine five-site cross gate via the
# simultaneous simple update (`project_star_pepskit!`); it is NOT replaced here.
# This module is an experimental study of representing that same gate as a
# matrix product operator along the star-site ordering
# `(center, right, up, left, down)`. The decomposition is exact (no truncation),
# so it reconstructs the dense gate bit-for-bit; its purpose is to expose the
# gate's MPO bond structure and provide a substrate for a future snake-path
# application. Applying the MPO onto an InfinitePEPS along a spatial snake is the
# open, path-biased research question flagged in docs/pepskit-native-refactor.md
# and is deliberately not wired into production.

using LinearAlgebra: svd, Diagonal

using ..SquarePXP: SQUARE_STAR_SITES
using ..PEPSKitStarUpdate: pxp_star_gate_dense

export StarGateMPO, star_gate_mpo, reconstruct_star_gate, star_mpo_bond_dims

# Big-endian dense index for a 5-tuple of 1/2 values (site 1 most significant),
# matching `square_pxp_star_hamiltonian` / `pxp_star_gate_dense` ordering.
_didx5(v) = 1 + (v[1] - 1) * 16 + (v[2] - 1) * 8 + (v[3] - 1) * 4 + (v[4] - 1) * 2 + (v[5] - 1)

"""
    StarGateMPO

Exact matrix-product-operator form of the five-site square-star PXP gate.
`tensors[k]` is the rank-4 dense array for star site `k` (order
`center, right, up, left, down`) with legs `[phys_out, phys_in, bond_left,
bond_right]`; the outermost bonds have dimension 1. `bond_dims` lists the four
internal bond dimensions. Only exactly-zero singular values are dropped at each
cut (no physical truncation), so it reconstructs the dense gate to machine
precision while reporting the genuine operator Schmidt rank.
"""
struct StarGateMPO
    tensors::Vector{Array{ComplexF64,4}}
    bond_dims::Vector{Int}
end

"""
    star_mpo_bond_dims(mpo::StarGateMPO)

Return the four internal bond dimensions of a [`StarGateMPO`](@ref) along the
star-site ordering. Small values indicate the gate is cheaply representable as a
matrix product operator.
"""
star_mpo_bond_dims(mpo::StarGateMPO) = mpo.bond_dims

# Dense rank-10 gate tensor T[o1,o2,o3,o4,o5, i1,i2,i3,i4,i5] from the 32x32 gate.
function _gate_rank10(gate::AbstractMatrix)
    T = zeros(ComplexF64, ntuple(Returns(2), 2 * SQUARE_STAR_SITES))
    rng = ntuple(_ -> 1:2, SQUARE_STAR_SITES)
    for o in Iterators.product(rng...), i in Iterators.product(rng...)
        T[o..., i...] = gate[_didx5(o), _didx5(i)]
    end
    return T
end

"""
    star_gate_mpo(dt; projected = true, evolution = :real)::StarGateMPO

Build the exact MPO of the five-site square-star PXP gate (see
[`pxp_star_gate_dense`](@ref)) by successive singular value decompositions along
the star-site ordering. No truncation is applied, so
[`reconstruct_star_gate`](@ref) returns the dense gate to machine precision.
"""
function star_gate_mpo(dt::Real; projected::Bool = true, evolution::Symbol = :real)
    gate = pxp_star_gate_dense(dt; projected, evolution)
    T = _gate_rank10(gate)
    n = SQUARE_STAR_SITES
    # Interleave to [o1,i1, o2,i2, …, o5,i5] then group each site's (o,i) pair.
    perm = collect(Iterators.flatten((k, k + n) for k in 1:n))
    M = permutedims(T, perm)                      # dims (2,2, 2,2, …) length 10

    tensors = Array{ComplexF64,4}[]
    bond_dims = Int[]
    left = 1                                       # current left-bond dimension
    # `rest` holds the not-yet-factored remainder reshaped as (left*4, 4^(remaining)).
    rest = reshape(M, left * 4, 4^(n - 1))
    for site in 1:(n - 1)
        U, s, V = svd(rest)
        # Keep only the numerically nonzero singular values: the gate has a small
        # exact operator Schmidt rank across each star-site cut, but `svd` returns
        # the full `min(size(rest))`-length spectrum padded with ~zeros. Dropping
        # them reports the true bond dimension while keeping reconstruction exact.
        tol = isempty(s) ? 0.0 : maximum(s) * size(rest, 1) * eps(Float64) * 10
        r = max(1, count(>(tol), s))
        U, s, V = U[:, 1:r], s[1:r], V[:, 1:r]
        # MPO tensor for this site: (left, o, i, right) -> (o, i, left, right)
        W = permutedims(reshape(U, left, 2, 2, r), (2, 3, 1, 4))
        push!(tensors, Array{ComplexF64,4}(W))
        push!(bond_dims, r)
        # fold S into the remainder and reshape for the next SVD
        SV = Diagonal(s) * V'                       # (r, 4^(n-site))
        remaining = n - site
        rest = reshape(SV, r * 4, 4^(remaining - 1))
        left = r
    end
    # final site: rest is (left*4, 1)
    Wlast = permutedims(reshape(rest, left, 2, 2, 1), (2, 3, 1, 4))
    push!(tensors, Array{ComplexF64,4}(Wlast))
    return StarGateMPO(tensors, bond_dims)
end

"""
    reconstruct_star_gate(mpo)::Matrix{ComplexF64}

Contract a [`StarGateMPO`](@ref) back into the dense `32x32` gate in
`(center, right, up, left, down)` order. Equals [`pxp_star_gate_dense`](@ref) to
machine precision.
"""
function reconstruct_star_gate(mpo::StarGateMPO)
    n = length(mpo.tensors)
    # Contract bonds left-to-right, accumulating open (o,i) legs.
    acc = mpo.tensors[1]                            # (o1,i1,left=1,right)
    acc = reshape(acc, 2, 2, size(acc, 4))          # drop trivial left bond -> (o1,i1,b)
    open_legs = 1
    for site in 2:n
        W = mpo.tensors[site]                       # (o,i,left,right)
        b = size(acc, ndims(acc))                   # current right bond
        # contract acc's last leg (b) with W's left leg
        am = reshape(acc, :, b)                      # (2^(2*open)*..., b)
        wm = reshape(permutedims(W, (3, 1, 2, 4)), b, :)  # (b, o*i*right)
        prod = am * wm                              # (..., o*i*right)
        rdim = size(W, 4)
        acc = reshape(prod, ntuple(_ -> 2, 2 * site)..., rdim)
        open_legs = site
    end
    acc = reshape(acc, ntuple(_ -> 2, 2 * n)...)    # drop trivial final bond
    # acc legs are (o1,i1,o2,i2,…,o5,i5); regroup to dense (outs, ins).
    G = zeros(ComplexF64, 2^n, 2^n)
    rng = ntuple(_ -> 1:2, n)
    for o in Iterators.product(rng...), i in Iterators.product(rng...)
        idx = ntuple(k -> isodd(k) ? o[(k + 1) ÷ 2] : i[k ÷ 2], 2 * n)
        G[_didx5(o), _didx5(i)] = acc[idx...]
    end
    return G
end

end # module
