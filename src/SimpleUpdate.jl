module SimpleUpdate

using LinearAlgebra
using ITensors
using ..Geometry: Coord, star_sites
using ..States: TriangularIPEPS, wrap_coord, opposite_direction

export SimpleUpdateDiagnostics, apply_star_gate_simple_update!

"""
    SimpleUpdateDiagnostics

Bookkeeping returned by `apply_star_gate_simple_update!`.

Fields:
- `discarded_weight`: total truncated weight across all bonds in this step
  (zero for the implemented site-product / identity paths).
- `affected_bonds`: list of `(rep, dir)` bond keys touched by the update.
- `output_bond_dims`: post-update dimensions of the affected bonds, in the
  same order as `affected_bonds`.
"""
struct SimpleUpdateDiagnostics
    discarded_weight::Float64
    affected_bonds::Vector{Tuple{Coord,Int}}
    output_bond_dims::Vector{Int}
end

const _STAR_NSITES = 7

"""
    _try_factorize_product_gate(G, n) -> Vector{Matrix{ComplexF64}} | nothing

Try to write a `2^n x 2^n` matrix `G` as a tensor product `u_1 ⊗ ... ⊗ u_n`
of `2x2` matrices. Returns the list of factors on success; `nothing`
otherwise. The factors satisfy `kron(u_1, kron(u_2, ..., u_n)) ≈ G`.
"""
function _try_factorize_product_gate(G::AbstractMatrix, n::Int; tol::Real = 1e-10)
    @assert size(G, 1) == size(G, 2) == 2^n
    factors = Matrix{ComplexF64}[]
    rest = Matrix{ComplexF64}(G)
    for _ in 1:(n - 1)
        m = size(rest, 1) ÷ 2
        rest_4d = reshape(rest, 2, m, 2, m)
        # Bring (out_first, in_first) to the front: (1,3,2,4)
        rest_perm = permutedims(rest_4d, (1, 3, 2, 4))
        rest_mat = reshape(rest_perm, 4, m * m)
        F = svd(rest_mat)
        if length(F.S) > 1 && F.S[2] / max(F.S[1], eps()) > tol
            return nothing
        end
        s = F.S[1]
        u_first = reshape(F.U[:, 1] * sqrt(s), 2, 2)
        rest = reshape(F.Vt[1, :] * sqrt(s), m, m)
        push!(factors, u_first)
    end
    push!(factors, rest)
    # Sanity check: reconstruct and compare.
    reconstructed = factors[end]
    for k in (n - 1):-1:1
        reconstructed = kron(factors[k], reconstructed)
    end
    if norm(reconstructed - G) > 1e-8 * max(norm(G), 1.0)
        return nothing
    end
    return factors
end

"""
    apply_star_gate_simple_update!(state, gate, center; cutoff=1e-12, maxdim=...) -> SimpleUpdateDiagnostics

Apply a 7-site star gate at `center`. The current implementation is
conservative and supports two correctness-validated paths:

1. **Identity gate**: detected numerically; no-op on the state.
2. **Site-product gate** `u_1 ⊗ ... ⊗ u_7`: each factor is applied to the
   physical index of the corresponding star site. For unit cells where
   multiple star positions wrap to the same representative, all wrapped
   factors must agree (translational invariance); otherwise an error is
   raised.

A general 7-site gate that does not factorize into a site-product is
*not yet implemented*: a clear error is raised. This is the documented
truncation gap until full Simple Update SVD bookkeeping lands.
"""
function apply_star_gate_simple_update!(state::TriangularIPEPS,
                                        gate::AbstractMatrix,
                                        center::Coord;
                                        cutoff::Real = 1e-12,
                                        maxdim::Union{Nothing,Integer} = nothing)
    size(gate) == (128, 128) || throw(ArgumentError("gate must be 128x128"))

    star = star_sites(center)
    affected = Tuple{Coord,Int}[]
    for d in 1:6
        push!(affected, (wrap_coord(state.unitcell, center), d))
    end

    G = Matrix{ComplexF64}(gate)
    Iref = Matrix{ComplexF64}(I, 128, 128)
    if norm(G - Iref) <= 1e-12 * max(norm(Iref), 1.0)
        dims = [dim(state.bond_inds[b]) for b in affected]
        return SimpleUpdateDiagnostics(0.0, affected, dims)
    end

    factors = _try_factorize_product_gate(G, _STAR_NSITES)
    if factors === nothing
        error("apply_star_gate_simple_update!: general non-product 7-site gates " *
              "are not yet implemented (full SU SVD bookkeeping is the documented " *
              "remaining truncation gap). Provide an identity or site-product gate, " *
              "or implement the missing path.")
    end

    # Group factors by representative; require all factors of a given rep to be
    # complex-scalar multiples of one another (translational invariance up to
    # gauge). Within each rep, recover a single canonical operator by averaging
    # the rescaled factors.
    rep_factors_list = Dict{Coord,Vector{Matrix{ComplexF64}}}()
    for (i, sc) in enumerate(star)
        rep = wrap_coord(state.unitcell, sc)
        push!(get!(rep_factors_list, rep, Matrix{ComplexF64}[]), factors[i])
    end

    rep_factor = Dict{Coord,Matrix{ComplexF64}}()
    for (rep, fs) in rep_factors_list
        ref = fs[1]
        ref_norm2 = sum(abs2, ref)
        ref_norm2 < 1e-28 && error(
            "apply_star_gate_simple_update!: vanishing factor at rep $rep")
        for f in fs
            c = sum(conj.(ref) .* f) / ref_norm2
            if norm(f - c * ref) > 1e-8 * max(norm(ref), 1.0)
                error("apply_star_gate_simple_update!: gate factors at star " *
                      "positions sharing unit-cell rep $rep are not scalar " *
                      "multiples of a common single-site operator " *
                      "(translational invariance violated).")
            end
        end
        # Pick a canonical representative scaled to unit Frobenius norm. The
        # absolute scale is absorbed into a state-wide normalization that
        # cancels in any expectation value.
        rep_factor[rep] = ref / sqrt(ref_norm2 / 2)
    end

    for (rep, u) in rep_factor
        T = state.tensors[rep]
        ph = state.phys_inds[rep]
        u_T = ITensor(u, prime(ph), ph)
        new_T = noprime(u_T * T)
        state.tensors[rep] = new_T
    end

    dims = [dim(state.bond_inds[b]) for b in affected]
    return SimpleUpdateDiagnostics(0.0, affected, dims)
end

end
