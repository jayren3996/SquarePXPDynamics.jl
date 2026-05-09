module SimpleUpdate

using LinearAlgebra
using ITensors
using ..Geometry: Coord, star_sites
using ..States: TriangularIPEPS, wrap_coord

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
    peeled = Matrix{ComplexF64}[]
    rest = Matrix{ComplexF64}(G)
    for _ in 1:(n - 1)
        m = size(rest, 1) ÷ 2
        rest_4d = reshape(rest, 2, m, 2, m)
        # Julia's kron indexing makes the rightmost two-state factor contiguous
        # in row/column storage, so this peels factors from right to left.
        rest_perm = permutedims(rest_4d, (1, 3, 2, 4))
        rest_mat = reshape(rest_perm, 4, m * m)
        F = svd(rest_mat)
        if length(F.S) > 1 && F.S[2] / max(F.S[1], eps()) > tol
            return nothing
        end
        s = F.S[1]
        u_right = reshape(F.U[:, 1] * sqrt(s), 2, 2)
        rest = reshape(F.Vt[1, :] * sqrt(s), m, m)
        push!(peeled, u_right)
    end
    push!(peeled, rest)
    factors = reverse(peeled)
    # Sanity check: reconstruct and compare.
    reconstructed = factors[1]
    for factor in factors[2:end]
        reconstructed = kron(reconstructed, factor)
    end
    if norm(reconstructed - G) > 1e-8 * max(norm(G), 1.0)
        return nothing
    end
    return factors
end

"""
    apply_star_gate_simple_update!(state, gate, center; cutoff=1e-12, maxdim=...) -> SimpleUpdateDiagnostics

Apply a 7-site star gate at `center`. The implementation supports these
correctness-validated paths:

1. **Identity gate**: detected numerically; no-op on the state.
2. **Site-product gate** `u_1 ⊗ ... ⊗ u_7`: each factor is applied to the
   physical index of the corresponding star site. For unit cells where
   multiple star positions wrap to the same representative, all wrapped
   factors must agree (translational invariance); otherwise an error is
   raised.

3. **Dense non-product gate**: `D=1` uses the exact dense product-state oracle.
   For `D>1`, the current Simple Update path uses a local product projection
   onto representative physical profiles while preserving the existing virtual
   bond dimensions up to `maxdim`.
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
    if factors === nothing && _is_d1_state(state)
        local_vectors = [_d1_site_vector(state, wrap_coord(state.unitcell, sc)) for sc in star]
        psi = local_vectors[1]
        for v in local_vectors[2:end]
            psi = kron(psi, v)
        end

        phi = G * psi
        product_factors = _try_factorize_product_state(phi, _STAR_NSITES)
        if product_factors !== nothing
            for (rep, positions) in _star_positions_by_rep(state, star)
                v = _common_factor_for_positions(product_factors, positions)
                _set_d1_site_vector!(state, rep, v)
            end
        else
            tensor_phi = reshape(phi, ntuple(_ -> 2, _STAR_NSITES)...)
            for (rep, positions) in _star_positions_by_rep(state, star)
                rho = _one_site_density_from_star_tensor(tensor_phi, positions[1])
                vals, vecs = eigen(Hermitian(rho))
                v = vecs[:, argmax(vals)]
                _set_d1_site_vector!(state, rep, ComplexF64.(v))
            end
        end

        dims = [dim(state.bond_inds[b]) for b in affected]
        return SimpleUpdateDiagnostics(0.0, affected, dims)
    end

    if factors === nothing
        return _apply_general_star_gate_simple_update!(
            state, G, center; cutoff = cutoff, maxdim = maxdim
        )
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

function _apply_general_star_gate_simple_update!(state::TriangularIPEPS,
                                                G::Matrix{ComplexF64},
                                                center::Coord;
                                                cutoff::Real,
                                                maxdim::Union{Nothing,Integer})
    size(G) == (128, 128) || throw(ArgumentError("gate must be 128x128"))
    maxdim === nothing && throw(ArgumentError("maxdim is required for general star updates"))
    maxdim >= 1 || throw(ArgumentError("maxdim must be >= 1"))
    cutoff >= 0 || throw(ArgumentError("cutoff must be nonnegative"))

    affected = _affected_star_bonds(state, center)
    for bond in affected
        current_dim = dim(state.bond_inds[bond])
        current_dim <= maxdim || throw(ArgumentError(
            "current bond dimension $current_dim exceeds requested maxdim $maxdim; " *
            "dimension-changing star writeback is not implemented yet"))
    end

    star = star_sites(center)
    reps = [wrap_coord(state.unitcell, sc) for sc in star]
    local_vectors = [_dominant_site_vector(state, rep) for rep in reps]
    psi = local_vectors[1]
    for v in local_vectors[2:end]
        psi = kron(psi, v)
    end

    phi = G * psi
    targets = _product_projection_targets(phi, _STAR_NSITES)
    projected = targets[1]
    for v in targets[2:end]
        projected = kron(projected, v)
    end
    residual = _relative_residual(phi, projected)

    for rep in unique(reps)
        old = _dominant_site_vector(state, rep)
        target = _representative_target(rep, reps, targets)
        op = _regularized_physical_map(old, target)
        _apply_physical_map!(state, rep, op)
        _normalize_site_tensor!(state, rep)
    end

    dims = [dim(state.bond_inds[b]) for b in affected]
    return SimpleUpdateDiagnostics(Float64(residual), affected, dims)
end

function _affected_star_bonds(state::TriangularIPEPS, center::Coord)
    rep = wrap_coord(state.unitcell, center)
    return [(rep, d) for d in 1:6]
end

function _is_d1_state(state::TriangularIPEPS)
    for rep in keys(state.tensors)
        ph = state.phys_inds[rep]
        for idx in inds(state.tensors[rep])
            idx == ph && continue
            dim(idx) == 1 || return false
        end
    end
    return true
end

function _d1_site_vector(state::TriangularIPEPS, rep::Coord)
    data = vec(array(state.tensors[rep]))
    length(data) == 2 || error("expected D=1 tensor to have exactly two entries")
    return ComplexF64[data[1], data[2]]
end

function _star_positions_by_rep(state::TriangularIPEPS, star)
    grouped = Dict{Coord,Vector{Int}}()
    for (i, sc) in enumerate(star)
        rep = wrap_coord(state.unitcell, sc)
        push!(get!(grouped, rep, Int[]), i)
    end
    return grouped
end

function _try_factorize_product_state(vec::AbstractVector, n::Int; tol::Real = 1e-10)
    length(vec) == 2^n || throw(ArgumentError("vec length must be 2^n"))
    factors = Vector{ComplexF64}[]
    rest = ComplexF64.(vec)
    for _ in 1:(n - 1)
        block = length(rest) ÷ 2
        rest_mat = Matrix{ComplexF64}(undef, 2, block)
        rest_mat[1, :] .= view(rest, 1:block)
        rest_mat[2, :] .= view(rest, (block + 1):(2 * block))
        F = svd(rest_mat)
        if length(F.S) > 1 && F.S[2] / max(F.S[1], eps()) > tol
            return nothing
        end
        s = F.S[1]
        push!(factors, ComplexF64.(F.U[:, 1] * sqrt(s)))
        rest = ComplexF64.(F.Vt[1, :] * sqrt(s))
    end
    push!(factors, rest)

    reconstructed = factors[1]
    for factor in factors[2:end]
        reconstructed = kron(reconstructed, factor)
    end
    if norm(reconstructed - vec) > 1e-8 * max(norm(vec), 1.0)
        return nothing
    end
    return factors
end

function _common_factor_for_positions(factors::Vector{Vector{ComplexF64}},
                                      positions::Vector{Int})
    ref = factors[positions[1]]
    ref_norm2 = sum(abs2, ref)
    ref_norm2 < 1e-28 && throw(ArgumentError("cannot use a zero product factor"))
    aligned = zeros(ComplexF64, 2)
    for pos in positions
        f = factors[pos]
        c = sum(conj.(ref) .* f) / ref_norm2
        if norm(f - c * ref) > 1e-8 * max(norm(ref), 1.0)
            return ref
        end
        aligned .+= f / c
    end
    return aligned / length(positions)
end

function _one_site_density_from_star_tensor(tensor_phi, position::Int)
    others = [i for i in 1:_STAR_NSITES if i != position]
    perm = (position, others...)
    psi = reshape(permutedims(tensor_phi, perm), 2, :)
    rho = psi * psi'
    tr = real(sum(diag(rho)))
    tr == 0 && return Matrix{ComplexF64}(I, 2, 2) / 2
    return rho / tr
end

function _set_d1_site_vector!(state::TriangularIPEPS, rep::Coord, v::Vector{ComplexF64})
    ph = state.phys_inds[rep]
    binds = Tuple(idx for idx in inds(state.tensors[rep]) if idx != ph)
    T = ITensor(ComplexF64, ph, binds...)
    nrm = norm(v)
    nrm == 0 && throw(ArgumentError("cannot set zero local vector"))
    v = v / nrm
    bind_assignments = [bind => 1 for bind in binds]
    for k in 1:2
        T[ph => k, bind_assignments...] = v[k]
    end
    state.tensors[rep] = T
    return nothing
end

function _dominant_site_vector(state::TriangularIPEPS, rep::Coord)
    T = state.tensors[rep]
    ph = state.phys_inds[rep]
    idxs = collect(inds(T))
    phpos = findfirst(i -> i == ph, idxs)
    phpos === nothing && error("physical index missing from tensor at $rep")
    data = array(T)
    perm = (phpos, (i for i in eachindex(idxs) if i != phpos)...)
    mat = reshape(permutedims(data, perm), 2, :)
    rho = mat * mat'
    trρ = real(tr(rho))
    if !(isfinite(trρ)) || trρ <= eps(Float64)
        return ComplexF64[1, 0]
    end
    vals, vecs = eigen(Hermitian(rho / trρ))
    v = ComplexF64.(vecs[:, argmax(vals)])
    return v / max(norm(v), eps(Float64))
end

function _product_projection_targets(vec::Vector{ComplexF64}, n::Int)
    tensor = reshape(vec, ntuple(_ -> 2, n)...)
    targets = Vector{ComplexF64}[]
    for pos in 1:n
        rho = _one_site_density_from_star_tensor(tensor, pos)
        vals, vecs = eigen(Hermitian(rho))
        v = ComplexF64.(vecs[:, argmax(vals)])
        push!(targets, v / max(norm(v), eps(Float64)))
    end
    return targets
end

function _representative_target(rep::Coord,
                                reps::Vector{Coord},
                                targets::Vector{Vector{ComplexF64}})
    positions = findall(==(rep), reps)
    ref = targets[positions[1]]
    aligned = zeros(ComplexF64, 2)
    for pos in positions
        v = targets[pos]
        phase = dot(ref, v)
        if abs(phase) > eps(Float64)
            v = v * exp(-im * angle(phase))
        end
        aligned .+= v
    end
    nrm = norm(aligned)
    return nrm <= eps(Float64) ? ref : aligned / nrm
end

function _regularized_physical_map(old::Vector{ComplexF64},
                                   target::Vector{ComplexF64};
                                   regularization::Float64 = 1e-10)
    oldn = old / max(norm(old), eps(Float64))
    targetn = target / max(norm(target), eps(Float64))
    projector = oldn * oldn'
    return targetn * oldn' + regularization * (Matrix{ComplexF64}(I, 2, 2) - projector)
end

function _apply_physical_map!(state::TriangularIPEPS, rep::Coord, op::Matrix{ComplexF64})
    T = state.tensors[rep]
    ph = state.phys_inds[rep]
    opT = ITensor(op, prime(ph), ph)
    state.tensors[rep] = noprime(opT * T)
    return nothing
end

function _normalize_site_tensor!(state::TriangularIPEPS, rep::Coord)
    nrm = sqrt(real(scalar(dag(state.tensors[rep]) * state.tensors[rep])))
    if isfinite(nrm) && nrm > eps(Float64)
        state.tensors[rep] ./= nrm
    end
    return nothing
end

function _relative_residual(target::Vector{ComplexF64}, projected::Vector{ComplexF64})
    denom = max(real(dot(target, target)), eps(Float64))
    scale = dot(projected, target) / max(sum(abs2, projected), eps(Float64))
    return Float64(real(sum(abs2, target - scale * projected) / denom))
end

end
