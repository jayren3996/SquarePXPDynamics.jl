module Observables

using LinearAlgebra
using ITensors
using ..Geometry: Coord, neighbor
using ..States: TriangularIPEPS, wrap_coord
using ..Models: blockade_projector
using ..SpinOps: projector_up

export local_expectation, tensor_norm, dense_blockade_violations
export local_blockade_violation, mean_blockade_violation

"""
    local_expectation(state, c, op) -> ComplexF64

Single-site expectation `<op>` at coordinate `c`. Treats the bond environment
as trivial — exact for `D=1` product states; for `D>1` it is a site-only
diagnostic, not the full PEPS expectation value.
"""
function local_expectation(state::TriangularIPEPS, c::Coord, op::AbstractMatrix)
    size(op) == (2, 2) || throw(ArgumentError("op must be 2x2"))
    rep = wrap_coord(state.unitcell, c)
    T = state.tensors[rep]
    ph = state.phys_inds[rep]
    ph_prime = prime(ph)
    op_T = ITensor(Matrix{ComplexF64}(op), ph_prime, ph)
    Tdag = prime(dag(T), ph)
    num = scalar(Tdag * op_T * T)
    denom = scalar(dag(T) * T)
    return num / denom
end

"""
    tensor_norm(state, c) -> Float64

Frobenius norm of the representative tensor at the unit-cell rep of `c`.
For a `D=1` product state with normalized local vectors, this equals 1.
"""
function tensor_norm(state::TriangularIPEPS, c::Coord)
    rep = wrap_coord(state.unitcell, c)
    T = state.tensors[rep]
    return sqrt(real(scalar(dag(T) * T)))
end

"""
    dense_blockade_violations(vec) -> Float64

Total weight of `vec` outside the 7-site star blockade subspace. Returns
`<vec| (I - P) |vec> / <vec|vec>` where `P` is the blockade projector.
Zero iff every nonzero amplitude is on a blockade-allowed configuration.
"""
function dense_blockade_violations(vec::AbstractVector)
    length(vec) == 128 || throw(ArgumentError("vec must have length 128"))
    P = blockade_projector()
    nrm2 = real(dot(vec, vec))
    nrm2 == 0 && return 0.0
    proj = P * vec
    allowed = real(dot(vec, proj))
    return real((nrm2 - allowed) / nrm2)
end

"""
    local_blockade_violation(state, c, d) -> Float64

Approximate nearest-neighbor blockade violation on bond `(c, neighbor(c, d))`
using the product of local `|up><up|` expectations. This is exact for `D=1`
product states and a local-environment diagnostic for `D>1`.
"""
function local_blockade_violation(state::TriangularIPEPS, c::Coord, d::Integer)
    1 <= d <= 6 || throw(ArgumentError("direction must be in 1:6"))
    pup_c = real(local_expectation(state, c, projector_up()))
    pup_n = real(local_expectation(state, neighbor(c, d), projector_up()))
    return pup_c * pup_n
end

function mean_blockade_violation(state::TriangularIPEPS, centers)
    vals = Float64[]
    for c in centers, d in 1:6
        push!(vals, local_blockade_violation(state, c, d))
    end
    return isempty(vals) ? 0.0 : sum(vals) / length(vals)
end

end
