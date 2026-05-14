module SpinOps

export pauli_x, pauli_y, pauli_z, identity2, projector_up, projector_down
export kron_all, embed_one_site

"""
    identity2()

Return the 2x2 identity matrix as `ComplexF64`.
"""
function identity2()
    return ComplexF64[1 0; 0 1]
end

"""
    pauli_x()

Return the spin-1/2 Pauli X matrix as `ComplexF64`.
"""
function pauli_x()
    return ComplexF64[0 1; 1 0]
end

"""
    pauli_y()

Return the spin-1/2 Pauli Y matrix as `ComplexF64`.
"""
function pauli_y()
    return ComplexF64[0 -im; im 0]
end

"""
    pauli_z()

Return the spin-1/2 Pauli Z matrix as `ComplexF64`.
"""
function pauli_z()
    return ComplexF64[1 0; 0 -1]
end

"""
    projector_up()

Return the single-site projector onto `|up> = |0>`.
"""
function projector_up()
    return ComplexF64[1 0; 0 0]
end

"""
    projector_down()

Return the single-site projector onto `|down> = |1>`.
"""
function projector_down()
    return ComplexF64[0 0; 0 1]
end

"""
    kron_all(ops)

Return the left-to-right Kronecker product of a nonempty vector of matrices.
All inputs are converted to `ComplexF64` matrices before multiplication.
"""
function kron_all(ops::AbstractVector{<:AbstractMatrix})
    isempty(ops) && throw(ArgumentError("ops must be nonempty"))
    out = Matrix{ComplexF64}(ops[1])
    for op in ops[2:end]
        out = kron(out, Matrix{ComplexF64}(op))
    end
    return out
end

"""
    embed_one_site(op, site, nsites)

Embed a 2x2 single-site operator into an `nsites` spin-1/2 Hilbert space.
Site numbering is one-based and follows the left-to-right tensor-product order.
"""
function embed_one_site(op::AbstractMatrix, site::Integer, nsites::Integer)
    size(op) == (2, 2) || throw(ArgumentError("op must be 2x2"))
    1 <= site <= nsites || throw(ArgumentError("site must be in 1:nsites"))
    return kron_all([i == site ? op : identity2() for i in 1:nsites])
end

end
