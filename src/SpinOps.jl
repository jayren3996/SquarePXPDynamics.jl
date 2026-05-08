module SpinOps

export pauli_x, pauli_y, pauli_z, identity2, projector_up, projector_down
export kron_all, embed_one_site

function identity2()
    return ComplexF64[1 0; 0 1]
end

function pauli_x()
    return ComplexF64[0 1; 1 0]
end

function pauli_y()
    return ComplexF64[0 -im; im 0]
end

function pauli_z()
    return ComplexF64[1 0; 0 -1]
end

function projector_up()
    return ComplexF64[1 0; 0 0]
end

function projector_down()
    return ComplexF64[0 0; 0 1]
end

function kron_all(ops::AbstractVector{<:AbstractMatrix})
    isempty(ops) && throw(ArgumentError("ops must be nonempty"))
    out = Matrix{ComplexF64}(ops[1])
    for op in ops[2:end]
        out = kron(out, Matrix{ComplexF64}(op))
    end
    return out
end

function embed_one_site(op::AbstractMatrix, site::Integer, nsites::Integer)
    1 <= site <= nsites || throw(ArgumentError("site must be in 1:nsites"))
    return kron_all([i == site ? op : identity2() for i in 1:nsites])
end

end
