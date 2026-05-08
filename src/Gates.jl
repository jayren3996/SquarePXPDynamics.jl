module Gates

using LinearAlgebra
using ..Models

export dense_gate, projected_gate

function validate_star_hamiltonian(H::AbstractMatrix)
    size(H) == (128, 128) || throw(ArgumentError("H must be 128x128"))
    return nothing
end

function dense_gate(H::AbstractMatrix, step::Real; evolution::Symbol = :real)
    validate_star_hamiltonian(H)

    if evolution == :real
        return exp(-im * step * Matrix{ComplexF64}(H))
    elseif evolution == :imaginary
        return exp(-step * Matrix{ComplexF64}(H))
    else
        throw(ArgumentError("evolution must be :real or :imaginary"))
    end
end

function projected_gate(H::AbstractMatrix, step::Real;
                        evolution::Symbol = :real,
                        projector::AbstractMatrix = blockade_projector())
    validate_star_hamiltonian(H)
    size(projector) == size(H) || throw(ArgumentError("projector must match H dimensions"))

    return Matrix{ComplexF64}(projector) * dense_gate(H, step; evolution)
end

end
