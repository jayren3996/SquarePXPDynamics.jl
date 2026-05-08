module Gates

using LinearAlgebra
using ..Models

export dense_gate, projected_gate

function dense_gate(H::AbstractMatrix, step::Real; evolution::Symbol = :real)
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
    return Matrix{ComplexF64}(projector) * dense_gate(H, step; evolution)
end

end
