module SolvableModels

using LinearAlgebra
using ..SpinOps

export cluster_center_z_expectation_exact, cluster_star_hamiltonian, dense_gate

const STAR_NEIGHBORS = 6

function cluster_star_hamiltonian()
    ops = Matrix{ComplexF64}[pauli_x()]
    append!(ops, [pauli_z() for _ in 1:STAR_NEIGHBORS])
    return kron_all(ops)
end

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

"""
    cluster_center_z_expectation_exact(t; initial = :z_plus)

Return the exact center-site `Z_c` expectation `cos(2t)` for the involutory
cluster-star Hamiltonian `K = X_c prod_j Z_j`, with the center and all
neighbors initialized in `+1` eigenstates of `Z`.

This helper is a small analytic regression target for gate signs and benchmark
plumbing; it is not a general triangular-lattice stabilizer solver.
"""
function cluster_center_z_expectation_exact(t::Real; initial::Symbol = :z_plus)
    if initial == :z_plus
        return cos(2t)
    else
        throw(ArgumentError("supported initial states: :z_plus"))
    end
end

end
