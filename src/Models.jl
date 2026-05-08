module Models

using LinearAlgebra
using ..SpinOps

export pxp_star_hamiltonian, blockade_projector, cluster_star_hamiltonian
export diagonal_star_hamiltonian, ising_bond_hamiltonian

const STAR_NSITES = 7
const CENTER_SITE = 1
const NEIGHBOR_SITES = 2:7

site_bit(state::Integer, site::Integer) = (state >> (STAR_NSITES - site)) & 1

function pxp_star_hamiltonian(projector::AbstractMatrix = projector_down(),
                              flip::AbstractMatrix = pauli_x())
    size(projector) == (2, 2) || throw(ArgumentError("projector must be 2x2"))
    size(flip) == (2, 2) || throw(ArgumentError("flip must be 2x2"))

    ops = Matrix{ComplexF64}[]
    push!(ops, Matrix{ComplexF64}(flip))
    append!(ops, [Matrix{ComplexF64}(projector) for _ in NEIGHBOR_SITES])
    return kron_all(ops)
end

function blockade_projector()
    diag = ones(ComplexF64, 2^STAR_NSITES)
    for state in 0:(2^STAR_NSITES - 1)
        forbidden = false
        for n in NEIGHBOR_SITES
            if site_bit(state, CENTER_SITE) == 0 && site_bit(state, n) == 0
                forbidden = true
            end
        end
        if forbidden
            diag[state + 1] = 0
        end
    end
    return Matrix(Diagonal(diag))
end

function cluster_star_hamiltonian()
    ops = Matrix{ComplexF64}[pauli_x()]
    append!(ops, [pauli_z() for _ in NEIGHBOR_SITES])
    return kron_all(ops)
end

function diagonal_star_hamiltonian()
    ops = Matrix{ComplexF64}[pauli_z()]
    append!(ops, [pauli_z() for _ in NEIGHBOR_SITES])
    return kron_all(ops)
end

function ising_bond_hamiltonian()
    return kron(pauli_z(), pauli_z())
end

end
