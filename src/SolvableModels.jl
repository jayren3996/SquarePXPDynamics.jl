module SolvableModels

export stabilizer_expectation_exact

"""
    stabilizer_expectation_exact(t; initial = :z_plus)

Return the exact single-star benchmark expectation `cos(2t)` for the involutory
cluster-star Hamiltonian `K = X_c prod_j Z_j`, with the center initialized in a
`+1` eigenstate of `Z` and all neighbors also in `+1` eigenstates of `Z`.

This helper is a small analytic regression target for gate signs and benchmark
plumbing; it is not a general triangular-lattice stabilizer solver.
"""
function stabilizer_expectation_exact(t::Real; initial::Symbol = :z_plus)
    if initial == :z_plus || initial == :plus
        return cos(2t)
    else
        throw(ArgumentError("supported initial states: :z_plus"))
    end
end

end
