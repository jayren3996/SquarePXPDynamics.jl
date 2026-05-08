module SolvableModels

export stabilizer_expectation_exact

"""
    stabilizer_expectation_exact(t::Real; initial::Symbol = :plus)

Exact single-star cluster evolution expectation for a center-site `|+>` state
with neighbors in `+1` Z eigenstates. For `initial = :plus`, this is `cos(2t)`.
"""
function stabilizer_expectation_exact(t::Real; initial::Symbol = :plus)
    if initial == :plus
        return cos(2t)
    else
        throw(ArgumentError("supported initial states: :plus"))
    end
end

end
