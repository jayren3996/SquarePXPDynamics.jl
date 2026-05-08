module Evolution

using ..Geometry: Coord
using ..States: TriangularIPEPS
using ..Schedules: first_order_colors, second_order_colors
using ..SimpleUpdate: apply_star_gate_simple_update!

export evolve_step!, color_canonical_center

"""
    color_canonical_center(color) -> Coord

Return a canonical center coordinate whose `star_color` equals `color`.
Used by `evolve_step!` to pick a representative center per color.
"""
function color_canonical_center(color::Integer)
    1 <= color <= 7 || throw(ArgumentError("color must be in 1:7"))
    # star_color(c) = (q + 3r) mod 7 + 1  →  pick (color-1, 0).
    return Coord(color - 1, 0)
end

"""
    evolve_step!(state, gate; order=:second, update=:simple) -> state

Apply one Trotter "step" to `state` by sweeping the star gate over all 7
colors of the triangular partition. `order` selects the schedule
(`:first` or `:second`); `update` selects the bond-update backend
(`:simple` is the only currently supported value).

The same `gate` is applied at each color; for translationally invariant
states the schedule reduces to repeated applications at color-canonical
centers.
"""
function evolve_step!(state::TriangularIPEPS, gate::AbstractMatrix;
                      order::Symbol = :second, update::Symbol = :simple)
    schedule = if order === :first
        first_order_colors()
    elseif order === :second
        second_order_colors()
    else
        throw(ArgumentError("order must be :first or :second"))
    end

    update === :simple || throw(ArgumentError("update must be :simple"))

    for color in schedule
        center = color_canonical_center(color)
        apply_star_gate_simple_update!(state, gate, center)
    end
    return state
end

end
