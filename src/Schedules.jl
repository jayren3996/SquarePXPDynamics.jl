module Schedules

export first_order_colors, second_order_colors, schedule_layers

function first_order_colors()
    return collect(1:7)
end

function second_order_colors()
    colors = first_order_colors()
    return vcat(colors, reverse(colors))
end

function schedule_layers(order::Symbol)
    if order === :first
        return [(color = c, scale = 1.0) for c in first_order_colors()]
    elseif order === :second
        return vcat(
            [(color = c, scale = 0.5) for c in first_order_colors()],
            [(color = c, scale = 0.5) for c in reverse(first_order_colors())],
        )
    else
        throw(ArgumentError("order must be :first or :second"))
    end
end

end
