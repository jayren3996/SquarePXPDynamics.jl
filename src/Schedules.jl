module Schedules

export first_order_colors, second_order_colors

function first_order_colors()
    return collect(1:7)
end

function second_order_colors()
    colors = first_order_colors()
    return vcat(colors, reverse(colors))
end

end
