module Internals

# Project-internal helpers shared across multiple submodules. None of these
# are intended as public API; they exist so each consumer module can `using
# ..Internals: foo` instead of redefining the same primitive.

export _csv_value
export _DIRECTIONS, _validate_direction, _opposite_direction

const _DIRECTIONS = (:right, :up, :left, :down)

"""
    _validate_direction(dir::Symbol) -> dir

Throw `ArgumentError` unless `dir` is one of `:right`, `:up`, `:left`, `:down`.
Returns `dir` on success.
"""
function _validate_direction(dir::Symbol)
    dir in _DIRECTIONS || throw(
        ArgumentError("direction must be :right, :up, :left, or :down (got :$dir)"),
    )
    return dir
end

"""
    _opposite_direction(dir::Symbol) -> Symbol

Return the four-way opposite of `dir` on the square lattice.
"""
function _opposite_direction(dir::Symbol)
    dir === :right && return :left
    dir === :up && return :down
    dir === :left && return :right
    dir === :down && return :up
    throw(ArgumentError("direction must be :right, :up, :left, or :down (got :$dir)"))
end

# CSV value formatters. Each consumer module used to redefine a similar
# fan-out by type; consolidating here lets the writers stay one-liners.

_csv_value(::Nothing) = ""
_csv_value(value::Bool) = string(value)
_csv_value(value::Real) = string(value)
_csv_value(value::Symbol) = String(value)
function _csv_value(value::AbstractString)
    escaped = replace(value, "\"" => "\"\"")
    return any(ch -> ch in escaped, (',', '"', '\n', '\r')) ? "\"$escaped\"" : escaped
end

end # module Internals
