module Internals

# Project-internal helpers shared across multiple submodules. None of these
# are intended as public API; they exist so each consumer module can `using
# ..Internals: foo` instead of redefining the same primitive.

export _csv_value
export _DIRECTIONS, _validate_direction, _opposite_direction
export _finite_positive, _finite_nonnegative, _positive_int, _nonnegative_int
export _nonnegative_float, _optional_finite_float, _optional_finite_nonnegative
export _finite_max_or_nothing

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

# Parameter validators shared across submodules. Each takes the value and a
# label used in the error message; the previously-scattered copies in
# PXPValidation, CTMTrust, ScarFinder, and PXPValidationCampaigns are gone.

function _finite_positive(value::Real, label::String)
    converted = Float64(value)
    isfinite(converted) && converted > 0 ||
        throw(ArgumentError("$label must be finite and positive"))
    return converted
end

function _finite_nonnegative(value::Real, label::String)
    converted = Float64(value)
    isfinite(converted) && converted >= 0 ||
        throw(ArgumentError("$label must be finite and nonnegative"))
    return converted
end

function _positive_int(value::Integer, label::String)
    converted = Int(value)
    converted >= 1 || throw(ArgumentError("$label must be at least 1"))
    return converted
end

function _nonnegative_int(value::Integer, label::String)
    converted = Int(value)
    converted >= 0 || throw(ArgumentError("$label must be nonnegative"))
    return converted
end

# Allows non-finite sentinels (Inf) for parameters used as "no limit" caps.
function _nonnegative_float(value::Real, label::String)
    converted = Float64(value)
    converted >= 0 || throw(ArgumentError("$label must be nonnegative"))
    return converted
end

_optional_finite_float(::Nothing, label::String) = nothing
function _optional_finite_float(value::Real, label::String)
    converted = Float64(value)
    isfinite(converted) || throw(ArgumentError("$label must be finite"))
    return converted
end

_optional_finite_nonnegative(::Nothing, label::String) = nothing
_optional_finite_nonnegative(value::Real, label::String) =
    _finite_nonnegative(value, label)

# Reduces a value collection to the maximum, or `nothing` if empty.
_finite_max_or_nothing(values) = isempty(values) ? nothing : maximum(values)

end # module Internals
