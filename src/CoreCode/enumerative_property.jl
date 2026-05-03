export
    EnumerativeProperty

"""
An `EnumerativeProperty` names an intrinsic property of an enumerative problem.

For a fixed enumerative problem, an enumerative property is expected to have one
mathematically correct value.
"""
struct EnumerativeProperty{T}
    name::String
end

get_type(::EnumerativeProperty{T}) where {T} = T
name(EProp::EnumerativeProperty) = EProp.name
Base.show(io::IO, EProp::EnumerativeProperty) = print(io, name(EProp))

const DEGREE = EnumerativeProperty{Int}("degree")
const SYSTEM = EnumerativeProperty{System}("system")
const INEQUATIONS = EnumerativeProperty{System}("inequations")
const BASE_FIBRE = EnumerativeProperty{Fibre}("base_fibre")
