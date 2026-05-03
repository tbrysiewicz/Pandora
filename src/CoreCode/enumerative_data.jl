export
    EnumerativeData,
    EnumerativeAttribute,
    enumerative_data,
    known_data,
    is_enumerative_property,
    is_enumerative_data

"""
An `EnumerativeData` object names computational data attached to an
enumerative problem.

Unlike an `EnumerativeProperty`, enumerative data is not expected to have a
unique mathematically correct value for a fixed problem. Repeated computations
may legitimately produce different values, such as different sampled monodromy
loops or cached fibres.
"""
struct EnumerativeData{T}
    name::String
end

const EnumerativeAttribute = Union{EnumerativeProperty, EnumerativeData}

get_type(::EnumerativeData{T}) where {T} = T
name(EData::EnumerativeData) = EData.name
Base.show(io::IO, EData::EnumerativeData) = print(io, name(EData))

is_enumerative_property(::EnumerativeProperty) = true
is_enumerative_property(::EnumerativeData) = false
is_enumerative_data(::EnumerativeProperty) = false
is_enumerative_data(::EnumerativeData) = true

