function convert(::Type{Vector{Vector{ComplexF64}}},result::Result)
    solutions(result)
end

function convert(::Type{Fibre},result::Tuple{Result,Vector{ComplexF64}})
    (solutions(result[1]),result[2])
end