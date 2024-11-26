const Fibre = Tuple{Vector{Vector{ComplexF64}},Vector{ComplexF64}}

function solutions(fibre::Fibre) :: Vector{Vector{ComplexF64}}
    return(fibre[1])
end

function parameters(fibre::Fibre) :: Vector{ComplexF64}
    return(fibre[2])
end
