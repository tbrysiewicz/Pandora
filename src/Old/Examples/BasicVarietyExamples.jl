export
    OrthogonalGroup,
    SpecialOrthogonalGroup

function OrthogonalGroup(n)
    @var x[1:n,1:n]
    M = x*x'- LinearAlgebra.I
    F = unique(vcat(collect(eachcol(M))...))
    Variety(System(F))
end

function SpecialOrthogonalGroup(n)
    @var x[1:n,1:n]
    M = x*x'- LinearAlgebra.I
    D = LinearAlgebra.det(x)-1
    F = unique(vcat(collect(eachcol(M))...))
    push!(F,D)
    Variety(System(F))
end
