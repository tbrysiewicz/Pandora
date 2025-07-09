

function sparse_system(Abullet::Vector{Matrix{Int}})
    n = length(Abullet)
    @assert all(A->size(A,1) == n, Abullet) # All supports must be for same number of variables and square
    m = maximum([size(A,2) for A in Abullet]) # Maximum number of monomials in any support
    @var a[1:n,1:m]
    @var x[1:n]
    equations = []
    for k in 1:n
        A = Abullet[k]
        eq = 0
        for i in 1:size(A,2)
            monomial = a[k,i] * prod([x[j]^A[j,i] for j in 1:n])
            eq += monomial
        end
        push!(equations, eq)
    end
    return EnumerativeProblem(System(equations, variables = vec(x), parameters = vec(a));torus_only = true)
end
