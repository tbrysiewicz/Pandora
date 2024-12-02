
function numerical_bijection(S1::Vector{T}, S2::Vector{T}) where T
    one_line = [findfirst(x->isapprox(x,s),S1) for s in S2]
    return(perm(one_line))
end

function monodromy_homomorphism(EP::EnumerativeProblem, loop::Vector{Vector{T}} where T)
    println("Base fibre: assuming loop begins and ends at base fibre")
    monodromy_homomorphism(EP,base_fibre(EP),loop)
end

function monodromy_homomorphism(EP::EnumerativeProblem, fibre::Fibre, loop::Vector{Vector{T}} where T)
    #Type forgiveness
    if typeof(loop) != Vector{Vector{ComplexF64}}
        loop = Vector{Vector{ComplexF64}}(loop)
        println("Parameter type: changing parameters in loop to complex floats")
    end

    #Loop agrees with fibre forgiveness
    if isapprox(parameters(fibre),first(loop))==false
        pushfirst!(loop,parameters(fibre))
        println("Inconsistent loop: starting at given fibre parameter")
    end

    #Incomplete loop forgiveness
    if isapprox(first(loop),last(loop))==false
        push!(loop,first(loop))
        println("Incompelte loop: using start parameter as end parameter")
    end

    (S,P) = fibre
    for l in loop[2:end]
        S = solve(EP,(S,P),l)
        P = l
    end
    
    g = numerical_bijection(solutions(fibre),S)
end


#When sampling many monodromy elements, one can make use of the parallelization in HC.jl by choosing loops all of the form
#   a->b[i]->c->a where a is the base parameter, and c is another fixed parameter.
#   This requires 2*N+1 fibres to be tracked. 1 from a->c, and N from a->b[i] and c->b[i]
function monodromy_sample(EP::EnumerativeProblem, N::Int; loop_scaling = 1.0)
    k = n_parameters(EP)
    c = loop_scaling.*randn(ComplexF64,k)
    b = [loop_scaling.*randn(ComplexF64,k) for i in 1:N]
    S1_bij = solve(EP, base_fibre(EP), b)
    S2 = solve(EP,c)
    S2_fibre = (S2,c)
    S2_bij = solve(EP,S2_fibre,b)
    perms = [numerical_bijection(S2_bij[i],S1_bij[i]) for i in 1:N]
end


function monodromy_group(EP::EnumerativeProblem)
    perms = monodromy_sample(EP,50)
    G = subgroup(unique(perms))
end

function galois_group(EP::EnumerativeProblem)
    monodromy_group(EP)
end