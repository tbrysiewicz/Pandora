
#p is interpretted as a function where p(i) = p[i]. But p(i) may be not injective or
#  surjective, in which case the following function should return false.
function is_valid_permutation(p::Union{Vector{Int},Vector{Union{Nothing,Int64}}},d::Int64)
    if in(nothing,p)
        return(false)
    end
    u=unique(p) #
    if length(u)<d || in(nothing,u) #nothing indicates not surjective, length indicates not injective
        return(false)
    else
        return(true)
    end
end

function numerical_bijection(S1::Vector{T} where T, S2::Vector{T} where T) 
    one_line = [findfirst(x->isapprox(x,s),S1) for s in S2]
    return(one_line)
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
function monodromy_sample(EP::EnumerativeProblem, N::Int; loop_scaling = 1.0, permutations_only = false)
    k = n_parameters(EP)
    c = loop_scaling.*randn(ComplexF64,k)
    b = [loop_scaling.*randn(ComplexF64,k) for i in 1:N]
    a = base_parameters(EP)
    S1_bij = solve(EP, base_fibre(EP), b)
    S2 = solve(EP,c)
    S2_fibre = (S2,c)
    S2_bij = solve(EP,S2_fibre,b)
    one_line_perms = [numerical_bijection(S2_bij[i],S1_bij[i]) for i in 1:N]
    indices_of_valid_permutations = findall(x->is_valid_permutation(x,degree(EP)),one_line_perms)
    println("# Loops computed:            ",N)
    println("# Valid permutations:        ",length(indices_of_valid_permutations))
    println("# Unique valid permutations: ",length(unique(one_line_perms[indices_of_valid_permutations])))
    if permutations_only
        return([perm(p) for p in one_line_perms[indices_of_valid_permutations]])
    else
        sampled_loops = [[a,bb,c,a] for bb in b]
        return(([perm(p) for p in one_line_perms[indices_of_valid_permutations]],sampled_loops[indices_of_valid_permutations]))
    end
end

function compute_monodromy_dictionary(EP::EnumerativeProblem; n_loops = 50) :: Dict{PermGroupElem,Vector{Vector{ComplexF64}}}
    monodromy_dictionary = Dict{PermGroupElem,Vector{Vector{ComplexF64}}}()
    (perms,loops) = monodromy_sample(EP,n_loops)
    U = unique(perms)
    unique_perm_locations = [findfirst(x->x==u,perms) for u in U]
    for upl in unique_perm_locations
        monodromy_dictionary[perms[upl]]=loops[upl]
    end
    return(monodromy_dictionary)

end

function compute_monodromy_group(EP::EnumerativeProblem)
    G = subgroup(unique(keys(monodromy_dictionary(EP))))
end




