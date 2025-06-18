export 
    branched_cover_decomposition,
    restrict

function restrict(G::PermGroup,O::Vector{Int64})
    S = sort(collect(O))
    relabel = Dict{Int,Int}()
    for i in 1:length(S)
        relabel[S[i]]=i
    end
    gg = gens(G)
    gens_on_O = [perm([relabel[g(o)] for o in O]) for g in gg]
    H,_ = sub(symmetric_group(length(S)),gens_on_O)
    return(H)
end


function branched_cover_decomposition(EP::EnumerativeProblem)
    G = galois_group(EP)
    if is_transitive(G)
        return([EP])
    end
    (sols,params) = base_fibre(EP)
    solution_partition = [sols[collect(o)] for o in orbits(G)]
    EP_List = [EnumerativeProblem(system(EP);populate=false) for sp in solution_partition]
    for i in eachindex(EP_List)
        know!(EP_List[i],BASE_FIBRE,(solution_partition[i],params))
        learn!(EP_List[i],DEGREE)
        o = orbits(G)[i]
        know!(EP_List[i],MONODROMY_GROUP,restrict(G,collect(o)))
    end
    return(EP_List)
end

#Given P = [p_1...p_k] parameters, this function considers the affine span of P 
#   as p_1+span(p_i-p_1) where span(p_i-p_1) = span(b_1...b_k-1) where the bi's are
#   an orthonormal basis for the span. 
#TODO: Which cache values can be inherited by the restriction?
function restrict(EP::EnumerativeProblem,P::Vector{Vector{Float64}})
	n = length(P)
    @var t[1:n-1]
    basis_vectors = gram_schmidt([(P[i]-P[1]) for i in 2:n])
    affine_span = P[1] + sum([t[i].*basis_vectors[i] for i in 1:n-1])
    new_expressions = [subs(f,parameters(EP)=>affine_span) for f in expressions(EP)]
    return(EnumerativeProblem(System(new_expressions,variables=variables(EP),parameters=t)))
end