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


function compute_components(EP::EnumerativeProblem)
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

