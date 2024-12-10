function compute_components(EP::EnumerativeProblem)
    G = galois_group(EP)
    if is_transitive(G)
        return([EP])
    end
    (sols,params) = base_fibre(EP)
    solution_partition = [sols[collect(o)] for o in orbits(G)]
    EP_List = [EnumerativeProblem(system(EP),(sp,params)) for sp in solution_partition]
    for i in 1:eachindex(EP_List)
        o = orbits(G)[i]
        E = EP_List[i]
        data(E)[:monodromy_group] = restrict(G,collect(o))
    end
    return(EP_List)
end

