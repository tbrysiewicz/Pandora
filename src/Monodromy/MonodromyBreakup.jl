function components(EP::EnumerativeProblem)
    G = galois_group(EP)
    if is_transitive(G)
        return([EP])
    end
    (sols,params) = base_fibre(EP)
    solution_partition = [sols[collect(o)] for o in orbits(G)]
    return([EnumerativeProblem(system(EP),(sp,params)) for sp in solution_partition])
end

