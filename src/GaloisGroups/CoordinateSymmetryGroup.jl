
function coordinate_symmetry_group(E::EnumerativeProblem; F=nothing)
    S = base_solutions(E)
    if F!=nothing
        S = [F(s) for s in S]
    end
    CoordValues = vcat(S...) #First extract the `unique` coordinate values (approximately)
    SC = [[findfirst(x->isapprox(x,c),CoordValues) for c in s] for s in S] #Then relabel coordinates by indices
    M = hcat(SC...)
    #Now we have a matrix of small integers and we need to find the nontrivial rows (corresponding to nontrivial blocks)
    CoordinateBlocks = filter(x->length(unique(x))!=length(x) && length(unique(x))!=1, eachrow(hcat(SC...)))
    if length(CoordinateBlocks)==0
        return(symmetric_group(degree(E)))
    end
    #For each of these nontrivial rows, write the incidence matrix describing the solutions w/ same coords
    CoordinateIncidences = [unique([findall(x->in(i,x),r) for i in r]) for r in CoordinateBlocks]
    IMs = [IncidenceMatrix(CI) for CI in CoordinateIncidences]

    #The automorphism groups of these incidences are the permutation groups which preserve that block
    #  so the coordinate symmetry group is the intersection of them
    AGlist = [automorphism_group(IM)[:on_cols] for IM in IMs]

    println(AGlist)
    AG,_ = intersect(AGlist...)
    return(AG)
end