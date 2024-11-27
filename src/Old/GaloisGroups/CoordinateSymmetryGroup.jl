export
    coordinate_symmetry_group,
    coordinate_symmetry_heatmap

function heatmap_matrix_cleanup(M)
    (n,m) = size(M)
    HM = []
    for i in 1:n
        r = M[i,:]
        if length(unique(r))==1
            push!(HM,[0 for j in 1:m])
        elseif length(unique(r)) == length(r)
            push!(HM,[-1 for j in 1:m])
        else
            push!(HM,r)
        end
    end
    HMU = []
    for h in HM
        U = sort(unique(vcat(h,[0,-1])))
        if h[1]==0 || h[1]==-1
            push!(HMU,h)
        else
            push!(HMU,[findfirst(x->x==c,U)-2 for c in h])
        end
    end
    H = transpose(hcat(HMU...))
    EC = collect(eachcol(H))
    EC = sort(EC)
    H = hcat(EC...)
    return(H)
end

function coordinate_symmetry_heatmap(E::EnumerativeProblem; F=nothing)
    HM = coordinate_symmetry_group(E;F=F,heatmatrix=true)
    MyColors = [:white,:black,:blue,:red,:orange,:purple]
    m = max(HM...)
    if m<0
        m=0
    end
    P = heatmap(HM;xaxis = :false, yaxis=:false, legend = :false, xticks = nothing, yticks = nothing, aspect_ratio = :equal, color = MyColors[1:m+2], yflip=:true)
end

function coordinate_symmetry_group(E::EnumerativeProblem; F=nothing, heatmatrix = false)
    S = base_solutions(E)
    if F!=nothing
        S = [F(s) for s in S]
    end
    CoordValues = vcat(S...) #First extract the `unique` coordinate values (approximately)
    SC = [[findfirst(x->isapprox(x,c),CoordValues) for c in s] for s in S] #Then relabel coordinates by indices
    M = hcat(SC...)
    if heatmatrix == true
        return(heatmap_matrix_cleanup(M))
    end
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