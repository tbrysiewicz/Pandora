
function generates_transitive_group(gg::Vector{PermGroupElem})
    H = subgroup(gg)
    is_transitive(H)
end


function transitivity_basis(G::PermGroup)
    id = perm(1:degree(G))
    candidate_basis = Vector{PermGroupElem}()
    if is_transitive(G)==false
        return("Error: given permutation group is not transitive")
    end
    while length(candidate_basis)==0 || generates_transitive_group(candidate_basis)==false
        newgen = rand(G)
        if newgen!=id
            push!(candidate_basis,newgen)
        end
    end
    return(candidate_basis)
end

function monodromy_basis(EP::EnumerativeProblem)
    G = monodromy_group(EP)
    N = order(G)
    id = perm(1:degree(G))
    gg = filter(x->x!=id,collect(keys(monodromy_dictionary(EP))))
    S = [1]
    current_group = subgroup(gg[S])
    for j in 2:length(gg)
        if in(gg[j],current_group)==false
            push!(S,j)
            current_group = subgroup(gg[S])
        end
    end
    return(gg[S])
end

function orbit_block_relabel!(EP::EnumerativeProblem)
    OBP = orbit_block_partition(EP)
    one_line_obp = vcat(vcat(OBP...)...)
    p = perm(one_line_obp)
    update_base_fibre!(EP,(base_solutions(EP)[one_line_obp],base_parameters(EP)))
    println("TODO: automatically recompute monodromy - for now you have to do it yourself with force_recompute")
end


#Group the solutions of an enumerative problem (identified implicitly as [1...d]) first by their orbits under monodromy
#  and then on the next level by their membership in the blocks of the largest block system
function orbit_block_partition(EP::EnumerativeProblem)
    G = monodromy_group(EP)
    O = sort([sort(collect(o)) for o in orbits(G)])
    orbit_block = []
    for o in O
        H = restrict(G,o) #This is the restriction of G to the orbit O (in particular, this group is transitive)
        mbrs = minimal_block_reps(H)
        max_block_size = max(map(x->length(x),mbrs)...)
        mbr = mbrs[findfirst(x->length(x)==max_block_size,mbrs)]
        block_system = collect(orbit(H,on_sets,mbr))
        push!(orbit_block,[sort(o[bs]) for bs in block_system])
    end
    orbit_block = sort(orbit_block, lt = (x,y)->length(x)<length(y))
    return(orbit_block)
end

#restrict a permutation group G<=S_n to a union S of orbits
#   of G; the output remains a subgroup of S_n, one that fixes
#   all elements outside of S
function restrict(G::PermGroup, S::Vector{Int64})
    gg = gens(G)
    d = length(S)
    Sd = symmetric_group(d)
    newgens = Vector{PermGroupElem}([])
    for g in gg
        c = cycles(g)
        fa = findall(x->issubset(x,S),c)
        relabeled_cycles = [[findfirst(x->x==a,S) for a in cyc] for cyc in c[fa]]
        newperm = cperm(Sd,relabeled_cycles)
        push!(newgens,newperm)
    end
    return(subgroup(newgens))
end