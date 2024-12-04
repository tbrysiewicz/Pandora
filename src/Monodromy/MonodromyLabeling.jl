


#The fibre of an enumerative problem constisting of d points must be ordered implicitly via its ordering
#   as a vector of points in julia. The orbit-block labeling of a fibre is as follows with respect to the mon group G
#  1) First it is labeled via G-orbits of weakly decreasing size
#       e.g. If [2,4,8],[1,3],[5,6,7] are orbits, it is relabeled
#               [1,2,3],[7,8],[4,5,6]
#  2) Now, restricting to an orbit O=[1,...,k], we write G_O for the monodromy group of the branched cover from that
#       component of the incidence variety to the parameter space. 
#       If G_O is imprimitive, with blocks B_1,...,B_r, then we reorder 1...k, so that 1...k/r = B_1, k/r+1...2k/r = B2, etc
#       Note: there could be multiple block systems - we only choose one
#  3) Finally, within blocks, we order the points by modulus. This may not be a total order now (there could be another block
#       system which tracks this modulus)
function orbit_block_relabeling(EP::EnumerativeProblem)
    G = monodromy_group(EP)
    O = sort([sort(collect(o)) for o in orbits(G)])
    G_restrictions = [restrict(G,o) for o in O]
    for o in o
    end
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