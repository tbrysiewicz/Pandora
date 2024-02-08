using Pandora
using HomotopyContinuation
using Oscar
using Combinatorics 


function Triangles(sets)
    sol = []
    for s in sets
        s = sort(s)
        T = combinations(s,3)
        [push!(sol,t) for t in T]
    end
    unique(sol)
end

function HeronFormula(a,b,c,A)
    (4*a*b-(a+b-c)^2-16*A)
end

function ArLenSystem(sets)
    trgs = Triangles(sets)
    Eqs = []
    Pars = []
    for tr in trgs
        lenvars = [(@var s[i...])[1][1] for i in collect(combinations(tr,2))]
        arvar = (@var t[tr...])[1][1]
        lst = vcat(lenvars,arvar)
        eq = push!(Eqs,HeronFormula(lst...))
        vs = push!(Pars,arvar)
    end
    System(Eqs, parameters=[i for i in Pars])
end


TwoSimplices = [[1,2,3,4,5],[1,2,3,4,6]];

E = AreaLengthSystem(TwoSimplices)
d = degree(E)
#G = galois_group(E)
#B = minimal_block_reps(G)
#order(G)

