


function EnumerativeProblem(X::Variety)
    #Insert check to see if the variety is actually set up (witness set computed)
        F = system(X)
        V = variables(F)
        N = length(V)
        D = dim(X)
        @var o[1:D,1:1+N]
        Equations = Vector{Expression}(expressions(F))
        for i in 1:D
            l = V'*o[i,1:end-1]-o[i,end]
            push!(Equations,l)
        end
        W = witness_set(X)
        LS = linear_subspace(W)
        A = (extrinsic(LS)).A
        b = (extrinsic(LS)).b
        Ab = hcat(A,b)
        P = vcat(Ab...)
        NewSystem = System(Equations,variables = V, parameters = vcat(o...))
        S = solutions(W)
        R = solve(NewSystem,S;start_parameters=P,target_parameters=P)
        E = EnumerativeProblem(NewSystem,(R,P))
    end
    