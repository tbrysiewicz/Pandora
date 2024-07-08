export
    solve_over_param,
    solve_over_params



#####Basic moving in parameter space

function solve_over_param(E::EnumerativeProblem,P; monodromy_recover=false)
    #Consider implementing special solvers when
    #  there is decomposability
    if is_populated(E)==false
        populate_base_fibre(E)
    end
    S = HomotopyContinuation.solve(system(E),base_solutions(E); start_parameters= base_parameters(E), target_parameters = P)
    if monodromy_recover==true && degree_check(E,S)==false
        println("Lost points during tracking...recovering via monodromy")
        M = monodromy_solve(system(E),solutions(S),P)
        S=solutions(M)
    end
    return S
end


function solve_over_params(E::EnumerativeProblem,P; start_fibre = nothing, monodromy_recover=false, checks=[degree_check, real_parity_check],retry=false, verbose=false)
    #Consider implementing special solvers when
    #  there is decomposability
    if is_populated(E)==false
        populate_base_fibre(E)
    end
    if start_fibre == nothing
        start_fibre = base_fibre(E)
    end
    S = HomotopyContinuation.solve(system(E),start_fibre[1]; start_parameters= start_fibre[2], target_parameters = P)
    verbose && println("Total number of fibres computed:",length(S))
    for f in checks
        B=filter(s->f(E,solutions(s[1]))==false,S)
        S = setdiff(S,B)

        #S=filter!(s->f(E,solutions(s[1])),S)
        verbose && println(string(length(S)),"/",length(P)," satisfies ",f)
        if retry==true && length(B)>0
            verbose && println("Retrying some failed runs")
            newS = solve_over_params(E,[b[2] for b in B]; checks=checks, retry=false)
            S=vcat(S,newS)
        end
    end
    if monodromy_recover==true && degree_check(E,S)==false
        println("Lost points during tracking...recovering via monodromy")
        M = monodromy_solve(E.F,solutions(S[1]),P)
        S=solutions(M)
    end
    return S
end
