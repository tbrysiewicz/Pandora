function solve(EP::EnumerativeProblem,fibre::Fibre, P::Vector{Vector{T}} where T <: Number)
    S = solve(system(EP),solutions(fibre); start_parameters= parameters(fibre), target_parameters = P)
    return(map(x->solutions(x[1]),S))
end

function solve(EP::EnumerativeProblem,fibre::Fibre, p::Vector{T} where T)
    S = solve(system(EP),solutions(fibre); start_parameters= parameters(fibre), target_parameters = p)
    return(solutions(S))
end

function solve(EP::EnumerativeProblem, p::Vector{T} where T)
    return(solve(EP,base_fibre(EP),p))
end

function solve(EP::EnumerativeProblem, P::Vector{Vector{T}} where T)
    return(solve(EP,base_fibre(EP),P))
end

function (EP::EnumerativeProblem)(fibre::Fibre, P::Vector{Vector{T}} where T <: Number)
    return(solve(EP,fibre,P))
end

function (EP::EnumerativeProblem)(fibre::Fibre, p::Vector{T} where T)
    return(solve(EP,fibre,p))    
end

function (EP::EnumerativeProblem)(P::Vector{Vector{T}} where T <:Number)   
    return(solve(EP,P))
end

function (EP::EnumerativeProblem)(p::Vector{T} where T <:Number)    
    return(solve(EP,p))
end



function (EP::EnumerativeProblem)()
    solve(EP)
end

function solve(EP::EnumerativeProblem)
    p = randn(ComplexF64,n_parameters(EP))
    if is_populated(EP)==false
        S = solve(system(EP);target_parameters=p)
        return(solutions(S))
    else
        return(solve(EP,p))
    end
end



#=
function rebase!(EP::EnumerativeProblem; new_param)
    if !is_populated(EP)
        return(populate!(EP))
    end
    if !isdefined(new_param)
        new_param = randn(ComplexF64,n_parameters(EP))
    end
    S = EP(new_param)
    if length(S) == degree(EP)
        EP.base_fibre = (S,new_param)
        return(nothing)
    else
        println("Failed to rebase")
        return(nothing)
    end
end

#####Basic moving in parameter space

function solve_over_param(E::EnumerativeProblem,P; monodromy_recover=false, show_progress=false)
    #Consider implementing special solvers when
    #  there is decomposability
    if is_populated(E)==false
        populate_base_fibre(E)
    end
    S = HomotopyContinuation.solve(system(E),base_solutions(E); start_parameters= base_parameters(E), target_parameters = P, show_progress=show_progress)
    if monodromy_recover==true && degree_check(E,S)==false
        println("Lost points during tracking...recovering via monodromy")
        M = monodromy_solve(system(E),solutions(S),P)
        S=solutions(M)
    end
    return S
end


function solve_over_params(E::EnumerativeProblem,P; start_fibre = nothing, monodromy_recover=false, checks=[degree_check],retry=false, show_progress=false, verbose=false)
    if is_populated(E)==false
        populate_base_fibre(E)
    end
    if start_fibre == nothing
        start_fibre = base_fibre(E)
    end
    S = HomotopyContinuation.solve(system(E),base_solutions(E); start_parameters= base_parameters(E), target_parameters = P,show_progress=show_progress)
    verbose && println("Total number of fibres computed:",length(S))
    for f in checks
        B=filter(s->f(E,solutions(s[1]))==false,S)
        S = setdiff(S,B)

        #S=filter!(s->f(E,solutions(s[1])),S)
        verbose && println(string(length(S)),"/",length(P)," satisfies ",f)
        verbose && println("Retry:",retry)
        if retry==true && length(B)>0
            verbose && println("Retrying some failed runs")
            rebase!(E)
            newS = solve_over_params(E,[b[2] for b in B]; checks=checks, retry=false, verbose=verbose)
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

=#
