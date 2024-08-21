using HomotopyContinuation

export
    monodromy_element,
    sample_monodromy_elements


#This function takes a system, solutions, and a sequence of parameters
#  and tracks the solutions along those parameters in sequence
function track_sequence(F::System,Sols::Vector{Vector{ComplexF64}},γ::Vector{Vector{ComplexF64}})
    S=Sols;
    P=γ[1]
    for NewP in γ[2:end]
        newS=solve(F,S; start_parameters=P, target_parameters=NewP)
        S=newS
        P=NewP
    end
    return S
end


function monodromy_element(E::EnumerativeProblem,γ::Vector{Vector{ComplexF64}})
    F = system(E)
    BP = γ[1] #This is the base parameters of the loop.
    if is_populated(E)==false
        populate_base_fibre(E)
    end
    if γ[1] != γ[end]
        println("Parameters do not describe a loop. Appending first parameter to last")
        push!(γ,γ[1])
    end
    key_loop = [base_parameters(E)]
    for v in γ
        push!(key_loop,v)
    end
    push!(key_loop,base_parameters(E))
    M = solutions(track_sequence(F,base_solutions(E),key_loop)) #Take the solutions after the loop
    ggamma= [findfirst(x->x≈base_solutions(E)[i],M) for i in 1:length(M)] #and see which solutions were permuted
    if nothing ∈ ggamma || length(unique(ggamma))!=length(ggamma)
        return(nothing)
    else
        return Oscar.perm(ggamma)
    end
end

function monodromy_element(E::EnumerativeProblem; radius = 1)
    if is_populated(E)==false
        populate_base_fibre(E)
    end
    γ = [base_parameters(E)] #starting at BaseParam
    for i in 1:3
        push!(γ,last(γ)+radius*randn(ComplexF64,n_parameters(E))) #with 3 other parameter values
    end
    push!(γ,base_parameters(E)) #and ending again at BaseParam
    return monodromy_element(E,γ)
end

function sample_monodromy_elements(E::EnumerativeProblem,n::Int64; radius = 1)
    println("Sampling ",n," elements from the monodromy group of ",E)
    S = [monodromy_element(E; radius=radius) for i in 1:n]
    S = Vector{PermGroupElem}(filter!(x->x!=nothing,S))
    println(length(S)," out of ", n, " are plausible permutations")
    return(S)
end
