export 
    MonodromyLoop,
    Path,
    monodromy

mutable struct MonodromyLoop
    F :: Fibre                      # Gives sols and ordering on them. (S,P)
    P :: Vector{Vector{ComplexF64}} # A list of parameters starting
    sigma :: Union{Perm,Nothing}
end


function MonodromyLoop(EP::EnumerativeProblem,P::Vector{Vector{ComplexF64}}) 
    F = base_fibre(EP)
    MonodromyLoop(F,P,nothing)
end


#p is interpretted as a function where p(i) = p[i]. But p(i) may be not injective or
#  surjective, in which case the following function should return false.
function is_valid_permutation(p::Union{Vector{Int},Vector{Union{Nothing,Int64}}},d::Int64) :: Bool
    u=unique(p) 
    return(!(length(u)<d || in(nothing,u))) #nothing indicates not surjective, length indicates not injective
end

function numerical_bijection(S1::Vector{T} where T, S2::Vector{T} where T; tol::Float64 = 1e-6) 
    one_line = [findfirst(x->isapprox(x,s),S1) for s in S2]
    return(one_line)
end

function check_and_fix_monodromy_input(fibre::Fibre, loop::Vector{Vector{T}})  where T
    if T != ComplexF64 #Make sure parameters are complex
        loop=Vector{Vector{ComplexF64}}(loop) #Parameter type: changing parameters in loop to complex floats
    end

    if isapprox(parameters(fibre),first(loop))==false #If input loop is doesn't start at fibre force it to
        pushfirst!(loop,parameters(fibre)) #Inconsistent loop: starting at given fibre parameter
    end

    if isapprox(first(loop),last(loop))==false #If input loop doesn't end at fibre, force it to
        push!(loop, first(loop)) #Incompelte loop: using start parameter as end parameter
    end
    return(loop)
end

function monodromy(EP::EnumerativeProblem, fibre::Fibre, loop::Vector{Vector{T}} where T)
    loop = check_and_fix_monodromy_input(fibre,loop)

    (S,P) = fibre
    for l in loop[2:end]
        (S,P) = (EP((S,P),l),l)
    end
    
    g = numerical_bijection(solutions(fibre),S)
    return is_valid_permutation(g,degree(EP)) ? perm(g) : g
end

monodromy(EP::EnumerativeProblem,loop::Vector{Vector{T}} where T) =monodromy(EP,base_fibre(EP),loop)
monodromy(EP::EnumerativeProblem,ML::MonodromyLoop) = monodromy(EP,ML.F,ML.P)
function monodromy!(EP::EnumerativeProblem,ML::MonodromyLoop)
    m = monodromy(EP,ML.F,ML.P)
    if typeof(m) <: Perm
        ML.sigma = m
    end
    return(m)
end

perm!(EP::EnumerativeProblem, ML::MonodromyLoop) = monodromy!(EP,ML)
