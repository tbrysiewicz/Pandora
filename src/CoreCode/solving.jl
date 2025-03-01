export 
    solve


#Solve over many parameters P of  F(EP) from fibre = (S1,P1) -> {(?,p)}_{p in P}
function solve(EP::EnumerativeProblem, fibre::Fibre, P::Vector{Vector{T}} where T <: Number)
    S = solve(system(EP),solutions(fibre); start_parameters= parameters(fibre), target_parameters = P, tracker_options = tracker_options(EP))
    return(map(x->solutions(x[1]),S))
end

#Solve over a single parameter p of  F(EP) from fibre = (S1,P1) -> (?,p)
function solve(EP::EnumerativeProblem,fibre::Fibre, p::Vector{T} where T)
    S = solve(system(EP),solutions(fibre); start_parameters= parameters(fibre), target_parameters = p, tracker_options = tracker_options(EP))
    return(solutions(S))
end

#solve when no fibre is given
#assume base-fibre is intended as start system
#  if base_fibre is unknown, provide warning that base fibre is being populated
#  (i.e. that EP is changing despite not calling solve!)
function solve(EP::EnumerativeProblem, p::Vector{T} where T)
    if knows(EP,base_fibre)
        return(solve(EP,base_fibre(EP),p))
    else
        println(AlterWarning)
        populate!(EP)
        return(solve(EP,p))
    end
end

#many solve when no fibre is given
function solve(EP::EnumerativeProblem, P::Vector{Vector{T}} where T)
    if knows(EP,base_fibre)
        return(solve(EP,base_fibre(EP),P))
    else
        println(AlterWarning)
        populate!(EP)
        return(solve(EP,P))
    end
end

(EP::EnumerativeProblem)(fibre::Fibre, P::Vector{Vector{T}} where T <: Number) = solve(EP,fibre,P)
(EP::EnumerativeProblem)(fibre::Fibre, p::Vector{T} where T) = solve(EP,fibre,p)    
(EP::EnumerativeProblem)(fibre1::Fibre, fibre2::Fibre) = solve(EP,fibre1,parameters(fibre2))    
(EP::EnumerativeProblem)(P::Vector{Vector{T}} where T <:Number) = solve(EP,P)
(EP::EnumerativeProblem)(p::Vector{T} where T <:Number) = solve(EP,p)
(EP::EnumerativeProblem)() = solve(EP)

function solve(EP::EnumerativeProblem)
    return(solve(EP,randn(ComplexF64,n_parameters(EP))))
end

function solve(EP::EnumerativeProblem,N::Int)
    return(solve(EP,[randn(ComplexF64,n_parameters(EP)) for i in 1:N]))
end

