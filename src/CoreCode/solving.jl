export 
    solve



#An EnumerativeSolver is an object ES which can be 
#  - called on a parameter ES(p) to return a Fibre
#  - called on a list of parameters ES(P) to return a Vector{Fibre}
#  - called on a Fibre and a parameter ES(F,p) to return a Fibre
#  - called on a Fibre and a list of parameters ES(F,P) to return a Vector{Fibre}
# This is more than just using HomotopyContinuation's solve functionality since the EnumerativeSolver 
#    can implement specialized algorithms which make solving easier, and it is abstract in the sense that
#    often solving is done with auxilliary variables, and if properly implemented, the EnumerativeProblem, 
#    although it uses the big system, will return the expected solution sizes. 
struct EnumerativeSolver
    solve_single :: Function
    solve_many   :: Function 
    solve_single_from_fibre   :: Function 
    solve_many_from_fibre   :: Function 
    function EnumerativeSolver(;solve_single = nothing, solve_many = nothing, 
                    solve_single_from_fibre = nothing, solve_many_from_fibre = nothing)
        if solve_single_from_fibre === nothing || solve_single === nothing
            error("Must be able to solve from Fibre to Fibre and solve over a single Fibre")
        else
            if solve_many === nothing
                function solve_many_from_single(P::Vector{Fibre})
                    return([solve_single(p) for p in P])
                end
                solve_many = solve_many_from_single
            end
            if solve_many_from_fibre === nothing
                function solve_many_fibre_from_single_fibre(P::Vector{Fibre})
                    return([solve_single_from_fibre(p) for p in P])
                end
                solve_many_from_fibre = solve_many_fibre_from_single_fibre
            end
        end
        return(EnumerativeSolver(solve_single,solve_many,solve_single_from_fibre,solve_many_from_fibre))
    end


    function EnumerativeSolver(S::System,BF::Fibre)
        function ssff(fibre::Fibre, p::Vector{T}  where T <: Number)
            sols = solutions(solve(S,solutions(fibre);start_parameters=parameters(fibre),target_parameters=p))
            return((sols,p))
        end
        function ss(p::Vector{T} where T)
            return(ssff(BF,p))
        end
        function smff(fibre::Fibre,P::Vector{Vector{T}} where T <: Number)
            sols = solutions(solve(S,solutions(fibre);start_parameters=parameters(fibre),target_parameters=P))
            return(map(x->(solutions(x[1]),x[2]),sols))
        end
        function sm(P::Vector{Vector{T}} where T <: Number)
            return(smff(BF,P))
        end
        return(EnumerativeSolver(ss,sm,ssff,smff))
    end
end

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
    if knows(EP,BASE_FIBRE)
        return(solve(EP,base_fibre(EP),p))
    else
        println(AlterWarning)
        populate!(EP)
        return(solve(EP,p))
    end
end

#many solve when no fibre is given
function solve(EP::EnumerativeProblem, P::Vector{Vector{T}} where T)
    if knows(EP,BASE_FIBRE)
        return(solve(EP,base_fibre(EP),P))
    else
        println(AlterWarning)
        populate!(EP)
        return(solve(EP,P))
    end
end

(EP::EnumerativeProblem)(fibre::Fibre, P::Vector{Vector{T}} where T <: Number) = solve(EP,fibre,P)
(EP::EnumerativeProblem)(fibre::Fibre, p::Vector{T} where T) = solve(EP,fibre,p)      
(EP::EnumerativeProblem)(P::Vector{Vector{T}} where T <:Number) = solve(EP,P)
(EP::EnumerativeProblem)(p::Vector{T} where T <:Number) = solve(EP,p)
(EP::EnumerativeProblem)() = solve(EP)

function solve(EP::EnumerativeProblem)
    return(solve(EP,randn(ComplexF64,n_parameters(EP))))
end

function solve(EP::EnumerativeProblem,N::Int)
    return(solve(EP,[randn(ComplexF64,n_parameters(EP)) for i in 1:N]))
end

