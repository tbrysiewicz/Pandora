export
    solve

# An EnumerativeSolver is an object ES which can be
#  - called on a parameter ES(p) to return a Fibre
#  - called on a list of parameters ES(P) to return a Vector{Fibre}
#  - called on a Fibre and a parameter ES(F, p) to return a Fibre
#  - called on a Fibre and a list of parameters ES(F, P) to return a Vector{Fibre}
# This is more than just using HomotopyContinuation's solve functionality since the EnumerativeSolver
#    can implement specialized algorithms which make solving easier, and it is abstract in the sense that
#    often solving is done with auxiliary variables, and if properly implemented, the EnumerativeProblem,
#    although it uses the big system, will return the expected solution sizes.
struct EnumerativeSolver
    solve_single
    solve_many
    solve_single_from_fibre
    solve_many_from_fibre
end

function EnumerativeSolver(;
    solve_single = nothing,
    solve_many = nothing,
    solve_single_from_fibre = nothing,
    solve_many_from_fibre = nothing,
)
    if solve_single_from_fibre === nothing || solve_single === nothing
        error("Must be able to solve from Fibre to Fibre and solve over a single Fibre")
    else
        if solve_many === nothing
            function solve_many_from_single(P::Vector{Fibre})
                return [solve_single(p) for p in P]
            end
            solve_many = solve_many_from_single
        end
        if solve_many_from_fibre === nothing
            function solve_many_fibre_from_single_fibre(P::Vector{Fibre})
                return [solve_single_from_fibre(p) for p in P]
            end
            solve_many_from_fibre = solve_many_fibre_from_single_fibre
        end
    end
    return EnumerativeSolver(solve_single, solve_many, solve_single_from_fibre, solve_many_from_fibre)
end

function EnumerativeSolver(S::System, BF::Fibre)
    function ssff(fibre::Fibre, p::Vector{T} where T <: Number)
        sols = solutions(solve(S, solutions(fibre); start_parameters = parameters(fibre), target_parameters = p))
        return (sols, p)
    end

    function ss(p::Vector{T} where T)
        return ssff(BF, p)
    end

    function smff(fibre::Fibre, P::Vector{Vector{T}} where T <: Number)
        sols = solutions(solve(S, solutions(fibre); start_parameters = parameters(fibre), target_parameters = P))
        return map(x -> (solutions(x[1]), x[2]), sols)
    end

    function sm(P::Vector{Vector{T}} where T <: Number)
        return smff(BF, P)
    end

    return EnumerativeSolver(ss, sm, ssff, smff)
end

# Solve over many parameters P of F(EP) from fibre = (S1, P1) -> {(?, p)}_{p in P}, tracking only solution at index idx.
# Warning: this forgets the parameters when returned. So when using, be sure to save parameters.
function solve(EP::EnumerativeProblem, fibre::Fibre, P::Vector{Vector{T}}, idx::Int) where T <: Number
    S = solve(system(EP), [solutions(fibre)[idx]]; start_parameters = parameters(fibre), target_parameters = P)
    return map(x -> x[1][1].solution, S)
end

# Solve over a single parameter p of F(EP) from fibre = (S1, P1) -> (?, p), tracking only solution at index idx.
function solve(EP::EnumerativeProblem, fibre::Fibre, p::Vector{T}, idx::Int) where T
    S = solve(system(EP), [solutions(fibre)[idx]]; start_parameters = parameters(fibre), target_parameters = p)
    return S[1].solution
end

# Solve when no fibre is given, tracking only solution at index idx.
function solve(EP::EnumerativeProblem, p::Vector{T}, idx::Int) where T
    if knows(EP, BASE_FIBRE)
        return solve(EP, base_fibre(EP), p, idx)
    else
        populate!(EP)
        return solve(EP, p, idx)
    end
end

# Many solve when no fibre is given, tracking only solution at index idx.
function solve(EP::EnumerativeProblem, P::Vector{Vector{T}}, idx::Int) where T
    if knows(EP, BASE_FIBRE)
        return solve(EP, base_fibre(EP), P, idx)
    else
        populate!(EP)
        return solve(EP, P, idx)
    end
end

# Solve over many parameters P of F(EP) from fibre = (S1, P1) -> {(?, p)}_{p in P}.
function solve(EP::EnumerativeProblem, fibre::Fibre, P::Vector{Vector{T}}) where T
    S = solve(system(EP), solutions(fibre); start_parameters = parameters(fibre), target_parameters = P)
    S = map(x -> solutions(x[1]), S)
    return S
end

# Solve over a single parameter p of F(EP) from fibre = (S1, P1) -> (?, p).
function solve(EP::EnumerativeProblem, fibre::Fibre, p::Vector{T} where T)
    S = solve(system(EP), solutions(fibre); start_parameters = parameters(fibre), target_parameters = p)
    return solutions(S)
end

# Solve when no fibre is given.
# Assume base-fibre is intended as start system.
# If base_fibre is unknown, provide warning that base fibre is being populated
# (i.e. that EP is changing despite not calling solve!).
function solve(EP::EnumerativeProblem, p::Vector{T} where T)
    if knows(EP, BASE_FIBRE)
        return solve(EP, base_fibre(EP), p)
    else
        populate!(EP)
        return solve(EP, p)
    end
end

# Many solve when no fibre is given.
function solve(EP::EnumerativeProblem, P::Vector{Vector{T}} where T)
    if knows(EP, BASE_FIBRE)
        return solve(EP, base_fibre(EP), P)
    else
        populate!(EP)
        return solve(EP, P)
    end
end

(EP::EnumerativeProblem)(fibre::Fibre, P::Vector{Vector{T}} where T <: Number) = solve(EP, fibre, P)
(EP::EnumerativeProblem)(fibre::Fibre, p::Vector{T} where T) = solve(EP, fibre, p)
(EP::EnumerativeProblem)(P::Vector{Vector{T}} where T <: Number) = solve(EP, P)
(EP::EnumerativeProblem)(p::Vector{T} where T <: Number) = solve(EP, p)
(EP::EnumerativeProblem)() = solve(EP)

(EP::EnumerativeProblem)(fibre::Fibre, P::Vector{Vector{T}}, idx::Int) where T <: Number = solve(EP, fibre, P, idx)
(EP::EnumerativeProblem)(fibre::Fibre, p::Vector{T}, idx::Int) where T = solve(EP, fibre, p, idx)
(EP::EnumerativeProblem)(P::Vector{Vector{T}}, idx::Int) where T <: Number = solve(EP, P, idx)
(EP::EnumerativeProblem)(p::Vector{T}, idx::Int) where T <: Number = solve(EP, p, idx)

function solve(EP::EnumerativeProblem)
    return solve(EP, randn(ComplexF64, n_parameters(EP)))
end

function solve(EP::EnumerativeProblem, N::Int)
    return solve(EP, [randn(ComplexF64, n_parameters(EP)) for i in 1:N])
end
