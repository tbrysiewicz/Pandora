export 
    is_real,
    sign,
    real_solutions,
    nonreal_solutions,
    dietmaier,
    dietmaier_pair,
    false_function,
    zero_function,
    n_real_solutions,
    valid_fibre_count,
    valid_real_fibre,
    valid_fibre,
    real_solutions

###solution, solution set, and fibre functions

is_real(s::Vector{ComplexF64}; tol::Float64=1e-6) = maximum(abs ∘ imag, s) < tol
sign(s::Vector{Float64}) = sign.(s)
real_solutions(S::Vector{Vector{ComplexF64}}; tol::Float64=1e-6) = filter(x->is_real(x; tol), S)
nonreal_solutions(S::Vector{Vector{ComplexF64}}; tol::Float64=1e-6) = filter(x->!is_real(x; tol), S)

function dietmaier(S::Vector{Vector{ComplexF64}})
    nrs = nonreal_solutions(S)
    return isempty(nrs) ? 0.0 : minimum(norm ∘ imag, nonreal_solutions(S)) 
end

dietmaier_pair(S::Vector{Vector{ComplexF64}}) = (n_real_solutions(S), -dietmaier(S))
false_function(s) = false
zero_function(s) = 0.0



n_real_solutions(fibre::Fibre; tol::Float64=1e-6)::Int = n_real_solutions(solutions(fibre); tol=tol)
n_real_solutions(S::Vector{Vector{ComplexF64}}; tol::Float64=1e-6)::Int = count(s -> is_real(s; tol=tol), S)
real_solutions(fibre::Fibre) = real_solutions(solutions(fibre))


valid_fibre_count(EP::EnumerativeProblem, fibre::Fibre) = length(solutions(fibre)) == degree(EP)
valid_real_fibre(EP::EnumerativeProblem, fibre::Fibre) =  n_real_solutions(fibre) % 2 == degree(EP) % 2
valid_real_solution_set(EP::EnumerativeProblem, S::Vector{Vector{ComplexF64}}) =  n_real_solutions(S) % 2 == degree(EP) % 2


function valid_fibre(EP::EnumerativeProblem, fibre::Fibre)
    # Check if the number of solutions matches the degree of the problem
    if valid_fibre_count(EP, fibre) == false
        return(false)
    end
    if is_real(fibre[2]) && valid_real_fibre(EP, fibre) == false
        return(false)
    end
    return(true)
end
