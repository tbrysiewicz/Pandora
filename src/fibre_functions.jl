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
valid_fibre(EP::EnumerativeProblem, fibre::Fibre) = length(solutions(fibre)) == degree(EP)
valid_real_fibre(EP::EnumerativeProblem, fibre::Fibre) = length(real_solutions(fibre)) % 2 == degree(EP) % 2
