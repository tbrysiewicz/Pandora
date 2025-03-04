
n_real_solutions(fibre::Fibre; tol::Float64=1e-6)::Int = n_real_solutions(solutions(fibre); tol=tol)
n_real_solutions(S::Vector{Vector{ComplexF64}}; tol::Float64=1e-6)::Int = count(s -> is_real(s; tol=tol), S)
real_solutions(fibre::Fibre) = real_solutions(solutions(fibre))
valid_fibre(EP::EnumerativeProblem, fibre::Fibre) = length(solutions(fibre)) == degree(EP)
valid_real_fibre(EP::EnumerativeProblem, fibre::Fibre) = length(real_solutions(fibre)) % 2 == degree(EP) % 2
