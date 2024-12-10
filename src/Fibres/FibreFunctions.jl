
function n_real_solutions(fibre::Fibre; tol::Float64=1e-6) :: Int
    n_real_solutions(solutions(fibre);tol=tol)
end
    
function n_real_solutions(S::Vector{Vector{ComplexF64}}; tol::Float64=1e-6) :: Int
    count(s->is_real(s;tol=tol), S)
end
    