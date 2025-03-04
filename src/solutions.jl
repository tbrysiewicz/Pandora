is_real(s::Vector{ComplexF64}; tol::Float64=1e-6) = maximum(abs ∘ imag, s) < tol
sign(s::Vector{Float64}; tol::Float64=1e-6) = sign.(s)
real_solutions(S::Vector{Vector{ComplexF64}}) = filter(is_real, S)
nonreal_solutions(S::Vector{Vector{ComplexF64}}) = filter(!is_real, S)

function dietmaier(S::Vector{Vector{ComplexF64}})
    nrs = nonreal_solutions(S)
    return isempty(nrs) ? 0.0 : minimum(norm ∘ imag, nonreal_solutions(S)) 
end

dietmaier_pair(S::Vector{Vector{ComplexF64}}) = (n_real_solutions(S), -dietmaier(S))
false_function(s) = false
zero_function(s) = 0.0
