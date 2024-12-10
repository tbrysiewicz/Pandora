function is_real(s::Vector{ComplexF64}; tol::Float64 = 1e-6)
    maximum(abs ∘ imag,s)<tol
end

function sign(s::Vector{Float64}; tol::Float64 = 1e-6)
    map(sign,s)
end

function real_solutions(S::Vector{Vector{ComplexF64}})
    filter(x->is_real(x),S)
end

function nonreal_solutions(S::Vector{Vector{ComplexF64}})
    filter(x->is_real(x)==false,S)
end

function dietmaier(S::Vector{Vector{ComplexF64}})
    minimum(norm ∘ imag, nonreal_solutions(S)) #is norm best here? or some 
end

function dietmaier_pair(S::Vector{Vector{ComplexF64}})
    (n_real_solutions(S),-minimum(norm ∘ imag, nonreal_solutions(S))) #is norm best here? or some 
end