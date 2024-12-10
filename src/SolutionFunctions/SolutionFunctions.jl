function is_real(s::Vector{ComplexF64}; tol::Float64 = 1e-6)
    maximum(abs ∘ imag,s)<tol
end

function sign(s::Vector{Float64}; tol::Float64 = 1e-6)
    map(sign,s)
end