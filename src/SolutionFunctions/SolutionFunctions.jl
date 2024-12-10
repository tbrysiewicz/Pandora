function is_real(s::Vector{ComplexF64}; tol::Float64 = 1e-6)
    maximum(imag|>a,s)<tol
end