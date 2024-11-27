import HomotopyContinuation.is_real
import HomotopyContinuation.nreal

export
    is_real,
    n_real


    
@doc raw"""
    is_real(s::Vector{ComplexF64};tol=1e-6)

 Determine if a vector of complex numbers is (approximately) real (i.e. the largest imaginary part is smaller than `tol`)
 """
function is_real(s::Vector{ComplexF64};tol=1e-6)
    for i in s
        if abs(imag(i))>tol
            return(false)
        end
    end
    return(true)
end


@doc raw"""
    n_real(S::Vector{Vector{ComplexF64}};tol=1e-6)

 Counts the number of complex vectors in `S` which are (approximately) real (i.e. the largest imaginary part is smaller than `tol`)
 """
function n_real(S::Vector{Vector{ComplexF64}};tol=1e-6)
    count(s->is_real(s),S)
end
