import HomotopyContinuation.is_real
import HomotopyContinuation.nreal

export
    solve_over_param,
    solve_over_params,
    degree_check,
    real_parity_check,
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




@doc raw"""
    degree_check(E::EnumerativeProblem, S)

 Returns `true` if the number of solutions in `S` equal the degree of `E` and returns `false` otherwise.
 """
function degree_check(E::EnumerativeProblem, S)
    if length(S)!=degree(E)
        return(false)
    else
        return(true)
    end
end




@doc raw"""
    real_parity_check(E::EnumerativeProblem, S)

 Returns `true` if the number of real solutions in `S` has the same parity as the degree of `E` and returns `false` otherwise.
    Note: This should return `true` whenever `S` describes a fibre of a real enumerative problem `E` over a real parameter.
 """
function real_parity_check(E::EnumerativeProblem, S)
    if iseven(n_real(S)-degree(E))==false
        return(false)
    else
        return(true)
    end
end



#####Basic moving in parameter space

function solve_over_param(E::EnumerativeProblem,P; monodromy_recover=false)
    #Consider implementing special solvers when
    #  there is decomposability
    if is_populated(E)==false
        populate_base_fibre(E)
    end
    S = solutions(HomotopyContinuation.solve(system(E),base_solutions(E); start_parameters= base_parameters(E), target_parameters = P))
    if monodromy_recover==true && degree_check(E,S)==false
        println("Lost points during tracking...recovering via monodromy")
        M = monodromy_solve(system(E),solutions(S),P)
        S=solutions(M)
    end
    return S
end


function solve_over_params(E::EnumerativeProblem,P; monodromy_recover=false, checks=[degree_check, real_parity_check])
    #Consider implementing special solvers when
    #  there is decomposability
    if is_populated(E)==false
        populate_base_fibre(E)
    end
    S = HomotopyContinuation.solve(system(E),base_solutions(E); start_parameters= base_parameters(E), target_parameters = P)
    println("Total number of fibres computed:",length(S))
    for f in checks
        S=filter!(s->f(E,solutions(s[1])),S)
        println(string(length(S)),"/",length(P)," satisfies ",f)
    end
    if monodromy_recover==true && degree_check(E,S)==false
        println("Lost points during tracking...recovering via monodromy")
        M = monodromy_solve(E.F,solutions(S),P)
        S=solutions(M)
    end
    return S
end
