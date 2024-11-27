export
    degree_check,
    real_parity_check


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