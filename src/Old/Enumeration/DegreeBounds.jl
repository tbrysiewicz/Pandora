export
    bezout,
    bkk


@doc raw"""
    bezout(E::EnumerativeProblem)

 Returns the Bezout number of the polynomial system underlying E (if that system is square) defined
  as the product of the degrees of the polynomials in that system.
 The Bezout number bounds the number of isolated solutions to the polynomial system, and hence, the degree
  of the enumerative problem.
 # Examples
 ```jldoctest
 julia> T = TwentySevenLines();

 julia> bezout(T)
 81
 ```
 """
function bezout(E::EnumerativeProblem)
    F = system(E)
    G = [HomotopyContinuation.evaluate(f,F.parameters=>randn(ComplexF64,length(F.parameters))) for f in F.expressions]
    degs= ([HomotopyContinuation.degree(g) for g in G])
    return(prod(degs))
end



@doc raw"""
  bkk(E::EnumerativeProblem)

 Returns the BKK number of the polynomial system underlying E (if that system is square) defined
  as the mixed volume of the Newton polytopes of the polynomials in that system.
  By default, the origin is included in each Newton polytope so the BKK bound applies to the number
  of solutions over affine space, instead of its usual bound of isolated solutions over the algebraic torus.
 # Examples
 ```jldoctest

 julia> T = TwentySevenLines();

 julia> bkk(T)
 45
 
 ```
 """
function bkk(E::EnumerativeProblem;only_torus=false)
    P = randn(Float64,n_parameters(E))
    HomotopyContinuation.paths_to_track(specialized_system(E);only_torus=only_torus)
end

function specialized_system(E::EnumerativeProblem; P = nothing)
    if P==nothing   
        P = randn(Float64,length(parameters(system(E))))
    end
    F = system(E)
    return(System(HomotopyContinuation.evaluate(expressions(F), parameters(F)=>P)))
end