

################################################################
################################################################
##############       Simple Getters         ####################
################################################################
################################################################


@doc raw"""
    base_fibre(EnumerativeProblem)

 Returns the base fibre of an enumerative problem, which is a pair consisting of 
   (1) A representation of the solutions to an enumerative problem over some (2) specific set of parameters
 (1) Can be either a `Result` or a vector of vectors of complex floats.
 (2) Is a vector of complex floats

 # Examples
 ```jldoctest

 julia> E


           X := V(f_1) ⊆ C^1 x C^3
           |
           |
           | π   2-to-1
           |
           V
          C^3

 An enumerative problem in 1 variable(s) cut out by 1 condition(s) over 3 parameters.


 julia> base_fibre(E)
 (Result with 2 solutions
 =======================
 • 2 paths tracked
 • 2 non-singular solutions (0 real)
 • random_seed: 0x92b57fe0
 • start_system: :polyhedral
 , ComplexF64[-0.2846411030482042 - 0.19493642803540295im, -0.7376884483151196 + 0.8241198115166922im, 0.7632749225104664 + 0.45579670461441946im])

```
"""
function base_fibre(E::EnumerativeProblem)
    E.BaseFibre
end


@doc raw"""
    is_populated(E::EnumerativeProblem)

 `true` if a base fibre of the enumerative problem has been computed and `false` otherwise

 # Examples
 ```jldoctest
 julia> T=TwentySevenLines()


           X := V(f_1..f_4) ⊆ C^4 x C^20
           |
           |
           | π    ???-to-1
           |
           V
          C^20

 An enumerative problem in 4 variable(s) cut out by 4 condition(s) over 20 parameters.

 julia> is_populated(T)
 false

 julia> degree(T)
 Populating a base fibre of the enumerative problem
 27

 julia> is_populated(T)
 true

 julia> T


           X := V(f_1..f_4) ⊆ C^4 x C^20
           |
           |
           | π   27-to-1
           |
           V
          C^20

 An enumerative problem in 4 variable(s) cut out by 4 condition(s) over 20 parameters.


```
"""
function is_populated(E::EnumerativeProblem)
    if typeof(E.BaseFibre)!=Nothing
        return(true)
    else
        return(false)
    end
end


@doc raw"""
    base_solutions(E::EnumerativeProblem)

 Extracts the solutions in the base fibre of an enumerative problem (if computed).

"""
function base_solutions(E::EnumerativeProblem)
    if is_populated(E)
        return(solutions(base_fibre(E)[1]))
    end
    return(nothing)
end


@doc raw"""
    base_parameters(E::EnumerativeProblem)

 Extracts the parameters of the base fibre of an enumerative problem (if computed)
"""
function base_parameters(E::EnumerativeProblem)
    if typeof(E.BaseFibre)!=Nothing
        return(base_fibre(E)[2])
    end
    return(nothing)
end


@doc raw"""
    system(E::EnumerativeProblem)

 Returns the parametrized polynomial system underlying the enumerative problem
"""
function system(E::EnumerativeProblem)
    E.F
end


@doc raw"""
    ambient_dimension(E::EnumerativeProblem)

 Returns the ambient dimension of the solution variety of an enumerative problem. Equivalently, this is the number of variables.
"""
function ambient_dimension(E::EnumerativeProblem)
    length(variables(system(E)))
end


@doc raw"""
    n_polynomials(E::EnumerativeProblem)

 Returns the number of polynomials which are in the underlying parametrized polynomial system associated to an enumerative problem.
"""
function n_polynomials(E::EnumerativeProblem)
    length(expressions(system(E)))
end


@doc raw"""
    n_parameters(E::EnumerativeProblem)

 Returns the dimension of the parameter space of an enumerative problem.

"""
function n_parameters(E::EnumerativeProblem)
    length(parameters(system(E)))
end



@doc raw"""
    degree(E::EnumerativeProblem; Retry=false, Method="Explicit")

 Computes the degree of an enumerative problem.
 If a base fibre has already been computed and `Retry=false` this returns the number of solutions in the base fibre.
 If a base fibre has already been computed and `Retry=true`, or if a base fibre has yet to be computed, this will 
  compute a base fibre via the method indicated by `Method`. Options for `Method` include
  -`Explicit`: use the polyhedral homotopy (useful if it is unknown whether the incidence variety is irreducible)
  -`Monodromy`: use monodromy solving to populate a fibre (use only if it is known that the incidence variety is irreducible)

# Examples
```jldoctest

julia> T = TwentySevenLines();

julia> degree(E)
2

julia> T = TwentySevenLines();

julia> degree(T)
Populating a base fibre of the enumerative problem
27

julia> degree(T)
27

julia> degree(T; Retry=true)
Populating a base fibre of the enumerative problem
27

julia> degree(T; Retry=true, Method="Monodromy")
Populating a base fibre of the enumerative problem
27
```
"""
function degree(E::EnumerativeProblem; Retry=false, Method="Explicit")
    if typeof(E.BaseFibre)==Nothing || Retry==true
        populate_base_fibre(E,Method=Method)
    end
    return(length(base_solutions(E)))
end


@doc raw"""
    variables(E::EnumerativeProblem)

 Returns the variables of the parametrized polynomial system underlying the enumerative problem `E`.

"""
function variables(E::EnumerativeProblem)
    variables(system(E))
end


@doc raw"""
    parameters(E::EnumerativeProblem)

 Returns the parameters of the parametrized polynomial system underlying the enumerative problem `E`.

"""
function parameters(E::EnumerativeProblem)
    parameters(system(E))
end

