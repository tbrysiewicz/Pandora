import HomotopyContinuation.degree
import HomotopyContinuation.parameters
import HomotopyContinuation.variables
import HomotopyContinuation.system

export
    EnumerativeProblem,
    ambient_dimension,
    n_polynomials,
    n_parameters,
    degree,
    base_parameters,
    base_solutions,
    base_fibre,
    populate_base_fibre,
    is_populated


################################################################
################################################################
############## Constructor and Show Function####################
################################################################
################################################################

mutable struct EnumerativeProblem
    F::System
    BaseFibre::Union{Tuple{Result,Vector{ComplexF64}},Nothing}
    Context::Dict{Any,Any}
end


@doc raw"""
    EnumerativeProblem(F)

 Constructs the enumerative problem associated to a zero-dimensional parametrized polynomial system $F$.
 Evaluation is `lazy` in that the degree of the enumerative problem will not be computed unless required or asked for.

 # Examples
 ```jldoctest
 julia> @var a,b,c,x;

 julia> F = System([a*x^2+b*x+c],variables=[x],parameters=[a,b,c]);


 julia> E = EnumerativeProblem(F)


           X := V(f_1) ⊆ C^1 x C^3
           |
           |
           | π    ???-to-1
           |
           V
          C^3

 An enumerative problem in 1 variable(s) cut out by 1 condition(s) over 3 parameters.


 julia> degree(E)
 Populating a base fibre of the enumerative problem
 2

 julia> E


           X := V(f_1) ⊆ C^1 x C^3
           |
           |
           | π   2-to-1
           |
           V
          C^3

 An enumerative problem in 1 variable(s) cut out by 1 condition(s) over 3 parameters.

```
"""
function EnumerativeProblem(F::System)
    EnumerativeProblem(F,
                nothing,
                Dict{Any,Any}())
end

function Base.show(io::IO, E::EnumerativeProblem)
    tenspaces="          "
    print(io,"\n\n")
    print(io,tenspaces," X := V(")
    if n_polynomials(E)==1
        print(io,"f_1")
    else
        print(io,"f_1..f_",n_polynomials(E),"")
    end
    print(io,") ⊆ C^",ambient_dimension(E)," x C^",n_parameters(E),"\n");
    println(io,tenspaces," |")
    println(io,tenspaces," |")
    print(io,tenspaces," | π ")
    if is_populated(E)
        println(io,"  ",degree(E),"-to-1")
    else
        println(io,"   ???-to-1")
    end
    println(io,tenspaces," |")
    println(io,tenspaces," V")
    println(io,tenspaces,"C^",n_parameters(E),"\n")
    println(io,"An enumerative problem in ",ambient_dimension(E)," variable(s) cut out by ", 
                n_polynomials(E)," condition(s) over ", n_parameters(E)," parameters.")
end


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


################################################################
################################################################
##############       Populators             ####################
################################################################
################################################################



@doc raw"""
    populate_base_fibre(E::EnumerativeProblem; Retry=false, Method="Explicit")

 Computes a base fibre via the method indicated by `Method`. Options for `Method` include
  -`Explicit`: use the polyhedral homotopy (useful if it is unknown whether the incidence variety is irreducible)
  -`Monodromy`: use monodromy solving to populate a fibre (use only if it is known that the incidence variety is irreducible)

 # Examples
 ```jldoctest

 julia> @var a,b,c,x
 (a, b, c, x)

 julia> F = System([(a*x^2+b)*(a*x^3+b*x+c)],parameters=[a,b,c],variables=[x]);

 julia> E = EnumerativeProblem(F)


           X := V(f_1) ⊆ C^1 x C^3
           |
           |
           | π    ???-to-1
           |
           V
          C^3

 An enumerative problem in 1 variable(s) cut out by 1 condition(s) over 3 parameters.


 julia> populate_base_fibre(E)
 Populating a base fibre of the enumerative problem
 Tracking 5 paths... 100%|████████████████████████████████████████| Time: 0:00:01
  # paths tracked:                  5
  # non-singular solutions (real):  5 (0)
  # singular endpoints (real):      0 (0)
  # total solutions (real):         5 (0)
 (Result with 5 solutions
 =======================
 • 5 paths tracked
 • 5 non-singular solutions (0 real)
 • random_seed: 0xc25fc23d
 • start_system: :polyhedral
 , ComplexF64[-0.17638270351070895 + 0.14941862158991764im, 0.19022741321532455 + 0.7336285101613741im, 0.40856380115222646 + 1.0425948255083233im]) 

 julia> E


           X := V(f_1) ⊆ C^1 x C^3
           |
           |
           | π   5-to-1
           |
           V
          C^3

 An enumerative problem in 1 variable(s) cut out by 1 condition(s) over 3 parameters.

 julia> populate_base_fibre(E; Method = "Monodromy") #This is the incorrect output since the incidence variety is not irreducible
 Populating a base fibre of the enumerative problem
 Using monodromy
 (Result with 2 solutions
 ======================= 
 • 2 paths tracked
 • 2 non-singular solutions (0 real)
 • random_seed: 0x524585f8
 , ComplexF64[-0.7963478582824728 - 0.3920460116495942im, -0.4650125911193124 + 0.12450279217525324im, 0.6104405028524421 - 0.44787385667455143im])

 julia> E #This output is incorrect since monodromy solving should not have been used


           X := V(f_1) ⊆ C^1 x C^3
           |
           |
           | π   2-to-1
           |
           V
          C^3

 An enumerative problem in 1 variable(s) cut out by 1 condition(s) over 3 parameters.


```
"""
function populate_base_fibre(E::EnumerativeProblem; Method="Explicit")
    println("Populating a base fibre of the enumerative problem")
    P = []
    S = []
    F = system(E)
    if Method == "Explicit"
        P = randn(ComplexF64,n_parameters(E))
        S = solve(F,target_parameters=P)
    elseif Method == "Monodromy"
        println("Using monodromy")
        M = monodromy_solve(F)
        P = parameters(M)
        S = solutions(M)
        S = solve(F,S,start_parameters=P,target_parameters=P)
    end
    E.BaseFibre=(S,P)
end



