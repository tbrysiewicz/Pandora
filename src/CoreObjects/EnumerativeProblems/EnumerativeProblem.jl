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
    is_populated,
    VarietyToEnumerativeProblem


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

function EnumerativeProblem(F::System,B::Union{Tuple{Result,Vector{ComplexF64}},Nothing})
    EnumerativeProblem(F,
                B,
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
        S = solutions(S)
    elseif Method == "Monodromy"
        println("Using monodromy")
        M = monodromy_solve(F)
        P = parameters(M)
        S = solutions(M)
    end
    S = solve(F,S,start_parameters=P,target_parameters=P)
    E.BaseFibre=(S,P)
end




function EnumerativeProblem(X::Variety)
#Insert check to see if the variety is actually set up (witness set computed)
	F = system(X)
	V = variables(F)
	N = length(V)
	D = dim(X)
	@var o[1:D,1:1+N]
	Equations = Vector{Expression}(expressions(F))
	for i in 1:D
		l = V'*o[i,1:end-1]-o[i,end]
		push!(Equations,l)
	end
	W = witness_set(X)
	LS = linear_subspace(W)
	A = (extrinsic(LS)).A
	b = (extrinsic(LS)).b
	Ab = hcat(A,b)
	P = vcat(Ab...)
	NewSystem = System(Equations,variables = V, parameters = vcat(o...))
	S = solutions(W)
	R = solve(NewSystem,S;start_parameters=P,target_parameters=P)
	E = EnumerativeProblem(NewSystem,(R,P))
end
