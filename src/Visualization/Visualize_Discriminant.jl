export
    restrict_enumerative_problem_to_plane







@doc raw"""
    restrict_enumerative_problem_to_plane(EP::EnumerativeProblem)

 Restrict an enumerative problem with >2 parameters to one with 2 parameters by slicing 
   the parameter space with a random 2-plane. This is often used to visualize portions
   of the discriminant of an enumerative problem. 

# Examples
```jldoctest

julia> T = TwentySevenLines()


           X := V(f_1..f_4) ⊆ C^4 x C^20
           |
           |
           | π    ???-to-1
           |
           V
          C^20

An enumerative problem in 4 variable(s) cut out by 4 condition(s) over 20 parameters.


julia> restrict_enumerative_problem_to_plane(T)

           X := V(f_1..f_4) ⊆ C^4 x C^2
           |
           |
           | π    ???-to-1
           |
           V
          C^2

An enumerative problem in 4 variable(s) cut out by 4 condition(s) over 2 parameters.

```
"""
function restrict_enumerative_problem_to_plane(EP::EnumerativeProblem)
	P = [randn(Float64,n_parameters(EP)) for i in 1:3]
	return(restrict_enumerative_problem(EP,P))
end

function restrict_enumerative_problem(EP::EnumerativeProblem,P::Vector{Vector{Float64}})
	F = system(EP)
	xv = variables(F)
	xp = parameters(F)
	n = length(P)
	@var t[1:n-1]
	affine_span = P[n]+sum([t[i].*(P[i]-P[n]) for i in 1:n-1])
	NewEquations = [subs(f,xp=>affine_span) for f in expressions(F)]
	return(EnumerativeProblem(System(NewEquations,variables=xv,parameters=t)))
end
#=
EP = TwentySevenLines()
MyPlots = visualizationWithRefinement(restrict_enumerative_problem(EP,[randn(Float64,20) for i in 1:3]),[-1,1],[-1,1],100,5)
=#

