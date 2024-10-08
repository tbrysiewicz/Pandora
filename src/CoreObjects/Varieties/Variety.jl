import HomotopyContinuation.dim
import HomotopyContinuation.degree
import HomotopyContinuation.witness_set
import HomotopyContinuation.nid
import HomotopyContinuation.system
import Oscar.dim
import Oscar.degree

#TODO: Design all of this so that it works smoothly for the reducible and nonreduced setting

export
    Variety,
    populate_witness!,
    dim,
    degree,
    system,
    ambient_dimension,
    witness_points,
    nid,
    assign_component!,
    sample


mutable struct Variety
    F::System
    W::Union{WitnessSet,Nothing}
    NID::Union{NumericalIrreducibleDecomposition,Nothing}
end

@doc raw"""
    Variety(F)

 Returns the (possibly reduced, possibly reducible) variety defined as the zero set of the polynomial system F.
 # Examples
 ```jldoctest
 julia> @var x,y,z;

 julia> F = System([x^2+y^2-1],variables=[x,y,z]);

 julia> Variety(F)
 A variety in C^3 cut out by 1 polynomials in 3 variables.
 
 ```
 """
function Variety(F::System)
    Variety(F,nothing,nothing)
end


function Base.show(io::IO, V::Variety)
    tenspaces="          "
    if is_populated(V)==false
    println("A variety in C^",ambient_dimension(V)," cut out by ", length(system(V).expressions), " polynomials in ", length(system(V).variables), " variables.")
    else
    println("A variety in C^",ambient_dimension(V)," of degree ",degree(V), " and dimension ",dim(V), ".")
    end
end
@doc raw"""
    system(V::Variety)

 Returns the polynomial system underlying the variety V.

 # Examples
 ```jldoctest
 julia> @var x,y,z;

 julia> F = System([x^2+y^2-1],variables=[x,y,z]);

 julia> V = Variety(F)
 A variety in C^3 cut out by 1 polynomials in 3 variables.

 julia> system(V)
System of length 1
 3 variables: x, y, z

 -1 + x^2 + y^2
 
 ```
 """
function system(V::Variety)
    return(V.F)
end


@doc raw"""
    is_populated(V::Variety)

 Returns `true` if a witness set for V has been computed, and `false` otherwise.

 # Examples
 ```jldoctest
 julia> @var x,y,z;

 julia> F = System([x^2+y^2-1],variables=[x,y,z]);

 julia> V = Variety(F)
 A variety in C^3 cut out by 1 polynomials in 3 variables.

 julia> is_populated(V)
 false

 julia> degree(V)
 2

 julia> is_populated(V)
 true
 
 julia> V
 A variety in C^3 of degree 2 and dimension 2.
 ```
 """
function is_populated(V::Variety)
    if typeof(V.W) != Nothing
        return(true)
    else
        return(false)
    end
end

@doc raw"""
    witness_set(V::Variety)

 Returns a witness set for the variety $V$.

 # Examples
 ```jldoctest
 julia> @var x,y,z;

 julia> F = System([x^2+y^2-1],variables=[x,y,z]);

 julia> V = Variety(F)
 A variety in C^3 cut out by 1 polynomials in 3 variables.

 julia> witness_set(V)
 Witness set for dimension 2 of degree 2
 ```
 """
function witness_set(V::Variety; dim = nothing, codim = nothing)
    if is_populated(V)
        return(V.W)
    else
        populate_witness!(V; dim = dim, codim = codim)
        return(witness_set(V))
    end
end

# this function was created to support numerical_bases in AlgebraicMatroid.jl 
# you must know the dim, d to use
## TODO: perform homotopy continuation on one sample n times instead of pulling n pts from a witness set
        #issue 1: degree of variety may not be at least n
        #issue 2: it is a bit of a mysterious induced by sampling points from the same witness set, ideally they are independent samples of each pt
function get_n_points(V::Variety, d::Int, n)

	MS = monodromy_solve(V.F,dim=d;target_solutions_count=n)

	W = WitnessSet(MixedSystem(V.F),MS.parameters,MS.results)
	
	return(solutions(W))

end

function monodromy_witness_populate!(V::Variety, d; ts = nothing)
	MS = nothing
	if ts == nothing
		MS = monodromy_solve(V.F,dim=d)
	else
		MS = monodromy_solve(V.F,dim=d;target_solutions_count=ts)
	end
	L = MS.parameters
	S = solutions(MS)
	W = WitnessSet(MixedSystem(V.F),MS.parameters,MS.results)
	V.W=W
end


@doc raw"""
    witness_points(V::Variety)

 Returns a set of witness points for the variety $V$.

 # Examples
 ```jldoctest
 julia> @var x,y,z;

 julia> F = System([x^2+y^2-1],variables=[x,y,z]);

 julia> V = Variety(F)
 A variety in C^3 cut out by 1 polynomials in 3 variables.

 julia> witness_points(V)
 2-element Vector{Vector{ComplexF64}}:
  [-0.41446267677331433 - 0.17687291471637545im, 0.9304365983676646 - 0.07878798169660943im, -0.9386659592171888 - 0.9578774541416526im]
  [1.0014636136988397 + 0.017532264384045294im, -0.12765502768352693 + 0.13754197672415844im, 3.964337543997791 - 1.034387642650246im] 
 ```
 """
function witness_points(V::Variety)
    if is_populated(V)
        return(solutions(witness_set(V)))
    else
        populate_witness!(V)
        return(witness_points(V))
    end
end


@doc raw"""
    degree(V::Variety)

 Returns the degree of a variety $V$.

 # Examples
 ```jldoctest
 julia> @var x,y,z;

 julia> F = System([x^2+y^2-1],variables=[x,y,z]);

 julia> V = Variety(F)
 A variety in C^3 cut out by 1 polynomials in 3 variables.

 julia> degree(V)
 2
 ```
 """
function degree(V::Variety; dim = nothing, codim=nothing)
    if is_populated(V)
        return(length(witness_points(V)))
    else
        populate_witness!(V; dim = dim, codim = codim)
        return(degree(V))
    end
end




@doc raw"""
    populate_witness!(V::Variety; dim = nothing, codim = nothing)

 Populates a witness superset for $V$ in dimension/codimension given.
    If no dimension or codimension are provided, the code defaults
    to HC.jl's expected dim/codim.

 """
function populate_witness!(V::Variety; dim = nothing, codim = nothing)
    if dim==nothing && codim == nothing
        println("Warning: no dimension or codimension was given...defaulting to HC.jl's expected dimension")
    end
    W = witness_set(system(V); dim = dim, codim = codim)
    V.W=W
end




@doc raw"""
    ambient_dimension(V::Variety)

 Returns the ambient dimension of the variety $V$.

 """
function ambient_dimension(V::Variety)
    return(length(HomotopyContinuation.variables(V.F)))
end




@doc raw"""
    dim(V::Variety)

 Returns the  dimension of the variety $V$.

 """
function dim(V::Variety; dimm = nothing, codimm = nothing)
    if is_populated(V)==false
        populate_witness!(V; dim = dimm, codim = codimm)
        return(dim(V))
    end
    d = HomotopyContinuation.dim(witness_set(V))
    return(d)
end

function nid(V::Variety)
    NID = HomotopyContinuation.nid(V.F)
    V.NID=NID
end

function assign_component!(V::Variety,dim::Int,num::Int)
    V.W = V.NID.Witness_Sets[dim][num]
    return(nothing)
end

function assign_component!(V::Variety,dim::Int,nums::Vector{Int})
    W = V.NID.Witness_Sets[dim][nums[1]]
    R = copy(W.R)
    for i in nums[2:end] #for each additional component of dimension dim
        Wother = V.NID.Witness_Sets[dim][i] #pull that witness set
        Rnew = solutions((witness_set(Wother,W.L)))  #and move it to the first one's linear space
        for r in Rnew
            push!(R,r) #and collect the witness superset points all together
        end
    end
    WNew = WitnessSet(W.F,W.L,R) #and create a new witness set with that data
    V.W = WNew
    return(nothing)
end


function sample(W::WitnessSet) #Takes witness set and samples random solution from it 
	return(rand(solutions(W)))
end 

#TODO: make it possible to sample n points at random - so that the witness set is randomly
# moved on each iteration
function sample(V::Variety; random = false) #Takes witness set and samples random solution from it 
    if is_populated(V)==true
    	return(rand(solutions(witness_set(V))))
    else
        println("You must first compute a witness set")
    end
end 
