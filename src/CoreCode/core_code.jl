##############################################################
#######################  Exports  ############################
##############################################################

export EnumerativeProperty,
       Fibre,
       AlgorithmDatum,
       KnowledgeNode,
       EnumerativeProblem,
       Citation,
       algorithms_which_return,
       base_parameters,
       base_fibre,
       system,
       degree,
       specialize,
       solutions,
       parameters,
       n_solutions,
       knowledge,
       ambient_dimension,
       n_polynomials,
       n_parameters,
       variables,
       expressions,
       property,
       value,
       input_knowledge,
       input_kwargs,
       algorithm,
       is_real

##############################################################
#############  Warning and Error Messages ####################
##############################################################

const NOALG = "No algorithm datum for this function"


##############################################################
#############  Verbose Printing           ####################
##############################################################
const VERBOSE = Ref(false)

macro vprintln(args...)
    return :(if VERBOSE[]
        println($(esc.(args)...))
    end)
end
macro vprint(args...)
    return :(if VERBOSE[]
        print($(esc.(args)...))
    end)
end
##############################################################
###################    Fibre          ########################
##############################################################

"""
A `Fibre` is a tuple `(S, P)` of two vectors:

- `S::Vector{Vector{ComplexF64}}`: A vector of solutions, each of which is a vector of complex numbers
- `P::Vector{ComplexF64}`: A vector of parameters, each of which is a complex number
"""
const Fibre = Tuple{Vector{Vector{ComplexF64}}, Vector{ComplexF64}}

"""Return the vector of solutions in the fibre."""
solutions(fibre::Fibre) = fibre[1]

"""Return the vector of parameters in the fibre."""
parameters(fibre::Fibre) = fibre[2]

"""Return the number of solutions in the fibre."""
n_solutions(fibre::Fibre) = length(solutions(fibre))

##############################################################
################### Enumerative Property #####################
##############################################################

"""
An `EnumerativeProperty` is a type that represents a property of an 
enumerative problem. It is parameterized by a type `T` and has a name.

Examples:

```julia
const DEGREE = EnumerativeProperty{Int}("degree")
const MONODROMY_GROUP = EnumerativeProperty{Group}("monodromy group")
```
"""
struct EnumerativeProperty{T}
    name::String
end

get_type(::EnumerativeProperty{T}) where {T} = T
name(EProp::EnumerativeProperty) = EProp.name
Base.show(io::IO, EProp::EnumerativeProperty) = print(io, name(EProp))

const DEGREE = EnumerativeProperty{Int}("degree")
const SYSTEM = EnumerativeProperty{System}("system")
const BASE_FIBRE = EnumerativeProperty{Fibre}("base_fibre")
const NULL_ENUMERATIVE_PROPERTY = EnumerativeProperty{Nothing}("null")


include("../citations.jl")
##############################################################
###################   AlgorithmDatum   #######################
##############################################################

"""
`AlgorithmDatum` is a type that represents an algorithm used in enumerative problems.
Its fields are 
- `name`: A string representing the name of the algorithm.
- `description`: A string describing the algorithm.
- `input_properties`: A vector of `EnumerativeProperty` instances that the algorithm takes as input.
- `default_kwargs`: A dictionary of keyword arguments that the algorithm can take, with default values.
- `output_property`: An `EnumerativeProperty` that the algorithm outputs.
- `citations`: A list of `Citation`s that provide reference to the algorithm.
- `reliability`: A symbol indicating the reliability of the algorithm.
"""
Base.@kwdef struct AlgorithmDatum
    name::String = "Unnamed Algorithm"
    description::String = " an algorithm"
    input_properties::Vector{EnumerativeProperty} = EnumerativeProperty[]
    default_kwargs::Dict{Symbol, Any} = Dict{Symbol, Any}()
    output_property::EnumerativeProperty = NULL_ENUMERATIVE_PROPERTY
    citations::Vector{Citation} = [NULL_CITATION]
    reliability::Symbol = :null
end

export name, description, input_properties, default_kwargs, output_property, reliability, citation

name(AD::AlgorithmDatum) = AD.name
description(AD::AlgorithmDatum) = AD.description
input_properties(AD::AlgorithmDatum) = AD.input_properties
default_kwargs(AD::AlgorithmDatum) = AD.default_kwargs
output_property(AD::AlgorithmDatum) = AD.output_property
reliability(AD::AlgorithmDatum) = AD.reliability
citations(AD::AlgorithmDatum) = AD.citations

function core_function(AD::AlgorithmDatum)
    K = collect(keys(ALGORITHM_DATA))
    FK = filter(x -> ALGORITHM_DATA[x] == AD, K)
    if length(FK) == 0
        error("No known functions correspond to this algorithm datum")
    elseif length(FK) > 1
        error("More than one function corresponds to this algorithm datum")
    else
        return FK[1]
    end
end

function Base.show(io::IO, AD::AlgorithmDatum)
    println(io, "Algorithm Datum:        ", name(AD))
    println(io, description(AD))
    println(io, "Input:            ")
    for i in input_properties(AD)
        println(io, "                  ", i)
    end
    println(io, "Output:\n                  ", output_property(AD))
    println(io, "Core Function:    ", core_function(AD))
    println(io, "Citations:         ", [short(c) for c in citations(AD)])
    println(io, "Reliability:      ", reliability(AD))
end

# Intentionally Empty
function user_given()
end

const ANY = EnumerativeProperty{Any}("any")

const user_given_datum = AlgorithmDatum(
    name = "User given",
    description = "The user declared this property",
    input_properties = Vector{EnumerativeProperty}([]),
    output_property = ANY,
    reliability = :user_given
)

global ALGORITHM_DATA = Dict{Function, AlgorithmDatum}()
ALGORITHM_DATA[user_given] = user_given_datum

name(F::Function) = haskey(ALGORITHM_DATA, F) ? name(ALGORITHM_DATA[F]) : error(NOALG)
description(F::Function) = haskey(ALGORITHM_DATA, F) ? description(ALGORITHM_DATA[F]) : error(NOALG)
input_properties(F::Function) = haskey(ALGORITHM_DATA, F) ? input_properties(ALGORITHM_DATA[F]) : error(NOALG)
default_kwargs(F::Function) = haskey(ALGORITHM_DATA, F) ? default_kwargs(ALGORITHM_DATA[F]) : error(NOALG)
output_property(F::Function) = haskey(ALGORITHM_DATA, F) ? output_property(ALGORITHM_DATA[F]) : error(NOALG)
citations(F::Function) = haskey(ALGORITHM_DATA, F) ? citations(ALGORITHM_DATA[F]) : error(NOALG)
reliability(F::Function) = haskey(ALGORITHM_DATA, F) ? reliability(ALGORITHM_DATA[F]) : error(NOALG)

##############################################################
###################   KnowledgeNode   ########################
##############################################################

"""
`KnowledgeNode` is a type that represents a piece of knowledge about an enumerative problem.
It contains the following fields:
- `property`: An `EnumerativeProperty{T}` that this knowledge node represents.
- `value`: The value of the property, of type `T`.
- `input_knowledge`: A vector of `KnowledgeNode` instances that are the inputs to the algorithm that computed this knowledge.
- `input_kwargs`: A dictionary of keyword arguments that were used when computing this knowledge.
- `algorithm`: The function that was used to compute this knowledge.
"""
mutable struct KnowledgeNode{T}
    property::EnumerativeProperty{T}
    value::T
    input_knowledge::Vector{KnowledgeNode}
    input_kwargs::Dict{Symbol, Any}
    algorithm::Function
end

const Knowledge = Vector{KnowledgeNode}

"""
 property(K::KnowledgeNode) returns the `EnumerativeProperty` of the knowledge node `K`.
"""
property(K::KnowledgeNode) = K.property
"""
 value(K::KnowledgeNode) returns the value of the knowledge node `K`.
"""
value(K::KnowledgeNode) = K.value
"""
 input_knowledge(K::KnowledgeNode) returns the input knowledge nodes that were used to compute `K`.
"""
input_knowledge(K::KnowledgeNode) = K.input_knowledge
"""
 input_kwargs(K::KnowledgeNode) returns the keyword arguments that were used to compute `K`.
"""
input_kwargs(K::KnowledgeNode) = K.input_kwargs
"""
 algorithm(K::KnowledgeNode) returns the algorithm that was used to compute `K`.
"""
algorithm(K::KnowledgeNode) = K.algorithm

known_properties(K::Knowledge) = unique([property(k) for k in K])
#TODO: reliability(K::KnowledgeNode) must track all the way back

export 
    property, 
    value, 
    input_knowledge, 
    input_kwargs, 
    algorithm

function Base.show(io::IO, K::Knowledge)
    for (i,k) in enumerate(K)
        print(io, i,") ", string(k), "\n")
    end
end


function Base.display(K::KnowledgeNode)
    print("Property:                ", (property(K)), "\n")
    print("Value:                   ", value(K), "\n")
    print("Algorithm:               ", name(algorithm(K)), "\n")
    print("Input Knowledge:         \n  ", input_knowledge(K), "")
    if length(input_kwargs(K))>0
        print("Input Keyword Arguments: ", input_kwargs(K), "\n")
    end
end
function Base.string(K::KnowledgeNode{T}) where T

    s = "["*string(property(K))
    s *= "] as computed by (" * string(name(algorithm(K))) * ") applied to "
    if length(input_knowledge(K)) == 0
        s *= "(nothing)."
    else
        s *= "[" * string(property(input_knowledge(K)[1]))
        for i in input_knowledge(K)[2:end]
            s *= ", " * string(property(i))
        end
        s *= "]"
    end
    return s
end

function Base.show(io::IO, K::KnowledgeNode{T}) where T
    s = string(K)
    print(io,s)
end
##############################################################
##############   ENUMERATIVE PROBLEM   #######################
##############################################################

"""
`AbstractEnumerativeProblem` is an abstract type that represents a problem in enumerative geometry.
It is the parent type for `EnumerativeProblem`, which contains a parametrized system of equations
for which there are finitely many solutions given a generic parameter. Homotopy continuation
is used to move solutions over one parameter to another. This ability to perform analytic continuation
over the implicit branched cover is what is necessary for an AbstractEnumerativeProblem.
"""
abstract type AbstractEnumerativeProblem end

"""
`EnumerativeProblem` is a concrete type that represents an enumerative problem.
It contains a `System` representing the equations of the problem and a `Knowledge` object
that stores the properties which are known about the problem, including how this knowledge
was obtained, and how reliable it is (e.g. whether it is known with probability 1, high
probability, is a proven result, or is some other form of "knowledge"). 

The `EnumerativeProblem` is the main concrete type used in Pandora.jl.
"""
mutable struct EnumerativeProblem <: AbstractEnumerativeProblem
    system::System
    knowledge::Knowledge

    function EnumerativeProblem(F::System; populate = true)
        EP = new()
        EP.system = F
        EP.knowledge = Knowledge([])
        know!(EP, SYSTEM, F)
        if populate
            populate!(EP)
        end
        return EP
    end
end



function populate!(EP::EnumerativeProblem; kwargs...)
    learn!(EP, BASE_FIBRE; algorithm = polyhedral_homotopy, kwargs...)
    learn!(EP, DEGREE; algorithm = n_solutions, kwargs...)
end

"""
    degree(EP::EnumerativeProblem; kwargs...)

Return the degree of the enumerative problem.
"""
degree(EP::EnumerativeProblem; kwargs...) = DEGREE(EP; kwargs...)

"""
    system(EP::EnumerativeProblem; kwargs...)

Return the system of equations defining the enumerative problem.
"""
system(EP::EnumerativeProblem; kwargs...) = SYSTEM(EP; kwargs...)

"""
    base_fibre(EP::EnumerativeProblem; kwargs...)

Return the base_fibre of the enumerative problem.
"""
function base_fibre(EP::EnumerativeProblem; kwargs...)
    BASE_FIBRE(EP; kwargs...)
end

"""
    knowledge(EP::EnumerativeProblem)
Return the vector of knowledge nodes associated with the enumerative problem.
"""
knowledge(EP::EnumerativeProblem) = EP.knowledge
"""
    variables(EP::EnumerativeProblem)
Return the variables of system of the enumerative problem.
"""
variables(EP::EnumerativeProblem) = variables(SYSTEM(EP))
"""
    parameters(EP::EnumerativeProblem)
Return the parameters of system of the enumerative problem.
"""
parameters(EP::EnumerativeProblem) = parameters(SYSTEM(EP))
"""
    expressions(EP::EnumerativeProblem)
Return the expressions of system of the enumerative problem.
"""
expressions(EP::EnumerativeProblem) = expressions(SYSTEM(EP))
"""
    ambient_dimension(EP::EnumerativeProblem)
Return the ambient dimension of the enumerative problem, which is the number of variables.
"""
ambient_dimension(EP::EnumerativeProblem) = length(variables(EP))
"""
    n_variables(EP::EnumerativeProblem)
Return the number of variables in the system of the enumerative problem.
"""
n_variables(EP::EnumerativeProblem) = ambient_dimension(EP)
"""
    n_polynomials(EP::EnumerativeProblem)
Return the number of polynomials in the system of the enumerative problem.
"""
n_polynomials(EP::EnumerativeProblem) = length(expressions(EP))
"""
    n_parameters(EP::EnumerativeProblem)
Return the number of parameters in the system of the enumerative problem.
"""
n_parameters(EP::EnumerativeProblem) = length(parameters(EP))
"""
    n_parameters(F::System)
Return the number of parameters in a system.
"""
n_parameters(F::System) = length(parameters(F))
"""
    base_parameters(EP::EnumerativeProblem)
Return the parameters of the base fibre of the enumerative problem.
"""
base_parameters(EP::EnumerativeProblem) = base_fibre(EP)[2]
"""
    base_solutions(EP::EnumerativeProblem)
Return the solutions of the base fibre of the enumerative problem.
"""
base_solutions(EP::EnumerativeProblem) = base_fibre(EP)[1]


function is_real(EP::EnumerativeProblem)
    T = typeof(map(E->coefficients(E,vcat(variables(EP),parameters(EP))), expressions(EP)))
    return promote_type(T, Vector{Vector{Float64}}) == Vector{Vector{Float64}}
end


"""
    specialize(F::System; P = nothing)
Return a specialized system for the given system `F` with parameters `P`.
If `P` is not provided, it will generate a random (complex) parameter vector.
"""
function specialize(F::System; P = nothing, real = false)
    if P === nothing
        if !real
            P = randn(ComplexF64, length(parameters(F)))
        else
            P = randn(Float64, length(parameters(F)))
        end
    end
    return System(evaluate(expressions(F), parameters(F) => P))
end

##############################################################
############# Knowledge Learning & Assignment ################
##############################################################

function knowledge_agrees_with_kwargs(k::KnowledgeNode; kwargs...)
    isempty(kwargs) && return true
    ikwargs = input_kwargs(k)
    isempty(ikwargs) && return true
    d = Dict(kwargs)
    for i in keys(ikwargs)
        if haskey(d, i)
            if d[i] != ikwargs[i]
                return false
            end
        end
    end
    return true
end

function knows(EP::EnumerativeProblem, EProp::EnumerativeProperty; kwargs...)
    #TODO: must check that the knowledge node has input kwargs which agree with kwargs...
    candidate_knowledge = filter(k->property(k) == EProp, knowledge(EP))
    kwarg_agreement  = filter(k->knowledge_agrees_with_kwargs(k;kwargs...),candidate_knowledge)
    if length(kwarg_agreement)>0
        return(true)
    else
        return(false)
    end
end


function find_algorithm(EProp::EnumerativeProperty, EP::EnumerativeProblem)
    potential_algorithms = algorithms_which_return(EProp)
    @vprintln("Pandora.jl is automatically finding an algorithm to compute ", EProp, ". To specify an algorithm, call again with algorithm=>[nameofalgorithm]")
    @vprintln("There is a total of ", length(potential_algorithms), " algorithm(s) in Pandora.jl which compute(s) ", name(EProp), ":")
    counter = 0
    length(potential_algorithms) == 0 && return nothing
    algorithm_to_use = potential_algorithms[1]
    for p in potential_algorithms
        counter += 1
        @vprint("      ", counter, ") ")
        if p == algorithm_to_use
            @vprintln("[USING] ", name(p))
        else
            @vprintln(name(p))
        end
    end
    return algorithm_to_use
end

function learn!(EP::EnumerativeProblem, EProp::EnumerativeProperty{T}; 
    algorithm = nothing, kwargs...)      where T

    #If the algorithm is not given, find an algorithm known to Pandora which
    #  computes the given enumerative property
    if algorithm === nothing
        algorithm = find_algorithm(EProp, EP)
    end
    f = algorithm
    @assert(ALGORITHM_DATA[f].output_property == EProp)
    input_knowledge = [get_knowledge!(i, EP; kwargs...) for i in input_properties(f)]
    input_knowledge_values = [value(kn) for kn in input_knowledge]
    kwargs_to_pass = copy(default_kwargs(f))
    if !isempty(kwargs)
        for kv in kwargs
            if haskey(kwargs_to_pass, kv[1])
                kwargs_to_pass[kv[1]] = kv[2]
            end
        end
    end
    o = f(input_knowledge_values...; kwargs_to_pass...)
    new_knowledge = KnowledgeNode{T}(EProp, o, input_knowledge, kwargs_to_pass, f)
    return know!(EP, new_knowledge)
end

function know!(EP::EnumerativeProblem, K::KnowledgeNode)
    push!(EP.knowledge, K)
    return K
end

function know!(EP::EnumerativeProblem, EProp::EnumerativeProperty{T}, value::T) where T
    K = KnowledgeNode(EProp, value, Vector{KnowledgeNode}(), Dict{Symbol, Any}(), user_given)
    know!(EP, K)
end

############ Getter System for EnumerativeProperty of Enumerative Problem #############

function get_knowledge(EProp::EnumerativeProperty, EP::EnumerativeProblem; kwargs...)
    if !knows(EP, EProp; kwargs...)
        return nothing
    end
    candidate_knowledge = filter(k -> property(k) == EProp, knowledge(EP))
    kwarg_agreement = filter(k -> knowledge_agrees_with_kwargs(k; kwargs...), candidate_knowledge)
    if length(kwarg_agreement) == 1
        return kwarg_agreement[1]
    else
        @vprintln("Warning: the property [", EProp, "] with given keyword arguments has more than one knowledge node.")
        return kwarg_agreement[1]
    end
end

function get_knowledge_value(EProp::EnumerativeProperty, EP::EnumerativeProblem; kwargs...)
    K = get_knowledge(EProp, EP; kwargs...)
    if K !== nothing
        return value(K)
    else
        return nothing
    end
end

function get_knowledge!(EProp::EnumerativeProperty, EP::EnumerativeProblem; kwargs...)
    gk = get_knowledge(EProp, EP; kwargs...)
    if gk !== nothing
        return gk
    else
        learn!(EP, EProp; kwargs...)
        return get_knowledge!(EProp, EP; kwargs...)
    end
end

get_knowledge_value!(EProp::EnumerativeProperty, EP::EnumerativeProblem; kwargs...) = value(get_knowledge!(EProp, EP; kwargs...))

#This is the core "getter" for enumerative properties of enumerative problems
function (EProp::EnumerativeProperty)(EP::EnumerativeProblem; kwargs...)
    !knows(EP, EProp; kwargs...)
    get_knowledge_value!(EProp, EP; kwargs...)
end

include("core_enumerative_algorithms.jl")

function Base.show(io::IO, EP::EnumerativeProblem)
    tenspaces="          "
    print(io,"\n\n")
    print(io,tenspaces," X := V(")
    if n_polynomials(EP)==1
        print(io,"f_1")
    else
        print(io,"f_1..f_",n_polynomials(EP),"")
    end
    print(io,") ⊆ C^",ambient_dimension(EP)," x C^",n_parameters(EP),"\n");
    println(io,tenspaces," |")
    println(io,tenspaces," |")
    print(io,tenspaces," | π ")
    if knows(EP,DEGREE)
        println(io,"  ",DEGREE(EP),"-to-1")
    else
        println(io,"   ???-to-1")
    end
    println(io,tenspaces," |")
    println(io,tenspaces," V")
    println(io,tenspaces,"C^",n_parameters(EP),"\n")
    println(io,"An enumerative problem in ",ambient_dimension(EP)," variable(s) cut out by ", 
                n_polynomials(EP)," condition(s) over ", n_parameters(EP)," parameter(s).")
    println(io,"The following information is known about this problem:")
    for k in known_properties(knowledge(EP))
        l = length(filter(x->property(x)==k,knowledge(EP)))
        if l == 1
         println(io,"-",k)
        else
         println(io,"-",k, " (# ways known: ",l,")")
        end
    end
end


include("enumerative_solver.jl")

function algorithms_which_return(EProp::EnumerativeProperty)
    filter(A -> output_property(ALGORITHM_DATA[A]) == EProp, collect(keys(ALGORITHM_DATA)))
end


#=
function knowledge_tree(K::KnowledgeNode; depth = 0)
    if depth == 0
        println("Knowledge Tree for ", property(K), ":")
    end
    println("  "*string(depth), "└── ", string(K))
    if length(input_knowledge(K)) > 0
        for i in input_knowledge(K)
            knowledge_tree(i; depth = depth + 1)
        end 
    else
        println("  "*string(depth), "    └── (no input knowledge)")
    end
end
export knowledge_tree
=#