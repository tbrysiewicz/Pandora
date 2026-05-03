##############################################################
#######################  Exports  ############################
##############################################################

export EnumerativeProperty,
       EnumerativeData,
       EnumerativeAttribute,
       Fibre,
       AlgorithmDatum,
       KnowledgeNode,
       EnumerativeProblem,
       Citation,
       algorithms_which_return,
       base_parameters,
       base_fibre,
       inequations,
       system,
       degree,
       specialize,
       solutions,
       parameters,
       n_solutions,
       knowledge,
       enumerative_data,
       properties,
       data,
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
       is_real,
       verbose

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

function verbose()
    if VERBOSE[]
        VERBOSE[] = false
        println("Verbose mode is now off.")
    else
        VERBOSE[] = true
        println("Verbose mode is now on.")
    end
end
##############################################################
###################    Fibre          ########################
##############################################################

@kwdef mutable struct FibreDatum
    parameters::Vector{ComplexF64} = ComplexF64[]
    solutions::Vector{Vector{ComplexF64}} = Vector{Vector{ComplexF64}}()
    function_values::Dict{Function, Any} = Dict{Function, Any}()
    certificates::Union{CertificationResult, Nothing} = nothing
    F::Union{System, Nothing} = nothing
end

# Essentially makes FibreDatum behave like a pair of solutions and parameters (which is a Fibre)
function Base.getindex(F::FibreDatum, i::Int)
    if i == 1
        return F.solutions
    elseif i == 2
        return F.parameters
    end
end

FibreDatum(FD::FibreDatum) = FD

"""
A `Fibre` is a tuple `(S, P)` of two vectors:

- `S::Vector{Vector{ComplexF64}}`: A vector of solutions, each of which is a vector of complex numbers
- `P::Vector{ComplexF64}`: A vector of parameters, each of which is a complex number
"""
const Fibre = Union{Tuple{Vector{Vector{ComplexF64}}, Vector{ComplexF64}}, FibreDatum}

Fibre(S::Vector{Vector{ComplexF64}}, P::Vector{Float64}) = (S, Vector{ComplexF64}(P))
Fibre(S::Vector{Vector{ComplexF64}}, P::Vector{ComplexF64}) = (S, Vector{ComplexF64}(P))
Fibre(F::Tuple{Vector{Vector{ComplexF64}}, Vector{Float64}}) = (F[1], Vector{ComplexF64}(F[2]))
Fibre(F::Tuple{Vector{Vector{ComplexF64}}, Vector{ComplexF64}}) = F
Fibre(F::Tuple{Result, Vector{Float64}}) = (solutions(F[1]), Vector{ComplexF64}(F[2]))
Fibre(F::Tuple{Result, Vector{ComplexF64}}) = (solutions(F[1]), Vector{ComplexF64}(F[2]))

FibreDatum(F::Tuple{Vector{Vector{ComplexF64}}, Vector{ComplexF64}}) = FibreDatum(parameters = F[2], solutions = F[1])
FibreDatum(F::Tuple{Vector{Vector{ComplexF64}}, Vector{Float64}}) = FibreDatum(parameters = Vector{ComplexF64}(F[2]), solutions = F[1])

"""Return the vector of solutions in the fibre."""
solutions(fibre::Fibre) = fibre[1]

"""Return the vector of parameters in the fibre."""
parameters(fibre::Fibre) = fibre[2]

"""Return the number of solutions in the fibre."""
n_solutions(fibre::Fibre) = length(solutions(fibre))

include("enumerative_property.jl")
include("enumerative_data.jl")

include("../citations.jl")
##############################################################
###################   AlgorithmDatum   #######################
##############################################################

"""
`AlgorithmDatum` is a type that represents an algorithm used in enumerative problems.
Its fields are
- `name`: A string representing the name of the algorithm.
- `description`: A string describing the algorithm.
- `input_properties`: A vector of `EnumerativeProperty` or `EnumerativeData` instances that the algorithm takes as input.
- `default_kwargs`: A dictionary of keyword arguments that the algorithm can take, with default values.
- `output_property`: An `EnumerativeProperty` or `EnumerativeData` instance that the algorithm outputs.
- `citations`: A list of `Citation`s that provide reference to the algorithm.
- `reliability`: A symbol indicating the reliability of the algorithm.
"""
Base.@kwdef struct AlgorithmDatum
    name::String = "Unnamed Algorithm"
    description::String = " an algorithm"
    input_properties::Vector{EnumerativeAttribute} = EnumerativeAttribute[]
    default_kwargs::Dict{Symbol, Any} = Dict{Symbol, Any}()
    output_property::EnumerativeAttribute = NULL_ENUMERATIVE_PROPERTY
    citations::Vector{Citation} = [NULL_CITATION]
    reliability::Symbol = :null
    automated::Bool = true
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
    name = "user_given",
    description = "The user declared this property",
    input_properties = EnumerativeAttribute[],
    output_property = ANY,
    reliability = :user_given,
    automated = false
)

function conjunction()
end

const conjunction_datum = AlgorithmDatum(
    name = "conjunction",
    description = "Combines multiple knowledge nodes into one",
    input_properties = EnumerativeAttribute[],
    output_property = ANY,
    reliability = :certified,
    automated = false
)

global ALGORITHM_DATA = Dict{Function, AlgorithmDatum}()
ALGORITHM_DATA[user_given] = user_given_datum
ALGORITHM_DATA[conjunction] = conjunction_datum
const DO_NOT_AUTOMATE = EnumerativeData{Vector{FibreDatum}}("fibre_data")

name(F::Function) = haskey(ALGORITHM_DATA, F) ? name(ALGORITHM_DATA[F]) : error(NOALG)
description(F::Function) = haskey(ALGORITHM_DATA, F) ? description(ALGORITHM_DATA[F]) : error(NOALG)
input_properties(F::Function) = haskey(ALGORITHM_DATA, F) ? input_properties(ALGORITHM_DATA[F]) : error(NOALG)
default_kwargs(F::Function) = haskey(ALGORITHM_DATA, F) ? default_kwargs(ALGORITHM_DATA[F]) : error(NOALG)
output_property(F::Function) = haskey(ALGORITHM_DATA, F) ? output_property(ALGORITHM_DATA[F]) : error(NOALG)
citations(F::Function) = haskey(ALGORITHM_DATA, F) ? citations(ALGORITHM_DATA[F]) : error(NOALG)
reliability(F::Function) = haskey(ALGORITHM_DATA, F) ? reliability(ALGORITHM_DATA[F]) : error(NOALG)

include("knowledge.jl")

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
It contains a `System` representing the equations of the problem, a list of
intrinsic properties known about the problem, and a list of computational data
attached to the problem.

The `EnumerativeProblem` is the main concrete type used in Pandora.jl.
"""
mutable struct EnumerativeProblem <: AbstractEnumerativeProblem
    properties::Knowledge
    data::Knowledge

    function EnumerativeProblem(F::System; inequations = System([]), populate = true, torus_only = false, certify = false, monodromy = false)
        EP = new()
        EP.properties = Knowledge([])
        EP.data = Knowledge([])
        know!(EP, SYSTEM, F)
        if torus_only
            inequations = System([prod(variables(F))])
        end
        know!(EP, INEQUATIONS, inequations)
        if populate
            populate!(EP; monodromy = monodromy)
        end
        if certify
            println("Certifying")
            FD = fibre_datum(EP)
            println(FD)
            remove_knowledge!(EP, BASE_FIBRE)
            remove_knowledge!(EP, DEGREE)
            know!(EP, BASE_FIBRE, (FD[1], FD[2]))
            learn!(EP, DEGREE; algorithm = n_solutions)
        end
        return EP
    end

    function EnumerativeProblem(F::System, inequations::System; populate = true, torus_only = false, monodromy = false)
        EP = new()
        EP.properties = Knowledge([])
        EP.data = Knowledge([])
        know!(EP, SYSTEM, F)
        if torus_only
            inequations = System(vcat(expressions(inequations), [prod(variables(F))]))
        end
        know!(EP, INEQUATIONS, inequations)

        if populate
            populate!(EP; monodromy = monodromy)
        end
        return EP
    end
 
    function EnumerativeProblem(EX::Vector{Expression}; variables::Vector{Variable} = [],
        parameters::Vector{Variable} = [], populate = true, inequations = System([]), torus_only = false, monodromy = false)

        F = System(EX, variables = variables, parameters = parameters)
        return EnumerativeProblem(F, inequations; torus_only = torus_only, populate = populate, monodromy = monodromy)
    end
end

function populate!(EP::EnumerativeProblem; kwargs...)
    # Check if monodromy is true in kwargs; if so, learn the base fibre using monodromy.
    if get(kwargs, :monodromy, false)
        learn!(EP, BASE_FIBRE; algorithm = monodromy, kwargs...)
    else
        learn!(EP, BASE_FIBRE; algorithm = polyhedral_homotopy, kwargs...)
    end
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
    inequations(EP::EnumerativeProblem; kwargs...)

Return the inequations imposed in the enumerative problem.
"""
inequations(EP::EnumerativeProblem; kwargs...) = INEQUATIONS(EP; kwargs...)

"""
    base_fibre(EP::EnumerativeProblem; kwargs...)

Return the base_fibre of the enumerative problem.
"""
function base_fibre(EP::EnumerativeProblem; kwargs...)
    BASE_FIBRE(EP; kwargs...)
end

"""
    knowledge(EP::EnumerativeProblem)
Return all property and data nodes associated with the enumerative problem.
"""
knowledge(EP::EnumerativeProblem) = vcat(properties(EP), data(EP))
"""
    properties(EP::EnumerativeProblem)
Return the vector of property nodes associated with the enumerative problem.
"""
properties(EP::EnumerativeProblem) = EP.properties
"""
    enumerative_data(EP::EnumerativeProblem)
Return the vector of data nodes associated with the enumerative problem.
"""
enumerative_data(EP::EnumerativeProblem) = data(EP)
"""
    data(EP::EnumerativeProblem)
Return the vector of data nodes associated with the enumerative problem.
"""
data(EP::EnumerativeProblem) = EP.data
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
    n_inequations(EP::EnumerativeProblem)
Return the number of inequations in the system of the enumerative problem.
"""
n_inequations(EP::EnumerativeProblem) = length(expressions(inequations(EP)))

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
    T = typeof(map(E -> coefficients(E, vcat(variables(EP), parameters(EP))), expressions(EP)))
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

function Base.show(io::IO, EP::EnumerativeProblem)
    tenspaces = "          "
    print(io, "\n\n")
    print(io, tenspaces, " X := V(")
    if n_polynomials(EP) == 1
        print(io, "f_1")
    else
        print(io, "f_1..f_", n_polynomials(EP), "")
    end
    print(io, ") ⊆ C^", ambient_dimension(EP), " x C^", n_parameters(EP), "\n")
    println(io, tenspaces, " |")
    println(io, tenspaces, " |")
    print(io, tenspaces, " | π ")
    if knows(EP, DEGREE)
        println(io, "  ", DEGREE(EP), "-to-1")
    else
        println(io, "   ???-to-1")
    end
    println(io, tenspaces, " |")
    println(io, tenspaces, " V")
    println(io, tenspaces, "C^", n_parameters(EP), "\n")
    println(
        io,
        "An enumerative problem in ", ambient_dimension(EP), " variable(s) cut out by ",
        n_polynomials(EP), " condition(s) over ", n_parameters(EP), " parameter(s).",
    )
    @vprintln(io, "The following information is known about this problem:")
    for k in known_properties(properties(EP))
        l = length(filter(x -> property(x) == k, properties(EP)))
        if l == 1
            @vprintln(io, "-", k)
        else
            @vprintln(io, "-", k, " (# ways known: ", l, ")")
        end
    end
    @vprintln(io, "The following data has been computed for this problem:")
    for k in known_data(data(EP))
        l = length(filter(x -> property(x) == k, data(EP)))
        if l == 1
            @vprintln(io, "-", k)
        else
            @vprintln(io, "-", k, " (# datasets: ", l, ")")
        end
    end
end

include("knowledge_management.jl")
include("core_enumerative_algorithms.jl")
include("enumerative_solver.jl")
