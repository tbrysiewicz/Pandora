export #main types
    EnumerativeProperty,
    AlgorithmDatum,
    KnowledgeNode,
    EnumerativeProblem


export 
    algorithms_which_return,
    base_parameters,
    base_fibre,
    system,
    degree

##############################################################
#############  Warning and Error Messages ####################
##############################################################
NOALG = "No algorithm datum for this function"

##############################################################
###################    Fibre          #######################
##############################################################
const Fibre = Tuple{Vector{Vector{ComplexF64}},Vector{ComplexF64}}
#const RealFibre = Tuple{Vector{Vector{Float64}},Vector{Float64}}

solutions(fibre::Fibre) = fibre[1]
parameters(fibre::Fibre) = fibre[2]
n_solutions(fibre::Fibre) = length(solutions(fibre))

export solutions, parameters, n_solutions

##############################################################
###################Enumerative Property#######################
##############################################################
struct EnumerativeProperty{T} 
    name :: String
end

get_type(::EnumerativeProperty{T}) where {T} = T
name(EProp::EnumerativeProperty) = EProp.name
Base.show(io::IO, EProp::EnumerativeProperty) =  print(io,name(EProp))

#We cannot define the enumerative property 'degree' since we
#   import functions from HC and Oscar called 'degree'
const DEGREE = EnumerativeProperty{Int}("degree")
const SYSTEM  = EnumerativeProperty{System}("system")
const BASE_FIBRE = EnumerativeProperty{Fibre}("base fibre")
const NULL_ENUMERATIVE_PROPERTY = EnumerativeProperty{Nothing}("null")

##############################################################
#############          Citation          #####################
##############################################################


struct Citation
    authors :: Vector{String}
    title :: String
    journal :: String
    year :: Int
end

global NULL_CITATION = Citation([""],"","",0)



##############################################################
###################      AlgorithmDatum     ##################
##############################################################

struct AlgorithmDatum
    name :: String
    description :: String
    input_properties :: Vector{EnumerativeProperty}
    default_kwargs :: Dict{Symbol,Any}
    output_property :: EnumerativeProperty
	citation :: Citation #the citation for this algorithm (if follows from definitions, give citation for definition)
	reliability :: Symbol
end

function AlgorithmDatum(;
    name :: String  = "Unnamed Algorithm",
    description :: String = " an algorithm",
    input_properties::Union{Vector{EnumerativeProperty},Vector{EnumerativeProperty{T}}} where T  = Vector{EnumerativeProperty}([]),
	default_kwargs :: Dict{Symbol,Any} = Dict{Symbol,Any}(), # a dictionary of non-EP field inputs (e.g. # loops, etc) and the default values
	output_property :: EnumerativeProperty = NULL_ENUMERATIVE_PROPERTY,  #the EP field which is the output of this algorithm 
	citation :: Citation = NULL_CITATION, #the citation for this algorithm (if follows from definitions, give citation for definition)
	reliability :: Symbol = :null)
    AlgorithmDatum(
        name,
        description,
        input_properties,
        default_kwargs,
        output_property,
        citation,
        reliability)
end

name(AD::AlgorithmDatum) = AD.name
description(AD::AlgorithmDatum) = AD.description
input_properties(AD::AlgorithmDatum) = AD.input_properties
default_kwargs(AD::AlgorithmDatum) = AD.default_kwargs
output_property(AD::AlgorithmDatum) = AD.output_property
reliability(AD::AlgorithmDatum) =  AD.reliability
citation(AD::AlgorithmDatum) =  AD.citation

function core_function(AD::AlgorithmDatum)
    K = collect(keys(ALGORITHM_DATA))
    FK = filter(x->ALGORITHM_DATA[x]==AD,K)
    if length(FK)==0
        error("No known functions correpsond to this algorithm datum")
    elseif length(FK)>1
        error("More than one function correpsonds to this algorithm datum")
    else
        return(FK[1])
    end
end


function Base.show(io::IO, AD::AlgorithmDatum)
    println(io,"Algorithm Datum:        ",name(AD))
    println(io,description(AD))
    println(io,"Input:            ")
    if reliability(AD)==:user_given
    print("None (the algorithm always returns the user-given value) \n")
    end
    for i in input_properties(AD)
        println(io,"                  ",i)
    end
    println(io,"Output:\n                  ",output_property(AD))
    println(io,"Core Function:    ",core_function(AD))
    println(io,"Citation:         ",citation(AD))
    println(io,"Reliability:      ",reliability(AD))
end


function user_given(EProp::EnumerativeProperty{T},value::T) where T 
    function c()
        return(value)
    end
    return(c)
end

const user_given_datum = AlgorithmDatum(;
name = "User given",
description = "The user declared this property",
input_properties = Vector{EnumerativeProperty}([]),
citation = NULL_CITATION,
reliability = :user_given
)


global ALGORITHM_DATA = Dict{Function,AlgorithmDatum}([])

name(F::Function) = haskey(ALGORITHM_DATA, F) ? name(ALGORITHM_DATA[F]) : error(NOALG)
description(F::Function) = haskey(ALGORITHM_DATA, F) ? description(ALGORITHM_DATA[F]) : error(NOALG)
input_properties(F::Function) = haskey(ALGORITHM_DATA, F) ? input_properties(ALGORITHM_DATA[F]) : error(NOALG)
default_kwargs(F::Function) = haskey(ALGORITHM_DATA, F) ? default_kwargs(ALGORITHM_DATA[F]) : error(NOALG)
output_property(F::Function) = haskey(ALGORITHM_DATA, F) ? output_property(ALGORITHM_DATA[F]) : error(NOALG)
citation(F::Function) = haskey(ALGORITHM_DATA, F) ? citation(ALGORITHM_DATA[F]) : error(NOALG)
reliability(F::Function) = haskey(ALGORITHM_DATA, F) ? reliability(ALGORITHM_DATA[F]) : error(NOALG)



##############################################################
###################      KnowledgeNode     ###################
##############################################################


mutable struct KnowledgeNode{T}
    property :: EnumerativeProperty{T}      #For example MonodromyGroup is of EnumerativeProperty{Group}
    value :: T                                  #The monodromy group G as a Group computed via 
                                                    # algorithm(input_knowledge...;input_parameters...)
    input_knowledge :: Vector{KnowledgeNode}      #The input knowledge nodes of the enumerative problem which form 
                                                    #the input to the corresponding algorithm
    input_kwargs :: Dict{Symbol,Any}        #Some algorithms have optional arguments, these are those
    algorithm :: Function           #The algorithm itself
end

const Knowledge = Vector{KnowledgeNode}


property(K::KnowledgeNode) = K.property
value(K::KnowledgeNode) = K.value
input_knowledge(K::KnowledgeNode) = K.input_knowledge
input_kwargs(K::KnowledgeNode) = K.input_kwargs
algorithm(K::KnowledgeNode) = K.algorithm

known_properties(K::Knowledge) = unique([property(k) for k in K])
    

#TODO: reliability(K::KnowledgeNode) must track all the way back

export property, value, input_knowledge, input_kwargs, algorithm

function Base.show(io::IO, K::KnowledgeNode)
    print(io,"[",property(K),"] is known as a consequence of [", name(algorithm(K)),"] ")
    if length(keys(input_kwargs(K)))>0
        print(io, "(with keyword arguments ",[input_kwargs(K)[i] for i in keys(input_kwargs(K))],") ")
    end
    print(io,"applied to knowledge of ")
    if length(input_knowledge(K))==0
        print(io,"(nothing).")
    else
        print(io,"[",property(input_knowledge(K)[1]))
        for i in input_knowledge(K)[2:end]
            print(io,", ",property(i))
        end
        print(io,"].")
    end
end



################################################################
################################################################
##############        ENUMERATIVE PROBLEM   ####################
################################################################
################################################################

abstract type AbstractEnumerativeProblem end


mutable struct EnumerativeProblem <: AbstractEnumerativeProblem
    system :: System
    knowledge :: Knowledge

    function EnumerativeProblem(F::System; populate = true)
        EP = new()
        EP.system = F
        EP.knowledge = Knowledge([])
        #EP.hc_options = Dict{Any,Any}()
        #EP.hc_options[:tracker_options]=TrackerOptions()

        know!(EP,SYSTEM,F)

        if populate
            populate!(EP)
        end
        return(EP)
    end
end

function populate!(EP::EnumerativeProblem; kwargs...)
    learn!(EP,BASE_FIBRE; algorithm = polyhedral_homotopy, kwargs...)
    learn!(EP,DEGREE; algorithm = n_solutions, kwargs...)
end


degree(EP::EnumerativeProblem; kwargs...) = DEGREE(EP; kwargs...)
system(EP::EnumerativeProblem; kwargs...) = SYSTEM(EP; kwargs...)
base_fibre(EP::EnumerativeProblem; kwargs...) = BASE_FIBRE(EP; kwargs...)



tracker_options(EP::EnumerativeProblem) = EP.hc_options[:tracker_options]
knowledge(EP::EnumerativeProblem) = EP.knowledge
variables(EP::EnumerativeProblem) = variables(SYSTEM(EP))
parameters(EP::EnumerativeProblem) = parameters(SYSTEM(EP))
expressions(EP::EnumerativeProblem) = expressions(SYSTEM(EP))
ambient_dimension(EP::EnumerativeProblem) = length(variables(EP))
n_variables(EP::EnumerativeProblem) = ambient_dimension(EP)
n_polynomials(EP::EnumerativeProblem) = length(expressions(EP))
n_parameters(EP::EnumerativeProblem) = length(parameters(EP))
base_parameters(EP::EnumerativeProblem) = base_fibre(EP)[2]




export knowledge, ambient_dimension, n_polynomials, n_parameters, variables, parameters, expressions, degree


#######################################################################################################################
###Knowledge Learning: learn!(EP,EProperty,algorithm) -> learn EProperty for EP via algorithm #########################
###                    learn!(EP,Eproperty) -> find an algorithm to learn EProperty and learn it ######################
#######################################################################################################################
###Knowledge Assignment:                                                                                ###############
###    know!(EP,KnowledgeNode) -> know the knowledge node                                                        ######
#######################################################################################################################
#######################################################################################################################
###User given knowledge:                                                                                ###############
### know!(EP,EProperty,val) -> Create user-given knowledge node that EP.Eproperty = val and discover that knowledge####
#######################################################################################################################

function knowledge_agrees_with_kwargs(k::KnowledgeNode; kwargs...)
    isempty(kwargs) && return(true)
    ikwargs = input_kwargs(k)
    isempty(ikwargs) && return(true)
    d = Dict(kwargs)
    for i in keys(ikwargs)
        if haskey(d,i)
            if d[i]!=ikwargs[i]
                return(false)
            end
        end
    end
    return(true)
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


function find_algorithm(EProp::EnumerativeProperty,EP::EnumerativeProblem)
    potential_algorithms = algorithms_which_return(EProp)
    println("Pandora.jl is automatically finding an algorithm to compute ",EProp,". To specify an algorithm, call again with algorithm=>[nameofalgorithm]")
    println("There is a total of ",length(potential_algorithms), " algorithm(s) in Pandora.jl which compute(s) ",name(EProp),":")
    counter=0
    length(potential_algorithms)==0 && return(nothing)
    algorithm_to_use = potential_algorithms[1]
    for p in potential_algorithms
        counter=counter+1

        print("      ",counter,") ")
        if p == algorithm_to_use
            println("[USING] ",name(p))
        else
            println(name(p))
        end
    end

    return(algorithm_to_use)
end

function learn!(EP::EnumerativeProblem, EProp::EnumerativeProperty{T}; 
    algorithm = nothing, kwargs...)      where T

    if algorithm === nothing
        algorithm = find_algorithm(EProp,EP)  #find and algorithm which will obtain it
    end
    
    f = algorithm
    input_knowledge = [get_knowledge!(i,EP) for i in input_properties(f)]
    input_knowledge_values = [value(kn) for kn in input_knowledge]
    kwargs_to_pass = copy(default_kwargs(f))
    if !isempty(kwargs)
       for kv in kwargs
           if haskey(kwargs_to_pass,kv[1])
               kwargs_to_pass[kv[1]] = kv[2]
           end
        end
    end
	o = f(input_knowledge_values...;kwargs_to_pass...)
    new_knowledge = KnowledgeNode{T}(EProp,o,input_knowledge,kwargs_to_pass,f)
	return(know!(EP,new_knowledge))
end

function know!(EP::EnumerativeProblem,K::KnowledgeNode)
    push!(EP.knowledge,K)
    return(K)
end

function know!(EP::EnumerativeProblem, EProp::EnumerativeProperty{T}, value::T) where T
    f = user_given(EProp,value)
    ALGORITHM_DATA[f]=user_given_datum
    learn!(EP,EProp; algorithm = f)
end



############Getter System for EnumerativeProperty of Enumerative Problem#############

function get_knowledge(EProp::EnumerativeProperty,EP::EnumerativeProblem; kwargs...)
    if knows(EP,EProp; kwargs...)==false
        return(nothing)
    end
    candidate_knowledge = filter(k->property(k) == EProp, knowledge(EP))
    kwarg_agreement  = filter(k->knowledge_agrees_with_kwargs(k;kwargs...),candidate_knowledge)
    if  length(kwarg_agreement)==1
        return(kwarg_agreement[1])
    else
        println("Warning: the property [",EProp,"] with given keyword arguments has more than one knowledge node.")
        return(kwarg_agreement[1])
    end
end

function get_knowledge_value(EProp::EnumerativeProperty,EP::EnumerativeProblem; kwargs...)
    K = get_knowledge(EProp,EP; kwargs...)
    if K !== nothing
        return(value(K))
    else
        return(nothing)
    end
end


function get_knowledge!(EProp::EnumerativeProperty,EP::EnumerativeProblem; kwargs...)
    gk = get_knowledge(EProp,EP; kwargs...)
    if gk !== nothing #Check if knowledge of this property is already known
                     #TODO: With input_kwargs which agree with kwargs!!!!
        return(gk)   #If so, return it
    else             #Otherwise,
        learn!(EP,EProp; kwargs...)
        return(get_knowledge!(EProp,EP;kwargs...))
    end
end

#Collects a knowledge node and then returns its value
get_knowledge_value!(EProp::EnumerativeProperty,EP::EnumerativeProblem; kwargs...)  =  value(get_knowledge!(EProp,EP; kwargs...))

#This is the core "getter" for enumerative properties of enumerative problems
function (EProp::EnumerativeProperty)(EP::EnumerativeProblem; kwargs...) 
    !knows(EP, EProp; kwargs...) 
    get_knowledge_value!(EProp,EP; kwargs...)
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


function update_base_fibre!(EP::EnumerativeProblem,F::Fibre)
    know!(EP,BASE_FIBRE,F)
    learn!(EP,DEGREE; algorithm = n_solutions)
end

include("solving.jl")

function algorithms_which_return(EProp::EnumerativeProperty)
    filter(A->output_property(ALGORITHM_DATA[A])==EProp,collect(keys(ALGORITHM_DATA)))
end

