export #main types
    EnumerativeProperty,
    EnumerativeAlgorithm,
    KnowledgeNode,
    EnumerativeProblem


export #callable on enumerative problems
    degree

export #enumerative properties (also callable on enumerative problems)
    system


##############################################################
#############  Warning and Error Messages ####################
##############################################################
const NeedToComputeWarning = "Computing...\n"

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
const enumerative_degree = EnumerativeProperty{Int}("degree")
const system  = EnumerativeProperty{System}("system")
const base_fibre = EnumerativeProperty{Fibre}("base fibre")

export enumerative_degree, system, base_fibre

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
#############  Enumerative Algorithm     #####################
##############################################################

struct EnumerativeAlgorithm
    name::String
	input_properties::Vector{EnumerativeProperty} 
	default_kwargs :: Dict{Symbol,Any}
	core_function :: Function
    output_property :: EnumerativeProperty
	citation :: Citation #the citation for this algorithm (if follows from definitions, give citation for definition)
	reliability :: Symbol
end

function EnumerativeAlgorithm(;
    name :: String  = "Unnamed Algorithm",
    input_properties::Union{Vector{EnumerativeProperty},Vector{EnumerativeProperty{T}}} where T  = Vector{EnumerativeProperty}([]),
	default_kwargs :: Dict{Symbol,Any} = Dict{Symbol,Any}(), # a dictionary of non-EP field inputs (e.g. # loops, etc) and the default values
	core_function :: Function = x->nothing,
	output_property :: EnumerativeProperty,  #the EP field which is the output of this algorithm 
	citation :: Citation = NULL_CITATION, #the citation for this algorithm (if follows from definitions, give citation for definition)
	reliability :: Symbol = :null)
    EnumerativeAlgorithm(
        name,
        input_properties,
        default_kwargs,
        core_function,
        output_property,
        citation,
        reliability)
end

name(EA::EnumerativeAlgorithm) = EA.name
input_properties(EA::EnumerativeAlgorithm) = EA.input_properties
default_kwargs(EA::EnumerativeAlgorithm) = EA.default_kwargs
core_function(EA::EnumerativeAlgorithm) = EA.core_function
output_property(EA::EnumerativeAlgorithm) = EA.output_property
reliability(EA::EnumerativeAlgorithm) =  EA.reliability
citation(EA::EnumerativeAlgorithm) =  EA.citation

function Base.show(io::IO, EA::EnumerativeAlgorithm)
    tenspaces="          "
    println(io,"Algorithm:       ",name(EA))
    println(io,"----------------------------------------------------")
    println(io,"Input:            ")
    if reliability(EA)==:user_given
    print("None (the algorithm always returns the user-given value) \n")
    end
    for i in input_properties(EA)
        println(io,"                    ",i)
    end
    println(io,"Output:           ",output_property(EA))
    println(io,"Citation:         ",citation(EA))
    println(io,"Reliability: ",reliability(EA))
end

function user_given(EProp::EnumerativeProperty{T},value::T) where T 
    alg_name = "User given "*name(EProp)
    input_properties = Vector{EnumerativeProperty}([])
    default_kwargs = Dict{Symbol,Any}()
    function c()
        return(value)
    end
    core_function = c
    output_property = EProp
    citation = NULL_CITATION
    reliability = :user_given
    return(EnumerativeAlgorithm(
        alg_name,
        input_properties,
        default_kwargs,
        core_function,
        output_property,
        citation,
        reliability))
end


global MAIN_ALGORITHMS = Vector{EnumerativeAlgorithm}([])



##############################################################
###################      KnowledgeNode     ###################
##############################################################


mutable struct KnowledgeNode{T}
    property :: EnumerativeProperty{T}      #For example MonodromyGroup is of EnumerativeProperty{Group}
    value :: T                                  #The monodromy group G as a Group computed via 
                                                    # algorithm(input_knowledge...;input_parameters...)
    input_knowledge :: Vector{KnowledgeNode}      #The input knowledge nodes of the enumerative problem which form 
                                                    #the input to the EnumerativeAlgorithm
    input_kwargs :: Dict{Symbol,Any}        #Some algorithms have optional arguments, these are those
    algorithm :: EnumerativeAlgorithm           #The algorithm itself
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
    println(io,"Knowledge node for [",property(K),"]:")
 #   println(io,value(K))
    println(io,"----------------------------")
    println(io,"Algorithm: [",name(algorithm(K)),"]")
    if length(input_knowledge(K))>0
        println(io,"Ran on knowledge nodes for ")
    end
    for i in input_knowledge(K)
        println(io,"     ",name(property(i)))
    end
    if length(keys(input_kwargs(K)))>0
        println(io,"Keyword arguments")
        for ip in keys(input_kwargs(K))
            println(io,"      ",ip,": ",(input_kwargs(K))[ip])
        end
    end
    print("")
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
    hc_options :: Dict{Any,Any}

    function EnumerativeProblem(F::System; populate = true)
        EP = new()
        EP.system = F
        EP.knowledge = Knowledge([])
        EP.hc_options = Dict{Any,Any}()
        EP.hc_options[:tracker_options]=TrackerOptions()

        know!(EP,system,F)

        if populate
            populate!(EP)
        end
        return(EP)
    end
end

function populate!(EP::EnumerativeProblem; kwargs...)
    learn!(EP,base_fibre; algorithm = solve_generic_via_start_system, kwargs...)
    learn!(EP,enumerative_degree; algorithm = degree_from_base_fibre, kwargs...)
end

#Getters for EnumerativeProblems should only be implemented if we do not expect to
#  hold a 'knowledge node' for them, since getting knowledge nodes is automatic. 

tracker_options(EP::EnumerativeProblem) = EP.hc_options[:tracker_options]
knowledge(EP::EnumerativeProblem) = EP.knowledge
ambient_dimension(EP::EnumerativeProblem) = length(variables(system(EP)))
n_polynomials(EP::EnumerativeProblem) = length(expressions(system(EP)))
n_parameters(EP::EnumerativeProblem) = length(parameters(system(EP)))
variables(EP::EnumerativeProblem) = variables(system(EP))
parameters(EP::EnumerativeProblem) = parameters(system(EP))
expressions(EP::EnumerativeProblem) = expressions(system(EP))
degree(EP::EnumerativeProblem; kwargs...) =  enumerative_degree(EP; kwargs...)

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
    potential_algorithms = filter(EA->output_property(EA)==EProp,MAIN_ALGORITHMS)
    println("There is a total of ",length(potential_algorithms), " algorithm(s) in Pandora.jl which compute(s) ",name(EProp),":")
    counter=0
    for p in potential_algorithms
        counter=counter+1
        println("      ",counter,") ",name(p))
    end
    if length(potential_algorithms)==0
        return(nothing)
    else
        #Extend to find the algorithm which requires the least amount of EP info which is unknown
        return(potential_algorithms[1])
    end
end

function learn!(EP::EnumerativeProblem, EProp::EnumerativeProperty{T}; 
    algorithm = nothing, kwargs...)      where T

    if algorithm == nothing
        algorithm = find_algorithm(EProp,EP)  #find and algorithm which will obtain it
    end
    
    EA = algorithm
    f = core_function(EA)
    #=
    for i in input_properties(EA)
        if knows(EP,i)==false
            error("Property \n    ["*string(i)*"]\nis required for algorithm \n    ["*name(EA)*"]\nbut is not known")
        end
    end
    =#
    input_knowledge = [get_knowledge!(i,EP) for i in input_properties(EA)]
    input_knowledge_values = [value(kn) for kn in input_knowledge]
    kwargs_to_pass = copy(default_kwargs(EA))
    if !isempty(kwargs)
       for kv in kwargs
           if haskey(kwargs_to_pass,kv[1])
               kwargs_to_pass[kv[1]] = kv[2]
           end
        end
    end
	o = f(input_knowledge_values...;kwargs_to_pass...)
    new_knowledge = KnowledgeNode{T}(EProp,o,input_knowledge,kwargs_to_pass,EA)
	return(know!(EP,new_knowledge))
end

function know!(EP::EnumerativeProblem,K::KnowledgeNode)
    push!(EP.knowledge,K)
    return(K)
end

function know!(EP::EnumerativeProblem, EProp::EnumerativeProperty{T}, value::T) where T
    EA = user_given(EProp,value)
    learn!(EP,EProp; algorithm = EA)
end



############Getter System for EnumerativeProperty of Enumerative Problem#############

function get_knowledge(EProp::EnumerativeProperty,EP::EnumerativeProblem; kwargs...)
    if knows(EP,EProp; kwargs...)==false
        return(nothing)
    end
    candidate_knowledge = filter(k->property(k) == EProp, knowledge(EP))
    kwarg_agreement  = filter(k->knowledge_agrees_with_kwargs(k;kwargs),candidate_knowledge)
    if  length(kwarg_agreement)==1
        return(kwarg_agreement[1])
    else
        println("Warning: the property [",EProp,"] with given keyword arguments has more than one knowledge node.")
        return(kwarg_agreement[1])
    end
end

function get_knowledge_value(EProp::EnumerativeProperty,EP::EnumerativeProblem; kwargs...)
    K = get_knowledge(EProp,EP; kwargs...)
    if K != nothing
        return(value(K))
    else
        return(nothing)
    end
end


function get_knowledge!(EProp::EnumerativeProperty,EP::EnumerativeProblem; kwargs...)
    gk = get_knowledge(EProp,EP; kwargs...)
    if gk != nothing #Check if knowledge of this property is already known
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
    !knows(EP, EProp; kwargs...) && print(NeedToComputeWarning)
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
    if knows(EP,enumerative_degree)
        println(io,"  ",degree(EP),"-to-1")
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
    know!(EP,base_fibre,F)
    learn!(EP,enumerative_degree; algorithm = degree_from_base_fibre)
end

include("solving.jl")