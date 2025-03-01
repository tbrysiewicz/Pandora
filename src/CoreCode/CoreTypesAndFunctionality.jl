export 
    EnumerativeProperty,
    EnumerativeAlgorithm,
    KnowledgeNode


##############################################################
###################    Core Aliases    #######################
##############################################################
const Fibre = Tuple{Vector{Vector{ComplexF64}},Vector{ComplexF64}}
#const RealFibre = Tuple{Vector{Vector{Float64}},Vector{Float64}}

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

function Base.show(io::IO, EA::EnumerativeAlgorithm)
    tenspaces="          "
    println(io,"Algorithm:       ",name(EA))
    println(io,"----------------------------------------------------")
    println(io,"Input:            ")
    if epistemic_status(EA)==:user_given
    print("None (the algorithm always returns the user-given value) \n")
    end
    for i in input_properties(EA)
        println(io,"                    ",i)
    end
    println(io,"Output:           ",output_property(EA))
    println(io,"Citation:         ",citation(EA))
    println(io,"Epistemic Status: ",epistemic_status(EA))
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
    epistemic_status = :user_given
    return(EnumerativeAlgorithm(
        alg_name,
        input_properties,
        default_kwargs,
        core_function,
        output_property,
        citation,
        epistemic_status))
end


global MAIN_ALGORITHMS = Vector{EnumerativeAlgorithm}([])



##############################################################
###################      KnowledgeNode     #######################
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

function Base.show(io::IO, K::KnowledgeNode)
    println(io,"Knowledge node for [",property(K),"]:")
    println(io,value(K))
    println(io,"----------------------------")
    println(io,"Algorithm used: [",name(algorithm(K)),"]")
    if length(input_knowledge(K))>0
        println(io,"Ran on knowledge nodes for ")
    end
    for i in input_knowledge(K)
        println(io,"     ",name(property(i)))
    end
    if length(keys(input_kwargs(K)))>0
        println(io," with input parameters")
        for ip in keys(input_kwargs(K))
            println(io," ",ip,": ",(input_kwargs(K))[ip])
        end
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
    hc_options :: Dict{Any,Any}

    function EnumerativeProblem(F::System; populate = true)
        EP = new()
        EP.system = F
        EP.knowledge = Knowledge([])
        EP.hc_options = Dict{Any,Any}()
        EP.hc_options[:tracker_options]=TrackerOptions()

        know!(EP,system,F)

        if populate
            learn!(EP,base_fibre; algorithm = polyhedral)
            learn!(EP,enumerative_degree; algorithm = degree_from_base_fibre)
        end
        return(EP)
    end
end

#Getters for EnumerativeProblems should only be implemented if we do not expect to
#  hold a 'knowledge node' for them, since getting knowledge nodes is automatic. 

knowledge(EP::EnumerativeProblem) = EP.knowledge
ambient_dimension(EP::EnumerativeProblem) = length(variables(system(EP)))
n_polynomials(EP::EnumerativeProblem) = length(expressions(system(EP)))
n_parameters(EP::EnumerativeProblem) = length(parameters(system(EP)))
variables(EP::EnumerativeProblem) = variables(system(EP))
parameters(EP::EnumerativeProblem) = parameters(system(EP))
expressions(EP::EnumerativeProblem) = expressions(system(EP))



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

function knows(EP::EnumerativeProblem, EProp::EnumerativeProperty)
    count(p -> property(p) == EProp, knowledge(EP)) == 0 ? false : true
end


function learn!(EP::EnumerativeProblem, EProp::EnumerativeProperty{T}; 
    algorithm = nothing, kwargs...)      where T

    if algorithm == nothing
        error("Taylor needs to code how to find best algorithm")
    end
    
    EA = algorithm
    f = core_function(EA)
    for i in default_kwargs(EA)
        if knows(EP,i)==false
            error("Property \n    ["*string(i)*"]\nis required for algorithm \n    ["*name(EA)*"]\nbut is not known")
        end
    end
    input_knowledge = [get_knowledge(i,EP) for i in input_properties(EA)]
    input_knowledge_values = [value(kn) for kn in input_knowledge]
    kwargs_to_pass = default_kwargs(EA)
    if !isempty(kwargs)
       if :k in keys(kwargs...)
           if haskey(:k,kwargs_to_pass)
               kwargs_to_pass[:k] = kwargs[:k]
           end
        end
    end
    flattened_dict = [kwargs_to_pass[k] for k in keys(kwargs_to_pass)]
	o = f(input_knowledge_values...;flattened_dict...)
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
    if knows(EP,EProp)==false
        return(nothing)
    end
    K = knowledge(EP)
    relevant_knowledge = filter(p->property(p)==EProp,K)
    if  length(relevant_knowledge)==1
        return(relevant_knowledge[1])
    else
        println("Warning: the property [",EProp,"] has more than one knowledge node.")
        return(relevant_knowledge[1])
    end
end

function get_knowledge_value(EProp::EnumerativeProperty,EP::EnumerativeProblem; user_parameters :: Dict{Symbol,Any} = Dict{Symbol,Any}())
    K = get_knowledge(EProp,EP; user_parameters=user_parameters)
    if K != nothing
        return(value(K))
    else
        return(nothing)
    end
end

function get_knowledge!(EProp::EnumerativeProperty,EP::EnumerativeProblem; kwargs...)
    gk = get_knowledge(EProp,EP; kwargs...)
    if gk != nothing #Check if knowledge of this property is already known
        return(gk)   #If so, return it
    else             #Otherwise,
        algorithm = find_algorithm(EProp,EP)  #find and algorithm which will obtain it
        if algorithm != nothing
            return(learn!(EP,EProp; algorithm = algorithm, kwargs...))
        end
    end
end

#Collects a knowledge node and then returns its value
get_knowledge_value!(EProp::EnumerativeProperty,EP::EnumerativeProblem; kwargs...)  =  value(get_knowledge!(EProp,EP; kwargs...))

#This is the core "getter" for enumerative properties of enumerative problems
(EProp::EnumerativeProperty)(EP::EnumerativeProblem; kwargs...) = get_knowledge_value!(EProp,EP; kwargs...)



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
    for k in knowledge(EP)
        println(io,"-",property(k))
    end
end