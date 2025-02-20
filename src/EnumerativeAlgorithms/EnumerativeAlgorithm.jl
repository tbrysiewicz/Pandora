struct EnumerativeAlgorithm
    name :: String 
	input_properties::Vector{EnumerativeProperty} 
	parameters :: Dict{Symbol,Any} # a dictionary of non-EP field inputs (e.g. # loops, etc) and the default values
	core_function :: Function
	output_property :: EnumerativeProperty #the EP field which is the output of this algorithm 
	citation :: Citation #the citation for this algorithm (if follows from definitions, give citation for definition)
	epistemic_status :: Symbol
end

function EnumerativeAlgorithm(;
    name :: String  = "unnamed", 	
    input_properties::Union{Vector{EnumerativeProperty},Vector{EnumerativeProperty{T}}} where T  = Vector{EnumerativeProperty}([]),
	parameters :: Dict{Symbol,Any} = Dict{Symbol,Any}(), # a dictionary of non-EP field inputs (e.g. # loops, etc) and the default values
	core_function :: Function = x->nothing,
	output_property :: EnumerativeProperty,  #the EP field which is the output of this algorithm 
	citation :: Citation = NULL_CITATION, #the citation for this algorithm (if follows from definitions, give citation for definition)
	epistemic_status :: Symbol = :null)
    EnumerativeAlgorithm(
        name,
        input_properties,
        parameters,
        core_function,
       output_property,
        citation,
        epistemic_status)
end


function core_function(EA::EnumerativeAlgorithm)
    EA.core_function
end

function input_properties(EA::EnumerativeAlgorithm)
    EA.input_properties
end

function output_property(EA::EnumerativeAlgorithm)
    EA.output_property
end

function parameters(EA::EnumerativeAlgorithm)
    EA.parameters
end

function name(EA::EnumerativeAlgorithm)
    EA.name
end

function epistemic_status(EA::EnumerativeAlgorithm)
    EA.epistemic_status
end

function citation(EA::EnumerativeAlgorithm)
    EA.citation
end

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
    parameters = Dict{Symbol,Any}()
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
        parameters,
        core_function,
        output_property,
        citation,
        epistemic_status))
end


global MAIN_ALGORITHMS = Vector{EnumerativeAlgorithm}([])

