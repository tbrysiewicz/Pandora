struct EnumerativeAlgorithm
    name :: String
	inputs::Vector{Function} #a list of functions (getters), which when applied to EP, produce the inputs to the algorithm
	parameters::Vector{Symbol} # a list of non-EP field inputs (e.g. # loops, etc)
	core_function :: Function
	output::Symbol #the EP field which is the output of this algorithm 
	citation::Citation #the citation for this algorithm (if follows from definitions, give citation for definition)
	epistemic_status::Symbol
end

function core_function(EA::EnumerativeAlgorithm)
    EA.core_function
end

function inputs(EA::EnumerativeAlgorithm)
    EA.inputs
end

function output(EA::EnumerativeAlgorithm)
    EA.output
end


function name(EA::EnumerativeAlgorithm)
    EA.name
end

function run_algorithm(EA::EnumerativeAlgorithm,EP::EnumerativeProblem, parameters::Dict{Symbol,Any})
	f = core_function(EA)
	o = f([i(EP) for i in inputs(EA)]...;[s=>parameters[s] for s in keys(parameters)]...)
	data(EP)[output(EA)] = o
    justify(EP,output(EA),"Computed via "*name(EA))
end


