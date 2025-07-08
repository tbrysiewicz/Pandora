#TODO: Put into Pandora format with algorithms etc. 

import Plots: histogram, scatter
export 
    SampleDatum, 
    histogram, 
    scatter

mutable struct SampleDatum
    f::Function
    sampler::Sampler
    function_values::Vector{Any}
end

values(SD::SampleDatum) = SD.function_values

function SampleDatum(EP::EnumerativeProblem, f::Function; sampler::Sampler = UniformSampler{Float64}(n_parameters(EP)), n_samples::Int = 1000, kwargs...)
    return SampleDatum(f, sampler, explore(EP, [f]; n_samples = n_samples, sampler = sampler, kwargs...)[1])
end
function SampleDatum(EP::EnumerativeProblem, F::Vector{Function}; sampler::Sampler = UniformSampler{Float64}(n_parameters(EP)), n_samples::Int = 1000, kwargs...)
    EX = explore(EP, F; n_samples = n_samples, sampler = sampler, kwargs...)
    return [SampleDatum(f, sampler, ex) for (f,ex) in zip(F, EX)]
end

function Base.show(io::IO, SD::SampleDatum)
    print(io, "SampleDatum with function: ", SD.f)
    print(io, ", sampler: ", SD.sampler)
    print(io, ", function values: ", tally(SD.function_values))
end


function histogram(SDVec::Vector{SampleDatum}; kwargs...)
    #Plot a histogram with all SampleDatum function values
    n_vals = [length(sd.function_values) for sd in SDVec]
    n_sample_string = ""
    if length(unique(n_vals)) == 1
        n_sample_string = string(n_vals[1])
    else
        n_sample_string = string(n_vals)
    end
    title = "Histogram for $(length(SDVec)) functions with "*string(n_sample_string)*" samples."
    mylabels = reshape([string(sd.f) for sd in SDVec], 1, :)
    P = Plots.histogram([values(sd) for sd in SDVec]; label = mylabels, kwargs...)    #add title with function name and sampler name
    P = Plots.title!(P, title)
    Plots.xlabel!(P, "Function Values")
    Plots.ylabel!(P, "Frequency")
    display(P)
    return P
end


#visualize sample datum via histogram
function histogram(SD::SampleDatum; kwargs...)
    vals = values(SD)
    Plots.histogram(vals; label = "$(SD.f)", kwargs...)
    #add title with function name and sampler name
    title = "Histogram for $(length(vals)) samples"
    P = title!(title)
    xlabel!("Function Value")
    ylabel!("Frequency")
    display(P)
    return P
end


#visualize sample datuam via scatter plot
function scatter(SD::SampleDatum; kwargs...)
    #check that the function values are 2D
    if length(SD.function_values) != 2
        error("SampleDatum must have exactly two function values for scatter plot.")
    end
    x_values = collect(SD.function_values.values)[1]
    y_values = collect(SD.function_values.values)[2]
    scatter(x_values, y_values; kwargs...)
end