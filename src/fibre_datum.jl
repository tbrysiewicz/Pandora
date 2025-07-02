export 
    FibreDatum,
    FibreData

@kwdef mutable struct FibreDatum
    parameters::Vector{ComplexF64} = ComplexF64[]
    solutions::Vector{Vector{ComplexF64}} = Vector{Vector{ComplexF64}}()
    function_values:: Dict{Symbol, Any} = Dict{Symbol, Any}()
    EP::EnumerativeProblem = nothing  # Reference to the EnumerativeProblem this FibreDatum belongs to. Relevant for certification. 
    certified::Bool
end

const FibreData = Vector{FibreDatum}




#convert a FibreDatum to a Fibre by making the pair (solutions, parameters)
function Base.convert(::Type{Fibre}, F::FibreDatum)
    return Fibre(F.solutions, F.parameters)
end

function Base.show(io::IO, F::FibreDatum)
    print(io, "FibreDatum with ", length(F.solutions), " solutions and ", length(F.parameters), " parameters.")
    if !isempty(F.function_values)
        print(io, " Function values: ", F.function_values)
    end
    print(io, " Certified: ", F.certified)
end


#####################################################################
#### Enumerative Properties and Algorithms for Fibre Data ####
#####################################################################
const FIBRE_DATA = EnumerativeProperty{FibreData}("fibre_data")

"""
    fibre_data(EP::EnumerativeProblem)

Return the cached and certified fibre data (a vector of FibreDatum) for the enumerative problem.
"""
function fibre_data(EP::EnumerativeProblem; kwargs...)
    FIBRE_DATA(EP; kwargs...)
end

function compute_fibre_data(F::System,BF::Fibre)::FibreData
    # Implement the actual computation here
    # Example: return [FibreDatum(...), ...]
    error("compute_fibre_data not yet implemented")
end

compute_fibre_data_datum = AlgorithmDatum(
    name = "Fibre Data",
    description = "Computes all fibre data for the enumerative problem.",
    input_properties = [SYSTEM, BASE_FIBRE],
    output_property = FIBRE_DATA,
    reliability = :certified
)

ALGORITHM_DATA[compute_fibre_data] = compute_fibre_data_datum


const FIBRE_DATUM = EnumerativeProperty{FibreDatum}("fibre_datum")

"""
    fibre_datum(EP::EnumerativeProblem)

Return a single fibre datum for the enumerative problem.
"""
function fibre_datum(EP::EnumerativeProblem; kwargs...)
    FIBRE_DATUM(EP; kwargs...)
end

function compute_fibre_datum(F::System)::FibreDatum
    # Implement the actual computation here
    # Example: return FibreDatum(...)
    error("compute_fibre_datum not yet implemented")
end

compute_fibre_datum_datum = AlgorithmDatum(
    name = "Fibre Datum",
    description = "Computes a single fibre datum for the enumerative problem.",
    input_properties = [SYSTEM],
    output_property = FIBRE_DATUM,
    reliability = :certified
)

ALGORITHM_DATA[compute_fibre_datum] = compute_fibre_datum_datum