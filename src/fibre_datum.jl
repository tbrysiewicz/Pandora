export 
    FibreDatum,
    FibreData

@kwdef mutable struct FibreDatum
    parameters::Vector{ComplexF64} = ComplexF64[]
    solutions::Vector{Vector{ComplexF64}} = Vector{Vector{ComplexF64}}()
    function_values:: Dict{Symbol, Any} = Dict{Symbol, Any}()
    certificates::Union{CertificationResult,Nothing} = nothing 
end

const FibreData = Vector{FibreDatum}

function is_certified(F::FibreDatum)
    count(x->x.certified, C.certificates)==length(solutions)
end



#convert a FibreDatum to a Fibre by making the pair (solutions, parameters)
function Base.convert(::Type{Fibre}, F::FibreDatum)
    return Fibre(F.solutions, F.parameters)
end

function Base.convert(::Type{FibreDatum}, F::Fibre)
    return FibreDatum(parameters = F[2], solutions = F[1])
end

function Base.show(io::IO, F::FibreDatum)
    print(io, "FibreDatum with ", length(F.solutions), " solutions and ", length(F.parameters), " parameters.")
    if !isempty(F.function_values)
        print(io, " Function values: ", F.function_values)
    end
    print(io, " Certified: ", is_certified(F)))
end



const FIBRE_DATUM = EnumerativeProperty{FibreDatum}("fibre_datum")

"""
    fibre_datum(EP::EnumerativeProblem; real = false, kwargs...)

Return a single fibre datum for the enumerative problem.
"""
function fibre_datum(EP::EnumerativeProblem; real = false, kwargs...)
    FIBRE_DATUM(EP; real = false, kwargs...)
end

function compute_fibre_datum(F::System, BF::Fibre; real = false, kwargs...)::FibreDatum
    field = real ? Float64 : ComplexF64
    P = randn(field, n_parameters(F))
    S = solve(F, BF[1]; start_parameters = BF[2], target_parameters = P)
    return(FibreDatum(
        parameters = P,
        solutions = S,
        certificates = certify(F, S; target_parameters = P)))
end

compute_fibre_datum_datum = AlgorithmDatum(
    name = "Fibre Datum",
    description = "Computes a single fibre datum for the enumerative problem.",
    input_properties = [SYSTEM, BASE_FIBRE],

    output_property = FIBRE_DATUM,
    reliability = :certified
)

ALGORITHM_DATA[compute_fibre_datum] = compute_fibre_datum_datum