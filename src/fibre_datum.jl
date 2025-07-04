export 
    FibreDatum,
    FibreData



const FibreData = Vector{FibreDatum}



#convert a FibreDatum to a Fibre by making the pair (solutions, parameters)
function Base.convert(::Type{Fibre}, F::FibreDatum)
    return Fibre(F.solutions, F.parameters)
end

function Base.convert(::Type{FibreDatum}, F::Fibre)
    return FibreDatum(parameters = F[2], solutions = F[1])
end

function Base.show(io::IO, F::FibreDatum)
    nreal = is_certified(F) ? count(x -> x.real, certificates(F)) : n_real_solutions(solutions(F))
    print(io, "FibreDatum with ", length(F.solutions), " complex solutions (",nreal," real) over ", length(F.parameters), " parameters")
    if is_certified(F) 
        print(io, " [certified].")
    else
        print(io,".")
    end
    if !isempty(F.function_values)
        for k in keys(F.function_values)
            print(io, "  ", k, " = ", F.function_values[k])
        end
    end
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
    S = solutions(solve(F, BF[1]; start_parameters = BF[2], target_parameters = P))
    return(FibreDatum(
        parameters = P,
        solutions = S,
        certificates = certify(F, S; target_parameters = P),
        F = F))
end

compute_fibre_datum_datum = AlgorithmDatum(
    name = "Fibre Datum",
    description = "Computes a single fibre datum for the enumerative problem.",
    input_properties = [SYSTEM, BASE_FIBRE],
    output_property = FIBRE_DATUM,
    reliability = :certified
)

ALGORITHM_DATA[compute_fibre_datum] = compute_fibre_datum_datum


fibre_data(EP::EnumerativeProblem) = [K.value for K in filter(k -> property(k) == FIBRE_DATUM, knowledge(EP))]
