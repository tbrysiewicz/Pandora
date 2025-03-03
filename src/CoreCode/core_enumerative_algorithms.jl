
##############################################################
######## Degree by counting points in the base fibre##########
##############################################################
const n_solutions_datum = AlgorithmDatum(
    name = "N_SOLUTIONS",
    input_properties = [BASE_FIBRE],
    output_property = DEGREE,
    reliability = :symbolic
)

ALGORITHM_DATA[n_solutions]=n_solutions_datum

##############################################################
######## Polyhedral homotopy to populate base fibre ##########
##############################################################
function polyhedral_homotopy(F::System) :: Fibre
    k = length(parameters(F))
    P = randn(ComplexF64,k)
    S = solutions(solve(F;target_parameters=P, start_system = :polyhedral))
    return((S,P))
end

const polyhedral_homotopy_datum = AlgorithmDatum(
    name = "polyhedral homotopy",
    input_properties = [SYSTEM],
    output_property = BASE_FIBRE,
    reliability = :numerical
)

ALGORITHM_DATA[polyhedral_homotopy]=polyhedral_homotopy_datum



##############################################################
####### Total degree homotopy to populate base fibre #########
##############################################################
function total_degree_homotopy(F::System) :: Fibre
    k = length(parameters(F))
    P = randn(ComplexF64,k)
    S = solutions(solve(F;target_parameters=P, start_system = :total_degree))
    return((S,P))
end

const total_degree_homotopy_datum = AlgorithmDatum(
    name = "total degree homotopy",
    input_properties = [SYSTEM],
    output_property = BASE_FIBRE,
    reliability = :numerical
)

ALGORITHM_DATA[total_degree_homotopy]=total_degree_homotopy_datum



