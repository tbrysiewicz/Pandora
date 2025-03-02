
##############################################################
######## Degree by counting points in the base fibre##########
##############################################################
const N_SOLUTIONS = EnumerativeAlgorithm(
    name = "count base fibre solutions",
    input_properties = [BASE_FIBRE],
    core_function = n_solutions,
    output_property = DEGREE,
    reliability = :symbolic
)

push!(MAIN_ALGORITHMS,N_SOLUTIONS)

##############################################################
######## Polyhedral homotopy to populate base fibre ##########
##############################################################
function hc_to_generic(F::System; start_system = :polyhedral) :: Fibre
    k = length(parameters(F))
    P = randn(ComplexF64,k)
    S = solutions(solve(F;target_parameters=P, start_system = start_system))
    return((S,P))
end

const HC_TO_GENERIC = EnumerativeAlgorithm(
    name = "start_system->generic_target_system via HC",
    input_properties = [SYSTEM],
    default_kwargs = Dict{Symbol,Any}(:start_system=>:polyhedral),
    core_function = hc_to_generic,
    output_property = BASE_FIBRE,
    reliability = :numerical
)

push!(MAIN_ALGORITHMS,HC_TO_GENERIC)