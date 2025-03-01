
##############################################################
######## Degree by counting points in the base fibre##########
##############################################################
const N_SOLUTIONS = EnumerativeAlgorithm(
    name = "definition: degree is the number of solutions in a generic base fibre",
    input_properties = [base_fibre],
    core_function = n_solutions,
    output_property = enumerative_degree,
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
    name = "solve generic system via start system using homotopy continuation",
    input_properties = [Pandora.system],
    default_kwargs = Dict{Symbol,Any}(:start_system=>:polyhedral),
    core_function = hc_to_generic,
    output_property = base_fibre,
    reliability = :numerical
)

push!(MAIN_ALGORITHMS,HC_TO_GENERIC)