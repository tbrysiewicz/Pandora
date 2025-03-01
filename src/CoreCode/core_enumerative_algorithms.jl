
##############################################################
######## Degree by counting points in the base fibre##########
##############################################################
const degree_from_base_fibre = EnumerativeAlgorithm(
    name = "definition: degree is the number of solutions in a generic base fibre",
    input_properties = [base_fibre],
    core_function = n_solutions,
    output_property = enumerative_degree,
    reliability = :symbolic
)

push!(MAIN_ALGORITHMS,degree_from_base_fibre)

##############################################################
######## Polyhedral homotopy to populate base fibre ##########
##############################################################
function homotopy_continuation_solve(F::System; start_system = :polyhedral) :: Fibre
    k = length(parameters(F))
    P = randn(ComplexF64,k)
    S = solutions(solve(F;target_parameters=P, start_system = start_system))
    return((S,P))
end

const solve_generic_via_start_system = EnumerativeAlgorithm(
    name = "solve generic system via start system using homotopy continuation",
    input_properties = [Pandora.system],
    default_kwargs = Dict{Symbol,Any}(:start_system=>:polyhedral),
    core_function = homotopy_continuation_solve,
    output_property = base_fibre,
    reliability = :numerical
)

push!(MAIN_ALGORITHMS,solve_generic_via_start_system)