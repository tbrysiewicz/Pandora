
const degree_from_base_fibre = EnumerativeAlgorithm(
    name = "Definition of degree",
    input_properties = [base_fibre],
    core_function = n_solutions,
    output_property = enumerative_degree,
    epistemic_status = :symbolic
)

push!(MAIN_ALGORITHMS,degree_from_base_fibre)