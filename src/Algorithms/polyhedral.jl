function base_fibre_via_polyhedral(F::System) :: Fibre
    k = length(parameters(F))
    P = randn(ComplexF64,k)
    S = solutions(solve(F;target_parameters=P))
    return((S,P))
end

const polyhedral = EnumerativeAlgorithm(
    name = "polyhedral homotopy",
    input_properties = [Pandora.system],
    core_function = base_fibre_via_polyhedral,
    output_property = base_fibre,
    epistemic_status = :numerical
)

push!(MAIN_ALGORITHMS,polyhedral)