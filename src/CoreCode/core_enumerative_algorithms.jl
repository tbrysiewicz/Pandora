##############################################################
#         Core Enumerative Algorithms for Pandora.jl         #
##############################################################

"""
    update_base_fibre!(EP::EnumerativeProblem, F::Fibre)

Update the base fibre of the enumerative problem and learn the degree using `n_solutions`.
"""
function update_base_fibre!(EP::EnumerativeProblem, F::Fibre)
    know!(EP, BASE_FIBRE, F)
    learn!(EP, DEGREE; algorithm = n_solutions)
end

##############################################################
#           n_solutions AlgorithmDatum Registration          #
##############################################################

const n_solutions_datum = AlgorithmDatum(
    name = "n_solutions",
    input_properties = [BASE_FIBRE],
    output_property = DEGREE,
    reliability = :certified
)

ALGORITHM_DATA[n_solutions] = n_solutions_datum

##############################################################
#   Implementation of BASE_FIBRE via Polyhedral Homotopy     #
##############################################################

"""
    polyhedral_homotopy(F::System) :: Fibre

Compute the base fibre using polyhedral homotopy.
"""
function polyhedral_homotopy(F::System, G::System) :: Fibre
    k = length(parameters(F))
    P = randn(ComplexF64, k)
    S = solutions(solve(F; target_parameters = P, start_system = :polyhedral))
    if length(expressions(G))>0
        filter!(s -> all(f -> norm(evaluate(f,variables(G)=>s)) <1e-6, expressions(G))==false, S)
    end
    return (S, P)
end

const polyhedral_homotopy_datum = AlgorithmDatum(
    name = "polyhedral_homotopy",
    input_properties = [SYSTEM, INEQUATIONS],
    output_property = BASE_FIBRE,
    reliability = :numerical
)

ALGORITHM_DATA[polyhedral_homotopy] = polyhedral_homotopy_datum

##############################################################
# Implementation of BASE_FIBRE via Total Degree Homotopy     #
##############################################################

"""
    total_degree_homotopy(F::System) :: Fibre

Compute the base fibre using total degree homotopy.
"""
function total_degree_homotopy(F::System, G::System) :: Fibre
    k = length(parameters(F))
    P = randn(ComplexF64, k)
    S = solutions(solve(F; target_parameters = P, start_system = :total_degree))
    if length(expressions(G))>0
        filter!(s -> all(f -> norm(evaluate(f,variables(G)=>s)) <1e-6, expressions(G))==false, S)
    end
    return (S, P)
end

const total_degree_homotopy_datum = AlgorithmDatum(
    name = "total degree homotopy",
    input_properties = [SYSTEM, INEQUATIONS],
    output_property = BASE_FIBRE,
    reliability = :numerical
)