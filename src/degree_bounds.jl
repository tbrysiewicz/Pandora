

export 
#    bkk_bound,
#    affine_bkk_bound,
#    degree_sequence,
#    newton_polytopes,
    bezout_bound




const BEZOUT_BOUND = EnumerativeProperty{Int}("bezout_bound")

"""
    bezout_bound(EP::EnumerativeProblem; kwargs...)

Return the Bézout bound of the enumerative problem — the product of the total degrees of the input polynomials.
"""
function bezout_bound(EP::EnumerativeProblem; kwargs...)
    BEZOUT_BOUND(EP; kwargs...)
end

# Core computation function takes only the input property System
function compute_bezout_bound(F::System)::Int
    G = specialize(F)
    Ms = support_coefficients(G)[1]
    deg_seq = [maximum(sum, eachcol(M)) for M in Ms]
    prod(deg_seq)
end

compute_bezout_bound_datum = AlgorithmDatum(
    name = "Bezout Bound",
    description = "Product of total degrees of the input polynomials",
    input_properties = [SYSTEM],
    output_property = BEZOUT_BOUND,
    reliability = :certified
)

ALGORITHM_DATA[compute_bezout_bound] = compute_bezout_bound_datum





#=



const BKK_BOUND = EnumerativeProperty{Int}("bkk (torus bound)")
bkk_bound(EP::EnumerativeProblem; kwargs...) = BKK_BOUND(EP; kwargs...)

const AFFINE_BKK_BOUND = EnumerativeProperty{Int}("bkk (affine bound)")
affine_bkk_bound(EP::EnumerativeProblem; kwargs...) = AFFINE_BKK_BOUND(EP; kwargs...)

const DEGREE_SEQUENCE = EnumerativeProperty{Vector{Int}}("degree sequence")
degree_sequence(EP::EnumerativeProblem; kwargs...) = DEGREE_SEQUENCE(EP; kwargs...)

const NEWTON_POLYTOPES = EnumerativeProperty{Vector{Polyhedron}}("newton polytopes")
newton_polytopes(EP::EnumerativeProblem; kwargs...) = NEWTON_POLYTOPES(EP; kwargs...)



#####################Compute Degree Sequence###########################
const compute_degree_sequence_datum = AlgorithmDatum(
name = "COMPUTE_DEGREE_SEQUENCE",
input_properties = [SYSTEM],
output_property = DEGREE_SEQUENCE,
reliability = :symbolic
)

ALGORITHM_DATA[compute_degree_sequence]=compute_degree_sequence_datum

function compute_degree_sequence(F::System) :: Vector{Int64}
  Ms = support_coefficients(F)[1]
  deg_seq = [max(map(x->sum(x),eachcol(M))...) for M in Ms]
  return(deg_seq)
end




#####################Compute Bezout Bound###########################
const compute_bezout_bound_datum = AlgorithmDatum(
name = "compute bezout bound",
input_properties = [SYSTEM],
core_function = compute_bezout_bound,
output_property = BEZOUT_BOUND,
reliability = :symbolic
)

ALGORITHM_DATA[compute_bezout_bound]=compute_bezout_bound_datum

function compute_bezout(F::System) :: Int64
  return(prod(degree_sequence(F)))
end


function compute_bkk(EP::EnumerativeProblem)
  paths_to_track(specialize(system(EP));only_torus=true)
end

function compute_affine_bkk(EP::EnumerativeProblem)
  paths_to_track(specialize(system(EP));only_torus=false)
end

function compute_bezout_bound(EP::EnumerativeProblem)
    return(prod(degree_sequence(EP)))
end

function degree_sequence(EP::EnumerativeProblem)
    Ms = support_coefficients(system(EP))[1]
    deg_seq = [max(map(x->sum(x),eachcol(M))...) for M in Ms]
    return(deg_seq)
end

=#