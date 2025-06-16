

export 
    bkk_bound,
    affine_bkk_bound,
    degree_sequence,
    newton_polytopes,
    bezout_bound,
    support



########## Implementation of Enumerative Property BEZOUT_BOUND#########

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



########## Implementation of Enumerative Property BKK_BOUND#########

const BKK_BOUND = EnumerativeProperty{Int}("bkk_bound")

"""
    bkk_bound(EP::EnumerativeProblem; kwargs...)

Return the bkk_bound of the enumerative problem.
"""
function bkk_bound(EP::EnumerativeProblem; kwargs...)
    BKK_BOUND(EP; kwargs...)
end

function compute_bkk_bound(F::System)::Int
    paths_to_track(specialize(F);only_torus=true)
end

compute_bkk_bound_datum = AlgorithmDatum(
    name = "bkk_bound",
    description = "Computes the BKK bound for the number of isolated solutions in the torus",
    input_properties = [SYSTEM],
    output_property = BKK_BOUND,
    reliability = :certified
)

ALGORITHM_DATA[compute_bkk_bound] = compute_bkk_bound_datum

########## Implementation of Enumerative Property AFFINE_BKK_BOUND#########

const AFFINE_BKK_BOUND = EnumerativeProperty{Int}("affine_bkk_bound")

"""
    affine_bkk_bound(EP::EnumerativeProblem; kwargs...)

Return the affine_bkk_bound of the enumerative problem.
"""
function affine_bkk_bound(EP::EnumerativeProblem; kwargs...)
  AFFINE_BKK_BOUND(EP; kwargs...)
end

function compute_affine_bkk_bound(F::System)::Int
  paths_to_track(specialize(F);only_torus=false)
end

compute_affine_bkk_bound_datum = AlgorithmDatum(
    name = "affine_bkk_bound",
    description = "Computes the bkk bound on affine space (bounding the number of isolated solutions in C^n)",
    input_properties = [SYSTEM],
    output_property = AFFINE_BKK_BOUND,
    reliability = :certified
)

ALGORITHM_DATA[compute_affine_bkk_bound] = compute_affine_bkk_bound_datum


########## Implementation of Enumerative Property DEGREE_SEQUENCE#########

const DEGREE_SEQUENCE = EnumerativeProperty{Vector{Int}}("degree_sequence")

"""
    degree_sequence(EP::EnumerativeProblem; kwargs...)

Return the degree_sequence of the enumerative problem.
"""
function degree_sequence(EP::EnumerativeProblem; kwargs...)
    DEGREE_SEQUENCE(EP; kwargs...)
end

function compute_degree_sequence(F::System)::Vector{Int}    
  G = specialize(F)
  return(degrees(G))
end

compute_degree_sequence_datum = AlgorithmDatum(
    name = "degree_sequence",
    description = "Computes the degree sequence of a system",
    input_properties = [SYSTEM],
    output_property = DEGREE_SEQUENCE,
    reliability = :certified
)

ALGORITHM_DATA[compute_degree_sequence] = compute_degree_sequence_datum


########## Implementation of Enumerative Property SUPPORT#########

const SUPPORT = EnumerativeProperty{Vector{Matrix{Int}}}("support")

"""
    support(EP::EnumerativeProblem; kwargs...)

Return the support of the enumerative problem.
"""
function support(EP::EnumerativeProblem; kwargs...)
    SUPPORT(EP; kwargs...)
end

function compute_support(F::System)::Vector{Matrix{Int}}
    G = specialize(F)
    return([Matrix{Int}(M') for M in support_coefficients(G)[1]])
end

compute_support_datum = AlgorithmDatum(
    name = "support",
    description = "Computes the supports of each expression in terms of matrices whose rows correspond to monomials",
    input_properties = [SYSTEM],
    output_property = SUPPORT,
    reliability = :certified
)

ALGORITHM_DATA[compute_support] = compute_support_datum

########## Implementation of Enumerative Property NEWTON_POLYTOPES#########

const NEWTON_POLYTOPES = EnumerativeProperty{Vector{Polyhedron}}("newton_polytopes")

"""
    newton_polytopes(EP::EnumerativeProblem; kwargs...)

Return the newton_polytopes of the enumerative problem.
"""
function newton_polytopes(EP::EnumerativeProblem; kwargs...)
    NEWTON_POLYTOPES(EP; kwargs...)
end

function compute_newton_polytopes(F::System)::Vector{Polyhedron}
  S = compute_support(F)
  return([convex_hull(s) for s in S])
end

compute_newton_polytopes_datum = AlgorithmDatum(
    name = "newton_polytopes",
    description = "Computes a vector of Newton polytopes of the system",
    input_properties = [SYSTEM],
    output_property = NEWTON_POLYTOPES,
    reliability = :certified
)

ALGORITHM_DATA[compute_newton_polytopes] = compute_newton_polytopes_datum
