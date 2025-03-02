#Make an EnumerativeProperty for the property you want to compute

#=
function compute_bkk(EP::EnumerativeProblem)
    paths_to_track(specialize(system(EP));only_torus=true)
  end
  
  function compute_affine_bkk(EP::EnumerativeProblem)
    paths_to_track(specialize(system(EP));only_torus=false)
  end
=#


const BKK_BOUND = EnumerativeProperty{Int}("bkk torus bound")
bkk_bound(EP::EnumerativeProblem; kwargs...) = BKK_BOUND(EP; kwargs...)
const AFFINE_BKK_BOUND = EnumerativeProperty{Int}("bkk affine bound")
affine_bkk_bound(EP::EnumerativeProblem; kwargs...) = AFFINE_BKK_BOUND(EP; kwargs...)
const DEGREE_SEQUENCE = EnumerativeProperty{Vector{Int}}("degree sequence")
degree_sequence(EP::EnumerativeProblem; kwargs...) = DEGREE_SEQUENCE(EP; kwargs...)
const BEZOUT_BOUND = EnumerativeProperty{Int}("bezout bound")
bezout_bound(EP::EnumerativeProblem; kwargs...) = BEZOUT_BOUND(EP; kwargs...)
const NEWTON_POLYTOPES = EnumerativeProperty{Vector{Polyhedron}}("newton polytopes")
newton_polytopes(EP::EnumerativeProblem; kwargs...) = NEWTON_POLYTOPES(EP; kwargs...)
  
function compute_bezout(F::System) :: Int64
  return(prod(degree_sequence(F)))
end

function compute_degree_sequence(F::System) :: Vector{Int64}
  Ms = support_coefficients(F)[1]
  deg_seq = [max(map(x->sum(x),eachcol(M))...) for M in Ms]
  return(deg_seq)
end



const COMPUTE_DEGREE_SEQUENCE = EnumerativeAlgorithm(
name = "degree sequence",
input_properties = [SYSTEM],
core_function = compute_degree_sequence,
output_property = DEGREE_SEQUENCE,
epistemic_status = :symbolic
)


push!(MAIN_ALGORITHMS,COMPUTE_DEGREE_SEQUENCE)

const BEZOUT_BOUND = EnumerativeAlgorithm(
name = "bezout bound",
input_properties = [SYSTEM],
core_function = compute_bezout,
output_property = BEZOUT_BOUND,
epistemic_status = :symbolic
)

push!(MAIN_ALGORITHMS,BEZOUT_BOUND)


function compute_bkk(EP::EnumerativeProblem)
  paths_to_track(specialize(system(EP));only_torus=true)
end

function compute_affine_bkk(EP::EnumerativeProblem)
  paths_to_track(specialize(system(EP));only_torus=false)
end

function compute_bezout(EP::EnumerativeProblem)
    return(prod(degree_sequence(EP)))
end

function degree_sequence(EP::EnumerativeProblem)
    Ms = support_coefficients(system(EP))[1]
    deg_seq = [max(map(x->sum(x),eachcol(M))...) for M in Ms]
    return(deg_seq)
end