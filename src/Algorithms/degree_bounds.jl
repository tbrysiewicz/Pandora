#Make an EnumerativeProperty for the property you want to compute

#=
function compute_bkk(EP::EnumerativeProblem)
    paths_to_track(specialize(system(EP));only_torus=true)
  end
  
  function compute_affine_bkk(EP::EnumerativeProblem)
    paths_to_track(specialize(system(EP));only_torus=false)
  end
=#


const bkk_bound = EnumerativeProperty{Int}("bkk torus bound")
const affine_bkk_bound = EnumerativeProperty{Int}("bkk affine bound")
const degree_sequence = EnumerativeProperty{Vector{Int}}("degree sequence")
const bezout_bound = EnumerativeProperty{Int}("bezout bound")
const newton_polytopes = EnumerativeProperty{Vector{Polyhedron}}("newton polytopes")
  
function compute_bezout(F::System) :: Int64
  return(prod(degree_sequence(F)))
end

function compute_degree_sequence(F::System) :: Vector{Int64}
  Ms = support_coefficients(F)[1]
  deg_seq = [max(map(x->sum(x),eachcol(M))...) for M in Ms]
  return(deg_seq)
end



const degree_sequence_algorithm = EnumerativeAlgorithm(
name = "degree sequence",
input_properties = [system],
core_function = compute_degree_sequence,
output_property = degree_sequence,
epistemic_status = :symbolic
)


push!(MAIN_ALGORITHMS,degree_sequence_algorithm)

const bezout_bound_algorithm = EnumerativeAlgorithm(
name = "bezout bound",
input_properties = [system],
core_function = compute_bezout,
output_property = bezout_bound,
epistemic_status = :symbolic
)

push!(MAIN_ALGORITHMS,bezout_bound_algorithm)


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