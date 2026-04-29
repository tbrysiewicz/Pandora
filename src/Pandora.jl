module Pandora

using Pkg

const PROJECT_TOML = Pkg.TOML.parsefile(joinpath(@__DIR__, "..", "Project.toml"))
const VERSION_NUMBER = VersionNumber(PROJECT_TOML["version"])

using Dates: today

# Use from HC only
using HomotopyContinuation: TrackerOptions, CertificationResult, Result, subs, coefficients, support_coefficients, paths_to_track, degrees, Expression, Variable, differentiate, monodromy_solve
# Use and extend from HC
import HomotopyContinuation: expressions, variables, parameters, solve, solutions, evaluate, certify, mixed_volume
# Use from HC with intent to export
using HomotopyContinuation: System, @var
# HC exports
export System, @var, mixed_volume

using LinearAlgebra: norm, isapprox, nullspace, I, det

# Use from Oscar only
using Oscar: Perm, PermGroupElem, PermGroup, symmetric_group, sub, Polyhedron, convex_hull, orbits, small_generating_set, minimal_generating_set
using Oscar: PointVector 
# Use and extend from Oscar
import Oscar: perm, degree
# Use from Oscar with intent to export
using Oscar: is_primitive, order, is_transitive, describe, gens, faces, volume, is_natural_alternating_group, is_natural_symmetric_group, minimal_block_reps
#imports as new names
import Oscar: vertices as oscar_vertices, dim as oscar_dim, ambient_dim as oscar_ambient_dim, orbits
# Oscar exports
export gens, order, is_primitive, is_transitive, describe, dim, vertices, faces, volume, orbits, is_natural_alternating_group, is_natural_symmetric_group, minimal_block_reps


using IntegerSmithNormalForm: snf, elementary_divisors
using Combinatorics: combinations

function __init__()
    print(raw"
        0ooo000oo oo oo
     0oo00o0o0o  0o0 o8oo
     0000o0o000o00o0000000
     ooo 0oo 0 oo00ooo00000
          oo00  /o o  ooooo
            \\\//  /////
               \\\////
                 |||/\
                 |||\/
                 |||||
           .....//||||\....
         _______/PANDO(RA)\_____
          ...................
                TCB
    ")
    print("Version")
    printstyled(" $VERSION_NUMBER ", color = :green)
end

include("CoreCode/core_code.jl")
include("degree_bounds.jl")
include("automation.jl")
include("monodromy.jl")
include("enumerative_problem_constructors.jl")
include("fibre_functions.jl")
include("fibre_visualization.jl")
include("samplers.jl")
include("optimization.jl")
include("Summarization/summarize.jl")
include("visualization.jl")
include("certification.jl")
include("fibre_datum.jl")
include("explore.jl")
include("sample_datum.jl")
include("lazy_brute.jl")
include("PandorasBox/alt_burmester.jl")
include("PandorasBox/benchmarks.jl")
include("PandorasBox/conic_lines.jl")
include("PandorasBox/cubic_bisecants.jl")
include("PandorasBox/fano_problems.jl")
include("PandorasBox/kuramoto.jl")
include("PandorasBox/planar_higher_order_contact.jl")
include("PandorasBox/planar_tangency.jl")
include("PandorasBox/rigid_graph_embeddings.jl")
include("PandorasBox/schubert.jl")
include("PandorasBox/sparse_polynomial_systems.jl")
include("PandorasBox/spectrahedral_nodes.jl")
include("PandorasBox/stewart_gough.jl")
include("PandorasBox/variety_projections.jl")
export ALGORITHM_DATA
end