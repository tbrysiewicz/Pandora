module Pandora

using Pkg

const PROJECT_TOML = Pkg.TOML.parsefile(joinpath(@__DIR__, "..", "Project.toml"))
const VERSION_NUMBER = VersionNumber(PROJECT_TOML["version"])

using Dates: today

#Use from HC only
using HomotopyContinuation: TrackerOptions, Result, subs, coefficients, support_coefficients, paths_to_track, degrees, Expression, Variable
#Use and extend from HC
import HomotopyContinuation: expressions, variables, parameters, solve, solutions, evaluate
#Use from HC with intent to export
using HomotopyContinuation: System,  @var
#HC exports
export System, @var

using LinearAlgebra: norm, isapprox, nullspace

using Oscar: Perm, PermGroupElem, PermGroup, symmetric_group, sub,  order, Polyhedron, describe, convex_hull
import Oscar: perm, degree
export gens, order, is_primitive, is_transitive, describe

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
				");

		print("Version")
		printstyled(" $VERSION_NUMBER ", color = :green)
end

include("CoreCode/core_code.jl")
include("degree_bounds.jl")

include("automation.jl")

include("Examples/named_examples.jl")

include("monodromy.jl")

include("solutions.jl")
include("fibres.jl")
include("fibre_visualization.jl")

#include("OldIntroduce/Samplers/Sampler.jl")

#include("OldIntroduce/Optimization/OptimizationGetters.jl")
#include("OldIntroduce/Optimization/Optimizer.jl")
#include("OldIntroduce/Optimization/DietmaierOptimization.jl")

include("samplers.jl")
include("optimization_structs.jl")
#include("optimization.jl")


include("Summarization/summarize.jl")

end # module Pandora


#=
Doc template

@doc raw"""
    b(n)

 Returns
 # Examples
 ```jldoctest

 ```
 """
=#
