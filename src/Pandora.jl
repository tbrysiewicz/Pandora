module Pandora

using Pkg

const PROJECT_TOML = Pkg.TOML.parsefile(joinpath(@__DIR__, "..", "Project.toml"))
const VERSION_NUMBER = VersionNumber(PROJECT_TOML["version"])

#Use from HC only
using HomotopyContinuation: TrackerOptions, Result, subs, coefficients
#Use and extend from HC
import HomotopyContinuation: expressions, variables, parameters, solve, solutions
#Use from HC with intent to export
using HomotopyContinuation: System,  @var
#HC exports
export System, @var

using Oscar: Perm, PermGroupElem, PermGroup, symmetric_group, sub,  order
import Oscar: perm, degree
export gens, order, is_primitive, is_transitive

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

include("automation.jl")

include("Examples/named_examples.jl")

include("monodromy.jl")




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
